// charge_flipping_worker.js
// version 138, 2 august 2026 -- GPU (WebGPU) path with CPU fallback,
// space group applied inside the iteration.
//
// Dual-space charge flipping (Oszlanyi & Suto, Acta Cryst A60 (2004) 134)
// run on Pawley-extracted integrated intensities.
//
// ---------------------------------------------------------------------------
//  What v137 got wrong
// ---------------------------------------------------------------------------
//  1. Neither path ran. The WGSL declared two different resources on
//     @binding(4) and @binding(5), which is a validation error, so the GPU
//     path always threw; and runChargeFlippingCPU still read prep.observedIdx
//     and prep.fObs, which buildFobsGrid had stopped returning, so the
//     fallback threw a TypeError on its first line.
//
//  2. The symmetry used for the orbit expansion was the HOLOHEDRY OF THE
//     CRYSTAL SYSTEM, not the Laue class of the space group. P4 (Laue 4/m,
//     multiplicity 8) was expanded with 4/mmm and got 16 grid points; the
//     measured intensity was then divided among eight positions that are not
//     equivalent to the other eight. Same for -3 vs -3m, 6/m vs 6/mmm and
//     m-3 vs m-3m. Here the orbit is generated from the actual symmetry
//     operators, so it is right by construction; the Laue tables below are
//     only a fallback for space-group JSON that predates the "symops" field.
//
//  3. Systematic absences were left free. Charge flipping is then allowed to
//     put density on reflections the lattice forbids -- for an F or I lattice
//     that is most of the reciprocal grid, and it lets the map break the
//     centring outright. They are now a hard zero every cycle.
//
//  4. The intensity repartition was inverted. v137 shared each observed
//     intensity freely among the members of one SYMMETRY ORBIT. Those members
//     are required by symmetry to have equal |F|; letting them float is
//     precisely what throws the space-group information away. The free
//     partition belongs between DISTINCT reflections that overlap in
//     2-theta, which v137 never handled at all. Both are now modelled:
//     an orbit gets equal amplitudes, a 2-theta cluster gets a Le Bail-style
//     split proportional to the current calculated intensities.
//
//  5. No Lorentz-polarisation correction. A Pawley intensity is a peak AREA,
//     which is proportional to m * LP * |F|^2. Without dividing by LP the
//     low-angle reflections are weighted several times too heavily and the
//     map is dominated by the first few peaks. The polarisation half of LP
//     is now whatever the host declares (none / lab / synchrotron), taken
//     per reflection from peak.lp when the host supplies it, so the map and
//     the refinement report cannot disagree about it.
//
//  6. Symmetry was applied only afterwards, by findOriginShift +
//     symmetryAverage. That is the orthodox choice and it is still available
//     (set symLambda = 0), but with symLambda > 0 the group is imposed on the
//     structure factors every cycle via F(hR) = F(h) exp(-2 pi i h.t), which
//     converges faster and fixes the origin as a side effect.
//
// ---------------------------------------------------------------------------
//  Everything below is shared by the GPU and CPU paths: the reflection model
//  is built once, in JavaScript, and the two backends consume the same flat
//  arrays. That is deliberate -- it is the only way to keep the fallback
//  numerically comparable to the GPU result.
// ---------------------------------------------------------------------------

'use strict';

// ---------------------------------------------------------------------------
//  Shared modules. constants.js and crystal.js are the single source of truth
//  for the physical constants and the lattice geometry; the copies of
//  metricTensor / fracDistance / cellVolume that used to live in this file
//  have been removed in favour of them.
//
//  This file is also loaded directly by node for its unit tests, where
//  importScripts does not exist and a required module's declarations are
//  module-scoped rather than global -- hence the two-branch loader.
// ---------------------------------------------------------------------------
if (typeof importScripts === 'function') {
    importScripts('constants.js', 'crystal.js');
    // Import comp-Harker modules for Wyckoff composition-driven search
    importScripts('symmetry_utils.js', 'observations.js', 'wyckoff_assign.js', 'scatterers.js', 'contacts.js', 'swarm_wyckoff.js', 'coord_refine.js');
} else if (typeof require === 'function' && typeof module !== 'undefined') {

    Object.assign(globalThis, require('./constants.js'), require('./crystal.js'));
}

// ---------------------------------------------------------------------------
//  1D complex FFT, radix-2, in place, on a strided view of a flat array.
//
//  ON THE TWIDDLE RECURRENCE (curRe/curIm updated by complex multiplication
//  rather than read from a precomputed table).
//
//  This is the classic textbook criticism of a radix-2 implementation, and it
//  does not apply here. It was measured, both ways, before being left alone.
//  Numbers from V8 on a 64^3 grid, against a reference DFT in double
//  precision (test/test_fft.js reproduces all of them):
//
//                          rel. error vs DFT        64^3 forward transform
//      n = 128    recurrence  2.9e-14                       --
//                 table       2.9e-14                       --
//      n = 1024   recurrence  2.7e-13                       --
//                 table       2.7e-13                       --
//      64^3       recurrence      --                     23.4 ms
//                 table           --                     26.8 ms
//                 table, per-direction (no sign multiply) 27.3 ms
//
//  ACCURACY: identical to two significant figures at every size tested. The
//  error is dominated by the O(log n) accumulation through the butterflies,
//  not by the twiddle recurrence, so removing the recurrence removes nothing
//  measurable. At 3e-14 relative it is twelve orders of magnitude below the
//  counting statistics on the intensities being transformed; claims of
//  "visible phase noise at 128^3" do not survive contact with a measurement.
//
//  SPEED: the table is 12-15% SLOWER. The recurrence keeps its two doubles in
//  registers across the inner loop, whereas the table costs two loads per
//  butterfly. That is the opposite of the usual C intuition and it is why this
//  was measured rather than assumed.
//
//  So: leave it. If you are here because a review flagged the recurrence, run
//  `node test/test_fft.js` before changing anything.
// ---------------------------------------------------------------------------

/**
 * In-place radix-2 complex FFT over a strided view of a flat array.
 * @param {Float64Array} re     Real parts, modified in place.
 * @param {Float64Array} im     Imaginary parts, modified in place.
 * @param {number} base         Index of element 0 of this transform.
 * @param {number} stride       Index step between consecutive elements.
 * @param {number} n            Transform length; MUST be a power of two.
 * @param {boolean} inverse     True for the inverse transform (UNSCALED; the
 *                              1/N^3 normalisation is applied once in fft3d).
 * @returns {void}
 */
function fft1d(re, im, base, stride, n, inverse) {
    for (let i = 1, j = 0; i < n; i++) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            const a = base + i * stride, b = base + j * stride;
            let t = re[a]; re[a] = re[b]; re[b] = t;
            t = im[a]; im[a] = im[b]; im[b] = t;
        }
    }

    const sign = inverse ? 1 : -1;
    for (let len = 2; len <= n; len <<= 1) {
        const ang = sign * 2 * Math.PI / len;
        const wRe = Math.cos(ang), wIm = Math.sin(ang);
        for (let i = 0; i < n; i += len) {
            let curRe = 1, curIm = 0;
            const half = len >> 1;
            for (let k = 0; k < half; k++) {
                const iu = base + (i + k) * stride;
                const iv = base + (i + k + half) * stride;
                const uRe = re[iu], uIm = im[iu];
                const vRe = re[iv] * curRe - im[iv] * curIm;
                const vIm = re[iv] * curIm + im[iv] * curRe;
                re[iu] = uRe + vRe; im[iu] = uIm + vIm;
                re[iv] = uRe - vRe; im[iv] = uIm - vIm;
                const nRe = curRe * wRe - curIm * wIm;
                curIm = curRe * wIm + curIm * wRe;
                curRe = nRe;
            }
        }
    }
}

// ---------------------------------------------------------------------------
//  3D FFT on an N x N x N grid stored as idx = h + k*N + l*N*N.
// ---------------------------------------------------------------------------
function fft3d(re, im, N, inverse) {
    const N2 = N * N;
    for (let l = 0; l < N; l++)
        for (let k = 0; k < N; k++)
            fft1d(re, im, k * N + l * N2, 1, N, inverse);
    for (let l = 0; l < N; l++)
        for (let h = 0; h < N; h++)
            fft1d(re, im, h + l * N2, N, N, inverse);
    for (let k = 0; k < N; k++)
        for (let h = 0; h < N; h++)
            fft1d(re, im, h + k * N, N2, N, inverse);

    if (inverse) {
        const s = 1 / (N * N2);
        for (let i = 0; i < re.length; i++) { re[i] *= s; im[i] *= s; }
    }
}

// ---------------------------------------------------------------------------
//  Reciprocal-space rotations for the 11 Laue classes.
//
//  ONLY A FALLBACK. When the space-group JSON carries the "symops" field the
//  orbit is generated from the real operators instead, which also supplies
//  the translations that the symmetrisation and the absence test need. These
//  tables are used when symops are missing, and they are keyed on the LAUE
//  CLASS, not on the crystal system -- that distinction is the whole point.
//
//  Each operator maps (h,k,l) -> (h',k',l'). Friedel mates are added
//  separately, so only proper rotations are listed. Note that every lower
//  class is a prefix of its holohedry, which is why the arrays are built by
//  slicing.
// ---------------------------------------------------------------------------
const OPS_4MMM = [
    (h, k, l) => [h, k, l],
    (h, k, l) => [-k, h, l],
    (h, k, l) => [-h, -k, l],
    (h, k, l) => [k, -h, l],
    (h, k, l) => [-h, k, -l],
    (h, k, l) => [h, -k, -l],
    (h, k, l) => [k, h, -l],
    (h, k, l) => [-k, -h, -l]
];

const OPS_6MMM = [
    (h, k, l) => [h, k, l],
    (h, k, l) => [-k, h + k, l],
    (h, k, l) => [-h - k, h, l],
    (h, k, l) => [-h, -k, l],
    (h, k, l) => [k, -h - k, l],
    (h, k, l) => [h + k, -h, l],
    (h, k, l) => [-k, -h, -l],
    (h, k, l) => [-h, h + k, -l],
    (h, k, l) => [h + k, -k, -l],
    (h, k, l) => [k, h, -l],
    (h, k, l) => [h, -h - k, -l],
    (h, k, l) => [-h - k, k, -l]
];

const OPS_M3M = [
    (h, k, l) => [h, k, l],     (h, k, l) => [-h, -k, l],
    (h, k, l) => [-h, k, -l],   (h, k, l) => [h, -k, -l],
    (h, k, l) => [l, h, k],     (h, k, l) => [l, -h, -k],
    (h, k, l) => [-l, -h, k],   (h, k, l) => [-l, h, -k],
    (h, k, l) => [k, l, h],     (h, k, l) => [-k, l, -h],
    (h, k, l) => [k, -l, -h],   (h, k, l) => [-k, -l, h],
    (h, k, l) => [k, h, -l],    (h, k, l) => [-k, -h, -l],
    (h, k, l) => [k, -h, l],    (h, k, l) => [-k, h, l],
    (h, k, l) => [h, l, -k],    (h, k, l) => [-h, l, k],
    (h, k, l) => [-h, -l, -k],  (h, k, l) => [h, -l, k],
    (h, k, l) => [l, k, -h],    (h, k, l) => [l, -k, h],
    (h, k, l) => [-l, k, h],    (h, k, l) => [-l, -k, -h]
];

const OPS_3BARM = [                       // 32, hexagonal axes
    (h, k, l) => [h, k, l],
    (h, k, l) => [-h - k, h, l],
    (h, k, l) => [k, -h - k, l],
    (h, k, l) => [-k, -h, -l],
    (h, k, l) => [-h, h + k, -l],
    (h, k, l) => [h + k, -k, -l]
];

const LAUE_OPS = {
    '-1':    OPS_M3M.slice(0, 1),
    '2/m':   [(h, k, l) => [h, k, l], (h, k, l) => [-h, k, -l]],   // unique axis b
    'mmm':   OPS_M3M.slice(0, 4),
    '4/m':   OPS_4MMM.slice(0, 4),
    '4/mmm': OPS_4MMM,
    '-3':    OPS_3BARM.slice(0, 3),
    '-3m':   OPS_3BARM,
    '6/m':   OPS_6MMM.slice(0, 6),
    '6/mmm': OPS_6MMM,
    'm-3':   OPS_M3M.slice(0, 12),
    'm-3m':  OPS_M3M
};

// Crystal system -> holohedry, retained only so that a job that supplies
// neither symops nor a Laue class still runs (with a warning).
const SYSTEM_HOLOHEDRY = {
    triclinic: '-1',
    monoclinic: '2/m',
    orthorhombic: 'mmm',
    tetragonal: '4/mmm',
    trigonal: '-3m',
    rhombohedral: '-3m',
    hexagonal: '6/mmm',
    cubic: 'm-3m'
};

// Centring translations, used for the absence test when symops are missing.
const CENTRING_VECTORS = {
    P: [],
    A: [[0, 0.5, 0.5]],
    B: [[0.5, 0, 0.5]],
    C: [[0.5, 0.5, 0]],
    I: [[0.5, 0.5, 0.5]],
    F: [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    R: [[2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]],
    H: [[2 / 3, 1 / 3, 0], [1 / 3, 2 / 3, 0]]
};

// ---------------------------------------------------------------------------
//  MOVED. metricTensor, fracDistance and cellVolume now live in crystal.js,
//  loaded above. They were duplicated here and in powder5.html / density3d.js,
//  each copy with its own idea of how to default a missing cell angle.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//  Deterministic PRNG, so a run can be reproduced from its reported seed.
// ---------------------------------------------------------------------------
function mulberry32(seed) {
    let t = seed >>> 0;
    return function () {
        t = (t + 0x6D2B79F5) | 0;
        let r = Math.imul(t ^ (t >>> 15), 1 | t);
        r = (r + Math.imul(r ^ (r >>> 7), 61 | r)) ^ r;
        return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
}

// ===========================================================================
//  SYMMETRY PLUMBING
// ===========================================================================

// (hR)_j = sum_i h_i R_ij, with r row-major (r[3*i+j] = R_ij).
function applyRowVector(h, r) {
    return [
        h[0] * r[0] + h[1] * r[3] + h[2] * r[6],
        h[0] * r[1] + h[1] * r[4] + h[2] * r[7],
        h[0] * r[2] + h[1] * r[5] + h[2] * r[8]
    ];
}

// Accept the operator list as it comes out of the space-group JSON, drop
// anything malformed, and remove exact duplicates.
function normalizeSymops(raw) {
    if (!Array.isArray(raw) || raw.length === 0) return null;
    const out = [];
    const seen = new Set();
    for (const op of raw) {
        if (!op || !Array.isArray(op.r) || op.r.length !== 9) continue;
        const r = op.r.map(Number);
        let t = Array.isArray(op.t) ? op.t.map(Number) : [0, 0, 0];
        if (op.t_num && op.t_den) t = op.t_num.map(n => Number(n) / Number(op.t_den));
        if (r.some(v => !Number.isFinite(v)) || t.some(v => !Number.isFinite(v))) continue;
        // Fold translations into [0,1) so the dedup key is canonical.
        const tw = t.map(v => { const f = v - Math.floor(v); return Math.abs(f - 1) < 1e-9 ? 0 : f; });
        const key = r.join(',') + '|' + tw.map(v => Math.round(v * 5040)).join(',');
        if (seen.has(key)) continue;
        seen.add(key);
        out.push({ r, t: tw, xyz: op.xyz || '' });
    }
    return out.length ? out : null;
}

// Turn either representation into a uniform list of
//   { apply(h,k,l) -> [h,k,l],  t: [tx,ty,tz] }
// The translation is what makes the symmetrisation phase factor meaningful;
// with the Laue fallback it is zero everywhere, which degrades symmetrize()
// to a plain orbit average.
function buildHklOps(symops, laueClass, system) {
    if (symops) {
        return {
            ops: symops.map(op => ({
                apply: (h, k, l) => applyRowVector([h, k, l], op.r),
                t: op.t
            })),
            source: 'symops',
            laue: laueClass || null
        };
    }
    const key = LAUE_OPS[laueClass] ? laueClass
              : (SYSTEM_HOLOHEDRY[String(system || '').toLowerCase()] || '-1');
    return {
        ops: LAUE_OPS[key].map(f => ({ apply: f, t: [0, 0, 0] })),
        source: 'laue',
        laue: key
    };
}

// ---------------------------------------------------------------------------
//  Systematic absences, straight from the operators.
//
//  h is extinct iff some operator (R, t) satisfies hR = h with h.t not an
//  integer. That single test covers centring (R = I, t = the centring
//  vector), screw axes and glide planes at once, so no reflection-condition
//  parser and no space-group database is needed inside the worker.
//
//  Cost is N^3 times the number of DISTINCT rotations, which is at most 48;
//  the translations attached to a rotation are only examined for the h that
//  the rotation actually fixes.
// ---------------------------------------------------------------------------
function findSystematicAbsences(N, symops, centring) {
    const N2 = N * N;
    const out = [];

    // Group translations by rotation.
    let rots;
    if (symops) {
        const byRot = new Map();
        for (const op of symops) {
            const key = op.r.join(',');
            if (!byRot.has(key)) byRot.set(key, { r: op.r, ts: [] });
            byRot.get(key).ts.push(op.t);
        }
        rots = [...byRot.values()];
    } else {
        const vecs = CENTRING_VECTORS[String(centring || 'P').toUpperCase()] || [];
        if (vecs.length === 0) return { absent: Uint32Array.from([]), tested: false };
        rots = [{ r: [1, 0, 0, 0, 1, 0, 0, 0, 1], ts: vecs }];
    }

    // Rotations that fix nothing but the origin can be skipped only after the
    // fixed-point test, which is cheap, so just run them all.
    const half = N >> 1;
    for (let ml = 0; ml < N; ml++) {
        const l = ml > half ? ml - N : ml;
        for (let mk = 0; mk < N; mk++) {
            const k = mk > half ? mk - N : mk;
            for (let mh = 0; mh < N; mh++) {
                const h = mh > half ? mh - N : mh;
                if (h === 0 && k === 0 && l === 0) continue;

                let extinct = false;
                for (let ri = 0; ri < rots.length && !extinct; ri++) {
                    const r = rots[ri].r;
                    // hR == h ?
                    if (h * r[0] + k * r[3] + l * r[6] !== h) continue;
                    if (h * r[1] + k * r[4] + l * r[7] !== k) continue;
                    if (h * r[2] + k * r[5] + l * r[8] !== l) continue;
                    const ts = rots[ri].ts;
                    for (let ti = 0; ti < ts.length; ti++) {
                        const t = ts[ti];
                        const ht = h * t[0] + k * t[1] + l * t[2];
                        if (Math.abs(ht - Math.round(ht)) > 1e-6) { extinct = true; break; }
                    }
                }
                if (extinct) out.push(mh + mk * N + ml * N2);
            }
        }
    }
    return { absent: Uint32Array.from(out), tested: true };
}

// ---------------------------------------------------------------------------
//  Lorentz-polarisation for a powder scan.
//
//  A Pawley intensity is the peak AREA, which carries m * LP * |F|^2. The
//  multiplicity m is absorbed by spreading the intensity over the m grid
//  points of the orbit, but LP has to be divided out explicitly or the
//  low-angle reflections dominate the map.
//
//  This used to hard-wire the unpolarised Bragg-Brentano case (with an
//  optional monochromator angle the host never actually sent), which meant a
//  synchrotron dataset was silently corrected as if it came off a sealed
//  tube. The host now chooses the model; this function speaks the same
//  parameterisation it does:
//
//      L(2th) = 1 / ( sin^2(th) . cos(th) )
//      P(2th) = ( 1 + K . cos^2(2th) ) / ( 1 + K )     P(0) = 1 for every K
//
//        K = 0                 no polarisation term (neutrons, or a beam
//                              polarised perpendicular to a vertical
//                              scattering plane)
//        K = 1                 unpolarised lab source
//        K = cos^2(2th_M)      lab source behind a monochromator
//        K = (1-f)/f           synchrotron, f = perpendicular-polarised
//                              fraction
//
//  `pol` may be:
//    * a descriptor { K } (or { mode, monoTth, polFrac }) from the host,
//    * a bare number, read as a monochromator 2-theta for backward
//      compatibility with the old two-argument call,
//    * absent, which reproduces the historical unpolarised behaviour.
// ---------------------------------------------------------------------------
function polarizationK(pol) {
    // Delegates to crystal.js, which every entry point loads. Kept as a name
    // here because callers in this file and in the host use it; the body is
    // not duplicated, because three copies of this arithmetic is how the
    // K-versus-descriptor confusion started.
    return sharkoPolarizationK(pol);
}

/**
 * Lp from 2-theta and a polarisation DESCRIPTOR (object, or a bare number read
 * as a legacy monochromator 2-theta).
 *
 * If you already hold the ratio K, call sharkoLorentzPolarization(tth, K)
 * instead. Passing K here would be wrong in a way that is easy to miss: K = 0,
 * which is what a neutron pattern needs, would be read as a monochromator
 * angle of zero degrees and fall back to K = 1.
 */
function lorentzPolarization(tthDeg, pol) {
    return sharkoLp(tthDeg, pol);
}

function describePolarization(pol) {
    const K = polarizationK(pol);
    const mode = (pol && typeof pol === 'object' && pol.mode) ? String(pol.mode) : 'lab';
    // Mode first, K second. A fully perpendicular-polarised synchrotron beam
    // also gives K = 0, and calling that "none" would misreport the geometry
    // even though the arithmetic is the same.
    if (mode === 'none') return 'none (Lorentz factor only)';
    if (mode === 'synchrotron') {
        const f = Number(pol.polFrac);
        return `synchrotron, f = ${Number.isFinite(f) ? f.toFixed(3) : '?'} (K = ${K.toFixed(4)})`;
    }
    if (K === 0) return 'none (Lorentz factor only)';
    const m = (pol && typeof pol === 'object') ? Number(pol.monoTth) : NaN;
    return (Number.isFinite(m) && m > 0)
        ? `lab tube, monochromator 2th_M = ${m.toFixed(2)} deg (K = ${K.toFixed(4)})`
        : `lab tube, unpolarised (K = ${K.toFixed(4)})`;
}

// ===========================================================================
//  THE REFLECTION MODEL
//
//  Produces the flat tables that both backends iterate on:
//
//    orbits[g]      = { start, count, targetI }   symmetry orbit
//    orbitIdx[p]    = grid index | (conjugate flag << 31)
//    phase[2p,2p+1] = exp(+2 pi i h0.t) for that member
//    clusters[c]    = { orbitStart, nOrbits, obsI }  2-theta overlap group
//    absent[]       = grid indices forced to zero
//
//  Orbits are emitted in cluster order so that a cluster is a contiguous
//  range, which is what lets one GPU thread handle a whole cluster.
// ===========================================================================
function buildReflectionModel(job, N) {
    const N2 = N * N, N3 = N2 * N;
    const wrap = v => ((v % N) + N) % N;

    const symops = normalizeSymops(job.symops);
    const built = buildHklOps(symops, job.laueClass, job.system);
    const hklOps = built.ops;

    const applyLP = job.applyLP !== false;
    // The host may send a full polarisation descriptor; `monochromatorTth`
    // stays supported so an older caller keeps working unchanged.
    const polModel = (job.polarization !== undefined && job.polarization !== null)
        ? job.polarization : job.monochromatorTth;
    let lpFromHost = 0, lpFromModel = 0;
    const tolTth = Number.isFinite(job.overlapTolTth) ? job.overlapTolTth : 0.05;

    let droppedZero = 0, droppedRange = 0, droppedDuplicate = 0;
    let multiplicityMismatch = 0, multiplicityExample = '';
    let maxIndexSeen = 0;
    const claimed = new Map();      // grid index -> orbit ordinal, for collisions

    // --- pass 1: build one orbit per unique reflection ----------------------
    const refs = [];
    for (const peak of job.hklList || []) {
        if (!peak) continue;
        const I = Number(peak.intensity);
        // A ZERO is kept. Dropping it would tell the algorithm the reflection
        // was never measured, which leaves it free to take any value the map
        // finds convenient; keeping it with a target of zero says what the
        // data actually say. Only a missing or negative value is refused --
        // the caller is expected to have applied a French-Wilson correction,
        // so anything still negative here is malformed input.
        if (!Number.isFinite(I) || I < 0) { droppedZero++; continue; }

        const h0 = peak.h_orig, k0 = peak.k_orig, l0 = peak.l_orig;
        if (!Number.isFinite(h0) || !Number.isFinite(k0) || !Number.isFinite(l0)) continue;
        if (h0 === 0 && k0 === 0 && l0 === 0) continue;

        // Distinct positions of this reflection, including Friedel mates.
        //
        // A position can be reached two ways: DIRECTLY, as hR for some
        // operator, or as the Friedel mate -hR of another. Both are valid --
        // the stored value differs only by a conjugation that the CONJ_BIT
        // records -- but they are not interchangeable when it comes to the
        // phase factor, because each carries the translation of the operator
        // that produced it.
        //
        // The old code took whichever arrived first, which made the result
        // depend on the order the space-group JSON happens to list its
        // operators. Two passes instead: every direct position first, then the
        // mates for whatever is still unclaimed. The direct entry is the one
        // the symmetrisation projection is actually derived for, so preferring
        // it is not just determinism for its own sake.
        const positions = new Map();
        const opPhases = [];
        for (const op of hklOps) {
            const [h, k, l] = op.apply(h0, k0, l0);
            const ht = h0 * op.t[0] + k0 * op.t[1] + l0 * op.t[2];
            const cr = Math.cos(2 * Math.PI * ht);
            const ci = Math.sin(2 * Math.PI * ht);
            opPhases.push({ h, k, l, cr, ci });
            const key = `${h},${k},${l}`;
            if (!positions.has(key)) positions.set(key, { h, k, l, cr, ci, conj: 0 });
        }
        for (const { h, k, l, cr, ci } of opPhases) {
            const key = `${-h},${-k},${-l}`;
            if (!positions.has(key)) positions.set(key, { h: -h, k: -k, l: -l, cr, ci, conj: 1 });
        }

        // The number of grid points this reflection occupies IS its powder
        // multiplicity: the Laue-group orbit including Friedel mates. The
        // observed intensity is m * Lp * |F|^2, and spreading it over these m
        // points is what divides the multiplicity out. If the host also sends
        // a multiplicity and the two disagree, the intensities are being
        // scaled by the wrong factor and every |F| in the map is wrong by
        // sqrt of the ratio -- silently. Count it and report it.
        const mTrue = positions.size;
        if (Number.isFinite(peak.multiplicity) && peak.multiplicity > 0
            && peak.multiplicity !== mTrue) {
            multiplicityMismatch++;
            if (!multiplicityExample) {
                multiplicityExample = `(${h0},${k0},${l0}): host says m = ${peak.multiplicity}, `
                                    + `the operators give ${mTrue}`;
            }
        }
        const members = [];
        let outOfRange = false;
        for (const p of positions.values()) {
            const m = Math.max(Math.abs(p.h), Math.abs(p.k), Math.abs(p.l));
            if (m > maxIndexSeen) maxIndexSeen = m;
            if (2 * m >= N) { outOfRange = true; continue; }
            members.push(p);
        }
        if (outOfRange) droppedRange++;
        if (members.length === 0) continue;

        // Prefer the Lp the host already computed for this reflection. It was
        // evaluated at the corrected 2-theta with the user's chosen model and
        // is the same number that appears in the refinement report, so taking
        // it verbatim is the only way the map, the report and the
        // French-Wilson prior can be guaranteed to agree. Recomputing here is
        // the fallback for callers that do not send it.
        let lp = 1;
        if (applyLP) {
            if (Number.isFinite(peak.lp) && peak.lp > 0) { lp = peak.lp; lpFromHost++; }
            else { lp = lorentzPolarization(peak.tth, polModel); lpFromModel++; }
        }
        // Truncated orbits keep only their share of the intensity; v137 spread
        // the whole of it over the survivors and inflated their amplitudes.
        const targetI = (I / lp) * (members.length / mTrue);

        refs.push({
            h0, k0, l0,
            tth: Number.isFinite(peak.tth) ? peak.tth : null,
            d: Number.isFinite(peak.d) ? peak.d : null,
            members, mTrue, targetI
        });
    }

    if (refs.length === 0) {
        return { error: 'no-reflections', maxIndexSeen, minGridForAll: 2 * maxIndexSeen + 2,
                 droppedZero, droppedRange };
    }

    // --- pass 2: sort by 2-theta and cluster overlapping reflections --------
    const sortKey = r => (r.tth !== null ? r.tth : (r.d !== null ? 1 / r.d : 0));
    refs.sort((a, b) => sortKey(a) - sortKey(b));

    const clusterOf = new Array(refs.length).fill(0);
    let nClusters = 0;
    for (let i = 0; i < refs.length; i++) {
        if (i === 0) { clusterOf[i] = nClusters; continue; }
        const prev = refs[i - 1], cur = refs[i];
        let same = false;
        if (prev.tth !== null && cur.tth !== null) {
            same = Math.abs(cur.tth - prev.tth) <= tolTth;
        } else if (prev.d !== null && cur.d !== null) {
            same = Math.abs(cur.d - prev.d) <= 1e-4 * Math.max(1e-9, prev.d);
        }
        if (!same) nClusters++;
        clusterOf[i] = nClusters;
    }
    nClusters++;

    // --- pass 3: flatten ----------------------------------------------------
    const orbitStart = [], orbitCount = [], orbitTargetI = [];
    const idxFlat = [], phaseFlat = [];
    const clusterOrbitStart = [], clusterNOrbits = [], clusterObsI = [];

    let cursor = 0;
    let curCluster = -1;
    let usedUnique = 0;

    for (let i = 0; i < refs.length; i++) {
        const r = refs[i];
        if (clusterOf[i] !== curCluster) {
            curCluster = clusterOf[i];
            clusterOrbitStart.push(orbitStart.length);
            clusterNOrbits.push(0);
            clusterObsI.push(0);
        }

        const start = cursor;
        let count = 0;
        for (const p of r.members) {
            const gi = wrap(p.h) + wrap(p.k) * N + wrap(p.l) * N2;
            // Two unique reflections should never land on the same grid point;
            // if they do the input list is redundant for this symmetry, and
            // silently letting the second one win would corrupt both orbits.
            if (claimed.has(gi)) { droppedDuplicate++; continue; }
            claimed.set(gi, orbitStart.length);
            idxFlat.push(p.conj ? (gi | 0x80000000) >>> 0 : gi >>> 0);
            phaseFlat.push(p.cr, p.ci);
            count++;
            cursor++;
        }
        if (count === 0) continue;

        // Rescale again if collisions cost us members.
        const scaled = r.targetI * (count / r.members.length);
        orbitStart.push(start);
        orbitCount.push(count);
        orbitTargetI.push(scaled);

        const c = clusterOrbitStart.length - 1;
        clusterNOrbits[c]++;
        clusterObsI[c] += scaled;
        usedUnique++;
    }

    const nOrbits = orbitStart.length;
    if (nOrbits === 0) {
        return { error: 'no-reflections', maxIndexSeen, minGridForAll: 2 * maxIndexSeen + 2,
                 droppedZero, droppedRange };
    }

    // Drop clusters that ended up empty (every orbit lost to collisions).
    const cOrbitStart = [], cNOrbits = [], cObsI = [];
    for (let c = 0; c < clusterOrbitStart.length; c++) {
        if (clusterNOrbits[c] > 0) {
            cOrbitStart.push(clusterOrbitStart[c]);
            cNOrbits.push(clusterNOrbits[c]);
            cObsI.push(clusterObsI[c]);
        }
    }

    // --- pass 4: systematic absences ---------------------------------------
    const abs = findSystematicAbsences(N, symops, job.centering);
    const absentList = [];
    let absentButObserved = 0;
    for (let i = 0; i < abs.absent.length; i++) {
        const gi = abs.absent[i];
        if (claimed.has(gi)) { absentButObserved++; continue; }
        absentList.push(gi);
    }

    // orbitIdx carries the absent list as a tail block, so the shader gets it
    // without a ninth storage binding.
    const nMembers = idxFlat.length;
    const orbitIdx = new Uint32Array(nMembers + absentList.length);
    orbitIdx.set(Uint32Array.from(idxFlat), 0);
    orbitIdx.set(Uint32Array.from(absentList), nMembers);

    const model = {
        nOrbits,
        nClusters: cOrbitStart.length,
        nMembers,
        orbitStart: Uint32Array.from(orbitStart),
        orbitCount: Uint32Array.from(orbitCount),
        orbitTargetI: Float32Array.from(orbitTargetI),
        orbitIdx,
        phase: Float32Array.from(phaseFlat),
        clusterOrbitStart: Uint32Array.from(cOrbitStart),
        clusterNOrbits: Uint32Array.from(cNOrbits),
        clusterObsI: Float32Array.from(cObsI),
        absentStart: nMembers,
        absentCount: absentList.length,
        absentButObserved,
        symmetrySource: built.source,
        laueClass: built.laue,
        nSymops: symops ? symops.length : hklOps.length,
        haveTranslations: !!symops,
        usedUnique, droppedZero, droppedRange, droppedDuplicate,
        multiplicityMismatch, multiplicityExample,
        lpApplied: applyLP,
        lpModel: applyLP ? describePolarization(polModel) : 'not applied',
        lpFromHost, lpFromModel,
        maxIndexSeen,
        minGridForAll: 2 * maxIndexSeen + 2,
        sumFobs: cObsI.reduce((s, v) => s + Math.sqrt(Math.max(0, v)), 0)
    };
    return model;
}

// Interleaved-complex packing helpers shared by both backends.
function modelInitialGrid(model, N, seed) {
    const init = new Float32Array(2 * N * N * N);
    const rand = mulberry32(seed);
    for (let g = 0; g < model.nOrbits; g++) {
        const start = model.orbitStart[g], count = model.orbitCount[g];
        const amp = Math.sqrt(Math.max(0, model.orbitTargetI[g]) / count);
        // One random phase per ORBIT, propagated through the operator phase
        // factors, so the starting point already obeys the space group.
        const ang = 2 * Math.PI * rand();
        const f0r = amp * Math.cos(ang), f0i = amp * Math.sin(ang);
        for (let i = 0; i < count; i++) {
            const p = start + i;
            const word = model.orbitIdx[p];
            const idx = word & 0x7fffffff;
            const cr = model.phase[2 * p], ci = model.phase[2 * p + 1];
            let a = f0r * cr + f0i * ci;
            let b = f0i * cr - f0r * ci;
            if (word & 0x80000000) b = -b;
            init[2 * idx] = a;
            init[2 * idx + 1] = b;
        }
    }
    return init;
}

// ---------------------------------------------------------------------------
//  Peak search on the final map.
// ---------------------------------------------------------------------------
/**
 * Radius, in angstrom, of the sphere integrated around each peak.
 *
 * WHY 0.65: it has to sit below half the shortest bond you expect to resolve,
 * or the integral swallows the neighbours. A sulfate S-O bond is 1.47 A, so
 * anything past ~0.73 A starts counting oxygen as sulfur. Going much smaller
 * is not free either -- at a grid spacing of ~0.2 A a 0.4 A sphere holds only
 * a handful of voxels and the integral becomes noise. 0.65 A is the widest
 * radius that is still safely inside the tightest common bond.
 * @type {number}
 */
const INTEGRATION_RADIUS = 0.65;

/**
 * Integrate the map over a sphere around a fractional position.
 *
 * WHY THIS EXISTS. The reported peak HEIGHT does not rank atoms by atomic
 * number, and on a real PbSO4 solution it ranked them wrongly: sulfur came out
 * at 1.000 and lead at 0.728, against an electron ratio of 16:82. Height is
 * the value of one voxel. It depends on where the atom happens to fall
 * relative to the grid, on its displacement parameters, and on the
 * series-termination ripples of its neighbours -- and lead, having the largest
 * ripples in the map, erodes its own maximum. The INTEGRAL over a fixed volume
 * is far less sensitive to all three, and is proportional to the electrons
 * enclosed.
 *
 * The number returned is NOT an absolute electron count. F(000) is held at
 * zero throughout charge flipping (zero_dc), so the map has zero mean and no
 * absolute scale, and the integral is "electrons above the mean" in arbitrary
 * units. Only RATIOS between sites are meaningful. Because the missing mean is
 * a constant times the sphere volume, it is the same deficit at every site --
 * a few electrons -- so it compresses the ratios slightly rather than
 * scrambling them. Restoring F(000) from an assumed composition would make
 * this an absolute count; that is a separate change.
 *
 * Negative density is summed as it stands rather than clamped to zero.
 * Clamping would bias every site upward by the depth of the ripples sitting in
 * its sphere, and those are deepest around the heavy atoms -- exactly where
 * the bias would do most damage.
 *
 * @param {Float64Array|Float32Array} rho  Map, N^3, index h + k*N + l*N^2.
 * @param {number} N            Grid size.
 * @param {number[][]} G        Real-space metric tensor.
 * @param {number} voxelVolume  V / N^3, angstrom^3.
 * @param {number} boxR         Half-width of the voxel box to scan, in voxels.
 * @param {number} radius       Sphere radius, angstrom.
 * @param {number} px @param {number} py @param {number} pz  Fractional centre.
 * @returns {number} Integrated density, arbitrary units.
 */
function integrateSphere(rho, N, G, voxelVolume, boxR, radius, px, py, pz) {
    const N2 = N * N;
    const wrap = (v) => ((v % N) + N) % N;
    const ch = px * N, ck = py * N, cl = pz * N;
    const h0 = Math.round(ch), k0 = Math.round(ck), l0 = Math.round(cl);
    // Two regions: the sphere itself, and a shell just outside it.
    //
    // THE SHELL IS NOT OPTIONAL. F(000) is pinned at zero, so the map's mean is
    // zero when the true mean is (total electrons)/V -- about 1.6 e/A^3 for a
    // structure like PbSO4. Every sphere therefore loses the same ~2 electrons,
    // and subtracting a constant from both members of a ratio does not cancel:
    // measured on a synthetic PbSO4 map it pushed S/Pb from the true 0.195 down
    // to 0.168 and O/Pb from 0.098 to 0.067, biasing the light atoms hardest --
    // precisely the ones that are hard to identify. The shell measures the
    // local baseline instead of assuming it, so the offset cancels, and it
    // absorbs slowly varying ripples at the same time.
    const rOut = radius * 1.6;
    let sum = 0, nIn = 0;
    let shellSum = 0, nShell = 0;
    const boxOut = Math.ceil(boxR * 1.6);
    for (let dl = -boxOut; dl <= boxOut; dl++) {
        for (let dk = -boxOut; dk <= boxOut; dk++) {
            for (let dh = -boxOut; dh <= boxOut; dh++) {
                const hh = h0 + dh, kk = k0 + dk, ll = l0 + dl;
                // Fractional offset of this voxel from the true peak centre,
                // which is NOT a grid point -- the centroid refinement moved
                // it. Measuring from the rounded index instead would bias the
                // sphere by up to half a voxel, systematically and per site.
                const fx = (hh - ch) / N, fy = (kk - ck) / N, fz = (ll - cl) / N;
                const d = fracDistance(G, fx, fy, fz);
                if (d > rOut) continue;
                const v = rho[wrap(hh) + wrap(kk) * N + wrap(ll) * N2];
                if (d <= radius) { sum += v; nIn++; }
                else { shellSum += v; nShell++; }
            }
        }
    }
    const baseline = nShell > 0 ? shellSum / nShell : 0;
    return (sum - baseline * nIn) * voxelVolume;
}

function findPeaks(rho, N, cell, opts) {
    const N2 = N * N;
    const minSeparation = opts.minSeparation;
    const maxPeaks = opts.maxPeaks;
    const G = metricTensor(cell);
    const voxelVolume = cellVolume(cell) / (N * N * N);
    // Voxels to scan on each side: the sphere radius in units of the COARSEST
    // axis spacing, so the box always encloses the sphere whatever the cell
    // shape. Capped at N/2 so a tiny cell cannot make the box wrap onto itself.
    const minSpacing = Math.min(cell.a / N, cell.b / N, cell.c / N);
    const boxR = Math.min(Math.floor(N / 2),
                          Math.max(1, Math.ceil(INTEGRATION_RADIUS / Math.max(1e-6, minSpacing))));

    let max = -Infinity;
    for (let i = 0; i < rho.length; i++) if (rho[i] > max) max = rho[i];
    if (!(max > 0)) return [];

    const cutoff = max * opts.heightCutoff;
    const wrap = (v) => ((v % N) + N) % N;
    const candidates = [];

    for (let l = 0; l < N; l++) {
        for (let k = 0; k < N; k++) {
            for (let h = 0; h < N; h++) {
                const v = rho[h + k * N + l * N2];
                if (v < cutoff) continue;
                let isMax = true;
                for (let dl = -1; dl <= 1 && isMax; dl++) {
                    for (let dk = -1; dk <= 1 && isMax; dk++) {
                        for (let dh = -1; dh <= 1; dh++) {
                            if (dh === 0 && dk === 0 && dl === 0) continue;
                            const n = wrap(h + dh) + wrap(k + dk) * N + wrap(l + dl) * N2;
                            if (rho[n] > v) { isMax = false; break; }
                        }
                    }
                }
                if (!isMax) continue;

                let sumW = 0, sumH = 0, sumK = 0, sumL = 0;
                for (let ddl = -2; ddl <= 2; ddl++) {
                    for (let ddk = -2; ddk <= 2; ddk++) {
                        for (let ddh = -2; ddh <= 2; ddh++) {
                            const n = wrap(h + ddh) + wrap(k + ddk) * N + wrap(l + ddl) * N2;
                            const w = Math.max(0, rho[n]);
                            sumW += w;
                            sumH += w * ddh;
                            sumK += w * ddk;
                            sumL += w * ddl;
                        }
                    }
                }
                const dh = sumW > 0 ? Math.max(-2.0, Math.min(2.0, sumH / sumW)) : 0;
                const dk = sumW > 0 ? Math.max(-2.0, Math.min(2.0, sumK / sumW)) : 0;
                const dl = sumW > 0 ? Math.max(-2.0, Math.min(2.0, sumL / sumW)) : 0;
                const px = ((h + dh) / N + 1) % 1;
                const py = ((k + dk) / N + 1) % 1;
                const pz = ((l + dl) / N + 1) % 1;
                candidates.push({
                    x: px, y: py, z: pz,
                    height: v,
                    charge: integrateSphere(rho, N, G, voxelVolume, boxR, INTEGRATION_RADIUS, px, py, pz)
                });
            }
        }
    }

    candidates.sort((p, q) => q.height - p.height);

    const kept = [];
    for (const c of candidates) {
        let tooClose = false;
        for (const p of kept) {
            if (fracDistance(G, c.x - p.x, c.y - p.y, c.z - p.z) < minSeparation) { tooClose = true; break; }
        }
        if (tooClose) continue;
        kept.push(c);
        if (kept.length >= maxPeaks) break;
    }

    const top = kept.length ? kept[0].height : 1;
    // The charge ratio is normalised to the LARGEST INTEGRAL, not to the
    // tallest peak. They are not the same site in general -- which is the
    // whole reason this exists.
    let topQ = 0;
    for (const p of kept) if (p.charge > topQ) topQ = p.charge;
    return kept.map((p, i) => ({
        rank: i + 1,
        x: p.x, y: p.y, z: p.z,
        height: p.height,
        relative: p.height / top,
        charge: p.charge,
        chargeRel: topQ > 0 ? p.charge / topQ : 0
    }));
}

// ---------------------------------------------------------------------------
//  Shared reporting block, so the GPU and CPU results are interchangeable.
// ---------------------------------------------------------------------------
function modelReport(model) {
    return {
        unique: model.usedUnique,
        gridPoints: model.nMembers,
        orbits: model.nOrbits,
        clusters: model.nClusters,
        overlapped: model.nOrbits - model.nClusters,
        absencesZeroed: model.absentCount,
        absentButObserved: model.absentButObserved,
        symmetrySource: model.symmetrySource,
        laueClass: model.laueClass,
        nSymops: model.nSymops,
        droppedOutOfRange: model.droppedRange,
        droppedZeroIntensity: model.droppedZero,
        droppedDuplicate: model.droppedDuplicate,
        lpApplied: model.lpApplied,
        lpModel: model.lpModel,
        lpFromHost: model.lpFromHost,
        lpFromModel: model.lpFromModel,
        maxIndex: model.maxIndexSeen,
        minGridForAll: model.minGridForAll,
        // Non-zero means the host's multiplicities and the space-group
        // operators disagree, i.e. the observed intensities are being divided
        // by the wrong m. Surfaced rather than swallowed: it makes every |F|
        // in the map wrong by a factor the user would otherwise never see.
        multiplicityMismatch: model.multiplicityMismatch,
        multiplicityExample: model.multiplicityExample
    };
}

// ===========================================================================
//  CPU PATH
//
//  Mirrors the GPU kernels one for one and in the same order, so a run with
//  backend: 'cpu' can be compared against the GPU result directly. It is the
//  reference implementation, not a simplified stand-in.
// ===========================================================================
function runChargeFlippingCPU(job) {
    const N = job.gridSize;
    const N3 = N * N * N;
    const maxIter = job.maxIterations;
    const deltaSigma = job.thresholdSigma;
    const cell = job.cell;
    const lambda = Math.min(1, Math.max(0, Number(job.symLambda) || 0));

    if (!Number.isInteger(N) || N < 8 || (N & (N - 1)) !== 0) {
        return { error: `Grid size must be a power of two (the FFT is radix-2); got ${N}.` };
    }

    const model = buildReflectionModel(job, N);
    if (model.error) {
        return { error: `No reflection fits a ${N}x${N}x${N} grid. The largest index in the fit is ${model.maxIndexSeen}, which needs a grid of at least ${model.minGridForAll}.` };
    }

    const rand = mulberry32(job.seed ^ 0x9e3779b9);
    const init = modelInitialGrid(model, N, job.seed);

    const re = new Float64Array(N3);
    const im = new Float64Array(N3);
    for (let i = 0; i < N3; i++) { re[i] = init[2 * i]; im[i] = init[2 * i + 1]; }
    re[0] = 0; im[0] = 0;

    const target = Float64Array.from(model.orbitTargetI);
    const rHistory = new Float32Array(maxIter);
    const bestFre = new Float64Array(N3);
    const bestFim = new Float64Array(N3);
    let haveBest = false, bestR = Infinity, bestIter = 0, lastReport = -1;

    const idx = model.orbitIdx, ph = model.phase;
    const oStart = model.orbitStart, oCount = model.orbitCount;
    const cStart = model.clusterOrbitStart, cN = model.clusterNOrbits, cObs = model.clusterObsI;

    for (let iter = 0; iter < maxIter; iter++) {
        // --- reciprocal -> real ---
        fft3d(re, im, N, true);

        // --- flip threshold in units of the map's own sigma ---
        let mean = 0;
        for (let i = 0; i < N3; i++) mean += re[i];
        mean /= N3;
        let varSum = 0;
        for (let i = 0; i < N3; i++) { const d = re[i] - mean; varSum += d * d; }
        const delta = deltaSigma * Math.sqrt(varSum / N3);

        for (let i = 0; i < N3; i++) {
            if (re[i] < delta) re[i] = -re[i];
            im[i] = 0;
        }

        // --- real -> reciprocal ---
        fft3d(re, im, N, false);

        // --- R factor on the OBSERVABLE quantity: cluster totals ---
        let sumDiff = 0;
        for (let c = 0; c < model.nClusters; c++) {
            let total = 0;
            for (let o = 0; o < cN[c]; o++) {
                const g = cStart[c] + o;
                for (let i = 0; i < oCount[g]; i++) {
                    const j = idx[oStart[g] + i] & 0x7fffffff;
                    total += re[j] * re[j] + im[j] * im[j];
                }
            }
            sumDiff += Math.abs(Math.sqrt(Math.max(0, total)) - Math.sqrt(Math.max(0, cObs[c])));
        }
        const R = model.sumFobs > 0 ? sumDiff / model.sumFobs : NaN;
        rHistory[iter] = R;

        // --- repartition the observed intensity inside each 2-theta cluster ---
        for (let c = 0; c < model.nClusters; c++) {
            let total = 0;
            for (let o = 0; o < cN[c]; o++) {
                const g = cStart[c] + o;
                let s = 0;
                for (let i = 0; i < oCount[g]; i++) {
                    const j = idx[oStart[g] + i] & 0x7fffffff;
                    s += re[j] * re[j] + im[j] * im[j];
                }
                target[g] = s;
                total += s;
            }
            if (total > 1e-20) {
                const k = cObs[c] / total;
                for (let o = 0; o < cN[c]; o++) target[cStart[c] + o] *= k;
            } else {
                let wsum = 0;
                for (let o = 0; o < cN[c]; o++) wsum += oCount[cStart[c] + o];
                if (wsum > 0) {
                    for (let o = 0; o < cN[c]; o++) {
                        target[cStart[c] + o] = cObs[c] * oCount[cStart[c] + o] / wsum;
                    }
                }
            }
        }

        // --- impose the space group on the structure factors ---
        if (lambda > 0) {
            for (let g = 0; g < model.nOrbits; g++) {
                const count = oCount[g];
                if (count < 2) continue;
                const start = oStart[g];
                let sr = 0, si = 0;
                for (let i = 0; i < count; i++) {
                    const p = start + i;
                    const word = idx[p];
                    const j = word & 0x7fffffff;
                    let a = re[j], b = im[j];
                    if (word & 0x80000000) b = -b;
                    const cr = ph[2 * p], ci = ph[2 * p + 1];
                    sr += a * cr - b * ci;
                    si += a * ci + b * cr;
                }
                sr /= count; si /= count;
                for (let i = 0; i < count; i++) {
                    const p = start + i;
                    const word = idx[p];
                    const j = word & 0x7fffffff;
                    const cr = ph[2 * p], ci = ph[2 * p + 1];
                    let a = sr * cr + si * ci;
                    let b = si * cr - sr * ci;
                    if (word & 0x80000000) b = -b;
                    re[j] += lambda * (a - re[j]);
                    im[j] += lambda * (b - im[j]);
                }
            }
        }

        // --- equal amplitudes across each symmetry orbit ---
        for (let g = 0; g < model.nOrbits; g++) {
            const start = oStart[g], count = oCount[g];
            if (count === 0) continue;
            const amp = Math.sqrt(Math.max(0, target[g]) / count);
            for (let i = 0; i < count; i++) {
                const j = idx[start + i] & 0x7fffffff;
                // PERF: Math.hypot is ~16x slower than the naive form (measured,
                // 4e6 iterations) because it guards against intermediate
                // overflow/underflow. Structure factors here are bounded well
                // inside double range, so the guard buys nothing. Identical
                // results to the last bit on every value this loop sees.
                const m = Math.sqrt(re[j] * re[j] + im[j] * im[j]);
                if (m > 1e-20) {
                    const s = amp / m;
                    re[j] *= s; im[j] *= s;
                } else {
                    const a = 2 * Math.PI * rand();
                    re[j] = amp * Math.cos(a);
                    im[j] = amp * Math.sin(a);
                }
            }
        }

        // --- systematic absences and F(000) ---
        for (let i = 0; i < model.absentCount; i++) {
            const j = idx[model.absentStart + i] & 0x7fffffff;
            re[j] = 0; im[j] = 0;
        }
        re[0] = 0; im[0] = 0;

        if (Number.isFinite(R) && R < bestR) {
            bestR = R; bestIter = iter + 1;
            bestFre.set(re); bestFim.set(im);
            haveBest = true;
        }

        const pct = Math.floor(((iter + 1) / maxIter) * 100);
        if (pct !== lastReport || iter === maxIter - 1) {
            lastReport = pct;
            postMessage({ type: 'cf-progress', iter: iter + 1, total: maxIter, R, bestR });
        }
    }

    if (haveBest) { re.set(bestFre); im.set(bestFim); }
    fft3d(re, im, N, true);
    const bestRho = Float32Array.from(re);

    let peak = 0;
    for (let i = 0; i < N3; i++) if (bestRho[i] > peak) peak = bestRho[i];
    if (peak > 0) for (let i = 0; i < N3; i++) bestRho[i] /= peak;

    const peaks = findPeaks(bestRho, N, cell, {
        minSeparation: job.minPeakSeparation,
        maxPeaks: job.maxPeaks,
        heightCutoff: job.peakHeightCutoff
    });

    return {
        map: bestRho, gridSize: N, peaks, rHistory, bestR, bestIter,
        finalR: rHistory[maxIter - 1], seed: job.seed, cell,
        volume: cellVolume(cell), backend: 'cpu', symLambda: lambda,
        reflections: modelReport(model)
    };
}

// ===========================================================================
//  DISPATCH
// ===========================================================================
async function runChargeFlipping(job) {
    const wantGPU = job.backend !== 'cpu';
    if (wantGPU && typeof runChargeFlippingGPU === 'function') {
        try {
            const out = await runChargeFlippingGPU(job);
            if (out && out.__noGPU) {
                postMessage({ type: 'cf-info', message: 'WebGPU unavailable; using CPU.' });
            } else {
                return out;
            }
        } catch (gpuErr) {
            console.warn('GPU charge flipping failed, falling back to CPU:', gpuErr);
            postMessage({ type: 'cf-info', message: 'GPU path failed; using CPU. (' + ((gpuErr && gpuErr.message) || gpuErr) + ')' });
        }
    }
    return runChargeFlippingCPU(job);
}

globalThis.onmessage = async function (e) {
    const job = e.data || {};
    try {
        if (job.type === 'run-charge-flipping') {
            const out = await runChargeFlipping(job);
            if (out.error) { postMessage({ type: 'cf-error', message: out.error }); return; }
            postMessage({ type: 'cf-result', results: out }, [out.map.buffer, out.rHistory.buffer]);
            return;
        }
if (job.type === 'build-structure') {
            const out = await buildStructure(job);
            if (out.error) { postMessage({ type: 'cf-error', message: out.error }); return; }
            postMessage({ type: 'cf-structure', results: out }, [out.map.buffer]);
            return;
        }

        // ------------------------------------------------------------------
        //  STANDALONE WYCKOFF SEARCH. NO MAP, NO CHARGE FLIPPING.
        //
        //  buildStructure exists to take a P1 density map and place it on a
        //  crystallographic origin: it searches origin shifts, decides a hand,
        //  and reports a symmetry correlation. Every one of those is an
        //  artefact of having solved a map. A Wyckoff search has none of those
        //  problems -- the projection operators put each atom on a
        //  crystallographic site by construction, so the origin is already
        //  fixed, and nothing is being inverted so there is no hand. This path
        //  therefore does not replace that stage, it simply has no use for it.
        //
        //  What it needs is what the Pawley fit already produced: a space
        //  group, a cell, a composition, and |Fobs|^2 with multiplicities and
        //  Lp. The map was never an input.
        // ------------------------------------------------------------------
        if (job.type === 'wyckoff-search') {
            const symops = job.symops || [];
            if (!symops.length) {
                postMessage({ type: 'cf-error', message: 'No symmetry operators for the chosen setting.' });
                return;
            }
            if (!job.targetComposition) {
                postMessage({ type: 'cf-error', message: 'A target composition is required for a Wyckoff search.' });
                return;
            }
            const N = job.gridSize || 32;
            const res = await runWyckoffDensityFit(null, N, job.cell, symops, job);
            if (res.error) { postMessage({ type: 'cf-error', message: res.error }); return; }
            if (!res.sites || !res.sites.length) {
                postMessage({ type: 'cf-error', message: 'The Wyckoff search produced no sites.' });
                return;
            }
            postMessage({
                type: 'cf-structure',
                results: {
                    sites: res.sites,
                    wyckoffResult: res,
                    refinement: res.refinement,
                    cell: job.cell,
                    nOps: symops.length,
                    symops: symops.map(o => o.xyz).filter(Boolean),
                    rawPeakCount: res.sites.length,
                    // No map, so nothing to report about one. Explicit nulls so
                    // the UI can tell "not applicable" from "zero".
                    map: null, gridSize: null,
                    originShift: null, hand: null,
                    symmetryCorrelation: null, originScore: null,
                    source: 'wyckoff'
                }
            });
            return;
        }
    } catch (err) {
        postMessage({ type: 'cf-error', message: (err && err.message) || String(err) });
    }
};

// ===========================================================================
//  GPU PATH
// ===========================================================================

let _gpuDevicePromise = null;
async function cfAcquireGPU() {
    if (_gpuDevicePromise) return _gpuDevicePromise;
    _gpuDevicePromise = (async () => {
        if (typeof navigator === 'undefined' || !navigator.gpu) return null;
        try {
            const adapter = await navigator.gpu.requestAdapter();
            if (!adapter) return null;
            const device = await adapter.requestDevice();
            // Record the loss so a RUNNING loop can see it. Clearing the
            // cache is not enough: when the driver resets the device mid-run
            // (a TDR, which a long 128^3 job can provoke) every pending
            // mapAsync simply never settles, so the loop sits on an await
            // forever -- busy cursor, live Stop button, no progress, no error.
            // A silent hang is the worst possible way to report this.
            device.__lostInfo = null;
            device.lost.then((info) => {
                device.__lostInfo = info || { reason: 'unknown', message: 'device lost' };
                _gpuDevicePromise = null;
                console.error('[charge flipping] WebGPU device lost:',
                              device.__lostInfo.reason, device.__lostInfo.message || '');
            });
            return device;
        } catch (e) {
            console.warn('WebGPU unavailable:', e);
            return null;
        }
    })();
    return _gpuDevicePromise;
}

let _cfShaderSrc = null;
async function cfLoadShader() {
    if (_cfShaderSrc) return _cfShaderSrc;
    const resp = await fetch('../charge_flipping.wgsl');
    if (!resp.ok) throw new Error('charge_flipping.wgsl failed to load (HTTP ' + resp.status + ')');
    _cfShaderSrc = await resp.text();
    return _cfShaderSrc;
}

// Two bind-group layouts. Splitting them keeps every pipeline inside the
// default limit of 8 storage buffers per shader stage: the FFT pipelines see
// three, the symmetry pipelines see eight.
async function cfGetPipelines(device) {
    if (device.__cfPipelines) return device.__cfPipelines;
    const code = await cfLoadShader();
    const module = device.createShaderModule({ code });

    const layout0 = device.createBindGroupLayout({
        entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
            { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
            { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }
        ]
    });
    const layout1 = device.createBindGroupLayout({
        entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
            { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
            { binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }
        ]
    });

    // Group 2 exists solely so that `best` does not have to live in group 0.
    // The per-stage storage-buffer limit is checked against the PIPELINE
    // LAYOUT, not against what the entry point reads: group 0 (three storage
    // buffers) plus group 1 (five) is already exactly the guaranteed maximum
    // of eight, so a fourth buffer in group 0 made every symmetry pipeline
    // nine and createPipelineLayout threw -- even though symmetrize, constrain
    // and the rest never touch `best`. Only keep_best binds group 2.
    //
    // Deliberately NOT solved by asking for a higher maxStorageBuffersPerShader
    // Stage in requestDevice(): many adapters report 16, but requiring it would
    // refuse the device outright on hardware that only offers the mandated 8.
    const layoutBest  = device.createBindGroupLayout({
        entries: [{ binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }]
    });
    // Group 1 is unused by keep_best, but the group INDICES have to line up
    // with the @group attributes in the WGSL, so slot 1 gets an empty layout.
    const layoutEmpty = device.createBindGroupLayout({ entries: [] });

    const plGrid = device.createPipelineLayout({ bindGroupLayouts: [layout0] });
    const plSym  = device.createPipelineLayout({ bindGroupLayouts: [layout0, layout1] });
    const plBest = device.createPipelineLayout({
        bindGroupLayouts: [layout0, layoutEmpty, layoutBest]
    });

    const mk = (entryPoint, pipelineLayout) => device.createComputePipeline({
        layout: pipelineLayout, compute: { module, entryPoint }
    });

    const pipelines = {
        layout0, layout1, layoutBest, layoutEmpty,
        fft_stage:        mk('fft_stage', plGrid),
        copy_buf:         mk('copy_buf', plGrid),
        scale_inverse:    mk('scale_inverse', plGrid),
        reduce_moments:   mk('reduce_moments', plGrid),
        finalize_moments: mk('finalize_moments', plGrid),
        flip:             mk('flip', plGrid),
        zero_dc:          mk('zero_dc', plGrid),
        keep_best:        mk('keep_best', plBest),
        commit_best:      mk('commit_best', plGrid),
        finalize_rfactor: mk('finalize_rfactor', plGrid),
        reduce_rfactor:   mk('reduce_rfactor', plSym),
        repartition:      mk('repartition', plSym),
        symmetrize:       mk('symmetrize', plSym),
        constrain:        mk('constrain', plSym),
        zero_absent:      mk('zero_absent', plSym)
    };

    device.__cfPipelines = pipelines;
    return pipelines;
}

const CF_WG = 64;
const CF_COPY_PER_THREAD = 4;   // must match COPY_PER_THREAD in the WGSL
// Iterations per R-factor readback. Large enough that the mapAsync stall stops
// dominating a small-grid run, small enough that the progress plot still moves.
// ---------------------------------------------------------------------------
//  ADAPTIVE BATCH SIZE
//
//  How many iterations may be queued before the host syncs. Too few and the
//  run is dominated by round-trip latency; too many and the GPU spends seconds
//  in one uninterrupted burst, which is how a driver watchdog decides the
//  device has stopped responding and resets it (a TDR). The reset is silent:
//  pending maps never settle and the loop hangs.
//
//  A table keyed on grid size was the first attempt and it is the wrong shape
//  of answer -- it guesses at hardware it cannot see. A 128^3 run is trivial
//  for a current discrete card and impossible for an old integrated one, and
//  the same table has to serve both. So the batch is MEASURED instead:
//  time each drain, divide by the iterations it covered, and pick the count
//  that lands the next burst near the target.
//
//  Starting at one iteration means the very first measurement is safe on any
//  device, however slow, before anything is queued in bulk.
// ---------------------------------------------------------------------------

/** Target GPU time per sync, ms. Windows' default TDR is 2000 ms; this leaves
 *  a factor of five, which also covers the burst being slower than the sample
 *  that sized it. */
const CF_BATCH_TARGET_MS = 400;

/** Above this, shrink immediately rather than easing down: the next burst at
 *  the current size would be uncomfortably close to a watchdog reset. */
const CF_BATCH_PANIC_MS = 900;

/** Never queue more than this, however fast the device looks. Past a point the
 *  latency saving is negligible and the exposure is not. */
const CF_R_BATCH_MAX = 25;

/**
 * Choose the next batch size from the last one and how long it took.
 *
 * @param {number} current   Iterations in the batch just drained.
 * @param {number} elapsedMs Wall time from first submit to drained.
 * @param {boolean} warmup   True for the first batch, whose time includes
 *                           pipeline compilation and is not representative.
 * @returns {number} Iterations for the next batch, at least 1.
 */
function cfNextBatch(current, elapsedMs, warmup, state) {
    if (!(elapsedMs > 0) || !(current > 0)) return current;
    const perIter = elapsedMs / current;

    // SIZE FROM THE WORST RECENT SAMPLE, NOT THE LAST ONE.
    //
    // A controller that trusts the most recent measurement is always one step
    // behind: if the device slows down sharply -- thermal throttling, another
    // application taking the GPU, a browser moving the tab to a different
    // adapter -- the batch it already chose runs at the new speed before
    // anything can react, and at a large batch that single burst is enough to
    // trip the watchdog. Simulating a 30x slowdown at batch 25 reaches 7.5 s
    // in one step, well past a 2 s reset.
    //
    // The estimator therefore keeps a running worst that decays slowly. It
    // reacts instantly to a device getting slower and forgets an old bad
    // sample over several batches, so a transient hiccup costs a little
    // throughput rather than a run.
    const st = state || {};
    st.worstPerIter = Math.max(perIter, (st.worstPerIter || 0) * 0.75);
    const estimate = st.worstPerIter;

    // A warm-up sample may only ever SHRINK the batch. Compilation cost is
    // charged once, and letting it set the steady-state size would peg a fast
    // device at one iteration for the whole run.
    if (warmup) {
        return (elapsedMs > CF_BATCH_PANIC_MS)
            ? 1
            : Math.max(1, Math.min(CF_R_BATCH_MAX, Math.floor(CF_BATCH_TARGET_MS / estimate)));
    }
    if (elapsedMs > CF_BATCH_PANIC_MS) return Math.max(1, Math.floor(current / 2));

    const want = CF_BATCH_TARGET_MS / estimate;
    // Damped: at most double per step. An undamped jump overshoots on the
    // first fast sample and then oscillates between 1 and the cap.
    return Math.max(1, Math.min(CF_R_BATCH_MAX, Math.floor(Math.min(want, current * 2 + 1))));
}
function cfCeilGroups(threads) { return Math.max(1, Math.ceil(threads / CF_WG)); }
function cfCopyGroups(N) { return cfCeilGroups((2 * N * N * N) / CF_COPY_PER_THREAD); }

// One 3D transform. Stages ping-pong between bufA and bufB through the two
// pre-baked bind-group variants; the single trailing copy runs only when the
// stage count is odd, instead of v137's copy after every stage.
function cfEncodeFFT(enc, pl, res, N, inverse) {
    const stageGroups = cfCeilGroups((N >> 1) * N * N);
    const copyGroups = cfCopyGroups(N);
    const bgs = inverse ? res.fftBindGroups.inverse : res.fftBindGroups.forward;

    let stage = 0;
    for (let axis = 0; axis < 3; axis++) {
        for (let L = 1; L < N; L <<= 1) {
            const pass = enc.beginComputePass();
            pass.setPipeline(pl.fft_stage);
            pass.setBindGroup(0, bgs[stage]);
            pass.dispatchWorkgroups(stageGroups);
            pass.end();
            stage++;
        }
    }

    // Odd number of stages leaves the answer in bufB; bring it home.
    if (stage & 1) {
        const pass = enc.beginComputePass();
        pass.setPipeline(pl.copy_buf);
        pass.setBindGroup(0, res.bgSwapped);
        pass.dispatchWorkgroups(copyGroups);
        pass.end();
    }

    if (inverse) {
        const pass = enc.beginComputePass();
        pass.setPipeline(pl.scale_inverse);
        pass.setBindGroup(0, res.bgMain);
        pass.dispatchWorkgroups(cfCeilGroups(N * N * N));
        pass.end();
    }
}

async function runChargeFlippingGPU(job) {
    const N = job.gridSize;
    const N3 = N * N * N;
    const maxIter = job.maxIterations;
    const cell = job.cell;
    const lambda = Math.min(1, Math.max(0, Number(job.symLambda) || 0));

    if (!Number.isInteger(N) || N < 8 || (N & (N - 1)) !== 0) {
        return { error: `Grid size must be a power of two; got ${N}.` };
    }

    const device = await cfAcquireGPU();
    if (!device) return { __noGPU: true };

    const model = buildReflectionModel(job, N);
    if (model.error) {
        return { error: `No reflection fits a ${N}x${N}x${N} grid. Largest index ${model.maxIndexSeen} needs at least ${model.minGridForAll}.` };
    }

    const pl = await cfGetPipelines(device);

    // -----------------------------------------------------------------------
    //  GPU allocations.
    //
    //  BUG FIX (leak on a failed allocation). Every buffer is registered with
    //  `track` AT THE MOMENT OF CREATION, and the try/finally opens BEFORE the
    //  first allocation rather than after the last one. Previously the eleven
    //  buffers below were created first and only then collected into an array
    //  inside the try, so a throw partway through allocation -- the adapter
    //  refusing a 50 MB buffer at 128^3, a device lost between two calls --
    //  leaked everything already allocated, with no reference left to destroy
    //  it. Repeated failed attempts would eventually take the adapter down,
    //  which is precisely the scenario in which allocation fails.
    //
    //  Nothing may be allocated outside `track`.
    // -----------------------------------------------------------------------
    /** @type {GPUBuffer[]} */
    const gpuBuffers = [];

    /**
     * Creates a device buffer and registers it for release in the finally
     * block below.
     * @param {GPUBufferDescriptor} desc
     * @returns {GPUBuffer}
     */
    const track = (desc) => {
        const b = device.createBuffer(desc);
        gpuBuffers.push(b);
        return b;
    };

    try {

    const cplxBytes = 2 * N3 * 4;
    const mkStorage = (bytes, extra = 0) => track({
        size: Math.max(4, bytes),
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC | extra
    });

    const bufA = mkStorage(cplxBytes);
    const bufB = mkStorage(cplxBytes);
    // The best iterate is now WRITTEN BY A KERNEL (keep_best), so it needs
    // STORAGE. It used to be a pure copy target driven from the host.
    const bestBuf = track({
        size: cplxBytes,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC
    });
    device.queue.writeBuffer(bufA, 0, modelInitialGrid(model, N, job.seed));

    const wgFull = cfCeilGroups(N3);
    const wgCluster = cfCeilGroups(model.nClusters);
    const partialsLen = 4 + 2 * Math.max(wgFull, wgCluster);
    const partialsBuf = mkStorage(partialsLen * 4);
    // partials[P_BESTR] is the record R that keep_best compares against. It has
    // to start at +inf so the first cycle always snapshots.
    device.queue.writeBuffer(partialsBuf, 0, Float32Array.of(0, 0, 0, Infinity));

    let rBatch = 1;   // grows once the first drain has been timed
    // R factors are read back in BATCHES of rBatch iterations instead of
    // one mapAsync per cycle. The device now picks the best iterate itself
    // (keep_best / commit_best), so nothing in the loop depends on the host
    // having seen R; the history is only wanted for the convergence plot, and
    // being up to CF_R_BATCH cycles behind costs nothing there.
    const headerRead = track({
        // Sized for the largest batch the controller may choose; a WebGPU
        // buffer cannot be resized, and only the first batchCount*16 bytes
        // are ever read.
        size: CF_R_BATCH_MAX * 16, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
    });

    // Orbit is now TWO u32 (stride 8): the third field, target_i, was written
    // here and read by no kernel. repartition fills orbitTarget[] from the
    // calculated intensities and constrain reads only that, so the copy in the
    // struct was a second source of truth for the same number -- the kind that
    // stays right until someone changes one of them.
    // Cluster is still three 4-byte scalars, stride 12.
    const orbitsData = new ArrayBuffer(Math.max(1, model.nOrbits) * 8);
    const ov = new DataView(orbitsData);
    for (let g = 0; g < model.nOrbits; g++) {
        ov.setUint32(g * 8 + 0, model.orbitStart[g], true);
        ov.setUint32(g * 8 + 4, model.orbitCount[g], true);
    }
    const clustersData = new ArrayBuffer(Math.max(1, model.nClusters) * 12);
    const cv = new DataView(clustersData);
    for (let c = 0; c < model.nClusters; c++) {
        cv.setUint32(c * 12 + 0, model.clusterOrbitStart[c], true);
        cv.setUint32(c * 12 + 4, model.clusterNOrbits[c], true);
        cv.setFloat32(c * 12 + 8, model.clusterObsI[c], true);
    }

    const orbitsBuf   = mkStorage(orbitsData.byteLength);
    const orbitIdxBuf = mkStorage(Math.max(1, model.orbitIdx.length) * 4);
    const phaseBuf    = mkStorage(Math.max(1, model.phase.length) * 4);
    const clustersBuf = mkStorage(clustersData.byteLength);
    const targetBuf   = mkStorage(Math.max(1, model.nOrbits) * 4);

    device.queue.writeBuffer(orbitsBuf, 0, orbitsData);
    device.queue.writeBuffer(orbitIdxBuf, 0, model.orbitIdx);
    device.queue.writeBuffer(phaseBuf, 0, model.phase.length ? model.phase : new Float32Array(1));
    device.queue.writeBuffer(clustersBuf, 0, clustersData);
    device.queue.writeBuffer(targetBuf, 0, model.orbitTargetI);

    // --- uniform blocks ----------------------------------------------------
    // Every field is static for the whole trial, so nothing is written to the
    // uniform buffer inside the loop. Only axis / stage_len / inverse differ
    // between blocks, and those are pre-baked, one block per FFT stage.
    const align = 256;
    const stagesPerFFT = 3 * Math.log2(N);
    const totalBlocks = 1 + stagesPerFFT * 2;
    const paramData = new ArrayBuffer(totalBlocks * align);
    const pv = new DataView(paramData);

    const writeP = (blk, axis, stageLen, inverse) => {
        const o = blk * align;
        pv.setUint32(o + 0, N, true);
        pv.setUint32(o + 4, axis, true);
        pv.setUint32(o + 8, stageLen, true);
        pv.setUint32(o + 12, inverse, true);
        pv.setFloat32(o + 16, job.thresholdSigma, true);
        pv.setUint32(o + 20, model.nOrbits, true);
        pv.setUint32(o + 24, model.nClusters, true);
        pv.setUint32(o + 28, model.absentStart, true);
        pv.setUint32(o + 32, model.absentCount, true);
        pv.setFloat32(o + 36, lambda, true);
        pv.setUint32(o + 40, job.seed >>> 0, true);
        pv.setUint32(o + 44, wgFull, true);
        pv.setUint32(o + 48, wgCluster, true);
        pv.setUint32(o + 52, 0, true);
        pv.setUint32(o + 56, 0, true);
        pv.setUint32(o + 60, 0, true);
    };
    writeP(0, 0, 1, 0);

    const paramBuf = track({
        size: totalBlocks * align,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    });

    const mkBg0 = (blk, swapped) => device.createBindGroup({
        layout: pl.layout0,
        entries: [
            { binding: 0, resource: { buffer: paramBuf, offset: blk * align, size: 64 } },
            { binding: 1, resource: { buffer: swapped ? bufB : bufA } },
            { binding: 2, resource: { buffer: swapped ? bufA : bufB } },
            { binding: 3, resource: { buffer: partialsBuf } }
        ]
    });

    const fftBindGroups = { forward: [], inverse: [] };
    let blk = 1;
    for (const inv of [0, 1]) {
        const arr = inv ? fftBindGroups.inverse : fftBindGroups.forward;
        let stage = 0;
        for (let axis = 0; axis < 3; axis++) {
            for (let L = 1; L < N; L <<= 1) {
                writeP(blk, axis, L, inv);
                arr.push(mkBg0(blk, (stage & 1) === 1));   // stage n reads A when n even
                blk++; stage++;
            }
        }
    }
    device.queue.writeBuffer(paramBuf, 0, paramData);

    const bgMain = mkBg0(0, false);
    const bgSwapped = mkBg0(0, true);
    const bg1 = device.createBindGroup({
        layout: pl.layout1,
        entries: [
            { binding: 0, resource: { buffer: orbitsBuf } },
            { binding: 1, resource: { buffer: orbitIdxBuf } },
            { binding: 2, resource: { buffer: phaseBuf } },
            { binding: 3, resource: { buffer: clustersBuf } },
            { binding: 4, resource: { buffer: targetBuf } }
        ]
    });

    // Group 2 carries the best-iterate grid, and group 1 gets an empty bind
    // group because keep_best's pipeline layout has a slot there that still
    // has to be set. See cfGetPipelines for why `best` is not in group 0.
    const bgBest = device.createBindGroup({
        layout: pl.layoutBest,
        entries: [{ binding: 0, resource: { buffer: bestBuf } }]
    });
    const bgEmpty = device.createBindGroup({ layout: pl.layoutEmpty, entries: [] });

    const res = { bufA, bufB, bgMain, bgSwapped, fftBindGroups };

    const pass1 = (enc, pipeline, groups, withSym) => {
        const p = enc.beginComputePass();
        p.setPipeline(pipeline);
        p.setBindGroup(0, bgMain);
        if (withSym) p.setBindGroup(1, bg1);
        p.dispatchWorkgroups(groups);
        p.end();
    };

    const wgOrbits = cfCeilGroups(model.nOrbits);
    const wgAbsent = cfCeilGroups(Math.max(1, model.absentCount));

    // Check every dispatch against the device limit BEFORE encoding anything.
    // WebGPU reports an oversized dispatch as an asynchronous validation
    // error: the command buffer is silently discarded, the grid is never
    // written, and the run completes in seconds with R = infinity instead of
    // failing. Throwing here instead sends the job to the CPU path with a
    // message that says which dispatch was too large.
    const maxWG = (device.limits && device.limits.maxComputeWorkgroupsPerDimension) || 65535;
    const dispatches = [
        ['FFT stage', cfCeilGroups((N >> 1) * N * N)],
        ['grid copy', cfCopyGroups(N)],
        ['full-grid pass', wgFull],
        ['cluster pass', wgCluster],
        ['orbit pass', wgOrbits],
        ['absence pass', wgAbsent]
    ];
    for (const [what, groups] of dispatches) {
        if (groups > maxWG) {
            throw new Error(
                `The ${what} would need ${groups} workgroups, but this device allows ` +
                `${maxWG} per dimension. Use a smaller grid.`);
        }
    }

    const rHistory = new Float32Array(maxIter).fill(NaN);
    let bestR = Infinity, bestIter = 0, lastReport = -1;

    // Everything from here on can throw and hand the job to the CPU path. The
    // allocations above are already registered with `track`, so the finally
    // block releases them whether the failure happens during allocation or
    // during the iteration. A 128^3 trial holds about 50 MB of device memory.

    /**
     * Map the staging buffer, copy `count` R-factor headers out of it and fold
     * them into rHistory. `firstIter` is the iteration the batch starts at.
     */
    const drainRBatch = async (firstIter, count) => {
        if (count <= 0) return;
        // Race the map against the device-lost promise and a watchdog. A lost
        // device never settles its pending maps, so awaiting one bare is an
        // unconditional hang; and a driver that is merely wedged reports
        // nothing at all, which the timeout turns into a real error instead of
        // a spinner that runs until the tab is closed.
        const WATCHDOG_MS = 30000;
        let timer = null;
        const guard = new Promise((_, reject) => {
            timer = setTimeout(() => reject(new Error(
                `The GPU stopped responding (no result for ${WATCHDOG_MS / 1000} s at ` +
                `grid ${N}). This usually means the driver reset the device; try a ` +
                `smaller grid or fewer iterations.`)), WATCHDOG_MS);
        });
        const lost = device.lost.then((info) => {
            throw new Error('WebGPU device lost during charge flipping: ' +
                            ((info && info.message) || (info && info.reason) || 'unknown') +
                            '. A smaller grid or fewer iterations usually avoids this.');
        });
        try {
            await Promise.race([headerRead.mapAsync(GPUMapMode.READ, 0, count * 16), lost, guard]);
        } finally {
            if (timer) clearTimeout(timer);
        }
        const hdr = new Float32Array(headerRead.getMappedRange(0, count * 16)).slice();
        headerRead.unmap();
        for (let k = 0; k < count; k++) {
            const num = hdr[4 * k + 1], den = hdr[4 * k + 2];
            const R = den > 0 ? num / den : NaN;
            rHistory[firstIter + k] = R;
            if (Number.isFinite(R) && R < bestR) { bestR = R; bestIter = firstIter + k + 1; }
        }
    };

    let batchStart = 0, batchCount = 0;
    let batchT0 = 0;              // when the first submit of this batch went out
    const batchState = {};        // running worst ms/iteration, for cfNextBatch
    let firstBatch = true;
    let lastBatchMs = 0;
    const nowMs = () => (typeof performance !== 'undefined' && performance.now)
        ? performance.now() : Date.now();

    for (let iter = 0; iter < maxIter; iter++) {
        // The first pass runs inside an error scope. Anything the driver
        // rejects -- an oversized dispatch, a buffer that is too big for this
        // adapter, a shader that failed to compile -- surfaces here as a
        // thrown error and sends the whole job to the CPU, rather than leaving
        // the loop to grind through maxIter cycles on a grid nothing ever
        // wrote to.
        if (iter === 0) device.pushErrorScope('validation');

        const enc = device.createCommandEncoder();

        cfEncodeFFT(enc, pl, res, N, true);              // F -> rho, in bufA
        pass1(enc, pl.reduce_moments, wgFull, false);
        pass1(enc, pl.finalize_moments, 1, false);       // writes partials[0] = delta
        pass1(enc, pl.flip, wgFull, false);
        cfEncodeFFT(enc, pl, res, N, false);             // rho -> F, in bufA

        pass1(enc, pl.reduce_rfactor, wgCluster, true);
        pass1(enc, pl.finalize_rfactor, 1, false);
        pass1(enc, pl.repartition, cfCeilGroups(model.nClusters), true);
        if (lambda > 0) pass1(enc, pl.symmetrize, wgOrbits, true);
        pass1(enc, pl.constrain, wgOrbits, true);
        if (model.absentCount > 0) pass1(enc, pl.zero_absent, wgAbsent, true);
        pass1(enc, pl.zero_dc, 1, false);

        // Snapshot the best iterate ON THE DEVICE. keep_best must run BEFORE
        // commit_best so every invocation compares against the old record.
        // This is what lets the R readback below be batched: the host no
        // longer has to see R in order to decide anything.
        //
        // Encoded by hand rather than through pass1(): keep_best is the only
        // kernel that binds three groups, because `best` sits in group 2.
        {
            const pB = enc.beginComputePass();
            pB.setPipeline(pl.keep_best);
            pB.setBindGroup(0, bgMain);
            pB.setBindGroup(1, bgEmpty);   // every group in the layout must be set
            pB.setBindGroup(2, bgBest);
            pB.dispatchWorkgroups(cfCopyGroups(N));
            pB.end();
        }
        pass1(enc, pl.commit_best, 1, false);

        enc.copyBufferToBuffer(partialsBuf, 0, headerRead, batchCount * 16, 16);
        if (batchCount === 0) batchT0 = nowMs();   // first submit of this batch
        batchCount++;
        device.queue.submit([enc.finish()]);

        if (iter === 0) {
            const gpuErr = await device.popErrorScope();
            if (gpuErr) throw new Error('WebGPU rejected the charge-flipping pass: ' + gpuErr.message);
        }

        const lastIteration = (iter === maxIter - 1);
        if (batchCount === rBatch || lastIteration) {
            const drained = batchCount;
            await drainRBatch(batchStart, batchCount);

            // Time the whole round trip -- submit to drained -- and resize the
            // next batch from it. This is the measurement the grid-size table
            // could not make: it reflects THIS device on THIS structure, so an
            // old integrated card converges to a small batch and a fast
            // discrete one to a large one, with no hardware list to maintain.
            lastBatchMs = nowMs() - batchT0;
            const next = cfNextBatch(drained, lastBatchMs, firstBatch, batchState);
            if (next !== rBatch) {
                console.info(`[charge flipping] sync batch ${rBatch} -> ${next} ` +
                             `(${(lastBatchMs / drained).toFixed(1)} ms per iteration)`);
            }
            rBatch = next;
            firstBatch = false;

            // A non-finite R in the very first batch means the GPU produced
            // nothing usable, whatever the driver did or did not report.
            if (batchStart === 0 && !Number.isFinite(rHistory[0])) {
                throw new Error('The first GPU cycle returned no usable R factor; the grid was not computed.');
            }

            batchStart = iter + 1;
            batchCount = 0;

            const pct = Math.floor(((iter + 1) / maxIter) * 100);
            if (pct !== lastReport || lastIteration) {
                lastReport = pct;
                postMessage({ type: 'cf-progress', iter: iter + 1, total: maxIter,
                              R: rHistory[iter], bestR });
            }
        }
    }
    const haveBest = Number.isFinite(bestR);

    // --- read the best iterate back exactly once ---------------------------
    const readbackBuf = track({
        size: cplxBytes, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
    });
    {
        const enc = device.createCommandEncoder();
        enc.copyBufferToBuffer(haveBest ? bestBuf : bufA, 0, readbackBuf, 0, cplxBytes);
        device.queue.submit([enc.finish()]);
    }
    await readbackBuf.mapAsync(GPUMapMode.READ);
    const factors = new Float32Array(readbackBuf.getMappedRange()).slice();
    readbackBuf.unmap();

    const re = new Float64Array(N3), im = new Float64Array(N3);
    for (let i = 0; i < N3; i++) { re[i] = factors[2 * i]; im[i] = factors[2 * i + 1]; }
    fft3d(re, im, N, true);
    const bestRho = Float32Array.from(re);

    let peak = 0;
    for (let i = 0; i < N3; i++) if (bestRho[i] > peak) peak = bestRho[i];
    if (peak > 0) for (let i = 0; i < N3; i++) bestRho[i] /= peak;

    const peaks = findPeaks(bestRho, N, cell, {
        minSeparation: job.minPeakSeparation,
        maxPeaks: job.maxPeaks,
        heightCutoff: job.peakHeightCutoff
    });

    return {
        map: bestRho, gridSize: N, peaks, rHistory, bestR, bestIter,
        finalR: rHistory[maxIter - 1], seed: job.seed, cell,
        volume: cellVolume(cell), backend: 'gpu', symLambda: lambda,
        syncBatch: rBatch, msPerIteration: lastBatchMs > 0 ? lastBatchMs / Math.max(1, rBatch) : null,
        reflections: modelReport(model)
    };

    } finally {
        gpuBuffers.forEach(b => { try { b.destroy(); } catch (e) { /* already gone */ } });
    }
}

// Exposed for the node-side unit test; harmless in a browser worker.
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        fft1d, fft3d, buildReflectionModel, findSystematicAbsences,
        lorentzPolarization, polarizationK, describePolarization,
        normalizeSymops, findPeaks,
        runChargeFlipping: runChargeFlippingCPU, runChargeFlippingCPU,
        metricTensor, fracDistance, buildStructure, findOriginShift,
        symmetryAverage, reduceToAsymmetricUnit, applyOp, LAUE_OPS
    };
}

// ===========================================================================
//  STRUCTURE BUILDING
//
//  With symLambda > 0 the map already obeys the space group and sits on a
//  standard origin, so the origin search below is a verification step rather
//  than a necessity. With symLambda = 0 it is doing the real work, exactly as
//  before. Either way:
//
//    1. find the origin shift that makes the map obey the symmetry operators
//    2. average the density over those operators
//    3. reduce the peaks to one representative per orbit, with multiplicity
//
//  Step 1 has a closed form. For x' = xR + t, a correctly positioned
//  structure satisfies
//
//        F(hR) = F(h) . exp(-2*pi*i * h.t)
//
//  Our map is displaced by an unknown s, so its phases carry an extra
//  -2*pi*h.s. Writing psi_g(h) = phi(hR) - phi(h) + 2*pi*h.t and
//  k_g(h) = hR - h (an integer vector), the agreement for a trial shift s is
//
//        A(s) = sum_h w(h) cos( psi_g(h) - 2*pi*k_g(h).s )
//             = Re{ sum_k C[k] exp(-2*pi*i*k.s) }     with C[k] = sum w e^{i psi}
//
//  which is a single forward FFT of C. One transform therefore scores EVERY
//  origin on the grid at once.
// ===========================================================================

function signedIndex(m, N) { return m > N / 2 ? m - N : m; }

function findOriginShift(re, im, N, symops, opts) {
    const N2 = N * N, N3 = N2 * N;
    const wrap = v => ((v % N) + N) % N;

    let maxF = 0;
    // PERF: see the note in the orbit loop -- naive magnitude, not Math.hypot.
    // This runs over the whole N^3 grid, so at 128^3 it is 2.1e6 calls.
    for (let i = 0; i < N3; i++) {
        const m = Math.sqrt(re[i] * re[i] + im[i] * im[i]);
        if (m > maxF) maxF = m;
    }
    const cut = maxF * (opts.amplitudeCutoff ?? 0.05);

    const cRe = new Float64Array(N3);
    const cIm = new Float64Array(N3);
    let contributions = 0;
    let sumW = 0;

    for (const op of symops) {
        const r = op.r, t = op.t;
        if (r[0] === 1 && r[4] === 1 && r[8] === 1 &&
            r[1] === 0 && r[2] === 0 && r[3] === 0 &&
            r[5] === 0 && r[6] === 0 && r[7] === 0 &&
            Math.abs(t[0]) < 1e-9 && Math.abs(t[1]) < 1e-9 && Math.abs(t[2]) < 1e-9) continue;

        for (let ml = 0; ml < N; ml++) {
            const hl = signedIndex(ml, N);
            for (let mk = 0; mk < N; mk++) {
                const hk = signedIndex(mk, N);
                for (let mh = 0; mh < N; mh++) {
                    const hh = signedIndex(mh, N);
                    if (!hh && !hk && !hl) continue;

                    const idx = mh + mk * N + ml * N2;
                    const a = re[idx], b = im[idx];
                    const mag = Math.sqrt(a * a + b * b);
                    if (mag < cut) continue;

                    const hp = applyRowVector([hh, hk, hl], r);
                    if (Math.abs(hp[0]) * 2 >= N || Math.abs(hp[1]) * 2 >= N || Math.abs(hp[2]) * 2 >= N) continue;

                    const jdx = wrap(hp[0]) + wrap(hp[1]) * N + wrap(hp[2]) * N2;
                    const c = re[jdx], d = im[jdx];
                    const mag2 = Math.sqrt(c * c + d * d);
                    if (mag2 < cut) continue;
                    
                    const delta = 2 * Math.PI * (hh * t[0] + hk * t[1] + hl * t[2]);
                    const cosD = Math.cos(delta);
                    const sinD = Math.sin(delta);
                    
                    const dot = a * c + b * d;
                    const cross = a * d - b * c;

                    sumW += mag * mag2;
                    const kIdx = wrap(hp[0] - hh) + wrap(hp[1] - hk) * N + wrap(hp[2] - hl) * N2;
                    cRe[kIdx] += dot * cosD - cross * sinD;
                    cIm[kIdx] += dot * sinD + cross * cosD;
                    contributions++;
                }
            }
        }
    }

    if (contributions === 0) return { candidates: [], sumW: 0, contributions: 0 };

    fft3d(cRe, cIm, N, false);

    // Rank the grid points and keep the strongest few as CANDIDATES rather
    // than trusting the single maximum: the maximum is genuinely degenerate
    // for groups whose rotations are all diagonal, and A(s) is maximised at
    // the NEGATIVE of the shift that symmetryAverage() needs.
    const order = [];
    for (let i = 0; i < N3; i++) order.push(i);
    order.sort((a, b) => cRe[b] - cRe[a]);

    const candidates = [];
    const keep = Math.max(1, opts.originCandidates ?? 12);
    for (const i of order) {
        const mh = i % N, mk = Math.floor(i / N) % N, ml = Math.floor(i / N2);
        candidates.push({
            shift: [((N - mh) % N) / N, ((N - mk) % N) / N, ((N - ml) % N) / N],
            score: sumW > 0 ? cRe[i] / sumW : 0
        });
        if (candidates.length >= keep) break;
    }

    return { candidates, sumW, contributions };
}

function sampleGrid(map, N, x, y, z) {
    const N2 = N * N;
    const fx = x * N, fy = y * N, fz = z * N;
    let ix = Math.floor(fx), iy = Math.floor(fy), iz = Math.floor(fz);
    const tx = fx - ix, ty = fy - iy, tz = fz - iz;

    ix = ((ix % N) + N) % N;
    iy = ((iy % N) + N) % N;
    iz = ((iz % N) + N) % N;

    const ix1 = (ix + 1) % N;
    const iy1 = (iy + 1) % N;
    const iz1 = (iz + 1) % N;

    const wx0 = 1 - tx, wx1 = tx;
    const wy0 = 1 - ty, wy1 = ty;
    const wz0 = 1 - tz, wz1 = tz;

    const z0 = iz * N2, z1 = iz1 * N2;
    const y0 = iy * N,  y1 = iy1 * N;

    return wx0 * wy0 * wz0 * map[ix  + y0 + z0] +
           wx1 * wy0 * wz0 * map[ix1 + y0 + z0] +
           wx0 * wy1 * wz0 * map[ix  + y1 + z0] +
           wx1 * wy1 * wz0 * map[ix1 + y1 + z0] +
           wx0 * wy0 * wz1 * map[ix  + y0 + z1] +
           wx1 * wy0 * wz1 * map[ix1 + y0 + z1] +
           wx0 * wy1 * wz1 * map[ix  + y1 + z1] +
           wx1 * wy1 * wz1 * map[ix1 + y1 + z1];
}

function symmetryAverage(map, N, symops, shift) {
    const N2 = N * N, N3 = N2 * N;
    const out = new Float32Array(N3);
    for (let l = 0; l < N; l++) {
        for (let k = 0; k < N; k++) {
            for (let hh = 0; hh < N; hh++) {
                const x = hh / N, y = k / N, z = l / N;
                let acc = 0;
                for (const op of symops) {
                    const r = op.r, t = op.t;
                    const px = r[0] * x + r[1] * y + r[2] * z + t[0];
                    const py = r[3] * x + r[4] * y + r[5] * z + t[1];
                    const pz = r[6] * x + r[7] * y + r[8] * z + t[2];
                    acc += sampleGrid(map, N, px + shift[0], py + shift[1], pz + shift[2]);
                }
                out[hh + k * N + l * N2] = acc / symops.length;
            }
        }
    }
    return out;
}

function symmetryCorrelation(map, sym, N, shift) {
    const N2 = N * N, N3 = N2 * N;
    let sa = 0, sb = 0;
    const a = new Float64Array(N3);
    for (let l = 0; l < N; l++) for (let k = 0; k < N; k++) for (let hh = 0; hh < N; hh++) {
        const i = hh + k * N + l * N2;
        a[i] = sampleGrid(map, N, hh / N + shift[0], k / N + shift[1], l / N + shift[2]);
        sa += a[i]; sb += sym[i];
    }
    sa /= N3; sb /= N3;
    let num = 0, da = 0, db = 0;
    for (let i = 0; i < N3; i++) {
        const u = a[i] - sa, v = sym[i] - sb;
        num += u * v; da += u * u; db += v * v;
    }
    return (da > 0 && db > 0) ? num / Math.sqrt(da * db) : 0;
}

function applyOp(op, p) {
    const r = op.r, t = op.t;
    return [
        r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + t[0],
        r[3] * p[0] + r[4] * p[1] + r[5] * p[2] + t[1],
        r[6] * p[0] + r[7] * p[1] + r[8] * p[2] + t[2]
    ];
}

function reduceToAsymmetricUnit(peaks, symops, cell, tolAngstrom) {
    const G = metricTensor(cell);
    const used = new Array(peaks.length).fill(false);
    const unique = [];

    for (let i = 0; i < peaks.length; i++) {
        if (used[i]) continue;
        const p = peaks[i];
        used[i] = true;

        const images = [];
        for (const op of symops) {
            const q = applyOp(op, [p.x, p.y, p.z]).map(v => ((v % 1) + 1) % 1);
            if (!images.some(im => fracDistance(G, im[0] - q[0], im[1] - q[1], im[2] - q[2]) < tolAngstrom)) {
                images.push(q);
            }
        }

        // Average over the symmetry images that were actually found, for the
        // integral as well as the height. Averaging matters more here: the
        // images of one site should carry identical charge, so the spread
        // between them is a direct measure of how noisy the integral is, and
        // the mean is the better estimate.
        let heightSum = p.height, heightN = 1;
        let chargeSum = Number.isFinite(p.charge) ? p.charge : 0;
        for (let j = i + 1; j < peaks.length; j++) {
            if (used[j]) continue;
            const q = peaks[j];
            if (images.some(im => fracDistance(G, im[0] - q.x, im[1] - q.y, im[2] - q.z) < tolAngstrom)) {
                used[j] = true;
                heightSum += q.height; heightN++;
                if (Number.isFinite(q.charge)) chargeSum += q.charge;
            }
        }

        unique.push({
            x: p.x, y: p.y, z: p.z,
            height: heightSum / heightN,
            charge: chargeSum / heightN,
            multiplicity: images.length,
            observedImages: heightN,
            special: images.length < symops.length
        });
    }

    unique.sort((a, b) => b.height - a.height);
    const top = unique.length ? unique[0].height : 1;
    let topQ = 0;
    for (const u of unique) if (u.charge > topQ) topQ = u.charge;
    return unique.map((u, i) => ({
        ...u, rank: i + 1,
        relative: u.height / top,
        chargeRel: topQ > 0 ? u.charge / topQ : 0
    }));
}

async function buildStructure(job) {
    const N = job.gridSize;
    const N3 = N * N * N;
    const symops = normalizeSymops(job.symops);
    const cell = job.cell;

    if (!symops) {
        return { error: 'No usable symmetry operators were supplied. The space-group JSON needs the "symops" field added by the v7 generator script.' };
    }

    const re = new Float64Array(N3);
    const im = new Float64Array(N3);
    for (let i = 0; i < N3; i++) re[i] = job.map[i];
    fft3d(re, im, N, false);

    const invertMap = (m) => {
        const wrap = v => ((v % N) + N) % N;
        const N2 = N * N;
        const inv = new Float32Array(N3);
        for (let l = 0; l < N; l++) for (let k = 0; k < N; k++) for (let hh = 0; hh < N; hh++) {
            inv[hh + k * N + l * N2] = m[wrap(-hh) + wrap(-k) * N + wrap(-l) * N2];
        }
        return inv;
    };

    // Charge flipping cannot distinguish rho(x) from rho(-x), so both hands
    // are searched. For each hand the origin search returns several candidate
    // shifts; every one is then symmetrised AND SCORED.
    let chosen = null;
    const tried = [];
    for (const hand of [1, -1]) {
        const workMap = hand === 1 ? job.map : invertMap(job.map);
        const r2 = Float64Array.from(re);
        const i2 = new Float64Array(N3);
        for (let i = 0; i < N3; i++) i2[i] = hand * im[i];
        const found = findOriginShift(r2, i2, N, symops, job);

        for (const cand of found.candidates) {
            const sym = symmetryAverage(workMap, N, symops, cand.shift);
            const corr = symmetryCorrelation(workMap, sym, N, cand.shift);
            tried.push({ hand, shift: cand.shift, score: cand.score, correlation: corr });
            if (!chosen || corr > chosen.correlation) {
                chosen = { hand, shift: cand.shift, score: cand.score, correlation: corr, sym, workMap };
            }
        }
    }

    if (!chosen) {
        return { error: 'The origin search found nothing to work with. The map may be empty, or the amplitude cutoff too high.' };
    }

    tried.sort((a, b) => b.correlation - a.correlation);
    const sym = chosen.sym;
    const correlation = chosen.correlation;

    let peak = 0;
    for (let i = 0; i < N3; i++) if (sym[i] > peak) peak = sym[i];
    if (peak > 0) for (let i = 0; i < N3; i++) sym[i] /= peak;

    let sites = [];
    let rawPeakCount = 0;
    let wyckoffResult = null;

    if (job.targetComposition) {
        // Execute the Replica-Exchange Monte Carlo swarm against the converged density map
        wyckoffResult = await runWyckoffDensityFit(sym, N, cell, symops, job);
        sites = wyckoffResult.sites || [];
        rawPeakCount = sites.length;
    } else {
        const rawPeaks = findPeaks(sym, N, cell, {
            minSeparation: job.minPeakSeparation,
            maxPeaks: job.maxPeaks ?? 60,
            heightCutoff: job.peakHeightCutoff ?? 0.04
        });
        rawPeakCount = rawPeaks.length;
        sites = reduceToAsymmetricUnit(rawPeaks, symops, cell, job.minPeakSeparation * 0.6);
    }

    return {
        map: sym,
        gridSize: N,
        originShift: chosen.shift,
        hand: chosen.hand,
        originScore: chosen.score,
        runnerUpCorrelation: tried.length > 1 ? tried[1].correlation : null,
        symmetryCorrelation: correlation,
        sites,
        rawPeakCount: rawPeakCount,
        nOps: symops.length,
        symops: symops.map(o => o.xyz).filter(Boolean),
        cell,
        wyckoffResult
    };
}

// ===========================================================================
//  WYCKOFF DENSITY FIT
//  Fits exact stoichiometry and Wyckoff constraints against the 
//  converged electron density map.
// ===========================================================================

// ===========================================================================
//  WYCKOFF DENSITY FIT ORCHESTRATOR
// ===========================================================================
const CF_ATOM_DATA = {
    "H":  { z: 1,  r: 1.20, m: 1.008, rc: 0.31   }, "HE": { z: 2,  r: 1.40, m: 4.003, rc: 0.28   }, "LI": { z: 3,  r: 1.82, m: 6.941, rc: 1.28   }, "BE": { z: 4,  r: 1.53, m: 9.012, rc: 0.96   },
    "B":  { z: 5,  r: 1.92, m: 10.811, rc: 0.84  }, "C":  { z: 6,  r: 1.70, m: 12.011, rc: 0.76  }, "N":  { z: 7,  r: 1.55, m: 14.007, rc: 0.71  }, "O":  { z: 8,  r: 1.52, m: 15.999, rc: 0.66  },
    "F":  { z: 9,  r: 1.47, m: 18.998, rc: 0.57  }, "NE": { z: 10, r: 1.54, m: 20.180, rc: 0.58  }, "NA": { z: 11, r: 2.27, m: 22.990, rc: 1.66  }, "MG": { z: 12, r: 1.73, m: 24.305, rc: 1.41  },
    "AL": { z: 13, r: 1.18, m: 26.982, rc: 1.21  }, "SI": { z: 14, r: 2.10, m: 28.085, rc: 1.11  }, "P":  { z: 15, r: 1.80, m: 30.974, rc: 1.07  }, "S":  { z: 16, r: 1.80, m: 32.06, rc: 1.05   },
    "CL": { z: 17, r: 1.75, m: 35.45, rc: 1.02   }, "AR": { z: 18, r: 1.88, m: 39.948, rc: 1.06  }, "K":  { z: 19, r: 2.75, m: 39.098, rc: 2.03  }, "CA": { z: 20, r: 2.31, m: 40.078, rc: 1.76  },
    "SC": { z: 21, r: 2.11, m: 44.956, rc: 1.70  }, "TI": { z: 22, r: 2.00, m: 47.867, rc: 1.60  }, "V":  { z: 23, r: 1.90, m: 50.942, rc: 1.53  }, "CR": { z: 24, r: 1.90, m: 51.996, rc: 1.39  },
    "MN": { z: 25, r: 1.90, m: 54.938, rc: 1.50  }, "FE": { z: 26, r: 1.26, m: 55.845, rc: 1.42  }, "CO": { z: 27, r: 1.90, m: 58.933, rc: 1.38  }, "NI": { z: 28, r: 1.63, m: 58.693, rc: 1.24  },
    "CU": { z: 29, r: 1.40, m: 63.546, rc: 1.32  }, "ZN": { z: 30, r: 1.39, m: 65.38, rc: 1.22   }, "GA": { z: 31, r: 1.87, m: 69.723, rc: 1.22  }, "GE": { z: 32, r: 2.11, m: 72.630, rc: 1.20  },
    "AS": { z: 33, r: 1.85, m: 74.922, rc: 1.19  }, "SE": { z: 34, r: 1.90, m: 78.971, rc: 1.20  }, "BR": { z: 35, r: 1.85, m: 79.904, rc: 1.20  }, "KR": { z: 36, r: 2.02, m: 83.798, rc: 1.16  },
    "RB": { z: 37, r: 3.03, m: 85.468, rc: 2.20  }, "SR": { z: 38, r: 2.49, m: 87.62, rc: 1.95   }, "Y":  { z: 39, r: 2.19, m: 88.906, rc: 1.90  }, "ZR": { z: 40, r: 2.06, m: 91.224, rc: 1.75  },
    "NB": { z: 41, r: 1.98, m: 92.906, rc: 1.64  }, "MO": { z: 42, r: 1.90, m: 95.95, rc: 1.54   }, "TC": { z: 43, r: 1.83, m: 98.0, rc: 1.47    }, "RU": { z: 44, r: 1.78, m: 101.07, rc: 1.46  },
    "RH": { z: 45, r: 1.73, m: 102.906, rc: 1.42 }, "PD": { z: 46, r: 1.63, m: 106.42, rc: 1.39  }, "AG": { z: 47, r: 1.72, m: 107.868, rc: 1.45 }, "CD": { z: 48, r: 1.58, m: 112.411, rc: 1.44 },
    "IN": { z: 49, r: 1.93, m: 114.818, rc: 1.42 }, "SN": { z: 50, r: 2.17, m: 118.71, rc: 1.39  }, "SB": { z: 51, r: 2.06, m: 121.760, rc: 1.39 }, "TE": { z: 52, r: 2.06, m: 127.6, rc: 1.38   },
    "I":  { z: 53, r: 1.98, m: 126.904, rc: 1.39 }, "XE": { z: 54, r: 2.16, m: 131.293, rc: 1.40 }, "CS": { z: 55, r: 3.43, m: 132.905, rc: 2.44 }, "BA": { z: 56, r: 2.68, m: 137.327, rc: 2.15 },
    "LA": { z: 57, r: 2.43, m: 138.905, rc: 2.07 }, "CE": { z: 58, r: 2.42, m: 140.116, rc: 2.04 }, "PR": { z: 59, r: 2.40, m: 140.908, rc: 2.03 }, "ND": { z: 60, r: 2.39, m: 144.242, rc: 2.01 },
    "PM": { z: 61, r: 2.38, m: 145.0, rc: 1.99   }, "SM": { z: 62, r: 2.36, m: 150.36, rc: 1.98  }, "EU": { z: 63, r: 2.35, m: 151.964, rc: 1.98 }, "GD": { z: 64, r: 2.34, m: 157.25, rc: 1.96  }, 
    "TB": { z: 65, r: 2.33, m: 158.925, rc: 1.94 }, "DY": { z: 66, r: 2.31, m: 162.500, rc: 1.92 }, "HO": { z: 67, r: 2.30, m: 164.930, rc: 1.92 }, "ER": { z: 68, r: 2.29, m: 167.259, rc: 1.89 }, 
    "TM": { z: 69, r: 2.27, m: 168.934, rc: 1.90 }, "YB": { z: 70, r: 2.26, m: 173.054, rc: 1.87 }, "LU": { z: 71, r: 2.24, m: 174.967, rc: 1.87 }, "HF": { z: 72, r: 2.12, m: 178.49, rc: 1.75  }, 
    "TA": { z: 73, r: 2.06, m: 180.948, rc: 1.70 }, "W":  { z: 74, r: 2.02, m: 183.84, rc: 1.62  }, "RE": { z: 75, r: 1.97, m: 186.207, rc: 1.51 }, "OS": { z: 76, r: 1.92, m: 190.23, rc: 1.44  }, 
    "IR": { z: 77, r: 1.87, m: 192.217, rc: 1.41 }, "PT": { z: 78, r: 1.75, m: 195.084, rc: 1.36 }, "AU": { z: 79, r: 1.66, m: 196.967, rc: 1.36 }, "HG": { z: 80, r: 1.55, m: 200.592, rc: 1.32 }, 
    "TL": { z: 81, r: 1.96, m: 204.38, rc: 1.45  }, "PB": { z: 82, r: 2.02, m: 207.2, rc: 1.46   }, "BI": { z: 83, r: 2.07, m: 208.980, rc: 1.48 }, "PO": { z: 84, r: 1.97, m: 209.0, rc: 1.40   }, 
    "AT": { z: 85, r: 2.02, m: 210.0, rc: 1.50   }, "RN": { z: 86, r: 2.20, m: 222.018, rc: 1.50 }, "FR": { z: 87, r: 3.48, m: 223.0, rc: 2.60   }, "RA": { z: 88, r: 2.83, m: 226.0, rc: 2.21   }, 
    "AC": { z: 89, r: 2.47, m: 227.0, rc: 2.15   }, "TH": { z: 90, r: 2.37, m: 232.038, rc: 2.06 }, "PA": { z: 91, r: 2.43, m: 231.036, rc: 2.00 }, "U":  { z: 92, r: 1.86, m: 238.029, rc: 1.96 }
};

/**
 * Pawley peaks -> the observation rows the coordinate refinement needs.
 *
 * DELEGATES to normaliseObservations() rather than doing the arithmetic here.
 * That is the point: it is the same function the reflection route calls, so
 * the coordinate refinement is guaranteed to be fitting the same |Fo|^2 the
 * rest of the program means by that symbol. Three things went wrong when this
 * was hand-rolled, and all three are silent - they produce finite numbers and
 * a refinement that converges somewhere slightly wrong:
 *
 *   Lp WHEN THE HOST DID NOT SUPPLY IT. Falling back to lp = 1 leaves the
 *   Lorentz-polarisation factor in the observation. Lp spans about a factor
 *   of 47 across a pattern, so this is not a small correction - it tilts the
 *   whole intensity distribution towards low angle. The map builder computes
 *   Lp from 2-theta in exactly this case, and so must this.
 *
 *   Lp WHEN IT HAS ALREADY BEEN REMOVED. job.applyLP === false means the
 *   host handed over intensities that are already Lp-free. Dividing again
 *   applies the correction twice, in the opposite direction to the first
 *   error and no less wrong.
 *
 *   MULTIPLICITY. The file's own m is authoritative, because it is the value
 *   the Pawley extraction assumed when it partitioned the pattern.
 *   Recomputing it from the operators gives a number that is usually the same
 *   and occasionally is not, and where it differs the extraction is right and
 *   the recomputation is wrong.
 *
 * Returns the whole normaliseObservations result, so the caller can report
 * which route was taken and surface any warnings instead of discarding them.
 */
/**
 * The raw observation rows, before normalisation.
 *
 * cfObservationsFromHkl builds these and then hands them to
 * normaliseObservations. runWyckoffSearch wants them UNNORMALISED, because it
 * runs the normalisation itself and derives the Wilson B from the result --
 * so this shares the row construction rather than duplicating it, and the map
 * and the search cannot end up disagreeing about Lp.
 *
 * @param {object} job
 * @returns {object[]} rows carrying h, k, l, Ihkl, lp, multiplicity, d, tth
 */
function cfRawObservationRows(job) {
    const applyLP = job.applyLP !== false;
    const polModel = (job.polarization !== undefined && job.polarization !== null)
        ? job.polarization : job.monochromatorTth;

    const rows = [];
    for (const peak of job.hklList || []) {
        if (!peak) continue;
        const I = Number(peak.intensity);
        const h = peak.h_orig, k = peak.k_orig, l = peak.l_orig;
        if (!Number.isFinite(I) || I < 0) continue;
        if (!Number.isFinite(h) || !Number.isFinite(k) || !Number.isFinite(l)) continue;
        if (h === 0 && k === 0 && l === 0) continue;
        // Ten bits per index in the kernel's packing, biased by 512.
        if (Math.abs(h) > 511 || Math.abs(k) > 511 || Math.abs(l) > 511) continue;

        rows.push({
            h, k, l, Ihkl: I,
            multiplicity: peak.multiplicity,
            d: peak.d, twoTheta: peak.tth,
            sigma: peak.sigma,
            lp: applyLP
                ? ((Number.isFinite(peak.lp) && peak.lp > 0) ? peak.lp
                                                            : lorentzPolarization(peak.tth, polModel))
                : 1
        });
    }
    return rows;
}

function cfObservationsFromHkl(job, symops) {
    const applyLP = job.applyLP !== false;
    const polModel = (job.polarization !== undefined && job.polarization !== null)
        ? job.polarization : job.monochromatorTth;

    const src = [];
    for (const peak of job.hklList || []) {
        if (!peak) continue;
        const I = Number(peak.intensity);
        const h = peak.h_orig, k = peak.k_orig, l = peak.l_orig;
        if (!Number.isFinite(I) || I < 0) continue;
        if (!Number.isFinite(h) || !Number.isFinite(k) || !Number.isFinite(l)) continue;
        if (h === 0 && k === 0 && l === 0) continue;

        const row = { h, k, l, Ihkl: I, multiplicity: peak.multiplicity,
                      d: peak.d, twoTheta: peak.tth };
        if (applyLP) {
            // Prefer the host's Lp - it was evaluated at the refined 2-theta
            // including the zero shift, and is the same number the map was
            // built from. Fall back to computing it HERE, with this file's own
            // lorentzPolarization(), rather than leaving the field absent and
            // letting normaliseObservations reconstruct it. Two reasons, and
            // both matter:
            //
            //   It is the identical call buildReflectionModel() makes a few
            //   hundred lines up, so the map and the refinement cannot end up
            //   with different Lp for the same reflection.
            //
            //   normaliseObservations' own fallback route calls
            //   sharkoLorentzPolarization(), which lives in crystal.js and is
            //   NOT in this worker's scope - it would throw. Supplying lp on
            //   every row keeps it on the m_lp route, which needs nothing
            //   beyond what is already here.
            row.lp = (Number.isFinite(peak.lp) && peak.lp > 0)
                ? peak.lp
                : lorentzPolarization(peak.tth, polModel);
        } else {
            // Already Lp-free: divide out the multiplicity and nothing else.
            row.lp = 1;
        }
        src.push(row);
    }
    if (!src.length) return { rows: [], route: null, warnings: [], errors: ['No usable Pawley intensities.'] };

    return normaliseObservations(src, {
        rotations: (symops || []).map(op => op.r),
        wavelength: job.wavelength,
        radiation: job.radiation || 'xray',
        polarisationK: polarizationK(polModel)
    });
}

async function runWyckoffDensityFit(densityMap, N, cell, symops, job) {

    const formula = job.targetComposition;
    if (!formula) return { sites: [] };

    let device;
    try {
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) throw new Error('WebGPU unavailable in worker.');
        device = await adapter.requestDevice();
    } catch (e) {
        console.error('Wyckoff Density Fit: WebGPU error:', e);
        return { sites: [] }; 
    }

    // WHICH OBJECTIVE. With a map, score by density overlap; without one,
    // score |Fcalc|^2 against the Pawley intensities. The Wyckoff search does
    // not NEED a map -- it needs a space group, a composition, and the
    // observations, all of which come from the Pawley fit -- so the mapless
    // path is the general one and the density path is the special case.
    const shaderFile = densityMap ? '../swarm_density.wgsl' : '../swarm_reflection.wgsl';
    const shaderResp = await fetch(shaderFile);
    if (!shaderResp.ok) throw new Error(shaderFile + ' could not be loaded.');
    const swarmShaderSrc = await shaderResp.text();

    // Raw observation rows for the reflection objective. runWyckoffSearch
    // normalises them itself and estimates the Wilson B from them, so they are
    // handed over unprocessed rather than pre-reduced.
    let rawRows = [];
    if (!densityMap) {
        rawRows = cfRawObservationRows(job);
        if (!rawRows.length) {
            return { sites: [], error: 'No usable Pawley intensities for the Wyckoff search.' };
        }
    }

    // Loaded for BOTH objectives, not just the mapless one, because the FINAL
    // REFINEMENT needs them too and could not get them: job.formFactor is a
    // function, and a function cannot cross a postMessage boundary, so
    // coord_refine's `o.formFactor || (() => null)` fell through to f = Z on
    // every single run.
    //
    // That is the discrepancy behind a search reporting CC = 0.99 and the
    // refinement then reporting R = 50%: THEY WERE NOT SCORING THE SAME MODEL.
    // The search used tabulated f(s) with its angular fall-off; the refinement
    // used a flat atomic number, which overweights every heavy atom at high
    // angle. No arrangement of atoms satisfies both descriptions at once, so
    // the refinement moved the structure away from the search's answer and
    // still could not fit.
    //
    // '../scatters', not the default 'scatters': a relative fetch in a worker
    // resolves against the WORKER script's URL, and this worker lives in js/.
    let scatterTables = null;
    try {
        scatterTables = await loadScatteringTables('../scatters');
    } catch (e) {
        console.warn('[CF Wyckoff] No scattering tables; falling back to f = Z.', e);
    }
    const refineFormFactor = (scatterTables && typeof makeFormFactor === 'function')
        ? makeFormFactor(scatterTables, { radiation: job.radiation || 'xray', ions: job.ions })
        : null;

    const stg = { sym_ops: symops, wyckoff: job.wyckoffTable || [] };
    
    let Z = job.wyckoffZ;
    if (!Z || isNaN(Z)) {
        const comp = parseFormula(formula, 1);
        const vol = cellVolume(cell);
        const perFormula = comp.reduce((n, c) => n + c.count, 0);
        const sug = suggestZ(vol, perFormula, symops.length);
        Z = sug ? sug.likely : 1;
    }


    postMessage({ type: 'cf-wyckoff-start', message: 'Compiling WebGPU Swarm...' });

    let out;
    try {
    out = await runWyckoffSearch({
            device, 
            ccShaderSource: swarmShaderSrc,
            // No groupStride here: swPackReflections decides it, and
            // hard-coding 3 would have overridden the weighted layout
            // it now emits, leaving the kernel reading Iobs where the
            // weight should be.

            setting: stg,
            cell,
            reflections: rawRows,
            wavelength: job.wavelength,
            radiation: job.radiation || 'xray',
            polarisationK: polarizationK((job.polarization !== undefined && job.polarization !== null)
                                          ? job.polarization : job.monochromatorTth),
            densityMap: densityMap,
            gridSize: N,
            formula, 
            Z: Z,
            maxSitesPerElement: job.wyckoffMaxSites || undefined,
            maxRepeat: job.wyckoffMaxRepeat || undefined,
            wyckoffCapCeiling: 24,
            atomData: CF_ATOM_DATA, 
            scatterTables: scatterTables,
            // Accept the pre-parsed array directly from the UI
            windows: job.distanceConstraints || [],
            harkerSites: [],
            numParticles: job.swarmParticles || 512,
            generations: job.swarmIterations || 1000,
            restarts: job.swarmRestarts || 4,
            minContact: job.minContact || 1.0,
onLog: m => console.log('[CF Wyckoff Swarm]', m),
            onProgress: p => {
                postMessage({ 
                    type: 'cf-wyckoff-progress', 
                    wave: p.wave, waves: p.waves, 
                    generation: p.generation, generations: p.generations, 
                    restart: p.restart, restarts: p.restarts, 
                    best: p.best 
                });
            },
            shouldStop: () => false
        });
} catch (e) {
        console.error('Wyckoff Density Fit failed during execution:', e);
        return { sites: [] };
    } finally {
        if (device) device.destroy();
    }

if (out && out.candidates && out.candidates.length > 0) {
        const G = metricTensor(cell);

        // ------------------------------------------------------------------
        //  A CONTACT FLOOR IS A CONSTRAINT, NOT A PREFERENCE.
        //
        //  The kernel charges a clash as a soft penalty, scaled by a ramp that
        //  starts gentle so the swarm can move through crowded configurations
        //  early. That is the right behaviour DURING the search. It is the
        //  wrong behaviour for the answer: a solution containing a 0.47 A Pb-O
        //  contact is not a slightly worse structure, it is not a structure,
        //  and no correlation coefficient should be able to buy it.
        //
        //  runWyckoffSearch returns `floors` precisely so the caller can reject
        //  on the same numbers the search enforced -- and nothing called it.
        //  Candidates are now walked in score order and the first one whose
        //  contacts are all legal is taken. Rejections are reported rather than
        //  silently dropped: if every candidate fails, that is a statement
        //  about the composition or the floor, and the user needs to see it.
        // ------------------------------------------------------------------
        let cfWorstContact = null;
        const floorOf = (typeof out.floors === 'function')
            ? out.floors
            : () => (Number.isFinite(job.minContact) ? job.minContact : 1.0);

        const worstContact = (sites) => {
            // Expand to the full cell, then take the shortest unlike-pair
            // distance as a fraction of that pair's floor. < 1 means illegal.
            const exp = [];
            for (const st of sites || []) {
                for (const op of symops) {
                    // applyOp, NOT op.apply. The latter is the RECIPROCAL
                    // transform used for hkl, which is the transpose of the
                    // rotation and carries no translation -- correct for
                    // indices, silently wrong for coordinates.
                    const p = applyOp(op, [st.x, st.y, st.z])
                        .map(v => ((v % 1) + 1) % 1);
                    if (!exp.some(q => q.el === st.element &&
                                       fracDistance(G, q.p[0] - p[0], q.p[1] - p[1], q.p[2] - p[2]) < 1e-3)) {
                        exp.push({ el: st.element, p });
                    }
                }
            }
            // fracDistance handles the minimum image itself, including the
            // 27-point search an oblique cell needs. An earlier version did
            // that search HERE, passing dx+1 into fracDistance -- which
            // immediately rounded it back, making the whole loop a no-op and
            // the test that 'verified' it meaningless, since both sides were
            // wrapping.
            let worst = Infinity, pair = null;
            for (let i = 0; i < exp.length; i++) {
                for (let j = i + 1; j < exp.length; j++) {
                    let dx = exp[i].p[0] - exp[j].p[0];
                    let dy = exp[i].p[1] - exp[j].p[1];
                    let dz = exp[i].p[2] - exp[j].p[2];
                    dx -= Math.round(dx); dy -= Math.round(dy); dz -= Math.round(dz);
                    const d = fracDistance(G, dx, dy, dz);
                    const fl = floorOf(exp[i].el, exp[j].el);
                    if (!(fl > 0)) continue;
                    const ratio = d / fl;
                    if (ratio < worst) {
                        worst = ratio;
                        pair = { a: exp[i].el, b: exp[j].el, d, floor: fl };
                    }
                }
            }
            return { worst, pair };
        };

        // Kept in scope: the refinement below has to be judged by the SAME
        // test the candidates were, or the filter only guards the seed.
        cfWorstContact = worstContact;

        let best = null;
        const contactRejects = [];
        for (const cand of out.candidates) {
            const chk = worstContact(cand.sites);
            if (chk.worst >= 1.0 || !chk.pair) { best = cand; break; }
            contactRejects.push(
                `${cand.assignment}: ${chk.pair.a}-${chk.pair.b} at ${chk.pair.d.toFixed(2)} A ` +
                `against a floor of ${chk.pair.floor.toFixed(2)} A`);
        }
        if (!best) {
            return {
                sites: [],
                error: `Every candidate breaks a minimum-contact distance. ` +
                       `Closest offenders: ${contactRejects.slice(0, 3).join('; ')}. ` +
                       `Either the composition or Z is wrong for this cell, or the contact ` +
                       `floor is set higher than the structure allows.`
            };
        }
        if (contactRejects.length) {
            console.warn(`[CF Wyckoff] ${contactRejects.length} higher-scoring candidate(s) ` +
                         `rejected on contact distance: ${contactRejects.slice(0, 3).join('; ')}`);
        }

        // ------------------------------------------------------------------
        // FINAL COORDINATE REFINEMENT, AGAINST THE DATA RATHER THAN THE MAP.
        //
        // Everything above scored structures by the density underneath them,
        // sampled from an N^3 grid by trilinear interpolation - and a
        // multilinear interpolant takes its maximum over a cell AT A VERTEX,
        // so an atom sitting on real density is pulled onto a grid node and
        // pinned there. The coordinates arriving here are quantised to 1/N,
        // about 0.26 A at N = 32 and 0.13 A at N = 64. Raising N halves the
        // error and costs eight times the memory; it does not remove the
        // bias, because the interpolant has the same shape at every N.
        //
        // The map was only ever a device for getting out of the phase
        // problem. The measurement is the Pawley intensity list, so the last
        // step drops the grid entirely and refines the free coordinates
        // against |Fo|^2 by least squares, with every site re-projected onto
        // its Wyckoff position at each step. Failure here is not fatal: the
        // unrefined coordinates are still the answer the map gave, so the
        // reason is recorded and the pipeline carries on.
        // ------------------------------------------------------------------
        // THE SEARCH RESULT, before the refinement is allowed to touch it.
        //
        // The refinement mutates best.sites in place, so without this copy the
        // Monte Carlo answer is unrecoverable by the time anything is
        // reported -- and the two are worth comparing: they optimise different
        // quantities, and when they disagree that disagreement is the
        // interesting part.
        const searchSites = best.sites.map(x => ({ ...x }));

        let refinement = null;
        try {
            const obs = cfObservationsFromHkl(job, symops);
            const obsRows = obs.rows || [];
            for (const w of (obs.warnings || [])) console.warn('[CF refine] ' + w);
            for (const e of (obs.errors || [])) console.warn('[CF refine] ' + e);
            if (obsRows.length) {
                const seed = best.sites.map(s => ({
                    element: s.element,
                    zn: Number.isFinite(s.zn)
                        ? s.zn
                        : (CF_ATOM_DATA[String(s.element).toUpperCase()] || {}).z,
                    x: s.x, y: s.y, z: s.z, w: s.w
                }));
                if (seed.every(s => Number.isFinite(s.zn) && s.w)) {
                    // out.overallB is the WILSON B the search itself used to
                    // build its scattering table. job.overallB is only set if
                    // the user typed one, so this was refining with B = 0 on
                    // every ordinary run -- fitting sharp point atoms to data
                    // that falls off, which loads the error onto the light
                    // atoms first. The search and the refinement must use the
                    // same B or they are not fitting the same model.
                    const refB = Number.isFinite(job.overallB) ? job.overallB
                               : (Number.isFinite(out.overallB) ? out.overallB : 0);
                    const r = refineCoordinatesAgainstPawley({
                        sites: seed, symOps: symops, obsRows,
                        overallB: refB, formFactor: job.formFactor || refineFormFactor
                    });

                    // ------------------------------------------------------
                    //  THE REFINEMENT IS SUBJECT TO THE CONTACT FLOOR TOO.
                    //
                    //  The candidate filter runs on the SEARCH output, and the
                    //  refinement then moves every free coordinate with no
                    //  contact term at all -- so it could and did walk a
                    //  legal seed into an illegal structure: a candidate that
                    //  passed at 1.4 A came out with S-O at 0.85 A against a
                    //  floor of 1.35, because shortening that contact happened
                    //  to lower R. Guarding only the input to a step that
                    //  moves the atoms guards nothing.
                    //
                    //  A refinement that breaks a floor is discarded rather
                    //  than clamped. Clamping would leave the coordinates at
                    //  an arbitrary point on the constraint surface and report
                    //  an R that belongs to neither position; the seed is at
                    //  least a structure the search actually scored.
                    // ------------------------------------------------------
                    let contactBreak = null;
                    if (r && r.sites && !r.error && cfWorstContact) {
                        const chk = cfWorstContact(r.sites);
                        if (chk.pair && chk.worst < 1.0) contactBreak = chk.pair;
                    }

                    if (contactBreak) {
                        // Not a `return`: this runs inside the try block of
                        // runWyckoffDensityFit, and returning here would skip
                        // the site formatting and hand back an empty result.
                        console.warn(`[CF refine] Refinement rejected: it produced ` +
                            `${contactBreak.a}-${contactBreak.b} at ${contactBreak.d.toFixed(2)} A ` +
                            `against a floor of ${contactBreak.floor.toFixed(2)} A. ` +
                            `Keeping the unrefined positions.`);
                        refinement = {
                            skipped: `the refined coordinates broke the ${contactBreak.a}-${contactBreak.b} ` +
                                     `contact floor (${contactBreak.d.toFixed(2)} A against ` +
                                     `${contactBreak.floor.toFixed(2)} A), so the search positions were kept`,
                            rejectedR: r.R
                        };
                    } else if (!r.error) {
                        r.sites.forEach((rs, i) => {
                            best.sites[i].x = rs.x;
                            best.sites[i].y = rs.y;
                            best.sites[i].z = rs.z;
                        });
                        refinement = {
                            // Which normalisation route the observations took,
                            // so a surprising R can be traced to how |Fo|^2 was
                            // obtained rather than guessed at.
                            route: obs.route, lpApplied: job.applyLP !== false,
                            R: r.R, Rstart: r.Rstart,
                            iterations: r.iterations, converged: r.converged,
                            nObs: r.nObs, nParams: r.nParams,
                            // How far each site actually moved, in fractional
                            // units. A shift of about 1/N is the grid bias
                            // being taken out; a much larger one means the map
                            // and the intensities disagree about that site and
                            // is worth looking at rather than accepting.
                            shifts: r.shifts
                        };
                    } else {
                        refinement = { skipped: r.error };
                    }
                } else {
                    refinement = { skipped: 'Sites lack an atomic number or a Wyckoff position.' };
                }
            } else {
                refinement = { skipped: (obs.errors && obs.errors[0]) ||
                                        'No usable Pawley intensities were supplied.' };
            }
        } catch (err) {
            refinement = { skipped: (err && err.message) || String(err) };
        }

        // Minimum image over the 27 neighbouring translations. The delta is
        // wrapped to the nearest cell first, which is enough on its own for a
        // near-orthogonal cell; the shell search covers the oblique ones,
        // where the nearest image is not always the wrapped one.
        const minImage = (dx, dy, dz) => {
            dx -= Math.round(dx); dy -= Math.round(dy); dz -= Math.round(dz);
            let best = Infinity;
            for (let i = -1; i <= 1; i++)
                for (let j = -1; j <= 1; j++)
                    for (let k = -1; k <= 1; k++) {
                        const d = fracDistance(G, dx + i, dy + j, dz + k);
                        if (d < best) best = d;
                    }
            return best;
        };

        // Expand the asymmetric unit to the whole cell once, then measure
        // every unique site against it.
        const expanded = [];
        best.sites.forEach((s, si) => {
            const imgs = [];
            for (const op of symops) {
                const q = applyOp(op, [s.x, s.y, s.z]).map(v => ((v % 1) + 1) % 1);
                if (!imgs.some(im => minImage(im[0] - q[0], im[1] - q[1], im[2] - q[2]) < 1e-3)) {
                    imgs.push(q);
                }
            }
            for (const q of imgs) {
                expanded.push({ element: s.element, siteIdx: si, x: q[0], y: q[1], z: q[2] });
            }
        });

        const formattedSites = best.sites.map((s, i) => {
            // THE DENSITY ACTUALLY UNDER THE SITE.
            //
            // This was hardcoded to 1.0, so every site claimed a full peak
            // whether it sat on one or on vacuum, and the one column that
            // would have shown a site fitted to empty map read the same as a
            // site fitted to the heaviest atom in the cell. The map arrives
            // scaled so its maximum is 1, so this is directly comparable
            // between sites and between runs.
            // No map, no map density. Reporting 0 would read as "sits on
            // nothing" rather than "was never measured".
            const mapDensity = densityMap ? sampleGrid(densityMap, N, s.x, s.y, s.z) : null;

            // Nearest neighbours, as measured - not as assumed. Nothing here
            // asks what coordination the site OUGHT to have; it reports what
            // the fitted structure has, so a missing bond or an impossible
            // contact is visible in the results table instead of having to be
            // computed by hand afterwards.
            const self = expanded.find(a => a.siteIdx === i);
            const near = [];
            for (const b of expanded) {
                if (b === self) continue;
                const d = minImage(self.x - b.x, self.y - b.y, self.z - b.z);
                if (d > 1e-3) near.push({ element: b.element, distance: d });
            }
            near.sort((p, q) => p.distance - q.distance);

            return {
                rank: i + 1,
                element: s.element,
                x: s.x,
                y: s.y,
                z: s.z,
                multiplicity: s.multiplicity,
                wyckoff: s.wyckoff,
                special: s.multiplicity < symops.length,
                // TWO DIFFERENT NUMBERS, BOTH USEFUL, NEITHER A PEAK SEARCH.
                //
                //   mapDensity  the density under this site as a fraction of
                //               the MAP's maximum. An absolute statement: 0.93
                //               means this site sits on something nearly as
                //               strong as the strongest feature in the cell,
                //               and a value at or below zero means it sits on
                //               nothing at all - the map had no opinion about
                //               it and it was placed by composition and
                //               penalties alone.
                //
                //   relative    the same sample against the strongest FITTED
                //               site. A ranking among the sites actually
                //               proposed, which is what the peak-picking path
                //               reports, so the two tables mean the same
                //               thing. Its top entry is 1.000 by construction
                //               and therefore says nothing on its own.
                //
                // `height` is kept as an alias of mapDensity so existing
                // callers keep working; new code should read mapDensity and
                // say which of the two it is showing.
                mapDensity,
                height: mapDensity,
                contacts: near.slice(0, 6),
                // The shortest contact of all, which is the single number that
                // says whether this site is physically possible.
                shortestContact: near.length ? near[0].distance : null
            };
        });

        // `relative` is each site's height against the tallest FITTED site,
        // matching what the peak-picking path reports, so the two routes'
        // tables mean the same thing.
        const topH = formattedSites.reduce(
            (m, s) => Number.isFinite(s.mapDensity) ? Math.max(m, s.mapDensity) : m, 0);
        for (const s of formattedSites) {
            s.relative = (topH > 0 && Number.isFinite(s.mapDensity)) ? s.mapDensity / topH : null;
        }

        return {
            sites: formattedSites,
            assignment: best.assignment,
            candidates: out.candidates,
            refinement,
            // The search's own numbers, kept separate from the refinement's.
            // searchCC is the FULL-RESOLUTION correlation from the quench, not
            // the ramped figure the progress readout shows.
            searchSites,
            searchCC: Number.isFinite(best.cc) ? best.cc : null,
            searchScore: Number.isFinite(best.score) ? best.score : null,
            z: Z
        };
    }
    
    return { sites: [] };
}