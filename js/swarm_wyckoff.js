/* ------------------------------------------------------------------
   swarm_wyckoff.js - the whole Wyckoff-constrained structure search.

   Harko.html calls runWyckoffSearch() once and gets ranked structures back.
   Everything between - normalising the observations, fixing the temperature
   factor, enumerating Wyckoff assignments, running the swarm and scoring the
   result - lives here rather than in the page.

   Requires, already loaded: symmetry_utils.js, wyckoff_assign.js,
   cc_fitness.js, observations.js, scatterers.js.

   THE SHAPE OF THE SEARCH

   Each particle carries an ASSIGNMENT - a choice of Wyckoff position for every
   independent site - so the cell contents are the requested composition by
   construction and cannot drift. Assignments are searched simultaneously in
   one dispatch, as independent sub-swarms: particles on different assignments
   must never attract each other, or the swarm tears itself between
   incompatible structures.

   Fitness runs in two phases. Phase 1 is the Patterson vector sum already in
   the retired swarm_multi.wgsl, which is smooth and gets particles into the right
   neighbourhood cheaply. Phase 2 is the multiplicity-weighted correlation of
   |F|^2, which IS the Patterson map correlation by Parseval and is what
   actually decides between near-solutions. Measured on real PbSO4 data the
   correct structure scores 0.96 against 0.30 for random ones - a 5.9 sigma
   separation over 200 trials.

   Assignments are raced rather than each given a full budget: most are
   visibly hopeless within a few tens of generations.
   ------------------------------------------------------------------ */

const SW_WG = 64;

// Bytes per particle in bufState / mcmcState. Was 8 float/u32 slots (32
// bytes): fit, cc, pen, stepSize, temp, seed, accepted, pad. Now 12 (48
// bytes): the same eight plus bestFit, bestCC, bestPen, bestPad - the
// per-particle best-ever cache that lets the CPU batch several generations
// of dispatches between readbacks without losing a chain's best discovery
// to a later, worse Metropolis accept. Must match ParticleState in
// swarm_reflection.wgsl exactly.
const SW_STATE_STRIDE = 48;

/* ------------------------------------------------------------------
   Search tunables, in one place.

   These were scattered through the file as `?? 0.002`, `|| 512`, `?? 24` and
   bare locals inside the update loop, which made "what does this default to"
   a question you answered by reading the whole file. Every one is still
   overridable per run through the options object; this is only where the
   fallback lives.
   ------------------------------------------------------------------ */
const SW_DEFAULTS = Object.freeze({
    // --- Observations ---
    overlapTol: 0.002,      // fractional d spacing within which powder lines are one observation

    // --- Wyckoff enumeration ---
    maxSitesPerElement: 4,  // raised automatically when the composition demands more
    maxRepeatPerPosition: 4,

    // --- Swarm size ---
    numParticles: 512,
    generations: 400,
    minParticlesPerAssignment: 24,

    // --- Resolution ramp ---
    rampStart: 0.25,        // fraction of the groups active at generation 0
    rampFull: 0.6,          // fraction of the run by which all groups are active

    // --- Penalties, in CC units ---
    penClash: 0.05,         // per clashing pair
    // Per Angstrom of unmet distance constraint. Raised from 0.02, which was
    // set when the only upper-bound rule was a loose "some O within 1.65 A"
    // nearest-neighbour test. As the restoring force of a coordination
    // constraint it was far too weak to argue with the correlation: it held
    // PbSO4's sulfate at 2.00 A against a window ending at 1.90 for an entire
    // run. The ramp keeps it gentle early - 0.2 x this at generation 0 - so
    // exploration is not what pays for it.
    penBond: 0.10,
    penCoord: 0.03,         // per neighbour missing from, or surplus to, a coordination number
    // The penalty ramp. Everything above is multiplied by a factor sweeping
    // from penRampStart to penRampEnd over the run: soft early, so the swarm
    // has a gradient out of the clashes every random start produces; decisive
    // late, so no impossible structure survives on correlation alone.
    penRampStart: 0.2,
    penRampEnd: 4.0,

    // --- Metropolis update ---
    // Replaced particle swarm. PSO's social term pulls every particle of an
    // assignment toward that assignment's best, so the population collapses
    // onto whichever basin was found first and the restarts exist to undo it.
    // On PbSO4 that showed as a run-to-run coin flip between structures with
    // Pb-O at 1.01 A and the true one - the correct answer scores BEST on both
    // CC and R when it is found, so the objective was never the problem and
    // the sampling was. Independent chains cannot collapse: each one keeps its
    // own state and its own step size, and nothing shares a direction.
    // Replica-exchange ladder, in correlation units. Temperatures are FIXED,
    // not annealed: a single cooling chain that anneals into a wrong basin has
    // no way back out, which is exactly what the PbSO4 logs showed - converged,
    // at the bottom of its basin, 0.016 in CC below a structure known to exist.
    // A ladder keeps hot replicas roaming for the whole run and swaps their
    // discoveries downward.
    tempHot: 0.05,          // top rung: moves freely, refines nothing
    tempCold: 5e-4,         // bottom rung: refines, barely explores
    rungs: 8,               // geometric between the two
    swapEvery: 10,          // generations between exchange sweeps
    stepInit: 0.08,         // proposal width, fractional coordinates
    stepMin: 2e-4,
    stepMax: 0.5,
    targetAccept: 0.3,      // Robbins-Monro target for the per-chain step size
    stepAdapt: 0.08,        // how hard the step size chases that target
    // The resolution ramp is held constant over a block of generations so a
    // proposal and the state it is tested against share one objective. More
    // steps means a smoother ramp and more re-measurements of the current
    // state: one extra dispatch each.
    rampSteps: 24,

    // --- Final quench ---
    // A greedy descent from each assignment's best, after the restarts. See
    // quench() for why the search cannot produce this itself.
    // Path length, not breadth. A greedy descent walks DOWN a valley one move
    // at a time, so 23 chains of 150 steps and 1 chain of 3450 are not the same
    // purchase: measured on the real PbSO4 data from a rank-1 structure at
    // R 11.9%, 23x150 reached 10.5%, 23x600 reached 10.0%, and 1x10000 reached
    // 9.9% on a third of the evaluations. Chains only help by trying different
    // directions from the same point; they cannot walk further for you.
    quenchSteps: 600,
    quenchStep0: 0.02,      // starting move size, fractional coordinates
    quenchStep1: 5e-5,      // and where it ends: well under the coordinate precision reported

    // --- Reporting ---
    topN: 20
});

/**
 * Largest symmetry-expanded atom count the device's workgroup memory allows.
 *
 * The CC kernel keeps five arrays per generated atom (x, y, z, type, and the
 * reduction scratch) rather than the six of the vector-sum kernel, so 20 bytes
 * per atom is the honest figure. The 0.85 factor leaves room for the compiler's
 * own workgroup allocations, which are not visible from here.
 */
function swMaxGenAtoms(device, floor = 128) {
    try {
        const BYTES_PER_ATOM = 20, REDUCTION_BYTES = SW_WG * 4 * 6;
        const budget = device?.limits?.maxComputeWorkgroupStorageSize || 16384;
        const cap = Math.floor((budget * 0.85 - REDUCTION_BYTES) / BYTES_PER_ATOM);
        return Math.max(floor, Math.min(4096, cap));
    } catch (e) {
        return floor;
    }
}

function swInject(src, maxGen, minImageShell = 0, groupStride = 3, countShell = 1) {
    let out = src.replace(
        /const\s+MAX_GEN_ATOMS\s*:\s*u32\s*=\s*\d+u\s*;\s*\/\/__MAX_GEN_ATOMS__/,
        `const MAX_GEN_ATOMS: u32 = ${maxGen}u; //__injected__`);
    if (out === src) throw new Error('MAX_GEN_ATOMS marker not found in the kernel source.');

    const before = out;
    out = out.replace(
        /const\s+MIN_IMAGE_SHELL\s*:\s*i32\s*=\s*-?\d+\s*;\s*\/\/__MIN_IMAGE_SHELL__/,
        `const MIN_IMAGE_SHELL: i32 = ${minImageShell}; //__injected__`);
    if (out === before) throw new Error('MIN_IMAGE_SHELL marker not found in the kernel source.');

    // COUNT_SHELL governs how many lattice translations the COORDINATION
    // COUNTING walks, which is a different question from MIN_IMAGE_SHELL's
    // "is the nearest-image distance trustworthy in this basis". Counting needs
    // the neighbouring cells whatever the basis, because two translates of one
    // atom can both be genuine neighbours; 1 matches expandForContacts() in
    // contacts.js, which is the reference these numbers have to agree with.
    const beforeCount = out;
    out = out.replace(
        /const\s+COUNT_SHELL\s*:\s*i32\s*=\s*-?\d+\s*;\s*\/\/__COUNT_SHELL__/,
        `const COUNT_SHELL: i32 = ${countShell}; //__injected__`);
    if (out === beforeCount) throw new Error('COUNT_SHELL marker not found in the kernel source.');

    // GROUP_STRIDE only exists in the reflection kernel, so a missing marker is
    // not an error here - the density kernel legitimately has none.
    out = out.replace(
        /const\s+GROUP_STRIDE\s*:\s*u32\s*=\s*\d+u\s*;\s*\/\/__GROUP_STRIDE__/,
        `const GROUP_STRIDE: u32 = ${groupStride}u; //__injected__`);
    return out;
}

/**
 * Is (-I | 0) one of the space group's operators?
 *
 * If it is, the cell contents the kernel generates are closed under inversion
 * through the origin: every atom at r has a partner at -r of the same element,
 * their sine terms cancel exactly, and F is real. The kernel then skips the
 * imaginary part outright - half the transcendentals in its hottest loop.
 *
 * TESTED, NOT ASSUMED, for two reasons. Being a centrosymmetric group is not
 * enough: the inversion centre has to be AT THE ORIGIN, and a setting that puts
 * it elsewhere carries (-I | t) with t nonzero, for which F is complex. And
 * powder data does not make this true by itself - Friedel's law makes the
 * PATTERN centrosymmetric, |F(h)| = |F(-h)|, which is a statement about
 * measured intensities and says nothing about whether F has an imaginary part.
 * A non-centrosymmetric structure measured on a powder still has one, and it
 * sits inside |F|^2 where it cannot be dropped.
 */
function hasInversionAtOrigin(symOps) {
    const nearInt = v => Math.abs(v - Math.round(v)) < 1e-6;
    return (symOps || []).some(op => {
        const r = op.r, t = op.t || [0, 0, 0];
        if (!r || r.length < 9) return false;
        for (let i = 0; i < 9; i++) {
            const want = (i % 4 === 0) ? -1 : 0;    // -I, row-major
            if (Math.abs(r[i] - want) > 1e-6) return false;
        }
        return nearInt(t[0]) && nearInt(t[1]) && nearInt(t[2]);
    });
}

/**
 * One standard normal deviate, Marsaglia polar.
 *
 * A Gaussian proposal rather than a uniform one because the acceptance test
 * assumes a symmetric kernel and a Gaussian keeps that property under the
 * projection onto a Wyckoff subspace, which a box does not - an axis-aligned
 * box projected onto a slanted subspace is no longer axis-aligned or uniform.
 *
 * WHY THIS FORM. This function is the single most expensive thing in the whole
 * search, which is not obvious and was measured rather than guessed: at 8192
 * chains the host-side half of a generation costs about 9 ms and 8.6 ms of that
 * is here. It is called once per coordinate per chain per step - 123,000 times
 * a generation at 8192 chains and five sites.
 *
 * Box-Muller was costing a log, a sqrt, a cos and two Math.random() calls per
 * deviate, and throwing away the sine half of the pair it computes. The polar
 * form has no trigonometry at all: it rejects about 21% of its samples for
 * landing outside the unit circle, and that is still cheaper than a cosine.
 * Both deviates are kept. Measured at 8192 chains: 8.6 ms -> 6.1 ms.
 *
 * The distribution is exactly normal, not an approximation. That matters: a
 * bounded proposal kernel would never make a long jump, and long jumps are how
 * a hot replica leaves a basin.
 */
let swSpareGaussian = 0, swHasSpare = false;
function swGaussian() {
    if (swHasSpare) { swHasSpare = false; return swSpareGaussian; }
    let u, v, q;
    do {
        u = 2 * Math.random() - 1;
        v = 2 * Math.random() - 1;
        q = u * u + v * v;
    } while (q === 0 || q >= 1);
    const f = Math.sqrt(-2 * Math.log(q) / q);
    swSpareGaussian = v * f;
    swHasSpare = true;
    return u * f;
}

// The two Robbins-Monro multipliers, one for an accepted proposal and one for
// a rejected one. See the adaptation step in the search loop.
const SW_STEP_UP = Math.exp(SW_DEFAULTS.stepAdapt * (1 - SW_DEFAULTS.targetAccept));
const SW_STEP_DOWN = Math.exp(SW_DEFAULTS.stepAdapt * (0 - SW_DEFAULTS.targetAccept));

function swBuffer(device, data, usage) {
    const buf = device.createBuffer({
        size: Math.max(4, Math.ceil(data.byteLength / 4) * 4),
        usage: usage | GPUBufferUsage.COPY_DST
    });
    device.queue.writeBuffer(buf, 0, data);
    return buf;
}

/* ------------------------------------------------------------------ */
/*  Grouping and packing                                               */
/* ------------------------------------------------------------------ */

/**
 * Groups reflections a powder pattern cannot resolve, and packs for the kernel.
 *
 * Pawley's split of intensity between overlapped lines is a fitting artefact,
 * not a measurement, so lines closer in d than the overlap tolerance are summed
 * and compared as one observation. Groups come out sorted by d* ascending,
 * which makes the resolution ramp a single bound on the group index.
 *
 * The observation is the multiplicity-WEIGHTED |Fo|^2. That weighting is not a
 * choice: Parseval sums over the full sphere, where each unique reflection
 * appears m times, so m-weighting is exactly what makes this the map
 * correlation. Weighting by m*Lp instead - comparing raw peak areas - looks
 * better (0.99 against 0.96) but discriminates worse (3.2 sigma against 5.9),
 * because Lp spans a factor of 47 across the pattern and a handful of strong
 * low-angle lines come to dominate both the true and the random structures.
 */
function swPackReflections(rows, options = {}) {
    const tol = options.overlapTol ?? SW_DEFAULTS.overlapTol;
    const sorted = [...rows].sort((a, b) => b.d - a.d);   // d descending = d* ascending

    const groups = [];
    let cur = null;
    for (const r of sorted) {
        if (cur && Math.abs(r.d - cur.d) / cur.d < tol) {
            cur.members.push(r);
        } else {
            cur = { members: [r], d: r.d };
            groups.push(cur);
        }
    }

    const nGroups = groups.length;
    const nRefl = sorted.length;
    const reflPack = new Uint32Array(nRefl * 2);
    // FOUR floats per group, not three: start, count, Iobs, WEIGHT.
    //
    // An intensity of 200 with an ESD of 500 is consistent with anything from
    // -300 to 700, and giving it the same say as one known to a few percent
    // over-trusts it exactly as badly as clamping it to zero would bias it.
    // The fourth slot is 1/sigma^2 on the group total, so a reflection that
    // measured nothing contributes almost nothing.
    const GROUP_STRIDE = 4;
    const groupMeta = new Float32Array(nGroups * GROUP_STRIDE);
    const groupD = new Float32Array(nGroups);
    const problems = [];
    let nUnweighted = 0;

    let w = 0;
    groups.forEach((g, gi) => {
        let io = 0, varSum = 0, haveSigma = true;
        groupMeta[gi * GROUP_STRIDE] = w;
        groupMeta[gi * GROUP_STRIDE + 1] = g.members.length;
        for (const r of g.members) {
            if (Math.abs(r.h) > 511 || Math.abs(r.k) > 511 || Math.abs(r.l) > 511) {
                problems.push(`Reflection ${r.h} ${r.k} ${r.l} exceeds the +/-511 packing range.`);
                continue;
            }
            reflPack[w * 2] = ((r.h + 512) & 0x3FF) | (((r.k + 512) & 0x3FF) << 10)
                            | (((r.l + 512) & 0x3FF) << 20);
            reflPack[w * 2 + 1] = r.mult;
            io += r.mult * r.Fo2;
            // The group is a SUM, so its variance is the sum of the members'.
            // Treating them as independent is an approximation -- reflections
            // this close in d are correlated in the Pawley fit, and the true
            // covariance is what the report's overlap clusters use -- but the
            // off-diagonals are not carried this far, and summing variances is
            // conservative rather than optimistic.
            if (Number.isFinite(r.sigma) && r.sigma > 0) {
                const sm = r.mult * r.sigma;
                varSum += sm * sm;
            } else {
                haveSigma = false;
            }
            w++;
        }
        groupMeta[gi * GROUP_STRIDE + 2] = io;

        // No sigma anywhere in the group means unit weight: an unweighted fit
        // is a defensible default, a fabricated weight is not.
        if (haveSigma && varSum > 0) {
            groupMeta[gi * GROUP_STRIDE + 3] = 1 / varSum;
        } else {
            groupMeta[gi * GROUP_STRIDE + 3] = 0;   // filled in below
            nUnweighted++;
        }
        groupD[gi] = 1 / g.d;
    });

    // Weights are only ever compared with each other, so their absolute scale
    // is free. Normalising to a mean of one keeps them in a range float32
    // represents comfortably: raw 1/sigma^2 on strong data can otherwise run
    // to 1e-8 or smaller and lose most of its precision in the kernel.
    let wSum = 0, wN = 0;
    for (let g = 0; g < nGroups; g++) {
        const v = groupMeta[g * GROUP_STRIDE + 3];
        if (v > 0) { wSum += v; wN++; }
    }
    if (wN > 0) {
        const mean = wSum / wN;
        for (let g = 0; g < nGroups; g++) {
            const v = groupMeta[g * GROUP_STRIDE + 3];
            groupMeta[g * GROUP_STRIDE + 3] = (v > 0) ? v / mean : 1.0;
        }
    } else {
        for (let g = 0; g < nGroups; g++) groupMeta[g * GROUP_STRIDE + 3] = 1.0;
    }
    if (nUnweighted && wN > 0) {
        problems.push(`${nUnweighted} of ${nGroups} groups carry no standard uncertainty ` +
                      `and were given unit weight.`);
    }

    return { reflPack, groupMeta, groupD, nRefl: w, nGroups, overlapTol: tol, groupStride: GROUP_STRIDE,
             weighted: wN > 0,
             overlapped: groups.filter(g => g.members.length > 1).length, problems };
}

/**
 * Scattering factors per (reflection, element type), folded with exp(-B s^2/4).
 *
 * On the MODEL side, never as a rescaling of the observations. Normalising the
 * data shell by shell to match point atoms was an earlier version of this and
 * it is wrong in a way that is easy to miss: it puts the true structure's
 * correlation at 0.88 instead of 1.000, because the data has been distorted and
 * the model has not.
 */
function swScatteringTable(refl, demand, options = {}) {
    const B = options.overallB ?? 0;
    const ff = options.formFactor || (() => null);
    const nElem = demand.length;
    const table = new Float32Array(refl.nRefl * nElem);

    const sOf = new Float32Array(refl.nRefl);
    for (let g = 0; g < refl.nGroups; g++) {
        const gs = refl.groupStride || 3;
        const st = refl.groupMeta[g * gs], c = refl.groupMeta[g * gs + 1];
        for (let m = 0; m < c; m++) sOf[st + m] = refl.groupD[g];
    }

    const missing = new Set();
    for (let r = 0; r < refl.nRefl; r++) {
        const s = sOf[r];
        const dw = B > 0 ? Math.exp(-B * s * s / 4) : 1;
        for (let e = 0; e < nElem; e++) {
            let f = ff(demand[e].element, s);
            if (!Number.isFinite(f)) { f = demand[e].z; missing.add(demand[e].element); }
            table[r * nElem + e] = f * dw;
        }
    }
    return { table, nElem, missing: [...missing] };
}

/**
 * How many reflection groups are active at a given point in the run.
 *
 * The |F|^2 landscape is far more oscillatory than the Patterson vector sum, so
 * a swarm turned loose on the full reflection set tends to stall in a local
 * maximum. Starting with the low-order data gives a smooth surface with few
 * maxima, and adding shells as the run proceeds sharpens it around whatever
 * basin the swarm has already found. Measured on real PbSO4 data at a 0.10
 * heavy-atom error, the full set reads 0.67 while the low-resolution quarter
 * reads 0.86 - the coarse surface still points uphill where the fine one has
 * begun to break up.
 *
 * Groups arrive sorted by d* ascending, so this is a single bound on the index.
 *
 * (cc_fitness.js has an identical copy. That module is the CPU reference
 * implementation used to validate the kernel; the app does not load it, and
 * duplicating thirty characters of arithmetic is better than making the page
 * fetch a second file for one function.)
 */
function rampedGroupCount(nGroups, generation, maxGen, startFrac = 0.25, fullBy = 0.6) {
    const t = Math.min(1, generation / Math.max(1, maxGen * fullBy));
    const frac = startFrac + (1 - startFrac) * t;
    return Math.max(8, Math.min(nGroups, Math.round(nGroups * frac)));
}

/* ------------------------------------------------------------------ */
/*  The search                                                         */
/* ------------------------------------------------------------------ */

/**
 * Runs the whole thing.
 *
 * @param {Object} o
 *   device            GPUDevice (already requested by the caller)
 *   ccShaderSource    text of swarm_reflection.wgsl
 *   setting           one entry from sg/<n>.json settings[]
 *   cell              {a,b,c,alpha,beta,gamma}
 *   reflections       parsed rows, whatever columns the file had
 *   wavelength        Angstrom, from the UI
 *   formula, Z        e.g. 'PbSO4', 4
 *   atomData          {PB:{z,r}, S:{...}, O:{...}}
 *   windows           [{a,b,dmin,dmax,bothWays}]
 *   harkerSites       consolidated candidate positions, optional
 *   scatterTables     from loadScatteringTables(), optional
 *   radiation         'xray' | 'neutron'
 *   generations, numParticles, topN
 *   onProgress(info)  called each generation
 *   shouldStop()      polled each generation
 */
async function runWyckoffSearch(o) {
    const log = [];
    const say = m => { log.push(m); if (o.onLog) o.onLog(m); };

    // Auto-Tuner: Probes hardware limits and prevents WebGPU/TDR crashes
    const maxExpectedAtoms = 150; // Safe upper ceiling for generated XYZ atoms per structure
    const bytesPerParticle = maxExpectedAtoms * 24; // 3 floats * 4 bytes * 2 arrays (current + best)
    const limitMem = Math.floor(o.device.limits.maxStorageBufferBindingSize / bytesPerParticle);
    const limitCompute = o.device.limits.maxComputeWorkgroupsPerDimension * 64;
    
    const safeMax = Math.floor(Math.min(limitMem, limitCompute) * 0.95); // 5% safety margin
    
    // o.numParticles, NOT o.swarmParticles.
    //
    // The caller passes `numParticles`, and every consumer downstream reads
    // `numParticles` -- so this block was testing an undefined property
    // (`undefined > safeMax` is false, so nothing was ever capped) and then
    // assigning the result to a property nobody reads. It also left
    // o.swarmParticles as NaN, since Math.floor(undefined / 64) is NaN. The
    // auto-tuner was a no-op in both directions: it could not cap a request
    // that would crash the device, and it could not snap one to the workgroup
    // size either.
    let want = Number.isFinite(o.numParticles) ? o.numParticles : SW_DEFAULTS.numParticles;
    if (want > safeMax) {
        say(`Hardware limit detected. Auto-capped chains from ${want} to ${safeMax} ` +
            `to stay inside this device's storage-buffer and dispatch limits.`);
        want = safeMax;
    }
    // Snap to a multiple of the workgroup size to prevent trailing-thread errors.
    o.numParticles = Math.max(SW_WG, Math.floor(want / SW_WG) * SW_WG);

    /* ---- 1. Composition ---- */
    const comp = parseFormula(o.formula, o.Z);


    const demand = comp.map(c => {
        const ad = (o.atomData || {})[c.element];
        if (!ad) throw new Error(`No atomic data for element "${c.element}".`);
        if (!Number.isInteger(c.count)) {
            throw new Error(`${c.element} needs ${c.count} atoms per cell, which is not a whole ` +
                            `number. Partial occupancy is not supported by this search.`);
        }
        // `r` is the covalent radius (van der Waals only where none is listed).
        // It no longer sets the clash floor - that is the minimum-contact
        // slider, applied to every pair alike, see buildRestraintTables - and
        // is carried here only as per-element metadata for anything that wants
        // a size.
        return { element: c.element, count: c.count, z: ad.z, r: ad.rc ?? ad.r };
    });
    const nTot = demand.reduce((n, d) => n + d.count, 0);
    say(`${o.formula} with Z=${o.Z}: ${demand.map(d => d.element + d.count).join(' ')} ` +
        `= ${nTot} atoms per cell.`);

   /* ---- 2. Observations ---- */
    let obs = { rows: [], errors: [], warnings: [], route: 'CF_Density', dMin: 0, dMax: 0 };
    let overallB = 0;
    let ff = null;

    {
        const rotations = (o.setting.sym_ops || []).map(op => op.r);
        obs = normaliseObservations(o.reflections, {
            rotations, wavelength: o.wavelength, radiation: o.radiation,
            polarisationK: o.polarisationK
        });
        if (obs.errors.length) throw new Error(obs.errors.join(' '));
        obs.warnings.forEach(say);
        say(`${obs.rows.length} reflections, d ${obs.dMin.toFixed(2)}-${obs.dMax.toFixed(2)} A, ` +
            `route "${obs.route}".`);

        /* ---- 3. Temperature factor ---- */
        // Structure independent, so it is fixed once here rather than refined
        // during the search. On real PbSO4 data B is worth 0.31 in correlation
        // (0.65 at B=0 against 0.97 at the optimum), so leaving it at zero is not
        // an option; but the peak is broad, so a good estimate is enough.
        ff = o.scatterTables
            ? makeFormFactor(o.scatterTables, { radiation: o.radiation || 'xray', ions: o.ions })
            : null;
        const wil = wilsonB(obs.rows, demand, ff);
        if (wil.note) say(wil.note);
        
        if (Number.isFinite(o.overallB)) {
            overallB = o.overallB;
            say(`Overall B = ${overallB.toFixed(2)} (User defined). Wilson plot estimates ${wil.B.toFixed(2)}.`);
        } else {
            overallB = wil.B;
            say(`Overall B = ${overallB.toFixed(2)}${wil.ok ? ` (Wilson, ${wil.shells} shells)` : ''}.`);
        }
    }

    /* ---- 4. Assignments ---- */
    if (!o.setting.wyckoff || !o.setting.wyckoff.length) {
        throw new Error(`Setting ${o.setting.symbol} carries no Wyckoff table. ` +
                        `Regenerate the database with cctbx_Harko_v1.py.`);
    }
    // What each element's coordination has to look like, read off the same
    // distance windows the kernel penalises against, so the combinatorial
    // filter and the energy cannot disagree about what the user asked for.
    // Only windows that state a COUNT say anything about the shape of a
    // coordination shell; a bare distance range does not.
    const coordination = {};
    for (const w of (o.windows || [])) {
        if (!w || !Number.isFinite(w.count) || w.count <= 0) continue;
        const key = String(w.a).toUpperCase();
        // Keep the most demanding statement if several name the same centre.
        if (!coordination[key] || w.count > coordination[key].count) {
            coordination[key] = { count: w.count, geometry: w.geometry || null };
        }
    }

    // The coordination filter is CORRECT and stays on. Checked against PbSO4 in
    // Pnma: without it the enumeration returns 35 assignments, with it 19, and
    // the true one (Pb 4c / S 4c / O 8d + O 4c + O 4c) is present either way.
    // All it removes is sulfur on 4a and 4b, both of which have an inversion
    // centre that no tetrahedron can sit on -- exactly what it claims to do.
    const en = enumerateAssignments(o.setting.wyckoff, demand, {
        maxSites: o.maxSitesPerElement ?? SW_DEFAULTS.maxSitesPerElement,
        maxRepeat: o.maxRepeat ?? SW_DEFAULTS.maxRepeatPerPosition,
        ceiling: o.wyckoffCapCeiling,
        coordination,
        symOps: o.setting.sym_ops
    });
    if (en.error) throw new Error(en.error);
    if (en.siteSymmetryExcluded && en.siteSymmetryExcluded.length) {
        // Report once per distinct reason rather than once per position: the
        // same rule usually rejects several positions and the repetition
        // buries anything else in the log.
        for (const m of [...new Set(en.siteSymmetryExcluded)]) say(`Site symmetry: ${m}`);
    }
    let assignments = en.assignments;
    if (!assignments.length) throw new Error('No Wyckoff assignment matches this composition.');
    // The enumeration drops assignments that need more independent sites than
    // the kernel can hold, and stops at maxTotal. Both are legitimate, and both
    // make the candidate list shorter than the combinatorics -- so both are
    // said out loud rather than leaving the count unexplained.
    if (en.droppedTooManySites) {
        say(`${en.droppedTooManySites} assignment(s) skipped: they need more than ` +
            `${WY_MAX_SITES} independent sites, which is the most the search kernel holds.`);
    }
    if (en.truncated) {
        say(`The enumeration hit its ${en.assignments.length}-assignment ceiling; ` +
            `plausible assignments beyond it were not searched.`);
    }

    // Reduced basis, so every distance in this run - the pre-search ranking
    // here, and the kernel's contact tests below - uses the same, correct,
    // nearest-image convention.
    const red = sharkoReducedCell(o.cell);
    rankAssignments(assignments, {
        orth: red, harkerSites: o.harkerSites || [],
        heavyZ: Math.max(...demand.map(d => d.z)) * 0.7
    });
    say(`${assignments.length} Wyckoff assignment(s) consistent with the composition.`);


    /* ---- 5. Pack ---- */
    const T = buildAssignmentTables(assignments, (o.setting.sym_ops || []).length, demand);
    T.warnings.forEach(say);
    
    let refl = { reflPack: new Uint32Array(0), groupMeta: new Float32Array(0), groupD: new Float32Array(0), nRefl: 0, nGroups: 0, groupStride: 3 };
    let ftab = { table: new Float32Array(0), missing: [] };
    
    {
        refl = swPackReflections(obs.rows, { overlapTol: o.overlapTol });
        refl.problems.forEach(say);
        say(`${refl.nGroups} reflection group(s), ${refl.overlapped} containing overlaps.`);

        ftab = swScatteringTable(refl, demand, { overallB, formFactor: ff });
        if (ftab.missing.length) {
            say(`No tabulated scattering factor for ${ftab.missing.join(', ')}; using f = Z. ` +
                `A missing element is indistinguishable from a wrong structure once it reaches ` +
                `the correlation, so check scatters/ is present.`);
        }
    }
    // o.minContact is the Minimum contact distance slider. It was passed in
    // here and never read - it reached a params slot the kernel's struct does
    // not have, so it landed in padding. The slider moved and nothing changed,
    // while the post-search filter DID use it, which is how a run could spend
    // itself on contacts that were then all rejected.
    const restraints = buildRestraintTables(demand, o.windows || [],
                                            { minContact: o.minContact });
    restraints.problems.forEach(say);
    say(`Contact floor ${Number.isFinite(o.minContact) ? o.minContact.toFixed(2) : '1.00'} A on ` +
        `every pair${(o.windows || []).some(w => Number.isFinite(w.dmin))
            ? ', except pairs given their own dmin' : ''}.`);

    /* ---- 6. Device ---- */
    const device = o.device;
    if (!device) throw new Error('No GPUDevice supplied.');
    const maxGen = swMaxGenAtoms(device);
    if (nTot > maxGen) {
        throw new Error(`${nTot} atoms per cell exceeds this GPU's workgroup budget of ${maxGen}. ` +
                        `Reduce Z, or use a device with more workgroup storage.`);
    }

    // Minimum image. The kernel works in a reduced basis for the same lattice,
    // where rounding each fractional component independently is provably the
    // nearest image out to red.safeRadius. That covers every threshold for
    // almost any cell; for a genuinely thin one - a layered structure with a
    // short in-plane repeat - it does not, and the kernel is compiled to search
    // the 27 neighbouring translations instead, which is exact at any distance
    // and costs 27x in the contact loop only when it is actually needed.
    const needExact = restraints.maxDistance > red.safeRadius;
    const minImageShell = needExact ? 1 : 0;
    if (red.changed) {
        say(`Contact distances use a reduced lattice basis (safe to ` +
            `${red.safeRadius.toFixed(2)} A); in the cell's own basis the ` +
            `minimum-image convention is not reliable for a skewed cell.`);
    }
    if (needExact) {
        say(`This cell is thin: the reduced basis is only exact to ` +
            `${red.safeRadius.toFixed(2)} A and the restraints reach ` +
            `${restraints.maxDistance.toFixed(2)} A. The kernel will search the ` +
            `neighbouring lattice translations as well - correct, but slower.`);
    }

    const module = device.createShaderModule({
        code: swInject(o.ccShaderSource, maxGen, minImageShell,
                       o.groupStride || refl.groupStride || 3)
    });
    const pipeline = device.createComputePipeline({ layout: 'auto', compute: { module, entryPoint: 'main' } });

    // Particle budget. Every assignment needs a floor of particles to behave as
    // a sub-swarm at all; that floor sets how many can run at once, and the
    // rest are queued into later waves in prior order.
    const numParticles = Math.max(SW_WG, o.numParticles || SW_DEFAULTS.numParticles);
    const minPer = o.minParticlesPerAssignment ?? SW_DEFAULTS.minParticlesPerAssignment;
    // Particles are shared out in proportion to the DIMENSIONALITY of each
    // assignment, not equally.
    //
    // An equal split systematically favours low-parameter models. A site on 4a
    // has no free parameter at all and a handful of particles finds its
    // optimum immediately; an assignment with eleven free parameters searching
    // the same budget is still wandering when the run ends, and reports a
    // correlation well below what it can actually reach. The ranking then
    // measures how hard each assignment was to search rather than how well it
    // fits - which is precisely backwards, since the extra parameters exist to
    // fit the data better.
    //
    // Weighting by parameter count is the cheapest correction that removes the
    // bias's direction. It does not remove it entirely - search difficulty
    // grows faster than linearly with dimension - so a close ranking between
    // assignments of very different size still deserves a longer run.
    const dims = assignments.map(A => Math.max(1, A.sites.reduce((n, s) => n + wyckoffFreedom(s.w), 0)));
    const waves = allocateParticles(assignments.length, numParticles, minPer, dims);
    const dimRange = `${Math.min(...dims)}-${Math.max(...dims)}`;
    say(`${numParticles} particles over ${waves.length} wave(s); assignments carry ` +
        `${dimRange} free parameters and receive particles in proportion.`);

/* ---- 7. Static buffers ---- */
    const symPacked = packSymOps(o.setting.sym_ops);
    // groupData and the scattering table share one binding: WebGPU guarantees
    // only eight storage buffers per stage and this kernel uses all eight.
    const groupData = new Float32Array(refl.groupMeta.length + ftab.table.length);
    groupData.set(refl.groupMeta, 0);
    groupData.set(ftab.table, refl.groupMeta.length);
    const fTabOff = refl.groupMeta.length;





// ------------------------------------------------------------------
    //  EVERY GPU BUFFER BELOW IS RELEASED WHETHER THIS FUNCTION RETURNS OR
    //  THROWS.
    //
    //  The ten destroy() calls used to sit on the normal exit path only, and
    //  there are two explicit `throw`s and ten `await`s between here and there
    //  -- any of which (a lost device, a failed mapAsync, a validation error)
    //  skipped the cleanup entirely.
    //
    //  That was invisible while a run did ONE search: the caller destroyed the
    //  whole device afterwards and the buffers went with it. It stopped being
    //  invisible when the caller began scanning Z, because the device is now
    //  created once and reused across every candidate, and a candidate that
    //  fails is CAUGHT AND SWALLOWED so the scan can continue. Each failed Z
    //  therefore stranded ten buffers in a context that stayed alive for the
    //  rest of the scan.
    //
    //  Declared with `let` before the try so the finally can see them; each is
    //  guarded because the throw may happen before the later ones exist.
    // ------------------------------------------------------------------
    const S = GPUBufferUsage.STORAGE;
    let bufGen, bufRefl, bufGroup, bufStatic, bufPos, bufAssign,
        bufState, bufPosRead, bufStateRead, bufParams;
    try {
    bufGen   = swBuffer(device, T.genPack, S);
    bufRefl  = swBuffer(device, refl.reflPack, S);
    
    bufGroup = swBuffer(device, groupData, S);
    
    const symLen = symPacked.length;
    const tabLen = restraints.tables.length;
    const projLen = T.siteProj.length;
    const staticFloats = new Float32Array(symLen + tabLen + projLen);
    staticFloats.set(symPacked, 0);
    staticFloats.set(restraints.tables, symLen);
    staticFloats.set(T.siteProj, symLen + tabLen);
    bufStatic = swBuffer(device, staticFloats, S);

    const coordsPerParticle = T.maxSites * 3;
    const totalCoordFloats = numParticles * coordsPerParticle;
    const positions = new Float32Array(totalCoordFloats);
    const proposals = new Float32Array(totalCoordFloats);
    // GPU-side best-ever coordinates for each chain, mirrored to the CPU
    // only at a Replica Exchange sync. See bufPos below: the GPU buffer is
    // twice this length, current half then best half.
    const bestPositions = new Float32Array(totalCoordFloats);
    const stepSize = new Float32Array(numParticles);
    const curFit = new Float32Array(numParticles);
    const curCC  = new Float32Array(numParticles).fill(NaN);
    const curPen = new Float32Array(numParticles);
    // CPU mirror of the GPU per-particle best cache, valid only right after
    // readStateFromGPU(). bestFit <= -1e29 means "nothing recorded yet",
    // matching the sentinel written in syncStateToGPU() and read in WGSL.
    const bestFit = new Float32Array(numParticles).fill(-1e30);
    const bestCC  = new Float32Array(numParticles);
    const bestPen = new Float32Array(numParticles);
    const tempOf = new Float32Array(numParticles);
    let acceptRate = NaN, swapRate = NaN;
    // MCMC dispatches since the last syncStateToGPU(). The kernel's per-chain
    // `accepted` counter is zeroed by that sync, so the denominator of the
    // acceptance rate is the number of proposals made in the same window --
    // which is exactly this. Without it acceptRate was declared, reported to
    // the UI on every progress message, and never assigned: NaN for the whole
    // run, so "are the chains moving at all" had no answer.
    let mcmcSinceSync = 0;

    // bufPos holds CURRENT coordinates in its first half and the GPU's
    // BEST-EVER coordinates per chain in its second half (see
    // swarm_reflection.wgsl). Doubling it is what lets the search go many
    // generations between CPU readbacks without silently losing a structure
    // a chain later moved away from.
    const bufPosBytes = totalCoordFloats * 2 * 4;
    if (bufPosBytes > device.limits.maxStorageBufferBindingSize) {
        throw new Error(`The current + best-ever position buffer needs ` +
            `${bufPosBytes} bytes per swarm (numParticles=${numParticles}, ` +
            `${coordsPerParticle} coordinate floats each, doubled for the ` +
            `best-ever cache), which exceeds this device's ` +
            `maxStorageBufferBindingSize of ${device.limits.maxStorageBufferBindingSize}. ` +
            `Reduce numParticles or the number of sites per assignment.`);
    }
    bufPos = device.createBuffer({
        size: bufPosBytes,
        usage: S | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST
    });
    bufAssign = swBuffer(device, new Uint32Array(numParticles), S);
    // Persistent readback buffers, reused for every synchronous readback
    // instead of being created and destroyed inside the generation loop -
    // creating a MAP_READ buffer per sync was another source of driver
    // overhead this design removes.
    bufPosRead = device.createBuffer({
        size: bufPosBytes,
        usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
    });
    bufStateRead = device.createBuffer({
        size: numParticles * SW_STATE_STRIDE,
        usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
    });

    const PARAM = Object.freeze({
        o0: 0, o1: 4, o2: 8,          
        r0: 12, r1: 16, r2: 20,       
        nTot: 24, maxSites: 25, numParticles: 26, nGroupsActive: 27,
        nElem: 28, nBondRules: 29, rMinOff: 30, ruleOff: 31,
        fTabOff: 32, nRefl: 33, penClash: 34, penBond: 35,
        penCoord: 36, penScale: 37, centro: 38, pad0: 39,
        symOff: 40, tablesOff: 41, projOff: 42, pad1: 43
    });
    const params = new Float32Array(44);
    
    bufParams = device.createBuffer({
        size: params.byteLength,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    });
    for (let i = 0; i < 3; i++) {
        const o4 = PARAM.o0 + i * 4, r4 = PARAM.r0 + i * 4;
        for (let k = 0; k < 3; k++) {
            params[o4 + k] = red.orth[i * 3 + k];
            params[r4 + k] = red.xform[i * 3 + k];
        }
    }
    params[PARAM.nTot] = nTot;
    params[PARAM.maxSites] = T.maxSites;
    params[PARAM.numParticles] = numParticles;
    params[PARAM.nElem] = restraints.nElem;
    params[PARAM.nBondRules] = restraints.nRules;
    params[PARAM.rMinOff] = restraints.rMinOff;
    params[PARAM.ruleOff] = restraints.ruleOff;
    params[PARAM.fTabOff] = fTabOff;
    params[PARAM.nRefl] = refl.nRefl;
    params[PARAM.penClash] = o.penClash ?? SW_DEFAULTS.penClash;

    // THE BOND PENALTY IS PER CHARGED NEIGHBOUR SLOT, NOT A RAW SUM.
    //
    // penBond and penCoord are documented as being "in CC units", and they were
    // not: the kernel accumulates them over every atom carrying a rule and every
    // neighbour slot that rule tracks, and hands the total to a score that is
    // CC minus the penalty, with CC confined to [0, 1]. On PbSO4 with the one
    // reasonable rule "S O 4 1.35/1.65" that is four sulfurs times four slots,
    // and at a random start with the nearest oxygens near 3 A the raw penalty
    // is 13.8 -- multiplied by the ramp, between 2.8 and 55. Against a
    // correlation that can only move by 1.0 in total, the diffraction data is
    // switched off: the sampler optimises geometry alone until a tetrahedron
    // exists, and only then does CC begin to matter, by which point the chain
    // is in whatever basin the geometry led it to.
    //
    // The symptom is a search that gets WORSE when given correct chemistry. The
    // run that produced this fix returned CC 0.699 with the constraint and CC
    // 0.993 without, and the structure it chose satisfied the constraint less
    // well than the one it rejected -- zero oxygens in the window against four,
    // penalty 0.87 against 0.00. It lost on both terms, which can only happen
    // if the objective it was actually descending was neither of them.
    //
    // Dividing by the number of charged slots restores the documented meaning:
    // one unit of penBond per Angstrom of error PER BOND. The same PbSO4 start
    // then costs 0.17 early and 3.46 late, so the data leads while the swarm is
    // choosing a basin and the geometry closes it out -- which is what the
    // ramp was written to do. A structure that satisfies its constraints is
    // charged nothing either way, so the answer's score is unchanged.
    //
    // Applied on the host, to penBond and penCoord only. Those two params
    // appear nowhere in the kernel except the rule branches, so the clash term
    // -- which is tuned, unnormalised, and works -- is untouched.
    const nOf = (el) => {
        const d = demand.find(x => String(x.element).toUpperCase() === String(el).toUpperCase());
        return d ? d.count : 0;
    };
    let bondSlots = 0, coordCentres = 0;
    for (const w of (o.windows || [])) {
        const nA = nOf(w.a);
        if (!nA) continue;
        const hasCount = Number.isFinite(w.count) && w.count > 0;
        const hasMax = Number.isFinite(w.dmax);
        if (hasCount) {
            bondSlots += nA * w.count;      // one slot per demanded neighbour
            coordCentres += nA;             // the surplus term charges once per centre
        } else if (hasMax) {
            bondSlots += nA;                // nearest-neighbour rule: one charge per centre
            // buildRestraintTables mirrors a bare upper bound, so B is charged too.
            if (w.bothWays !== false && String(w.a).toUpperCase() !== String(w.b).toUpperCase()) {
                bondSlots += nOf(w.b);
            }
        }
    }
    const bondNorm = Math.max(1, bondSlots);
    const coordNorm = Math.max(1, coordCentres);
    params[PARAM.penBond] = (o.penBond ?? SW_DEFAULTS.penBond) / bondNorm;
    params[PARAM.penCoord] = (o.penCoord ?? SW_DEFAULTS.penCoord) / coordNorm;
    if (bondSlots > 0) {
        say(`Distance restraints charge ${bondSlots} neighbour slot(s) across the cell; ` +
            `penBond is divided by that, so the penalty is per bond and stays commensurate ` +
            `with the correlation instead of swamping it.`);
    }
    params[PARAM.penScale] = SW_DEFAULTS.penRampStart;
    params[PARAM.nGroupsActive] = refl.nGroups;
    
    params[PARAM.centro] = hasInversionAtOrigin(o.setting.sym_ops) ? 1 : 0;
    if (params[PARAM.centro]) {
        say('Centrosymmetric group: F is real, so the kernel skips the imaginary part.');
    }

    params[PARAM.symOff] = 0;
    params[PARAM.tablesOff] = symLen;
    params[PARAM.projOff] = symLen + tabLen;

    bufState = device.createBuffer({
        size: numParticles * SW_STATE_STRIDE,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST
    });

    const bindEntries = [
        { binding: 0, resource: { buffer: bufPos } },
        { binding: 1, resource: { buffer: bufAssign } },
        { binding: 2, resource: { buffer: bufGen } },
        { binding: 3, resource: { buffer: bufStatic } },
        { binding: 5, resource: { buffer: bufGroup } },
        { binding: 6, resource: { buffer: bufState } },
        { binding: 7, resource: { buffer: bufParams } }
    ];
    
    bindEntries.push({ binding: 4, resource: { buffer: bufRefl } });

    const bind = device.createBindGroup({
        layout: pipeline.getBindGroupLayout(0),
        entries: bindEntries
    });
    // ...

        /**
         * Enqueues one dispatch - proposal, fitness, and (if isMCMC) the
         * accept/reject test - and returns immediately. No readback, no
         * await: this is what lets the CPU queue up many generations
         * without forcing the GPU to stop and report back on every one.
         * The WebGPU queue is ordered, so a run of these plus the
         * params/bufPos writes ahead of each stays correctly sequenced.
         */
        function dispatchCoords(coords, isMCMC = false) {
            if (isMCMC) mcmcSinceSync++;
            params[PARAM.pad0] = isMCMC ? 1.0 : 0.0;
            device.queue.writeBuffer(bufParams, 0, params);
            if (coords) device.queue.writeBuffer(bufPos, 0, coords);

            const enc = device.createCommandEncoder();
            const pass = enc.beginComputePass();
            pass.setPipeline(pipeline);
            pass.setBindGroup(0, bind);
            pass.dispatchWorkgroups(numParticles);
            pass.end();
            device.queue.submit([enc.finish()]);
        }



        /**
         * Safely awaits mapAsync against TDR hangs and device loss,
         * without leaking floating promise rejections.
         */
        async function safeMapAsync(buffers) {
            const WATCHDOG_MS = 30000;
            let timer = null;
            
            const guard = new Promise((_, reject) => {
                timer = setTimeout(() => reject(new Error(
                    `The GPU stopped responding (no result for ${WATCHDOG_MS / 1000} s). ` +
                    `This usually means the driver reset the device.`
                )), WATCHDOG_MS);
            });

            // Catch floating rejections in case the race resolves via guard or device.lost
            const mapPromise = Promise.all(buffers.map(b => b.mapAsync(GPUMapMode.READ)));
            mapPromise.catch(() => {}); 

            // Convert device.lost into a throwable rejection, caught safely
            const lostPromise = device.lost.then(info => {
                throw new Error('WebGPU device lost: ' + (info?.message || info?.reason || 'unknown'));
            });
            lostPromise.catch(() => {});

            try {
                await Promise.race([mapPromise, lostPromise, guard]);
            } finally {
                if (timer) clearTimeout(timer);
            }
        }




        /**
         * Dispatch, then read back only fit/cc, immediately. This is the
         * old evaluateCoords()'s exact behaviour and return shape, kept for
         * the two callers that need a synchronous answer for one dispatch:
         * quench(), which is a CPU-driven greedy descent and has to see
         * each step's result before proposing the next, and the rare
         * resolution-ramp / reassignment re-score below, which needs a real
         * number to run checkSpread() and recordBest() against right away.
         * Neither is the per-generation hot path this file exists to fix.
         */
        async function evaluateCoordsSync(coords, isMCMC = false) {
            dispatchCoords(coords, isMCMC);

            const enc = device.createCommandEncoder();
            enc.copyBufferToBuffer(bufState, 0, bufStateRead, 0, bufState.size);
            device.queue.submit([enc.finish()]);

            await safeMapAsync([bufStateRead]);
            const stateF32 = new Float32Array(bufStateRead.getMappedRange());
            const S32 = SW_STATE_STRIDE / 4;
            const out = new Float32Array(numParticles * 2);
            for (let i = 0; i < numParticles; i++) {
                out[i] = stateF32[i * S32 + 0];               // fit
                out[numParticles + i] = stateF32[i * S32 + 1]; // cc
            }
            bufStateRead.unmap();
            return out;
        }

        /**
         * The ONE synchronization point for the normal generation loop.
         * Pulls back current AND best-ever positions plus the full particle
         * state - everything CPU-side code (Replica Exchange, best-tracking)
         * needs - in one copy + one pair of mapAsync calls, using the
         * persistent readback buffers rather than allocating new ones.
         */
        async function readStateFromGPU() {
            const enc = device.createCommandEncoder();
            enc.copyBufferToBuffer(bufPos, 0, bufPosRead, 0, bufPos.size);
            enc.copyBufferToBuffer(bufState, 0, bufStateRead, 0, bufState.size);
            device.queue.submit([enc.finish()]);

            await safeMapAsync([bufPosRead, bufStateRead]);

            const posF32 = new Float32Array(bufPosRead.getMappedRange());
            positions.set(posF32.subarray(0, totalCoordFloats));
            bestPositions.set(posF32.subarray(totalCoordFloats, totalCoordFloats * 2));

            const range = bufStateRead.getMappedRange();
            const stateF32 = new Float32Array(range);
            const stateU32r = new Uint32Array(range);
            const S32 = SW_STATE_STRIDE / 4;
            let accepted = 0;
            for (let i = 0; i < numParticles; i++) {
                const b = i * S32;
                curFit[i]   = stateF32[b + 0];
                curCC[i]    = stateF32[b + 1];
                curPen[i]   = stateF32[b + 2];
                stepSize[i] = stateF32[b + 3];
                accepted   += stateU32r[b + 6];
                bestFit[i]  = stateF32[b + 8];
                bestCC[i]   = stateF32[b + 9];
                bestPen[i]  = stateF32[b + 10];
            }
            if (mcmcSinceSync > 0) {
                acceptRate = accepted / (numParticles * mcmcSinceSync);
            }

            bufPosRead.unmap();
            bufStateRead.unmap();
        }

    /* ---- 8. Run ---- */
    // Heavy-atom seeding lived here: it started a fraction of the particles
    // with the heaviest site on a consolidated Harker peak. Removed with the
    // checkbox that drove it. It existed to collapse the hardest part of a
    // GENERAL-POSITION search - finding the heavy atom's three free
    // coordinates - and the Wyckoff route does not have that problem, because
    // the assignment is what fixes the position. `harkerSites` is still
    // accepted and still used, by rankAssignments, to order the assignments
    // before the search; that use is unaffected.


    const generations = o.generations || SW_DEFAULTS.generations;
    const restarts = Math.max(1, o.restarts || 1);
    const results = [];
    let stopped = false;

    // Best over the WHOLE run, not the current wave. Reporting the wave's own
    // best made the fitness trace saw-tooth downwards: wave 1 holds the
    // most plausible assignments and scores highest, wave 4 holds the least
    // and starts from nothing, so the chart appeared to show the search
    // getting worse when it was simply starting a new problem.
    let runBestFit = -Infinity, runBestCC = NaN, runBestAssign = -1;

    for (let wi = 0; wi < waves.length && !stopped; wi++) {
        let ids = waves[wi].ids.slice();
        const assignOf = waves[wi].assignOf;

        // Per-assignment global bests. Particles on different assignments are
        // solving different problems, so a single shared best would pull each
        // sub-swarm toward a structure it cannot express.
        const gBestFit = new Float32Array(assignments.length).fill(-Infinity);
        const gBestCC  = new Float32Array(assignments.length).fill(NaN);
        const gBestPen = new Float32Array(assignments.length);   // raw, at weight 1
        const gBestPos = new Float32Array(assignments.length * coordsPerParticle);

        for (let i = 0; i < numParticles; i++) {
            const base = i * coordsPerParticle;
            for (let k = 0; k < coordsPerParticle; k++) positions[base + k] = Math.random();
            projectSites(positions, i, assignOf[i], T);
            curFit[i] = -Infinity; curCC[i] = NaN; curPen[i] = 0;
            stepSize[i] = SW_DEFAULTS.stepInit;
        }
        device.queue.writeBuffer(bufAssign, 0, assignOf);
        let ladders = [];

        /**
         * Builds the replica ladders and hands every chain its temperature.
         *
         * A ladder is `rungs` chains OF THE SAME ASSIGNMENT holding the
         * geometric temperature sequence from tempCold to tempHot. Chains of
         * different assignments are solving different problems in different
         * subspaces and their coordinates are not interchangeable, so an
         * exchange between them would be meaningless - the grouping by
         * assignment is what makes a swap legal.
         *
         * An assignment with a partial ladder left over puts those chains on the
         * COLD end. A partial ladder cannot exchange usefully, so the useful
         * thing for it to do is refine.
         */
        function buildLadders() {
            const byAssign = new Map();
            for (let i = 0; i < numParticles; i++) {
                const A = assignOf[i];
                if (!byAssign.has(A)) byAssign.set(A, []);
                byAssign.get(A).push(i);
            }
            const R = Math.max(2, SW_DEFAULTS.rungs);
            const ratio = Math.pow(SW_DEFAULTS.tempHot / SW_DEFAULTS.tempCold, 1 / (R - 1));
            const rungTemp = Array.from({ length: R },
                (_, r) => SW_DEFAULTS.tempCold * Math.pow(ratio, r));

            const out = [];
            for (const chains of byAssign.values()) {
                let c = 0;
                while (c + R <= chains.length) {
                    const lad = chains.slice(c, c + R);
                    lad.forEach((ci, r) => { tempOf[ci] = rungTemp[r]; });
                    out.push(lad);
                    c += R;
                }
                for (; c < chains.length; c++) tempOf[chains[c]] = rungTemp[0];
            }
            return out;
        }

        /**
         * One exchange sweep over every ladder.
         *
         * Adjacent rungs only, alternating which pairs are tried, which is the
         * standard way to keep the sweep a valid sequence of independent
         * two-body moves rather than a scramble.
         *
         *     P(swap) = min(1, exp[ (beta_cold - beta_hot) (f_hot - f_cold) ])
         *
         * So a hot replica that has found something better than the cold one
         * below it hands the structure down, and the cold replica's worse
         * structure goes up to be knocked about further. Nothing is lost and
         * nothing is duplicated: the two configurations trade places, the
         * temperatures stay with the rungs.
         */
        function exchangeSweep(parity, scale) {
            let tried = 0, taken = 0;
            const energy = i => Number.isFinite(curCC[i]) ? curCC[i] - curPen[i] * scale : curFit[i];
            for (const lad of ladders) {
                for (let r = parity; r + 1 < lad.length; r += 2) {
                    const a = lad[r], b = lad[r + 1];          // a colder than b
                    const fa = energy(a), fb = energy(b);
                    if (!Number.isFinite(fa) || !Number.isFinite(fb)) continue;
                    tried++;
                    const d = (1 / tempOf[a] - 1 / tempOf[b]) * (fb - fa);
                    if (d < 0 && Math.random() >= Math.exp(d)) continue;
                    taken++;
                    const ba = a * coordsPerParticle, bb = b * coordsPerParticle;
                    for (let k = 0; k < coordsPerParticle; k++) {
                        const t = positions[ba + k];
                        positions[ba + k] = positions[bb + k];
                        positions[bb + k] = t;
                    }
                    let t = curFit[a]; curFit[a] = curFit[b]; curFit[b] = t;
                    t = curCC[a];  curCC[a]  = curCC[b];  curCC[b]  = t;
                    t = curPen[a]; curPen[a] = curPen[b]; curPen[b] = t;
                    // The step size belongs to the RUNG, not to the structure:
                    // it is the move scale that rung's temperature accepts.
                }
            }
            if (tried) swapRate = taken / tried;
        }

        /**
         * One dispatch: coordinates in, fitness and CC out.
         *
         * `duringWait` runs after the work is submitted and before the readback
         * is awaited - i.e. in the window where the GPU is busy and the main
         * thread would otherwise be idle. See fillDeviates().
         */
  

        /**
         * Records anything evaluated against the per-assignment best.
         *
         * Separate from the accept/reject test on purpose. A chain may refuse a
         * proposal and that proposal still be the best structure the run has
         * produced; throwing it away because one Markov chain declined to move
         * there would lose answers the search had already found.
         */
        function recordBest(fitArr, coords, scale) {
            for (let i = 0; i < numParticles; i++) {
                const A = assignOf[i], base = i * coordsPerParticle;
                const f = fitArr[i];
                if (!Number.isFinite(f)) continue;
                const cc = fitArr[numParticles + i];
                const rawPen = Number.isFinite(cc) ? (cc - f) / scale : 0;
                // Re-score the archived entry at the current weight before
                // comparing: what is stored came from an older, gentler ramp.
                const gRef = Number.isFinite(gBestCC[A]) ? gBestCC[A] - gBestPen[A] * scale
                                                        : gBestFit[A];
                if (f > gRef) {
                    gBestFit[A] = f; gBestCC[A] = cc; gBestPen[A] = rawPen;
                    gBestPos.set(coords.subarray(base, base + coordsPerParticle),
                                 A * coordsPerParticle);
                } else if (Number.isFinite(gRef)) {
                    gBestFit[A] = gRef;
                }
            }
        }

        /**
         * Same job as recordBest(), but reading the GPU's per-particle best
         * cache (bestFit/bestCC/bestPen/bestPositions, just filled in by
         * readStateFromGPU()) instead of a single generation's current fit.
         * That cache holds, for each chain, the best structure the GPU ever
         * accepted since the last sync - including ones a later Metropolis
         * step accepted away from - so calling this once per sync recovers
         * exactly what per-generation recordBest() used to catch every
         * generation.
         *
         * bestFit/bestCC/bestPen are stored RAW (weight 1), same as
         * gBestPen, so both sides are re-scored to the CURRENT penScale
         * before comparing - never compare the GPU's raw bestFit straight
         * across, it was written at whatever scale was active when that
         * particular structure was found.
         */
        function updateGlobalBestFromGPU(scale) {
            for (let i = 0; i < numParticles; i++) {
                if (bestFit[i] <= -1e29) continue;   // nothing recorded for this chain yet

                const A = assignOf[i], base = i * coordsPerParticle;
                const cc = bestCC[i], rawPen = bestPen[i];
                const fCandidate = cc - rawPen * scale;
                const gRef = Number.isFinite(gBestCC[A]) ? gBestCC[A] - gBestPen[A] * scale
                                                          : gBestFit[A];
                if (fCandidate > gRef) {
                    gBestFit[A] = fCandidate; gBestCC[A] = cc; gBestPen[A] = rawPen;
                    gBestPos.set(bestPositions.subarray(base, base + coordsPerParticle),
                                 A * coordsPerParticle);
                } else if (Number.isFinite(gRef)) {
                    gBestFit[A] = gRef;
                }
            }
        }

        /**
         * Greedy descent from every assignment's best, once the restarts are done.
         *
         * WHY THE SEARCH CANNOT DO THIS ITSELF. What an assignment stores is the
         * best structure ever PROPOSED for it - and a proposal is accepted or
         * rejected by one chain, at the temperature of the moment, under a
         * partial reflection set. Nothing ever goes back and refines that
         * particular structure. Each restart anneals a fresh set of chains and
         * the winner is whichever anneal happened to finish nearest a minimum,
         * which is why a run could return the right assignment with the sulfate
         * stretched to 1.65 A: the right basin, not the bottom of it.
         *
         * Three things are different here from the search proper:
         *
         *   - T = 0. Only improvements are kept. This is not sampling and makes
         *     no pretence of being; it is the last descent.
         *   - The FINAL objective: every reflection group, full penalty weight.
         *     That is what the results are reported at, so it is what the last
         *     refinement should be against. The stored correlation came from
         *     somewhere along the ramp and is re-measured before anything is
         *     compared to it.
         *   - The step shrinks geometrically to well below the precision the
         *     coordinates are printed at, so the descent actually finishes
         *     rather than rattling around the minimum.
         *
         * Every assignment in the wave is quenched, not only those that survived
         * pruning. The candidate table ranks on R, and quenching some structures
         * and not others would put "did it survive the prune" into that ranking -
         * which is the same bias against high-dimensional assignments that the
         * particle weighting and the late prune point exist to remove.
         */
        async function quench(candidateIds) {
            const live = candidateIds.filter(A => Number.isFinite(gBestFit[A]) &&
                                                  Number.isFinite(gBestCC[A]));
            if (!live.length) return;

// Per assignment, not just the leader. Reporting only the best CC
            // hides the case that matters: the leader is already converged and
            // eleven others move. Those others are what the R ranking compares.
            const before = new Float32Array(assignments.length);
            live.forEach(A => { before[A] = gBestCC[A]; });

            // The quench scores every candidate at FULL resolution, so the
            // ramp is over and the whole group list is in play.
            params[PARAM.nGroupsActive] = refl.nGroups;
            const scale = SW_DEFAULTS.penRampEnd;
            params[PARAM.penScale] = scale;

            // Spread the chains over the assignments and start every one of them
            // on its own assignment's best. Chains sharing a start diverge
            // immediately, since each draws its own moves.
            for (let i = 0; i < numParticles; i++) assignOf[i] = live[i % live.length];
            device.queue.writeBuffer(bufAssign, 0, assignOf);
            for (let i = 0; i < numParticles; i++) {
                const base = i * coordsPerParticle, gb = assignOf[i] * coordsPerParticle;
                positions.set(gBestPos.subarray(gb, gb + coordsPerParticle), base);
            }

            // Re-measure at the final objective before descending from it.
            const f0 = await evaluateCoordsSync(positions);
            recordBest(f0, positions, scale);
            for (let i = 0; i < numParticles; i++) curFit[i] = f0[i];

            const steps = SW_DEFAULTS.quenchSteps;
            for (let q = 0; q < steps; q++) {
                if (o.shouldStop && o.shouldStop()) { stopped = true; break; }
                const t = steps > 1 ? q / (steps - 1) : 1;
                const sd = SW_DEFAULTS.quenchStep0 *
                    Math.pow(SW_DEFAULTS.quenchStep1 / SW_DEFAULTS.quenchStep0, t);

                for (let i = 0; i < numParticles; i++) {
                    const base = i * coordsPerParticle;
                    for (let k = 0; k < coordsPerParticle; k++) {
                        const x = positions[base + k] + swGaussian() * sd;
                        proposals[base + k] = x - Math.floor(x);
                    }
                    projectSites(proposals, i, assignOf[i], T);
                }

                const fq = await evaluateCoordsSync(proposals);
                recordBest(fq, proposals, scale);
                for (let i = 0; i < numParticles; i++) {
                    const f = fq[i];
                    if (!Number.isFinite(f) || !(f > curFit[i])) continue;
                    const base = i * coordsPerParticle;
                    positions.set(proposals.subarray(base, base + coordsPerParticle), base);
                    curFit[i] = f;
                }

                if (o.onProgress && q % 25 === 0) await new Promise(r => setTimeout(r, 0));
            }

            let moved = 0, biggest = 0, leadBefore = -Infinity, leadAfter = -Infinity;
            for (const A of live) {
                const g = gBestCC[A] - before[A];
                if (g > 5e-5) { moved++; if (g > biggest) biggest = g; }
                if (before[A] > leadBefore) leadBefore = before[A];
                if (gBestCC[A] > leadAfter) leadAfter = gBestCC[A];
            }
            say(`Quench: ${live.length} assignment(s) at full resolution; ` +
                (moved
                    ? `${moved} improved, largest +${biggest.toFixed(4)} CC; `
                    : 'none improved; ') +
                (leadAfter - leadBefore > 5e-5
                    ? `best CC ${leadBefore.toFixed(4)} -> ${leadAfter.toFixed(4)}.`
                    : `best CC unchanged at ${leadAfter.toFixed(4)}. The coordinates are at the ` +
                      `bottom of their basins, so anything still wrong is the basin, not the polish.`));
        }

        /**
         * A kernel that compiles, dispatches and writes nothing leaves the
         * fitness buffer at zero for every chain. That is not an obvious
         * failure: the ranking still returns candidates, the assignments are
         * still valid Wyckoff choices, and the list looks entirely reasonable -
         * it is just the enumeration order with a correlation of zero attached.
         */
        function checkSpread(fitArr) {
            let lo = Infinity, hi = -Infinity, finite = 0;
            for (let i = 0; i < numParticles; i++) {
                const f = fitArr[i];
                if (!Number.isFinite(f)) continue;
                finite++; if (f < lo) lo = f; if (f > hi) hi = f;
            }
            if (finite === 0) {
                throw new Error('The GPU returned no finite fitness values. The kernel ran but ' +
                                'produced NaN or Inf everywhere - check the browser console for ' +
                                'WGSL validation errors.');
            }
            if (hi - lo < 1e-9) {
                throw new Error(`Every chain scored exactly ${hi.toFixed(6)}. The kernel is not ` +
                            `discriminating between structures, so any ranking would be the ` +
                            `enumeration order rather than a result. Check the console for ` +
                            `WGSL errors and that the reflection buffers are non-empty ` +
                            `(${refl.nGroups} groups, ${refl.nRefl} reflections).`);
        }
        say(`Generation 0 fitness spread ${lo.toFixed(4)} to ${hi.toFixed(4)} ` +
            `across ${finite} chain(s).`);
    }

      for (let restart = 0; restart < restarts && !stopped; restart++) {
        if (restart > 0) {
            // A restart re-seeds the chains but keeps every assignment's
            // global best, so the run accumulates rather than throwing away
            // what it found. The chain states themselves are discarded: they
            // are the end of a cooling schedule and would not move again.
            for (let i = 0; i < numParticles; i++) {
                const base = i * coordsPerParticle;
                for (let k = 0; k < coordsPerParticle; k++) positions[base + k] = Math.random();
                projectSites(positions, i, assignOf[i], T);
                curFit[i] = -Infinity; curCC[i] = NaN; curPen[i] = 0;
                stepSize[i] = SW_DEFAULTS.stepInit;
            }
            say(`wave ${wi + 1}, restart ${restart + 1} of ${restarts}` +
                (Number.isFinite(acceptRate) ? `; acceptance ${(acceptRate * 100).toFixed(0)}%` : '') +
                (Number.isFinite(swapRate) ? `, exchange ${(swapRate * 100).toFixed(0)}%.` : '.'));
        }

        // Each generation is one dispatch of PROPOSALS, so the GPU cost per
        // generation is identical to the swarm it replaced. The extra
        // evaluations of the current state are the few marked below.
let needCurrentEval = true, firstEval = (wi === 0);
        let lastRampLevel = -1;
        // Rebuilt per restart because a prune may have moved chains between
        // assignments, and a ladder that spans two assignments would be
        // exchanging structures between incompatible subspaces.
        ladders = buildLadders();
        
        const S32 = SW_STATE_STRIDE / 4;
        const stateData = new ArrayBuffer(numParticles * SW_STATE_STRIDE);
        const stateF32 = new Float32Array(stateData);
        const stateU32 = new Uint32Array(stateData);

        const syncStateToGPU = () => {
            for (let i = 0; i < numParticles; i++) {
                const b = i * S32;
                stateF32[b + 0] = curFit[i];
                stateF32[b + 1] = curCC[i];
                stateF32[b + 2] = curPen[i];
                stateF32[b + 3] = stepSize[i];
                stateF32[b + 4] = tempOf[i];
                stateU32[b + 5] = Math.floor(Math.random() * 4294967296);
                stateU32[b + 6] = 0;
                // Reset the GPU's per-particle best cache. Safe to reset here
                // because every call site harvests it (updateGlobalBestFromGPU)
                // before calling syncStateToGPU() - see readStateFromGPU()'s
                // callers below. A fresh block of generations starts a fresh
                // cache.
                stateF32[b + 8]  = -1e30;  // bestFit: "nothing recorded yet"
                stateF32[b + 9]  = 0;      // bestCC
                stateF32[b + 10] = 0;      // bestPen
                stateF32[b + 11] = 0;      // bestPad
            }
            device.queue.writeBuffer(bufState, 0, stateData);
            device.queue.writeBuffer(bufPos, 0, positions);
            // The kernel's accept counters were just zeroed, so the window the
            // acceptance rate is measured over restarts here too.
            mcmcSinceSync = 0;
        };
        syncStateToGPU();

        if (wi === 0 && restart === 0) {
            const R = Math.max(2, SW_DEFAULTS.rungs);
            say(`Replica exchange: ${ladders.length} ladder(s) of ${R} rungs, ` +
                `T ${SW_DEFAULTS.tempCold} to ${SW_DEFAULTS.tempHot}, ` +
                `swapping every ${SW_DEFAULTS.swapEvery} step(s).`);
        }

        for (let gen = 0; gen < generations; gen++) {
            if (o.shouldStop && o.shouldStop()) { stopped = true; break; }

            const t = generations > 1 ? gen / (generations - 1) : 1;

            const rampLevel = Math.min(SW_DEFAULTS.rampSteps - 1,
                Math.floor(t * SW_DEFAULTS.rampSteps));
                
            params[PARAM.nGroupsActive] = rampedGroupCount(
                refl.nGroups, rampLevel, SW_DEFAULTS.rampSteps - 1,
                o.rampStart ?? SW_DEFAULTS.rampStart, o.rampFull ?? SW_DEFAULTS.rampFull);
            if (rampLevel !== lastRampLevel) { lastRampLevel = rampLevel; needCurrentEval = true; }

            const scale = SW_DEFAULTS.penRampStart +
                t * (SW_DEFAULTS.penRampEnd - SW_DEFAULTS.penRampStart);
            params[PARAM.penScale] = scale;

            if (needCurrentEval) {
                // Rare - a resolution-ramp step or a post-prune reassignment,
                // not every generation - so a synchronous round trip here is
                // not the cost this file exists to remove. dispatchCoords()
                // re-scores the CURRENT GPU state under the objective that
                // just changed; readStateFromGPU() pulls that back AND
                // harvests whatever the per-particle best cache accumulated
                // during the generations since the last sync, so nothing
                // found in that window is lost to the reset syncStateToGPU()
                // is about to do.
                dispatchCoords(null, false);
                await readStateFromGPU();
                if (firstEval) { checkSpread(curFit); firstEval = false; }
                updateGlobalBestFromGPU(scale);
                // The freshly re-measured CURRENT state, packed the way
                // recordBest() expects, so it too is tested against gBestPos:
                // the GPU cache only tracks MCMC accepts, and this state was
                // just measured outside that path.
                const f0 = new Float32Array(numParticles * 2);
                f0.set(curFit, 0);
                f0.set(curCC, numParticles);
                recordBest(f0, positions, scale);
                syncStateToGPU();
                needCurrentEval = false;
            }

            /* ---- GPU-resident MCMC dispatch: enqueue and move on ---- */
            dispatchCoords(null, true);

            const atSwap = ladders.length && (gen + 1) % SW_DEFAULTS.swapEvery === 0;

            if (atSwap) {
                // The ONE synchronization point for a normal generation. The
                // GPU has been proposing, accepting, and caching its own
                // per-particle best for up to swapEvery generations; pull
                // that back now, fold it into gBestPos, do the CPU-side
                // Replica Exchange swap, then push the (possibly swapped)
                // state back down with a freshly reset best cache.
                await readStateFromGPU();
                updateGlobalBestFromGPU(scale);
                exchangeSweep(((gen + 1) / SW_DEFAULTS.swapEvery) % 2, scale);
                syncStateToGPU();
            }

            // Successive halving. Most assignments are visibly hopeless within
            // a few tens of generations, and the budget is far better spent on
            // the few that are not.
            //
            // Pruned once per WAVE, not once per restart. `ids` persists across
            // restarts, so a per-restart prune compounds: on PbSO4 it went
            // 35 -> 18 -> 9 -> 5 -> 4 over the first four restarts of ten, and
            // the field was decided before any restart had a chance to
            // converge. Worse, it is biased - an assignment with eleven free
            // parameters needs longer to look good than one with six, so the
            // flexible ones, which include the right answer for anything but
            // the simplest structure, are retired first.
            const pruneAt = Math.max(20, Math.round(generations * 0.45));
            if (restart === 0 && gen > 0 && gen === pruneAt && ids.length > 4) {
                // gBestFit is only as fresh as the last sync, and pruneAt has
                // no reason to land on a swapEvery boundary. Force one now
                // rather than prune on stale numbers - up to swapEvery - 1
                // generations of "which assignments are hopeless" would
                // otherwise be missing from the decision.
                if (!atSwap) {
                    await readStateFromGPU();
                    updateGlobalBestFromGPU(scale);
                }
                const keep = pruneAssignments(ids, gBestFit, 0.5, 4);
                if (keep.length < ids.length) {
                    ids = keep;
                    // Reallocate by dimension, not round-robin. This is where the
                    // budget finally becomes generous - half the assignments have
                    // just retired - so distributing it uniformly here would undo
                    // the weighting at the moment it matters most.
                    const fresh = weightedAssign(ids, numParticles, minPer, dims);
                    for (let i = 0; i < numParticles; i++) {
                        const A = fresh[i];
                        if (assignOf[i] === A) continue;
                        assignOf[i] = A;
                        const base = i * coordsPerParticle;
                        // A reassigned chain starts fresh: its state belongs to a
                        // different subspace and is not a point of the new one.
                        for (let k = 0; k < coordsPerParticle; k++) {
                            positions[base + k] = Math.random();
                        }
                        projectSites(positions, i, A, T);
                        stepSize[i] = SW_DEFAULTS.stepInit;
                        curFit[i] = -Infinity; curCC[i] = NaN; curPen[i] = 0;
                    }
                    device.queue.writeBuffer(bufAssign, 0, assignOf);
                    // Push the reassigned chains' fresh random coordinates to
                    // the GPU immediately. needCurrentEval alone is not
                    // enough: the next needCurrentEval dispatch passes
                    // coords=null and would otherwise re-score whatever
                    // stale, pre-reassignment coordinates are still sitting
                    // in bufPos instead of these.
                    syncStateToGPU();
                    // Those chains have no energy for their new subspace, so the
                    // next generation has to measure one before it can test
                    // anything against it.
                    needCurrentEval = true;
                    ladders = buildLadders();
                    say(`gen ${gen}: ${ids.length} assignment(s) still in contention.`);
                }
            }

            if (o.onProgress && (gen % 10 === 0 || gen === generations - 1)) {
                let best = -Infinity, bestA = -1;
                for (const a of ids) if (gBestFit[a] > best) { best = gBestFit[a]; bestA = a; }
                if (best > runBestFit) {
                    runBestFit = best; runBestCC = gBestCC[bestA]; runBestAssign = bestA;
                }

                // Hand out the running best as sites, so the caller can rebuild
                // the calculated and difference maps as the search proceeds.
                // Watching those two panes converge is how a wrong atom shows
                // itself - a correlation climbing on its own says only that
                // something is improving, not what.
                // zn is the atomic number; x, y, z are fractional coordinates.
                // Naming both "z" is how the two get silently swapped.
                let bestSites = null;
                if (bestA >= 0) {
                    const A = assignments[bestA];
                    const gb = bestA * coordsPerParticle;
                    bestSites = A.sites.map((st, si) => ({
                        element: st.element,
                        zn: st.z,
                        multiplicity: st.w.multiplicity,
                        wyckoff: `${st.w.multiplicity}${st.w.letter}`,
                        siteSymmetry: st.w.site_symmetry,
                        x: gBestPos[gb + si * 3],
                        y: gBestPos[gb + si * 3 + 1],
                        z: gBestPos[gb + si * 3 + 2]
                    }));
                }
                o.onProgress({ wave: wi + 1, waves: waves.length, generation: gen,
                               generations, restart: restart + 1, restarts,
                               // Whether the chains are moving at all. Without
                               // it a stuck run and a working one look the same
                               // from outside: both report a fitness that is not
                               // going up.
                               acceptRate, swapRate,
                               // `best` is the run's best so far, so the trace
                               // only ever climbs; `waveBest` is this wave's own.
                               best: runBestFit, waveBest: best,
                               cc: runBestCC, active: ids.length,
                               assignment: bestA >= 0
                                   ? assignments[bestA].sites.map(x => `${x.element} ${x.w.multiplicity}${x.w.letter}`).join(', ')
                                   : null,
                               bestSites });
                await new Promise(r => setTimeout(r, 0));
            }
        }

      }

        if (!stopped) await quench(waves[wi].ids);

        for (const A of waves[wi].ids) {
            if (!Number.isFinite(gBestFit[A]) || gBestFit[A] <= -Infinity) continue;
            // The reported penalty is the RAW one, at weight 1. Reporting it at
            // whatever the ramp happened to reach makes two candidates recorded
            // at different moments incomparable, and makes the number depend on
            // a tuning constant rather than on the structure.
            results.push({ assignIdx: A, score: gBestCC[A] - gBestPen[A],
                           cc: Number.isFinite(gBestCC[A]) ? gBestCC[A] : gBestFit[A],
                           penalty: Number.isFinite(gBestCC[A]) ? gBestPen[A] : 0,
                           coords: gBestPos.slice(A * coordsPerParticle,
                                                  (A + 1) * coordsPerParticle) });
        }
    }

    /* ---- 9. Report ---- */
    // Ranked by CC alone - the agreement between the observed and calculated
    // Patterson maps, which is the quantity that actually means something.
    //
    // The penalty weights are arbitrary: 0.05 per clash, 0.02 per Angstrom, no
    // more principled than any other pair of numbers. They earn their place
    // inside the search, where they keep particles away from nonsense and cost
    // nothing if slightly wrong, but ranking on them would let an arbitrary
    // constant decide which structure is reported. Distance windows are
    // enforced afterwards as a filter, which needs no weight at all.
    results.sort((a, b) => b.cc - a.cc);
    // More than the UI shows. Chemistry filtering happens in the caller, where
    // the contact code lives, and filtering a list of five can leave one.
    const topN = o.topN ?? SW_DEFAULTS.topN;
    const top = results.slice(0, topN).map(r => {
        const A = assignments[r.assignIdx];
        const sites = A.sites.map((s, si) => ({
            element: s.element,
            wyckoff: `${s.w.multiplicity}${s.w.letter}`,
            multiplicity: s.w.multiplicity,
            siteSymmetry: s.w.site_symmetry,
            // Occupancy is 1.0 and the multiplicity comes from the Wyckoff
            // letter. Scaling occupancy by mult/order - which a general-position
            // model needs - would double-count here.
            occupancy: 1.0,
            // Carried out so a downstream coordinate refinement can rebuild
            // the structure factor and stay on the subspace without having to
            // re-derive either. `zn` is the ATOMIC NUMBER and x, y, z are
            // fractional coordinates - naming both "z" is how the two get
            // silently swapped, so they are never both called z here.
            zn: s.z,
            w: s.w,
            x: r.coords[si * 3], y: r.coords[si * 3 + 1], z: r.coords[si * 3 + 2]
        }));
        return { cc: r.cc, score: r.score, penalty: r.penalty, assignment: A.sites.map(s => `${s.element} ${s.w.multiplicity}${s.w.letter}`).join(', '),
                 harkerResidual: A.harkerResidual, sites, nFreeParams:
                     A.sites.reduce((n, s) => n + wyckoffFreedom(s.w), 0) };
    });

    say(top.length ? `Best score ${top[0].score.toFixed(4)} (CC ${top[0].cc.toFixed(4)}) for ${top[0].assignment}.`
                   : 'No candidate produced a finite score.');
    return { candidates: top, all: results.length, assignments: assignments.length,
             // How many the solver kept, so the caller can say so rather than
             // leaving the user to wonder why the table is shorter than the
             // number of assignments searched.
             kept: top.length, topN,
             // The clash floors the search actually enforced, so the caller's
             // filter can reject on the same numbers instead of its own.
             floors: restraints.floors,
             overallB, route: obs.route, reflections: obs.rows.length,
             groups: refl.nGroups, stopped, log,
             // The caller has to group the reflections the same way to count
             // independent observations for a significance test; handing back
             // the tolerance keeps the two from drifting apart.
             overlapTol: refl.overlapTol,
             // Kept so the caller can compute an R factor for whichever
             // candidate is selected, without re-normalising the file.
             obsRows: obs.rows, demand };
    } finally {
        // Guarded individually: a throw before the later buffers were created
        // leaves those undefined, and a cleanup that itself throws would mask
        // the error that caused it.
        for (const b of [bufGen, bufStatic, bufRefl, bufGroup, bufPos, bufAssign,
                         bufState, bufPosRead, bufStateRead, bufParams]) {
            if (b) { try { b.destroy(); } catch (e) {} }
        }
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { SHARKO_SWARM_WYCKOFF_VERSION, runWyckoffSearch, swPackReflections, swScatteringTable,
                       swMaxGenAtoms, rampedGroupCount };
}
