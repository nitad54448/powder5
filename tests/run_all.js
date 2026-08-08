// test/test_refactor_equivalence.js
// ---------------------------------------------------------------------------
//  Proves that constants.js + crystal.js + profile.js reproduce the
//  PRE-REFACTOR implementation BIT FOR BIT, over a wide sweep of profile
//  types, crystal systems, angles and reflections.
//
//  Run:  node test/test_refactor_equivalence.js
//
//  A refactor of numerical code is only safe if it is provably a no-op. This
//  is that proof. If it ever fails, the shared modules have drifted from the
//  behaviour the fits were validated against.
// ---------------------------------------------------------------------------

'use strict';

const ref = require('./reference_original.js');

// Load the shared modules the way a worker does: into the global scope, in
// order, so profile.js can see constants.js's declarations.
const C = require('../constants.js');
Object.assign(globalThis, C);
// MIN_PROFILE_FWHM_DEG is a live getter on the exports object; profile.js
// reads it as a bare global, so mirror it explicitly and keep it in sync.
function syncMin() { globalThis.MIN_PROFILE_FWHM_DEG = C.MIN_PROFILE_FWHM_DEG; }
syncMin();

const P = require('../profile.js');
const X = require('../crystal.js');

let failures = 0, checks = 0;

/** Bit-exact comparison, with NaN === NaN and +-0 distinguished. */
function same(a, b) {
    if (typeof a !== typeof b) return false;
    if (typeof a !== 'number') return a === b;
    if (Number.isNaN(a) && Number.isNaN(b)) return true;
    return Object.is(a, b);
}

function check(label, got, want) {
    checks++;
    if (!same(got, want)) {
        failures++;
        if (failures <= 20) {
            console.error(`  FAIL ${label}\n       got  ${got}\n       want ${want}`);
        }
    }
}

// ---------------------------------------------------------------------------
//  Deterministic parameter sweep
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
const rnd = mulberry32(20260808);
const pick = arr => arr[Math.floor(rnd() * arr.length)];
const span = (lo, hi) => lo + rnd() * (hi - lo);

function makeParams(profileType) {
    const p = {
        profileType,
        lambda: 1.5405981,
        zeroShift: span(-0.2, 0.2),
        SL: rnd() < 0.5 ? 0 : span(0, 0.05),
        HL: rnd() < 0.5 ? 0 : span(0, 0.05)
    };
    if (profileType === 'simple_pvoigt') {
        Object.assign(p, {
            GU: span(0, 20), GV: span(-8, 8), GW: span(-5, 30),
            GP: span(0, 5), LX: span(0, 15), eta: rnd(),
            shft: span(-40, 40), trns: span(-0.02, 0.02)
        });
    } else if (profileType === 'split_pvoigt') {
        Object.assign(p, {
            GU_L: span(0, 20), GV_L: span(-8, 8), GW_L: span(-5, 30), LX_L: span(0, 15),
            GU_R: span(0, 20), GV_R: span(-8, 8), GW_R: span(-5, 30), LX_R: span(0, 15),
            eta_split: rnd(), shft_split: span(-40, 40), trns_split: span(-0.02, 0.02)
        });
    } else {
        Object.assign(p, {
            U: span(0, 20), V: span(-8, 8), W: span(-5, 30), X: span(0, 12), Y: span(0, 12),
            S400: span(0, 400), S040: span(0, 400), S004: span(0, 400),
            S220: span(-200, 200), S202: span(-200, 200), S022: span(-200, 200)
        });
        if (rnd() < 0.5) p.eta_aniso = rnd();
    }
    return p;
}

// ---------------------------------------------------------------------------
console.log('1. profile: calculatePeakShift / calculateProfileWidths / getPeakFWHM');
// ---------------------------------------------------------------------------
for (let trial = 0; trial < 400; trial++) {
    const stepDeg = span(0.005, 0.05);
    const axis = [];
    for (let i = 0; i < 200; i++) axis.push(5 + i * stepDeg);
    ref.setMinProfileFwhmFromAxis(axis);
    C.setMinProfileFwhmFromAxis(axis);
    syncMin();

    const pt = pick(['simple_pvoigt', 'split_pvoigt', 'tch_aniso']);
    const params = makeParams(pt);
    const hkl = {
        h: Math.round(span(-6, 6)), k: Math.round(span(-6, 6)),
        l: Math.round(span(-6, 6)), d: span(0.7, 12)
    };

    for (const tth of [0.5, 3.2, 17.9, 45.0, 88.7, 120.4, 165.3, 179.5]) {
        check(`peakShift ${pt} @${tth}`,
              P.calculatePeakShift(tth, params), ref.calculatePeakShift(tth, params));

        for (const side of ['left', 'right', 'center']) {
            const g = P.calculateProfileWidths(tth, hkl, params, side);
            const w = ref.calculateProfileWidths(tth, hkl, params, side);
            check(`widths.G ${pt} ${side} @${tth}`, g.gamma_G, w.gamma_G);
            check(`widths.L ${pt} ${side} @${tth}`, g.gamma_L, w.gamma_L);
            check(`peakFWHM ${pt} ${side} @${tth}`,
                  P.getPeakFWHM(g.gamma_G, g.gamma_L),
                  ref.getPeakFWHM(w.gamma_G, w.gamma_L));
        }
    }
}

// ---------------------------------------------------------------------------
console.log('2. profile: prepareVoigt scalars and evalVoigt across the window');
// ---------------------------------------------------------------------------
const PREP_KEYS = ['type', 'x0', 'asym_param', 'eta', 'max_d', 'pedestal',
                   'fwhm_total', 'analyticArea', 'Cg', 'cL', 'cG',
                   'H_G', 'H_L', 'fwhm',
                   'H_G_L', 'H_L_L', 'H_G_R', 'H_L_R',
                   'cL_L', 'cG_L', 'cL_R', 'cG_R'];

for (let trial = 0; trial < 300; trial++) {
    const stepDeg = span(0.005, 0.05);
    const axis = [];
    for (let i = 0; i < 200; i++) axis.push(5 + i * stepDeg);
    ref.setMinProfileFwhmFromAxis(axis);
    C.setMinProfileFwhmFromAxis(axis);
    syncMin();

    const pt = pick(['simple_pvoigt', 'split_pvoigt', 'tch_aniso']);
    const params = makeParams(pt);
    const hkl = {
        h: Math.round(span(-6, 6)), k: Math.round(span(-6, 6)),
        l: Math.round(span(-6, 6)), d: span(0.7, 12)
    };
    const basePos = span(3, 175);
    const x0 = basePos + P.calculatePeakShift(basePos, params);

    const gp = P.prepareVoigt(basePos, x0, hkl, params);
    const wp = ref.prepareVoigt(basePos, x0, hkl, params);

    for (const k of PREP_KEYS) {
        if (wp[k] === undefined && gp[k] === undefined) continue;
        check(`prep.${k} ${pt}`, gp[k], wp[k]);
    }

    // Sample the profile densely inside the window and just outside it.
    for (let i = 0; i <= 60; i++) {
        const frac = -1.15 + (2.30 * i) / 60;
        const x = x0 + frac * wp.max_d;
        check(`evalVoigt ${pt} f=${frac.toFixed(3)}`,
              P.evalVoigt(x, gp), ref.evalVoigt(x, wp));
    }
    // Exactly on the centre and exactly on both edges.
    for (const x of [x0, x0 + wp.max_d, x0 - wp.max_d]) {
        check(`evalVoigt ${pt} edge`, P.evalVoigt(x, gp), ref.evalVoigt(x, wp));
    }
}

// ---------------------------------------------------------------------------
console.log('3. crystal: updateHklPositions, every system, every monoclinic axis');
// ---------------------------------------------------------------------------
const SYSTEMS = ['cubic', 'tetragonal', 'orthorhombic', 'hexagonal', 'trigonal',
                 'rhombohedral', 'monoclinic', 'triclinic'];

const hklTemplate = [];
for (let h = -5; h <= 5; h++)
    for (let k = -5; k <= 5; k++)
        for (let l = -5; l <= 5; l++)
            hklTemplate.push({ h_orig: h, k_orig: k, l_orig: l });
// A few malformed entries: the null and missing-index paths must match too.
hklTemplate.push(null, {}, { h_orig: 1 }, { h_orig: 0, k_orig: 0, l_orig: 0 });

for (const system of SYSTEMS) {
    for (let trial = 0; trial < 12; trial++) {
        const params = {
            a: span(3, 18), b: span(3, 18), c: span(3, 25),
            alpha: 90, beta: 90, gamma: 90,
            lambda: pick([1.5405981, 0.7093, 2.2897, 0.4])
        };
        if (system === 'triclinic') {
            params.alpha = span(70, 110);
            params.beta  = span(70, 110);
            params.gamma = span(70, 110);
        } else if (system === 'monoclinic') {
            // Exercise all three unique-axis settings.
            const axis = ['a', 'b', 'c'][trial % 3];
            params[X.MONOCLINIC_ANGLE_FOR_AXIS[axis]] = span(75, 125);
        } else if (system === 'hexagonal' || system === 'trigonal'
                || system === 'rhombohedral') {
            params.gamma = 120;
        }

        const A = hklTemplate.map(p => (p ? { ...p } : null));
        const B = hklTemplate.map(p => (p ? { ...p } : null));
        X.updateHklPositions(A, params, system);
        ref.updateHklPositions(B, params, system);

        for (let i = 0; i < A.length; i++) {
            if (A[i] === null) { check(`${system} null preserved`, B[i], null); continue; }
            check(`${system} tth[${i}]`, A[i].tth, B[i].tth);
            check(`${system} d[${i}]`,   A[i].d,   B[i].d);
        }
    }
}

// Degenerate and hostile inputs.
for (const bad of [
    { a: 0, lambda: 1.54 },
    { a: 5, lambda: 0 },
    { a: 5, b: 5, c: 5, alpha: 179.99, beta: 179.99, gamma: 179.99, lambda: 1.54 },
    { a: 5, b: 5, c: 5, alpha: 0.001, beta: 90, gamma: 90, lambda: 1.54 },
    { a: 1e-9, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90, lambda: 1.54 }
]) {
    for (const system of SYSTEMS) {
        const A = hklTemplate.map(p => (p ? { ...p } : null));
        const B = hklTemplate.map(p => (p ? { ...p } : null));
        X.updateHklPositions(A, bad, system);
        ref.updateHklPositions(B, bad, system);
        for (let i = 0; i < A.length; i++) {
            if (A[i] === null) continue;
            check(`degenerate ${system} tth[${i}]`, A[i].tth, B[i].tth);
            check(`degenerate ${system} d[${i}]`,   A[i].d,   B[i].d);
        }
    }
}

// ---------------------------------------------------------------------------
console.log('4. crystal: buildInvDsqEvaluator agrees with updateHklPositions');
// ---------------------------------------------------------------------------
// These two used to be separate implementations in separate files, with a
// comment promising they agreed. Now one is built from the other, and this
// asserts it. A generator that disagrees with the 2-theta filter silently
// drops or invents reflections.
let worstEval = 0, nEval = 0;
for (const system of SYSTEMS) {
    for (const axis of ['a', 'b', 'c']) {
        const params = {
            a: 7.31, b: 9.02, c: 11.47,
            alpha: 90, beta: 90, gamma: 90, lambda: 1.5405981
        };
        if (system === 'triclinic') {
            params.alpha = 87.2; params.beta = 103.4; params.gamma = 95.9;
        } else if (system === 'monoclinic') {
            params[X.MONOCLINIC_ANGLE_FOR_AXIS[axis]] = 104.3;
        } else if (system === 'hexagonal' || system === 'trigonal'
                || system === 'rhombohedral') {
            params.gamma = 120;
        }

        const ev = X.buildInvDsqEvaluator(system, params);
        checks++;
        if (!ev) { failures++; console.error(`  FAIL no evaluator for ${system}`); continue; }

        const list = [];
        for (let h = -4; h <= 4; h++)
            for (let k = -4; k <= 4; k++)
                for (let l = -4; l <= 4; l++)
                    if (h || k || l) list.push({ h_orig: h, k_orig: k, l_orig: l });
        X.updateHklPositions(list, params, system);

        for (const q of list) {
            if (q.d === null) continue;
            const fromPositions = 1 / (q.d * q.d);
            const fromEvaluator = ev.f(q.h_orig, q.k_orig, q.l_orig);
            worstEval = Math.max(worstEval,
                Math.abs(fromPositions - fromEvaluator) / fromPositions);
            nEval++;
        }

        // lSeparable must be honest: with no h*l or k*l cross term, 1/d^2 is
        // monotonic in |l| at fixed h,k, which is what lets a generator break
        // out of the l loop early. If the flag lies, reflections go missing.
        if (ev.lSeparable) {
            for (let h = -3; h <= 3; h++) {
                for (let k = -3; k <= 3; k++) {
                    let prev = -Infinity;
                    for (let l = 0; l <= 6; l++) {
                        const v = ev.f(h, k, l);
                        checks++;
                        if (v < prev - 1e-12) {
                            failures++;
                            console.error(`  FAIL ${system}: lSeparable but 1/d^2 ` +
                                          `not monotonic in l at h=${h} k=${k}`);
                        }
                        prev = v;
                    }
                }
            }
        }
        if (system !== 'monoclinic') break;
    }
}
checks++;
if (!(worstEval < 1e-13)) {
    failures++;
    console.error(`  FAIL evaluator/positions disagreement: ${worstEval}`);
}
console.log(`   ${nEval} reflections, max rel err ${worstEval.toExponential(2)}`);

// ---------------------------------------------------------------------------
console.log('5. crystal: the general metric reproduces the closed forms exactly');
// ---------------------------------------------------------------------------
function closedMonoclinic(a, b, c, ang, axis, h, k, l) {
    const r = ang * Math.PI / 180, s = Math.sin(r), cs = Math.cos(r);
    if (axis === 'b') return (h * h / (a * a) + l * l / (c * c) - 2 * h * l * cs / (a * c)) / (s * s) + k * k / (b * b);
    if (axis === 'a') return (k * k / (b * b) + l * l / (c * c) - 2 * k * l * cs / (b * c)) / (s * s) + h * h / (a * a);
    return (h * h / (a * a) + k * k / (b * b) - 2 * h * k * cs / (a * b)) / (s * s) + l * l / (c * c);
}
let worstMono = 0;
for (const axis of ['a', 'b', 'c']) {
    for (const ang of [72.5, 95.3, 102.7, 118.4]) {
        const cell = { a: 7.31, b: 9.02, c: 11.47, alpha: 90, beta: 90, gamma: 90 };
        cell[X.MONOCLINIC_ANGLE_FOR_AXIS[axis]] = ang;
        const f = X.buildInvDsq(cell);
        for (let h = -5; h <= 5; h++)
            for (let k = -5; k <= 5; k++)
                for (let l = -5; l <= 5; l++) {
                    if (!h && !k && !l) continue;
                    const g = X.invDsq(f, h, k, l);
                    const e = closedMonoclinic(cell.a, cell.b, cell.c, ang, axis, h, k, l);
                    worstMono = Math.max(worstMono, Math.abs(g - e) / Math.abs(e));
                }
    }
}
checks++;
if (!(worstMono < 1e-13)) {
    failures++;
    console.error(`  FAIL monoclinic closed-form agreement: ${worstMono}`);
}
console.log(`   monoclinic a/b/c vs closed form: max rel err ${worstMono.toExponential(2)}`);

// ---------------------------------------------------------------------------
console.log(`\n${checks} checks, ${failures} failures`);
if (failures) { console.error('REGRESSION: the refactor is NOT a no-op.'); process.exit(1); }
console.log('The shared modules are bit-for-bit identical to the originals.');
