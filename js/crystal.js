// crystal.js
// ---------------------------------------------------------------------------
//  SINGLE SOURCE OF TRUTH for lattice geometry: metric tensors, cell volumes,
//  fractional -> Cartesian conversion, d-spacings and the monoclinic
//  unique-axis convention.
//
//  Why this file exists
//  --------------------
//  Before this module the same geometry was implemented four times:
//
//    powder5.html          updateHklPositions  (general triclinic metric)
//    refinement_worker.js  updateHklPositions  (byte-identical copy)
//    charge_flipping_worker.js  metricTensor + cellVolume + fracDistance
//    density3d.js          cellBasis + fracToCart
//
//  and the MONOCLINIC UNIQUE AXIS was handled inconsistently across them:
//  powder5.html keyed off MONOCLINIC_ANGLE_FOR_AXIS, the workers spoke only of
//  beta, and updateHklPositions quietly did the right thing for all three axes
//  by using the general triclinic metric. That last point is worth stating
//  plainly, because it is easy to misread the code and "fix" it into being
//  wrong: for monoclinic, the unique axis is carried IMPLICITLY by whichever
//  of alpha/beta/gamma is not 90. Nothing needs to be told which axis it is,
//  as long as the caller has put the oblique angle in the right slot. The
//  helpers below make that contract explicit rather than implicit.
//
//  Loading
//  -------
//    browser :  script src="crystal.js"
//                 (written without angle brackets on purpose: a literal
//                  closing script tag inside this file would terminate the
//                  page early if anyone ever inlined it into powder5.html)   (after constants.js)
//    worker  :  importScripts('crystal.js')
//    node    :  require('./crystal.js')
// ---------------------------------------------------------------------------

'use strict';

/** @type {number} */
const DEG2RAD = Math.PI / 180;

/**
 * A unit cell.
 * @typedef {object} Cell
 * @property {number} a     Edge a, angstrom.
 * @property {number} b     Edge b, angstrom.
 * @property {number} c     Edge c, angstrom.
 * @property {number} [alpha] Angle alpha, degrees. Default 90.
 * @property {number} [beta]  Angle beta,  degrees. Default 90.
 * @property {number} [gamma] Angle gamma, degrees. Default 90.
 */

// ===========================================================================
//  1. Monoclinic unique-axis convention
// ===========================================================================

/**
 * Which cell ANGLE is the oblique one for each monoclinic unique axis.
 * Unique axis b (the common setting) => beta != 90.
 * @type {{a: string, b: string, c: string}}
 */
const MONOCLINIC_ANGLE_FOR_AXIS = { a: 'alpha', b: 'beta', c: 'gamma' };

/**
 * The two angles that are constrained to 90 for a given monoclinic unique
 * axis. Useful for enforcing the constraint on a parameter set.
 * @type {{a: string[], b: string[], c: string[]}}
 */
const MONOCLINIC_FIXED_ANGLES_FOR_AXIS = {
    a: ['beta', 'gamma'],
    b: ['alpha', 'gamma'],
    c: ['alpha', 'beta']
};

/**
 * Normalises a parameter object into a full {a,b,c,alpha,beta,gamma} cell,
 * applying the symmetry constraints of the crystal system.
 *
 * This is the ONE place the monoclinic convention is applied to numbers. Call
 * it before any geometry helper below when the parameters come from a UI or a
 * refinement vector, where the constrained angles may be stale or absent.
 *
 * @param {object} p        Parameters carrying at least `a`.
 * @param {string} system   Crystal system, lower case.
 * @param {('a'|'b'|'c')} [uniqueAxis='b'] Monoclinic unique axis. Ignored for
 *        every other system.
 * @returns {Cell|null} A full cell, or null if `p` has no usable `a`.
 */
function normaliseCell(p, system, uniqueAxis = 'b') {
    if (!p || !(p.a > 0)) return null;
    const s = String(system || '').toLowerCase();
    const a = p.a;
    const b = (p.b > 0) ? p.b : a;
    const c = (p.c > 0) ? p.c : a;
    const num = (v, d) => (typeof v === 'number' && isFinite(v)) ? v : d;

    switch (s) {
        case 'cubic':
            return { a, b: a, c: a, alpha: 90, beta: 90, gamma: 90 };
        case 'tetragonal':
            return { a, b: a, c, alpha: 90, beta: 90, gamma: 90 };
        case 'hexagonal':
        case 'trigonal':
            return { a, b: a, c, alpha: 90, beta: 90, gamma: 120 };
        case 'rhombohedral':
            // Hexagonal axes throughout this program (see the note in
            // powder5.html's setting resolver), so same as trigonal.
            return { a, b: a, c, alpha: 90, beta: 90, gamma: 120 };
        case 'orthorhombic':
            return { a, b, c, alpha: 90, beta: 90, gamma: 90 };
        case 'monoclinic': {
            const axis = (uniqueAxis === 'a' || uniqueAxis === 'c') ? uniqueAxis : 'b';
            const oblique = MONOCLINIC_ANGLE_FOR_AXIS[axis];
            const cell = { a, b, c, alpha: 90, beta: 90, gamma: 90 };
            cell[oblique] = num(p[oblique], 90);
            return cell;
        }
        case 'triclinic':
        default:
            return {
                a, b, c,
                alpha: num(p.alpha, 90),
                beta:  num(p.beta,  90),
                gamma: num(p.gamma, 90)
            };
    }
}

// ===========================================================================
//  2. Metric tensors
// ===========================================================================

/**
 * REAL-space metric tensor G, such that
 *     |v|^2 = x^T G x   for a fractional vector x.
 *
 * @param {Cell} cell
 * @returns {number[][]} Symmetric 3x3, angstrom^2.
 */
function metricTensor(cell) {
    const a = cell.a, b = cell.b, c = cell.c;
    const ca = Math.cos((cell.alpha ?? 90) * DEG2RAD);
    const cb = Math.cos((cell.beta  ?? 90) * DEG2RAD);
    const cg = Math.cos((cell.gamma ?? 90) * DEG2RAD);
    return [
        [a * a,      a * b * cg, a * c * cb],
        [a * b * cg, b * b,      b * c * ca],
        [a * c * cb, b * c * ca, c * c]
    ];
}

/**
 * Distance in angstrom between two points given their FRACTIONAL separation,
 * reduced to the shortest lattice image on each axis.
 *
 * @param {number[][]} G Real-space metric tensor from metricTensor().
 * @param {number} dx Fractional separation along a.
 * @param {number} dy Fractional separation along b.
 * @param {number} dz Fractional separation along c.
 * @returns {number} Distance, angstrom.
 */
function fracDistance(G, dx, dy, dz) {
    dx -= Math.round(dx); dy -= Math.round(dy); dz -= Math.round(dz);

    const q = (x, y, z) => x * (G[0][0] * x + G[0][1] * y + G[0][2] * z)
                         + y * (G[1][0] * x + G[1][1] * y + G[1][2] * z)
                         + z * (G[2][0] * x + G[2][1] * y + G[2][2] * z);

    // ROUNDING EACH COMPONENT INTO [-1/2, 1/2) IS NOT THE MINIMUM IMAGE.
    //
    // It is exact for an orthogonal cell and only for an orthogonal cell. When
    // the axes are oblique the shortest lattice vector between two points can
    // lie in a neighbouring translation that componentwise rounding never
    // visits, and the error is always in the same direction: the distance
    // comes out TOO LONG. So a real clash reads as a comfortable separation, a
    // symmetry image reads as a distinct atom, and two peaks that are the same
    // peak survive deduplication -- all silently, and all only in cells that
    // are not orthorhombic or higher.
    //
    // The off-diagonal terms of the metric tensor are exactly zero for an
    // orthogonal cell, so the cheap path is taken whenever they are, which is
    // most of the time. Only an oblique cell pays for the 27-point search,
    // and it is the only one that needs it.
    // Tested RELATIVE to the diagonal: G entries scale as a*b, so on a 100 A
    // cell the residual off-diagonal from floating-point rounding is far
    // larger than any fixed absolute epsilon would allow.
    const scale = Math.abs(G[0][0]) + Math.abs(G[1][1]) + Math.abs(G[2][2]);
    const eps = 1e-12 * (scale > 0 ? scale : 1);
    const oblique = Math.abs(G[0][1]) > eps ||
                    Math.abs(G[0][2]) > eps ||
                    Math.abs(G[1][2]) > eps;
    if (!oblique) return Math.sqrt(Math.max(0, q(dx, dy, dz)));

    // +/-2, not +/-1. One shell still missed about one oblique case in 200,
    // and the ones it missed were the badly wrong ones -- a shortest vector
    // two translations away turns up exactly when the cell is most skewed.
    // 125 metric evaluations sounds like a lot but only an oblique cell ever
    // reaches this branch, and the alternative is a Niggli reduction that
    // would have to be threaded through every caller.
    let best = Infinity;
    for (let a = -2; a <= 2; a++) {
        for (let b = -2; b <= 2; b++) {
            for (let c = -2; c <= 2; c++) {
                const s = q(dx + a, dy + b, dz + c);
                if (s < best) best = s;
            }
        }
    }
    return Math.sqrt(Math.max(0, best));
}

/**
 * Unit-cell volume.
 * @param {Cell} cell
 * @returns {number} Volume, angstrom^3.
 */
function cellVolume(cell) {
    const ca = Math.cos((cell.alpha ?? 90) * DEG2RAD);
    const cb = Math.cos((cell.beta  ?? 90) * DEG2RAD);
    const cg = Math.cos((cell.gamma ?? 90) * DEG2RAD);
    const t = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg;
    return cell.a * cell.b * cell.c * Math.sqrt(Math.max(1e-12, t));
}

/**
 * Fractional -> Cartesian basis, in the standard setting: a along x, b in the
 * xy plane. Returned as three basis vectors (the columns of the conversion
 * matrix), which is what a vertex transform and a cell wireframe both want.
 *
 * @param {Cell} cell
 * @returns {{av: number[], bv: number[], cv: number[]}}
 */
function cellBasis(cell) {
    const a = cell.a, b = cell.b, c = cell.c;
    const al = (cell.alpha ?? 90) * DEG2RAD;
    const be = (cell.beta  ?? 90) * DEG2RAD;
    const ga = (cell.gamma ?? 90) * DEG2RAD;
    const ca = Math.cos(al), cb = Math.cos(be), cg = Math.cos(ga), sg = Math.sin(ga);
    const vol = Math.sqrt(Math.max(1e-12,
        1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg));
    return {
        av: [a, 0, 0],
        bv: [b * cg, b * sg, 0],
        cv: [c * cb, c * (ca - cb * cg) / sg, c * vol / sg]
    };
}

/**
 * Applies a cellBasis() to a fractional coordinate.
 * @param {{av: number[], bv: number[], cv: number[]}} basis
 * @param {number} f0
 * @param {number} f1
 * @param {number} f2
 * @returns {number[]} Cartesian [x, y, z], angstrom.
 */
function fracToCart(basis, f0, f1, f2) {
    return [
        basis.av[0] * f0 + basis.bv[0] * f1 + basis.cv[0] * f2,
        basis.av[1] * f0 + basis.bv[1] * f1 + basis.cv[1] * f2,
        basis.av[2] * f0 + basis.bv[2] * f1 + basis.cv[2] * f2
    ];
}

// ===========================================================================
//  3. Reciprocal metric and d-spacings
// ===========================================================================

/**
 * Coefficients of the quadratic form 1/d^2 = f(h, k, l).
 * @typedef {object} InvDsqForm
 * @property {number} aa Coefficient of h^2.
 * @property {number} bb Coefficient of k^2.
 * @property {number} cc Coefficient of l^2.
 * @property {number} ab Coefficient of h*k (already includes the factor 2).
 * @property {number} bc Coefficient of k*l (already includes the factor 2).
 * @property {number} ac Coefficient of h*l (already includes the factor 2).
 */

/**
 * Builds the reciprocal quadratic form for ANY crystal system from a full
 * cell, via the general triclinic reciprocal metric.
 *
 * A NOTE ON THE MONOCLINIC CASE, because it has been "fixed" into being wrong
 * before: this routine is not told which axis is unique and does not need to
 * be. With the two constrained angles at 90 the general form reduces to the
 * closed monoclinic expression EXACTLY (agreement to 1 part in 1e16, verified
 * numerically), for unique axis a, b or c alike. The unique axis is carried by
 * which angle differs from 90, which is the caller's job -- use normaliseCell()
 * to guarantee it.
 *
 * Returns null for a degenerate cell (non-positive volume factor, or an angle
 * at 0 or 180), which callers must treat as "no d-spacing exists".
 *
 * @param {Cell} cell
 * @returns {InvDsqForm|null}
 */
function buildInvDsq(cell) {
    const a = cell.a;
    const b = (cell.b > 0) ? cell.b : a;
    const c = (cell.c > 0) ? cell.c : a;
    if (!(a > 0) || !(b > 0) || !(c > 0)) return null;

    const al = (cell.alpha ?? 90) * DEG2RAD;
    const be = (cell.beta  ?? 90) * DEG2RAD;
    const ga = (cell.gamma ?? 90) * DEG2RAD;

    const cosA = Math.cos(al), cosB = Math.cos(be), cosG = Math.cos(ga);
    const sinA = Math.sin(al), sinB = Math.sin(be), sinG = Math.sin(ga);

    const volume_factor = 1 - cosA * cosA - cosB * cosB - cosG * cosG
                            + 2 * cosA * cosB * cosG;
    if (volume_factor <= 1e-9 ||
        Math.abs(sinA) < 1e-9 || Math.abs(sinB) < 1e-9 || Math.abs(sinG) < 1e-9) {
        return null;
    }

    const V = a * b * c * Math.sqrt(volume_factor);

    const a_star = b * c * sinA / V;
    const b_star = a * c * sinB / V;
    const c_star = a * b * sinG / V;

    const cosA_star = (cosB * cosG - cosA) / (sinB * sinG);
    const cosB_star = (cosA * cosG - cosB) / (sinA * sinG);
    const cosG_star = (cosA * cosB - cosG) / (sinA * sinB);

    return {
        aa: a_star * a_star,
        bb: b_star * b_star,
        cc: c_star * c_star,
        ab: 2 * a_star * b_star * cosG_star,
        bc: 2 * b_star * c_star * cosA_star,
        ac: 2 * a_star * c_star * cosB_star
    };
}

/**
 * Evaluates 1/d^2 for one reflection.
 * @param {InvDsqForm} f
 * @param {number} h
 * @param {number} k
 * @param {number} l
 * @returns {number} 1/d^2, angstrom^-2.
 */
function invDsq(f, h, k, l) {
    return h * h * f.aa + k * k * f.bb + l * l * f.cc
         + h * k * f.ab + k * l * f.bc + h * l * f.ac;
}

/**
 * d-spacing for one reflection, or null if it does not exist.
 * @param {Cell} cell
 * @param {number} h
 * @param {number} k
 * @param {number} l
 * @returns {number|null} d, angstrom.
 */
function dSpacing(cell, h, k, l) {
    const f = buildInvDsq(cell);
    if (!f) return null;
    const q = invDsq(f, h, k, l);
    return (isFinite(q) && q > 1e-12) ? 1 / Math.sqrt(q) : null;
}

/**
 * An evaluator for 1/d^2, plus a flag saying whether the l loop in an HKL
 * generator may break early.
 * @typedef {object} InvDsqEvaluator
 * @property {function(number, number, number): number} f  (h,k,l) -> 1/d^2.
 * @property {boolean} lSeparable True when 1/d^2 has no h*l or k*l cross term,
 *           so 1/d^2 is monotonic in |l| at fixed h,k and a generator can stop
 *           at the first l that falls outside the 2-theta range.
 */

/**
 * Builds a closure evaluating 1/d^2 for a whole crystal system.
 *
 * This exists alongside buildInvDsq because an HKL generator wants one
 * hoisted closure and a monotonicity flag, whereas updateHklPositions wants
 * the raw coefficients. Both are built from the SAME metric here, which is
 * the point: a copy of this used to live in powder5.html carrying a comment
 * promising it used "exactly the same metric as updateHklPositions() so that
 * generation and the downstream 2-theta filter can never disagree". That
 * promise was maintained by hand across two files. Now it is structural.
 *
 * Returns null for a degenerate cell or an unknown system.
 *
 * @param {string} system Crystal system, lower case.
 * @param {object} params Lattice parameters; a is required.
 * @returns {InvDsqEvaluator|null}
 */
function buildInvDsqEvaluator(system, params) {
    const { a, b, c } = params;
    if (!(a > 0)) return null;
    const a_sq = a * a;

    switch (system) {
        case 'cubic':
            return { f: (h, k, l) => (h * h + k * k + l * l) / a_sq, lSeparable: true };

        case 'tetragonal': {
            const c_sq = (c > 0) ? c * c : a_sq;
            return { f: (h, k, l) => (h * h + k * k) / a_sq + (l * l) / c_sq, lSeparable: true };
        }

        case 'orthorhombic': {
            const b_sq = (b > 0) ? b * b : a_sq;
            const c_sq = (c > 0) ? c * c : a_sq;
            return {
                f: (h, k, l) => (h * h) / a_sq + (k * k) / b_sq + (l * l) / c_sq,
                lSeparable: true
            };
        }

        case 'hexagonal':
        case 'trigonal':
        case 'rhombohedral': {
            const c_sq = (c > 0) ? c * c : a_sq;
            return {
                f: (h, k, l) => 4 * (h * h + h * k + k * k) / (3 * a_sq) + (l * l) / c_sq,
                lSeparable: true
            };
        }

        // Monoclinic is not special-cased: the general metric reduces to the
        // monoclinic closed form exactly once the two fixed angles are 90, and
        // it does so for unique axis a, b and c alike.
        case 'monoclinic':
        case 'triclinic': {
            const form = buildInvDsq({
                a,
                b: (b > 0) ? b : a,
                c: (c > 0) ? c : a,
                alpha: params.alpha,
                beta: params.beta,
                gamma: params.gamma
            });
            if (!form) return null;
            const { aa, bb, cc, ab, bc, ac } = form;
            return {
                f: (h, k, l) => h * h * aa + k * k * bb + l * l * cc
                              + h * k * ab + k * l * bc + h * l * ac,
                lSeparable: false
            };
        }

        default:
            return null;
    }
}

// ===========================================================================
//  4. Peak positions
// ===========================================================================

/**
 * Recomputes `tth` (degrees) and `d` (angstrom) in place for every reflection
 * in a list, from the current cell parameters.
 *
 * Reflections whose indices are missing, whose cell is degenerate, or whose
 * Bragg condition cannot be satisfied at this wavelength get tth = d = null.
 * Callers everywhere treat null as "skip this reflection".
 *
 * The closed-form branches for the high-symmetry systems are kept rather than
 * routed through buildInvDsq because this runs once per LM iteration over
 * every reflection, and they avoid ~20 trig calls per call.
 *
 * @param {Array<object|null>} hklList Reflections carrying h_orig/k_orig/l_orig.
 * @param {object} params Must carry a, lambda; b/c/angles as the system needs.
 * @param {string} system Crystal system, lower case.
 * @returns {void}
 */
function updateHklPositions(hklList, params, system) {
    if (!hklList || hklList.length === 0) return;
    const clearAll = () => hklList.forEach(peak => {
        if (peak) { peak.tth = null; peak.d = null; }
    });

    if (!params) { clearAll(); return; }
    const { a, b, c, lambda } = params;
    if (!lambda || lambda <= 0 || !a || a <= 0) { clearAll(); return; }

    const lambda_sq_over_4 = (lambda * lambda) / 4.0;
    const a_sq = a * a;

    let b_sq, c_sq;
    /** @type {InvDsqForm|null} */
    let form = null;

    if (system === 'triclinic' || system === 'monoclinic') {
        form = buildInvDsq({
            a,
            b: (b > 0) ? b : a,
            c: (c > 0) ? c : a,
            alpha: params.alpha,
            beta: params.beta,
            gamma: params.gamma
        });
        if (!form) { clearAll(); return; }
    } else if (system === 'tetragonal' || system === 'hexagonal'
            || system === 'rhombohedral' || system === 'trigonal') {
        c_sq = (c && c > 0) ? (c * c) : a_sq;
    } else if (system === 'orthorhombic') {
        b_sq = (b && b > 0) ? (b * b) : a_sq;
        c_sq = (c && c > 0) ? (c * c) : a_sq;
    }

    hklList.forEach(peak => {
        if (!peak || peak.h_orig === undefined
                  || peak.k_orig === undefined
                  || peak.l_orig === undefined) {
            if (peak) { peak.tth = null; peak.d = null; }
            return;
        }
        const h = peak.h_orig, k = peak.k_orig, l = peak.l_orig;
        const h2 = h * h, k2 = k * k, l2 = l * l;

        let inv_d_sq = 0;
        try {
            switch (system) {
                case 'cubic':
                    inv_d_sq = (h2 + k2 + l2) / a_sq;
                    break;
                case 'tetragonal':
                    inv_d_sq = (h2 + k2) / a_sq + l2 / c_sq;
                    break;
                case 'orthorhombic':
                    inv_d_sq = h2 / a_sq + k2 / b_sq + l2 / c_sq;
                    break;
                case 'hexagonal':
                case 'rhombohedral':
                case 'trigonal':
                    if (a_sq <= 0 || c_sq <= 0) throw new Error('Invalid lattice param');
                    inv_d_sq = 4 * (h2 + h * k + k2) / (3 * a_sq) + l2 / c_sq;
                    break;
                case 'monoclinic':
                case 'triclinic':
                    if (!form) throw new Error('Invalid lattice param');
                    inv_d_sq = invDsq(form, h, k, l);
                    break;
                default:
                    throw new Error(`Unknown system: ${system}`);
            }

            if (!isFinite(inv_d_sq) || inv_d_sq <= 1e-12) {
                peak.tth = null;
                peak.d = null;
            } else {
                const sinThetaSq = lambda_sq_over_4 * inv_d_sq;
                if (sinThetaSq <= 1 && sinThetaSq > 0) {
                    const thetaRad = Math.asin(Math.sqrt(sinThetaSq));
                    peak.tth = 2 * thetaRad / DEG2RAD;
                    peak.d = 1 / Math.sqrt(inv_d_sq);
                } else {
                    peak.tth = null;
                    peak.d = null;
                }
            }
        } catch (error) {
            peak.tth = null;
            peak.d = null;
        }
    });
}

// ---------------------------------------------------------------------------
//  LORENTZ-POLARISATION
//
//  ONE implementation, here, because crystal.js is the only module loaded by
//  all three entry points - powder5.html, charge_flipping_worker.js and
//  refinement_worker.js - and because the three copies that preceded it did
//  not agree.
//
//  The split that caused the trouble is between the POLARISATION RATIO K,
//  which is a number, and a polarisation DESCRIPTOR, which is an object (or,
//  for backward compatibility, a bare monochromator 2-theta). Two functions
//  called lorentzPolarization took a descriptor; observations.js called a
//  third by a different name and passed it K. Aliasing one to the other looks
//  like the obvious repair and is a trap:
//
//      lorentzPolarization(tth, 0)
//
//  reads the 0 as a monochromator angle, fails the "> 0" test and falls back
//  to K = 1. But K = 0 is exactly what a NEUTRON pattern needs - an
//  unpolarised beam with no polarisation factor at all - so a neutron dataset
//  would be corrected as though it came off an X-ray tube. The error reaches
//  50% at 2-theta = 90 degrees and varies smoothly with angle, which in a
//  difference plot is nearly indistinguishable from a temperature factor.
//
//  So the two conventions are kept visibly apart: sharkoPolarizationK()
//  resolves a descriptor to K, once, and everything downstream takes K.
// ---------------------------------------------------------------------------

/** Default polarisation model, matching the UI controls. */
const SHARKO_POLARIZATION_DEFAULTS = Object.freeze({ mode: 'lab', monoTth: 0, polFrac: 0.95 });

/**
 * Polarisation descriptor -> the ratio K.
 *
 * Accepts, in order of precedence:
 *   { K }                        already resolved, used verbatim
 *   { mode, monoTth, polFrac }   'none' -> 0, 'synchrotron' -> (1-f)/f,
 *                                'lab'  -> cos^2(2theta_M), or 1 unmonochromated
 *   a bare number                legacy monochromator 2-theta, in degrees
 *   absent                       1, the unpolarised tube
 *
 * A BARE NUMBER IS A 2-THETA, NEVER A K. That is the historical meaning and
 * changing it silently would rescale every intensity in any caller still using
 * the old two-argument form. Pass { K: value } to supply a ratio.
 */
function sharkoPolarizationK(pol) {
    if (typeof pol === 'number') {
        return (Number.isFinite(pol) && pol > 0 && pol < 180)
            ? Math.cos(pol * DEG2RAD) ** 2 : 1;
    }
    if (pol && typeof pol === 'object') {
        if (Number.isFinite(pol.K) && pol.K >= 0) return pol.K;
        const mode = String(pol.mode || SHARKO_POLARIZATION_DEFAULTS.mode).toLowerCase();
        if (mode === 'none') return 0;
        if (mode === 'synchrotron') {
            const f = Number(pol.polFrac);
            if (!Number.isFinite(f) || f >= 1) return 0;
            // f -> 0 is a beam polarised IN the diffraction plane; K diverges
            // and P collapses to cos^2(2theta). Capped so nothing overflows.
            return f <= 1e-6 ? 1e6 : (1 - f) / f;
        }
        const m = Number(pol.monoTth);
        return (Number.isFinite(m) && m > 0 && m < 180)
            ? Math.cos(m * DEG2RAD) ** 2 : 1;
    }
    return 1;
}

/**
 * Lorentz factor for Bragg-Brentano powder geometry: 1 / (sin^2(th) cos(th)).
 * NaN where the geometry is degenerate, so callers can decide what that means.
 */
function sharkoLorentzFactor(tthDeg) {
    if (!Number.isFinite(tthDeg) || tthDeg <= 0 || tthDeg >= 180) return NaN;
    const th = 0.5 * tthDeg * DEG2RAD;
    const s = Math.sin(th), c = Math.cos(th);
    if (s < 1e-6 || c < 1e-6) return NaN;
    return 1 / (s * s * c);
}

/** Polarisation factor from K, normalised so P(0) = 1. */
function sharkoPolarizationFactor(tthDeg, K) {
    if (!Number.isFinite(tthDeg)) return NaN;
    const k = (Number.isFinite(K) && K >= 0) ? K : 1;
    return (1 + k * Math.cos(tthDeg * DEG2RAD) ** 2) / (1 + k);
}

/**
 * Lp from 2-theta and the RATIO K. The canonical form; everything else wraps it.
 *
 * Returns 1 where the geometry is degenerate rather than NaN, so a bad 2-theta
 * leaves an intensity alone instead of turning it into NaN and silently
 * dropping the reflection from a normalisation.
 */
function sharkoLorentzPolarization(tthDeg, K) {
    const L = sharkoLorentzFactor(tthDeg);
    if (!Number.isFinite(L)) return 1;
    const P = sharkoPolarizationFactor(tthDeg, K);
    return Number.isFinite(P) ? L * P : 1;
}

/** Lp from 2-theta and a DESCRIPTOR. Convenience wrapper; resolves K first. */
function sharkoLp(tthDeg, pol) {
    return sharkoLorentzPolarization(tthDeg, sharkoPolarizationK(pol));
}

// ---------------------------------------------------------------------------
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        DEG2RAD,
        MONOCLINIC_ANGLE_FOR_AXIS, MONOCLINIC_FIXED_ANGLES_FOR_AXIS,
        normaliseCell,
        metricTensor, fracDistance, cellVolume, cellBasis, fracToCart,
        buildInvDsq, invDsq, buildInvDsqEvaluator, dSpacing, updateHklPositions,
        SHARKO_POLARIZATION_DEFAULTS, sharkoPolarizationK, sharkoLorentzFactor,
        sharkoPolarizationFactor, sharkoLorentzPolarization, sharkoLp
    };
}
