// profile.js
// ---------------------------------------------------------------------------
//  SINGLE SOURCE OF TRUTH for the peak-shape model: peak shift, profile
//  widths, pseudo-Voigt preparation and evaluation.
//
//  Why this file exists
//  --------------------
//  This code used to exist twice, in powder5.html and in
//  refinement_worker.js, ~200 lines each. The two copies had ALREADY DIVERGED
//  by the time they were merged here:
//
//    * the worker had a reciprocal-precompute optimisation in evalVoigt
//      (delta * inv_H_G instead of delta / H_G) that the main thread never
//      received;
//    * the main thread had getPseudoVoigtArea, which the worker never had;
//    * the two had drifted apart in comments and in the arrangement of the
//      getPeakFWHM polynomial.
//
//  None of that changed the numbers -- this time. The point is that nothing in
//  the program could have told you if it had: the preview would have drawn one
//  peak shape and the refinement fitted another, and the report would have
//  been perfectly self-consistent about the wrong answer.
//
//  This file is the merge: the worker's evaluation path (it is the faster one
//  and it is the one the fit actually used) plus the main thread's
//  getPseudoVoigtArea.
//
//  Loading
//  -------
//    browser :  script src="profile.js"
//                 (written without angle brackets on purpose: a literal
//                  closing script tag inside this file would terminate the
//                  page early if anyone ever inlined it into powder5.html)   (after constants.js)
//    worker  :  importScripts('profile.js')
//    node    :  require('./profile.js')
//
//  Depends on constants.js for MIN_PROFILE_FWHM_DEG, softPositive,
//  pvMixCoeffs, the PV_* areas, the GSAS_* conversions, the STEPHENS_*
//  constants and the CALCULATION_WINDOW_* / PROFILE_WINDOW_MAX_DEG window.
//
//  All width helpers return widths in DEGREES 2-theta.
// ---------------------------------------------------------------------------

'use strict';

/**
 * A refinement parameter set. Only the fields each profile type consumes need
 * to be present; every read below is defensive (`params.GU || 0`).
 * @typedef {object} ProfileParams
 * @property {('simple_pvoigt'|'split_pvoigt'|'tch_aniso')} profileType
 * @property {number} [GU] @property {number} [GV] @property {number} [GW]
 * @property {number} [GP] @property {number} [LX] @property {number} [eta]
 * @property {number} [shft]  Sample displacement, micrometre-ish GSAS units.
 * @property {number} [trns]  Sample transparency.
 * @property {number} [GU_L] @property {number} [GV_L] @property {number} [GW_L]
 * @property {number} [LX_L] @property {number} [GU_R] @property {number} [GV_R]
 * @property {number} [GW_R] @property {number} [LX_R]
 * @property {number} [eta_split] @property {number} [shft_split]
 * @property {number} [trns_split]
 * @property {number} [U] @property {number} [V] @property {number} [W]
 * @property {number} [X] @property {number} [Y]
 * @property {number} [S400] @property {number} [S040] @property {number} [S004]
 * @property {number} [S220] @property {number} [S202] @property {number} [S022]
 * @property {number} [eta_aniso] Stephens G/L split; STEPHENS_DEFAULT_ETA if absent.
 * @property {number} [SL] Axial divergence, low side.
 * @property {number} [HL] Axial divergence, high side.
 */

/**
 * One reflection. h/k/l are the indices used for the anisotropic-strain sum;
 * `d` is the d-spacing in angstrom.
 * @typedef {object} Reflection
 * @property {number} [h] @property {number} [k] @property {number} [l]
 * @property {number} [d]
 */

/**
 * Gaussian and Lorentzian FWHM, both in DEGREES 2-theta.
 * @typedef {{gamma_G: number, gamma_L: number}} Widths
 */

/**
 * Everything evalVoigt needs for one peak, precomputed once per peak per
 * pattern evaluation. `type` selects the branch: 1 = simple, 2 = split,
 * 3 = TCH. Reciprocals (inv_*) are precomputed because the innermost loop runs
 * once per data point per reflection.
 * @typedef {object} VoigtPrep
 * @property {1|2|3} type
 * @property {number} x0          Peak centre, degrees 2-theta.
 * @property {number} asym_param  Axial-divergence asymmetry.
 * @property {number} eta         Lorentzian AREA fraction.
 * @property {number} max_d       Truncation half-width, degrees.
 * @property {number} pedestal    Profile value at max_d, subtracted off.
 * @property {number} fwhm_total  Convoluted FWHM, degrees.
 * @property {number} analyticArea Exact area of the UNTRUNCATED shape.
 * @property {number} Cg          4*ln2, the Gaussian exponent coefficient.
 */

/**
 * Peak position shift from sample displacement and sample transparency.
 *
 * @param {number} tth 2-theta, degrees.
 * @param {ProfileParams} params
 * @returns {number} Shift to ADD to the Bragg position, degrees 2-theta.
 *          Always 0 for tch_aniso, which models no shift.
 */
function calculatePeakShift(tth, params) {
     if (!params || !params.profileType) return 0;
    const profileType = params.profileType;
    
    const calcShift = (tth, shft, trns) => {
        const thetaRad = tth * (Math.PI / 180) / 2;
        if (Math.abs(thetaRad - Math.PI / 2.0) < 1e-6) return 0;
        const cosTheta = Math.cos(thetaRad);
        const sin2Theta = Math.sin(2 * thetaRad);
        const displacementShift = -(shft / 1000) * cosTheta * (180 / Math.PI);
        const transparencyShift = trns * sin2Theta * (180 / Math.PI);
        const totalShift = displacementShift + transparencyShift;
        return isFinite(totalShift) ? totalShift : 0;
    };

    switch (profileType) {
        case "simple_pvoigt":
            return calcShift(tth, params.shft || 0, params.trns || 0);
        case "split_pvoigt":
            return calcShift(tth, params.shft_split || 0, params.trns_split || 0);
        case "tch_aniso":
        default:
            return 0; // No shift for TCH profile
    }
}

// The GSAS / Stephens unit conventions and constants
// (GSAS_GAUSSIAN_TO_DEG, GSAS_LORENTZIAN_TO_DEG, STEPHENS_PREFACTOR,
//  STEPHENS_DEFAULT_ETA) now live in constants.js, together with the full
// explanation of the conventions and of every numeric value.

/**
 * Caglioti sigma^2 (GSAS centidegree^2 / 8ln2) -> Gaussian FWHM in degrees.
 * Component floors are a QUARTER of the minimum resolvable FWHM, so neither
 * component alone pins the total; the combined floor in calculateProfileWidths
 * does that, and it preserves the Gaussian/Lorentzian ratio (hence eta).
 *
 * @param {number} sigma_sq_cdeg2 Caglioti sigma^2 in GSAS units.
 * @returns {number} Gaussian FWHM, degrees 2-theta.
 */
function gaussianFromSigmaSq(sigma_sq_cdeg2) {
    const floorSq = Math.pow(0.25 * MIN_PROFILE_FWHM_DEG / GSAS_GAUSSIAN_TO_DEG, 2);
    return Math.sqrt(softPositive(sigma_sq_cdeg2, floorSq)) * GSAS_GAUSSIAN_TO_DEG;
}
/**
 * Lorentzian gamma (GSAS centidegrees) -> FWHM in degrees.
 * @param {number} gL_cdeg
 * @returns {number} Lorentzian FWHM, degrees 2-theta.
 */
function lorentzianFromGamma(gL_cdeg) {
    const floor = 0.25 * MIN_PROFILE_FWHM_DEG / GSAS_LORENTZIAN_TO_DEG;
    return softPositive(gL_cdeg, floor) * GSAS_LORENTZIAN_TO_DEG;
}

/**
 * Smooth lower bound at `f`: ~f for x <= f, within 0.6% of x by x = 2.5f, with a
 * nonzero derivative throughout.
 *
 * The exponent is a deliberate compromise. sqrt(x^2 + f^2) is the smoothest
 * choice but inflates a legitimate 2.5f width by 8%, which is not acceptable in
 * a quantity being refined; a very sharp knee keeps widths honest but makes the
 * derivative BELOW the floor so small that a parameter driven deep past it is
 * effectively frozen again -- the exact failure this floor exists to avoid.
 * Fourth order keeps the inflation under 1% at 2.5f while leaving a usable
 * gradient pointing back out.
 *
 * @param {number} x Value to floor.
 * @param {number} f Floor, strictly positive.
 * @returns {number}
 */
function softFloor(x, f) {
    if (!isFinite(x) || !(x > 0)) return f;
    if (!(f > 0)) return x;
    const r = x / f;
    if (r > 20) return x;                       // floor irrelevant; avoid pow overflow
    return f * Math.pow(Math.pow(r, 4) + 1, 0.25);
}

/**
 * Floors the profile widths at the resolvable FWHM.
 *
 * FIX: the floor used to be applied to GW and LX individually, via hard bounds
 * in getParameterMapping. But sigma^2 = GU tan^2(th) + GV tan(th) + GW +
 * GP/cos^2(th), so bounding GW alone is simultaneously too weak -- GV < 0 can
 * still collapse the peak at low angle -- and too strong, since it ignores what
 * the other terms already contribute. And a HARD bound is what killed the
 * Jacobian column: once a width sat on its floor the finite difference returned
 * exactly zero and the parameter was frozen for the rest of the fit.
 *
 * Both COMPONENTS are floored, not just the convoluted total: simple_pvoigt and
 * split_pvoigt give the Gaussian and the Lorentzian independent widths, so one
 * of them alone can be an unresolvable sub-step spike while the combined FWHM
 * still looks healthy -- and with the area-correct eta mixing that spike then
 * carries its full share of the peak area. Flooring the components implies the
 * floor on the total, so no second pass is needed (which also keeps the five
 * Math.pow calls in getPeakFWHM out of this hot path).
 *
 * @param {number} gamma_G Gaussian FWHM, degrees.
 * @param {number} gamma_L Lorentzian FWHM, degrees.
 * @returns {Widths}
 */
function applyFwhmFloor(gamma_G, gamma_L) {
    return {
        gamma_G: softFloor(gamma_G, MIN_PROFILE_FWHM_DEG),
        gamma_L: softFloor(gamma_L, MIN_PROFILE_FWHM_DEG)
    };
}

/**
 * Calculates widths for simple_pvoigt (GSAS units → degrees 2θ).
 */
/**
 * @param {number} tth   2-theta, degrees.
 * @param {Reflection} hkl
 * @param {ProfileParams} params
 * @param {number} safeThetaRad     theta in radians, clamped below pi/2.
 * @param {number} tanTheta
 * @param {number} cosTheta_safe    cos(theta), floored at 1e-9.
 * @param {number} cosThetaSq_safe  cos^2(theta), floored at 1e-9.
 * @returns {Widths} Unfloored widths for simple_pvoigt.
 */
function _calculateWidths_Simple(tth, hkl, params, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe) {
    let gamma_G = MIN_PROFILE_FWHM_DEG;
    let gamma_L = MIN_PROFILE_FWHM_DEG;

    const GU = params.GU || 0;
    const GV = params.GV || 0;
    const GW = params.GW || 0;
    const GP = params.GP || 0;
    const LX = params.LX || 0;

    //   Caglioti in cdeg²/(8 ln2). Convert √σ² to FWHM in degrees.
    const sigma_sq_cdeg2 = GU * tanTheta * tanTheta + GV * tanTheta + GW + GP / cosThetaSq_safe;
    gamma_G = gaussianFromSigmaSq(sigma_sq_cdeg2);
    //   Lorentzian (size only) in centidegrees → degrees.
    const gL_cdeg = LX / cosTheta_safe;
    gamma_L = lorentzianFromGamma(gL_cdeg);

    return { gamma_G, gamma_L };
}

/**
 * Calculates widths for split_pvoigt (GSAS units → degrees 2θ).
 */
/**
 * @param {number} tth   2-theta, degrees.
 * @param {Reflection} hkl
 * @param {ProfileParams} params
 * @param {('left'|'right'|'center')} side Which flank's parameters to use.
 * @param {number} safeThetaRad
 * @param {number} tanTheta
 * @param {number} cosTheta_safe
 * @param {number} cosThetaSq_safe
 * @returns {Widths} Unfloored widths for split_pvoigt.
 */
function _calculateWidths_Split(tth, hkl, params, side, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe) {
    let gamma_G = MIN_PROFILE_FWHM_DEG;
    let gamma_L = MIN_PROFILE_FWHM_DEG;
    let GU, GV, GW, LX;

    if (side === 'left') {
        GU = params.GU_L || 0;
        GV = params.GV_L || 0;
        GW = params.GW_L || 0;
        LX = params.LX_L || 0;
    } else { // 'right' or 'center'
        GU = params.GU_R || 0;
        GV = params.GV_R || 0;
        GW = params.GW_R || 0;
        LX = params.LX_R || 0;
    }

    const sigma_sq_cdeg2 = GU * tanTheta * tanTheta + GV * tanTheta + GW;  //  no GP for split
    gamma_G = gaussianFromSigmaSq(sigma_sq_cdeg2);
    const gL_cdeg = LX / cosTheta_safe;
    gamma_L = lorentzianFromGamma(gL_cdeg);

    return { gamma_G, gamma_L };
}

/**
 * Calculates widths for tch_aniso (GSAS units → degrees 2θ).
 *   U,V,W   are GSAS GU,GV,GW.
 *   X       is GSAS LY  (strain · tanθ).
 *   Y       is GSAS LX  (size   / cosθ).
 *   Stephens contribution is added in canonical TOPAS / GSAS-II form.
 */
/**
 * @param {number} tth   2-theta, degrees.
 * @param {Reflection} hkl
 * @param {ProfileParams} params
 * @param {number} safeThetaRad     theta in radians, clamped below pi/2.
 * @param {number} tanTheta
 * @param {number} cosTheta_safe    cos(theta), floored at 1e-9.
 * @param {number} cosThetaSq_safe  cos^2(theta), floored at 1e-9.
 * @returns {Widths} Unfloored widths for tch_aniso.
 */
function _calculateWidths_TCH(tth, hkl, params, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe) {
    let gamma_G = MIN_PROFILE_FWHM_DEG;
    let gamma_L = MIN_PROFILE_FWHM_DEG;

    const U = params.U || 0;
    const V = params.V || 0;
    const W = params.W || 0;
    const X = params.X || 0;     //  GSAS LY  (strain · tanθ)
    const Y = params.Y || 0;     //  GSAS LX  (size   / cosθ)

    const sigma_sq_cdeg2 = U * tanTheta * tanTheta + V * tanTheta + W;
    gamma_G = gaussianFromSigmaSq(sigma_sq_cdeg2);
    const gL_cdeg = X * tanTheta + Y / cosTheta_safe;
    gamma_L = lorentzianFromGamma(gL_cdeg);

    //   Stephens anisotropic strain (orthorhombic-and-up form).
    //   Higher-symmetry systems collapse via the symmetry constraints
    //   applied on the main thread before parameters arrive here.
    if (hkl && hkl.d && hkl.h_orig !== undefined) {
        const d_sq = hkl.d * hkl.d;
        if (d_sq > 1e-9) {
            const h_val = hkl.h_orig, k_val = hkl.k_orig, l_val = hkl.l_orig;
            const h2 = h_val*h_val, k2 = k_val*k_val, l2 = l_val*l_val;
            const h4 = h2*h2,        k4 = k2*k2,        l4 = l2*l2;

            const S400 = params.S400 || 0, S040 = params.S040 || 0, S004 = params.S004 || 0;
            const S220 = params.S220 || 0, S202 = params.S202 || 0, S022 = params.S022 || 0;


            let M_HKL = 0;
            if (params.system === 'hexagonal' || params.system === 'trigonal' || params.system === 'rhombohedral') {
                const hk_term = h2 + h_val * k_val + k2;
                M_HKL = S400 * hk_term * hk_term + S004 * l4 + S202 * hk_term * l2;
            } else {
                M_HKL = S400*h4 + S040*k4 + S004*l4 + S220*h2*k2 + S202*h2*l2 + S022*k2*l2;
            }

            if (M_HKL > 0 && isFinite(M_HKL)) {
                const pp = d_sq * Math.sqrt(M_HKL) / 1000.0;
                const aniso_total = STEPHENS_PREFACTOR * pp * tanTheta;
                const eta_a = (typeof params.eta_aniso === 'number')
                              ? Math.max(0, Math.min(1, params.eta_aniso))
                              : STEPHENS_DEFAULT_ETA;
                const aniso_L = aniso_total * eta_a;
                const aniso_G = aniso_total * (1 - eta_a);
                if (isFinite(aniso_L) && aniso_L > 0) gamma_L += aniso_L;
                if (isFinite(aniso_G) && aniso_G > 0) {
                    gamma_G = Math.sqrt(gamma_G*gamma_G + aniso_G*aniso_G);
                }
            }
        }
    }

    return { gamma_G, gamma_L };
}


/**
 * Calculates Gaussian and Lorentzian width components (gamma_G, gamma_L)
 * This is now a DISPATCHER function.
 */
/**
 * Gaussian and Lorentzian FWHM for one reflection, floored at the resolvable
 * width. The single entry point for widths -- nothing else should call the
 * per-type helpers.
 *
 * @param {number} tth 2-theta, degrees.
 * @param {Reflection} hkl
 * @param {ProfileParams} params
 * @param {('left'|'right'|'center')} [side='center'] Only meaningful for
 *        split_pvoigt.
 * @returns {Widths} Both in DEGREES 2-theta, both >= MIN_PROFILE_FWHM_DEG.
 */
function calculateProfileWidths(tth, hkl, params, side = 'center') {
    if (!params || !params.profileType) {
        return { gamma_G: MIN_PROFILE_FWHM_DEG, gamma_L: MIN_PROFILE_FWHM_DEG };
    }
    
    const profileType = params.profileType;
    const thetaRad = tth * (Math.PI / 180) / 2;

    const MAX_ANGLE_RAD = Math.PI / 2.0 - 1e-6;
    const safeThetaRad = Math.min(thetaRad, MAX_ANGLE_RAD);
     if (safeThetaRad < 1e-6) {
         return { gamma_G: MIN_PROFILE_FWHM_DEG, gamma_L: MIN_PROFILE_FWHM_DEG };
     }

    const tanTheta = Math.tan(safeThetaRad);
    const cosTheta = Math.cos(safeThetaRad);
    const cosTheta_safe = Math.max(cosTheta, 1e-9);
    const cosThetaSq_safe = Math.max(cosTheta * cosTheta, 1e-9);

    let widths = { gamma_G: MIN_PROFILE_FWHM_DEG, gamma_L: MIN_PROFILE_FWHM_DEG };

    //   Dispatch to specific width calculators  
    switch (profileType) {
        case "simple_pvoigt":
            widths = _calculateWidths_Simple(tth, hkl, params, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe);
            break;
        case "split_pvoigt":
            widths = _calculateWidths_Split(tth, hkl, params, side, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe);
            break;
        case "tch_aniso":
            widths = _calculateWidths_TCH(tth, hkl, params, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe);
            break;
    }

    return applyFwhmFloor(widths.gamma_G, widths.gamma_L);
}


/**
 * Calculates the total FWHM of a TCH/Split peak from its Gaussian and Lorentzian components.
 */
/**
 * Convoluted pseudo-Voigt FWHM from the two component widths, via the
 * Thompson-Cox-Hastings (1987) fifth-order polynomial approximation to the
 * Voigt width. Accurate to ~0.02% over the full Gaussian-to-Lorentzian range.
 *
 * @param {number} gamma_G Gaussian FWHM, degrees.
 * @param {number} gamma_L Lorentzian FWHM, degrees.
 * @returns {number} Total FWHM, degrees, floored at 1e-6.
 */
function getPeakFWHM(gamma_G, gamma_L) {
    const gG = Math.max(1e-9, gamma_G || 1e-9);
    const gL = Math.max(1e-9, gamma_L || 1e-9);
    const fwhm_g_5 = Math.pow(gG, 5);
    const fwhm_l_5 = Math.pow(gL, 5);
    const fwhm_g_4_l = 2.69269 * Math.pow(gG, 4) * gL;
    const fwhm_g_3_l_2 = 2.42843 * Math.pow(gG, 3) * Math.pow(gL, 2);
    const fwhm_g_2_l_3 = 4.47163 * Math.pow(gG, 2) * Math.pow(gL, 3);
    const fwhm_g_l_4 = 0.07842 * gG * Math.pow(gL, 4);
    const fwhm_pow5 = fwhm_g_5 + fwhm_g_4_l + fwhm_g_3_l_2 + fwhm_g_2_l_3 + fwhm_g_l_4 + fwhm_l_5;
     if (fwhm_pow5 < 0 || !isFinite(fwhm_pow5)) return Math.max(gG, gL, 1e-6);
     const fwhm = Math.pow(fwhm_pow5, 0.2);
     return Math.max(1e-6, fwhm);
}

/**
 * Precomputes everything evalVoigt needs for one peak. Called once per peak
 * per pattern evaluation; evalVoigt is then called once per data point in the
 * peak's window, so anything that can be hoisted belongs here.
 *
 * @param {number} tth_peak Angle at which to EVALUATE the widths (the Bragg
 *        position plus zero-shift, before the displacement shift).
 * @param {number} x0 Peak CENTRE, degrees 2-theta (after every shift).
 * @param {Reflection} hkl
 * @param {ProfileParams} params
 * @returns {VoigtPrep}
 */
function prepareVoigt(tth_peak, x0, hkl, params) {
    if (!params) return null;
    const profileType = params.profileType;
    const Cg = 2.772588722239781; // 4 * ln(2)
    
    // Precompute asymmetry coefficient once per peak (for TCH)
    let asym_param = 0;
    if (profileType === "tch_aniso" && (params.SL || params.HL) && tth_peak >= 0.1 && tth_peak < 180) {
        const theta_rad = tth_peak * (Math.PI / 180) / 2.0;
        const safe_theta_rad = Math.max(1e-6, Math.min(Math.PI / 2.0 - 1e-6, theta_rad));
        const tan_theta = Math.tan(safe_theta_rad);
        if (Math.abs(tan_theta) >= 1e-9) {
            asym_param = (params.SL || 0) / tan_theta + (params.HL || 0);
        }
    }

    let prep;
    let fwhm_total;

    if (profileType === "simple_pvoigt") {
        const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
        const H_G = Math.max(1e-9, gamma_G), H_L = Math.max(1e-9, gamma_L);
        const eta = Math.max(0, Math.min(1, params.eta || 0.5));
        fwhm_total = getPeakFWHM(H_G, H_L);
        // eta is the AREA mixing fraction; pvMixCoeffs converts it to weights
        // for the unit-height components evalVoigt actually evaluates.
        const mix1 = pvMixCoeffs(eta, H_G, H_L);

prep = { type: 1, x0, asym_param, H_G, H_L, inv_H_G: 1/H_G, inv_H_L: 1/H_L, eta,
                 cL: mix1.cL, cG: mix1.cG, Cg };
        prep.analyticArea = mix1.cL * PV_LORENTZ_AREA * H_L
                          + mix1.cG * PV_GAUSS_AREA   * H_G;
    } else if (profileType === "split_pvoigt") {

        const wL = calculateProfileWidths(tth_peak, hkl, params, 'left');
        const wR = calculateProfileWidths(tth_peak, hkl, params, 'right');
        const H_G_L = Math.max(1e-9, wL.gamma_G), H_L_L = Math.max(1e-9, wL.gamma_L);
        const H_G_R = Math.max(1e-9, wR.gamma_G), H_L_R = Math.max(1e-9, wR.gamma_L);
        const eta = Math.max(0, Math.min(1, params.eta_split || 0.5));
        fwhm_total = Math.max(getPeakFWHM(H_G_L, H_L_L), getPeakFWHM(H_G_R, H_L_R));
        // Each flank gets its own mixing weights (the widths differ), and each
        // is renormalised to 1 at delta = 0, so the profile stays continuous
        // across the centre.
        const mixL = pvMixCoeffs(eta, H_G_L, H_L_L);
        const mixR = pvMixCoeffs(eta, H_G_R, H_L_R);


prep = { type: 2, x0, asym_param, H_G_L, H_L_L, H_G_R, H_L_R, 
                 inv_H_G_L: 1/H_G_L, inv_H_L_L: 1/H_L_L, inv_H_G_R: 1/H_G_R, inv_H_L_R: 1/H_L_R, 
                 eta,
                 cL_L: mixL.cL, cG_L: mixL.cG, cL_R: mixR.cL, cG_R: mixR.cG, Cg };
        prep.analyticArea = 0.5 * (mixL.cL * PV_LORENTZ_AREA * H_L_L
                                 + mixL.cG * PV_GAUSS_AREA   * H_G_L)
                          + 0.5 * (mixR.cL * PV_LORENTZ_AREA * H_L_R
                                 + mixR.cG * PV_GAUSS_AREA   * H_G_R);
    } else { // tch_aniso
        const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
        const H_G = Math.max(1e-9, gamma_G), H_L = Math.max(1e-9, gamma_L);
        const fwhm = Math.max(1e-9, getPeakFWHM(H_G, H_L));
        let eta = 0;
        if (fwhm > 1e-9) {
            const ratio = H_L / fwhm;
            eta = Math.max(0, Math.min(1, 1.36603 * ratio - 0.47719 * (ratio * ratio) + 0.11116 * Math.pow(ratio, 3)));
        }
        fwhm_total = fwhm;
        // The TCH eta polynomial above is defined for AREA-normalised
        // components, so it must be converted before use -- see pvMixCoeffs.
        const mix3 = pvMixCoeffs(eta, fwhm, fwhm);

        prep = { type: 3, x0, asym_param, fwhm, inv_fwhm: 1/fwhm, eta,
                 cL: mix3.cL, cG: mix3.cG, Cg };
        prep.analyticArea = (mix3.cL * PV_LORENTZ_AREA + mix3.cG * PV_GAUSS_AREA) * fwhm;
    }

    // FIX: `max_d` is now the SINGLE authoritative truncation radius, and it
    // matches the window that calculatePattern uses to pick its loop bounds.
    // Previously the two disagreed (10*(H_G+H_L) here vs 8*FWHM there) and a
    // separate PEAK_HEIGHT_CUTOFF chopped the profile at 0.2% of peak height,
    // which put a step of 0.002*I_peak into the model. That step is large
    // against 1/I weights and it made the one-sided finite-difference
    // Jacobian discontinuous.
    //
    // Asymmetry broadens one flank by up to 1/(1-0.95) ... in practice the
    // clamp limits it to ~2x, so widen the window when asymmetry is active.
    prep.fwhm_total = Math.max(1e-6, fwhm_total);
    // Window scaled by Lorentzian content: 8 FWHM for a pure Gaussian up to
    // CALCULATION_WINDOW_MULTIPLIER_L for a pure Lorentzian, capped in degrees
    // so a pathologically broad peak cannot make the window span the pattern.
    const etaWin = Math.max(0, Math.min(1, prep.eta || 0));
    const windowMult = CALCULATION_WINDOW_MULTIPLIER_G
                     + (CALCULATION_WINDOW_MULTIPLIER_L - CALCULATION_WINDOW_MULTIPLIER_G) * etaWin;
    const truncationRadius = Math.min(PROFILE_WINDOW_MAX_DEG,
        windowMult * Math.max(0.01, prep.fwhm_total) * (asym_param !== 0 ? 2.0 : 1.0));

    // Pedestal: the profile value at the truncation radius. Subtracting it
    // (clamped at zero) is what makes the profile reach zero continuously at
    // the window edge instead of dropping off a cliff.
    //
    // Measured with max_d temporarily Infinity: probing at exactly
    // `x0 + truncationRadius` otherwise rounds a fraction of an ulp past the
    // cutoff, evalVoigt early-returns 0, and the pedestal silently comes out
    // as 0. An asymmetric shape has two different edge values, so take the
    // larger -- then neither flank can go negative. Twice per peak, not per
    // data point.
    prep.pedestal = 0;
    prep.max_d = Infinity;
    prep.pedestal = Math.max(evalVoigt(x0 + truncationRadius, prep),
                             evalVoigt(x0 - truncationRadius, prep));
    prep.max_d = truncationRadius;
    return prep;
}

/**
 * Profile value at x, truncated at p.max_d and pedestal-subtracted so it
 * reaches zero continuously at the window edge.
 * This is the innermost hot loop of the whole program -- it runs once per data
 * point per reflection per function evaluation -- so it is deliberately kept
 * as one self-contained function with no helper call.
 */

/**
 * @param {number} x Where to evaluate, degrees 2-theta.
 * @param {VoigtPrep} p From prepareVoigt.
 * @returns {number} Unit-height profile value, 0 outside p.max_d.
 */
function evalVoigt(x, p) {
    if (!p) return 0.0;
    let delta = x - p.x0;
    if (delta > p.max_d || delta < -p.max_d) return 0.0;
    if (p.asym_param !== 0 && Math.abs(delta) >= 1e-9) {
        const t = Math.max(-0.95, Math.min(0.95, p.asym_param * delta));
        delta /= (1.0 - t);
    }

    let inv_hg, inv_hl, cL, cG;
    if (p.type === 1) {
        inv_hg = p.inv_H_G; inv_hl = p.inv_H_L; cL = p.cL; cG = p.cG;
    } else if (p.type === 2) {
        if (delta < 0) { inv_hg = p.inv_H_G_L; inv_hl = p.inv_H_L_L; cL = p.cL_L; cG = p.cG_L; }
        else           { inv_hg = p.inv_H_G_R; inv_hl = p.inv_H_L_R; cL = p.cL_R; cG = p.cG_R; }
    } else {
        inv_hg = inv_hl = p.inv_fwhm; cL = p.cL; cG = p.cG;
    }

    const dg = (delta * inv_hg) * (delta * inv_hg);
    const dl = (delta * inv_hl) * (delta * inv_hl);
    const g = Math.exp(-p.Cg * dg);
    const l = 1.0 / (1.0 + 4.0 * dl);
    
    const v = cL * l + cG * g - p.pedestal;
    return v > 0 ? v : 0.0;
}

/**
 * Integrated area under a UNIT-HEIGHT pseudo-Voigt, free of the truncation
 * bias in the drawn profile. Used to convert peak heights to integrated
 * intensities for reporting.
 *
 * @param {number} tth_peak 2-theta, degrees.
 * @param {Reflection} hkl
 * @param {ProfileParams} params
 * @returns {number} Area in degrees 2-theta per unit height; 1.0 on failure.
 */
function getPseudoVoigtArea(tth_peak, hkl, params) {
    const GAUSS_AREA_CONST = PV_GAUSS_AREA;
    const LORENTZ_AREA_CONST = PV_LORENTZ_AREA;

    if (!params || !params.profileType) return 1.0;

    // FIX: delegate to prepareVoigt so this can never drift out of the
    // eta convention evalVoigt uses. prep.analyticArea is the exact
    // area of the untruncated unit-height mixture; the switch below is
    // kept only as a fallback if prepareVoigt is unavailable.
    try {
        const prep = prepareVoigt(tth_peak, tth_peak, hkl, params);
        if (prep && isFinite(prep.analyticArea) && prep.analyticArea > 0) {
            return prep.analyticArea;
        }
    } catch (e) { /* fall through to the analytic switch */ }

    const profileType = params.profileType;

    switch (profileType) {
        case "simple_pvoigt": {
            const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
            const gG = Math.max(1e-9, gamma_G);
            const gL = Math.max(1e-9, gamma_L);
            const area_G = gG * GAUSS_AREA_CONST;
            const area_L = gL * LORENTZ_AREA_CONST;
            const currentEta = Math.max(0, Math.min(1, params.eta || 0.5));
            const totalArea = currentEta * area_L + (1 - currentEta) * area_G;
            return isFinite(totalArea) && totalArea > 0 ? totalArea : 1.0;
        }
        
        case "split_pvoigt": {
            // We must approximate the area of the split function.
            // We'll average the areas calculated from the left and right widths.
            const { gamma_G: gG_L, gamma_L: gL_L } = calculateProfileWidths(tth_peak, hkl, params, 'left');
            const { gamma_G: gG_R, gamma_L: gL_R } = calculateProfileWidths(tth_peak, hkl, params, 'right');
            
            const currentEta = Math.max(0, Math.min(1, params.eta_split || 0.5));
            
            const area_G_L = Math.max(1e-9, gG_L) * GAUSS_AREA_CONST;
            const area_L_L = Math.max(1e-9, gL_L) * LORENTZ_AREA_CONST;
            const totalArea_L = currentEta * area_L_L + (1 - currentEta) * area_G_L;

            const area_G_R = Math.max(1e-9, gG_R) * GAUSS_AREA_CONST;
            const area_L_R = Math.max(1e-9, gL_R) * LORENTZ_AREA_CONST;
            const totalArea_R = currentEta * area_L_R + (1 - currentEta) * area_G_R;
            
            const totalArea = (totalArea_L + totalArea_R) / 2.0;
            return isFinite(totalArea) && totalArea > 0 ? totalArea : 1.0;
        }

        case "tch_aniso": {
            const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
            const gG = Math.max(1e-9, gamma_G);
            const gL = Math.max(1e-9, gamma_L);
            
            const fwhm = getPeakFWHM(gG, gL);
            const ratio = (fwhm > 1e-9) ? gL / fwhm : 0;
            const eta_calc = 1.36603 * ratio - 0.47719 * (ratio * ratio) + 0.11116 * Math.pow(ratio, 3);
            const currentEta = Math.max(0, Math.min(1, eta_calc));
            const area_G_combined = fwhm * GAUSS_AREA_CONST;
            const area_L_combined = fwhm * LORENTZ_AREA_CONST;
            const totalArea = currentEta * area_L_combined + (1 - currentEta) * area_G_combined;
            return isFinite(totalArea) && totalArea > 0 ? totalArea : 1.0;
        }
    }
    return 1.0; // Default
}

// ---------------------------------------------------------------------------
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        calculatePeakShift,
        gaussianFromSigmaSq, lorentzianFromGamma, softFloor, applyFwhmFloor,
        calculateProfileWidths, getPeakFWHM, getPseudoVoigtArea,
        prepareVoigt, evalVoigt
    };
}
