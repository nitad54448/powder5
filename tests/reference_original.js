// AUTO-EXTRACTED reference copy of the PRE-REFACTOR worker code.
// Used only to prove the shared modules are numerically identical.
'use strict';
//   2. Define Constants & Global Worker State  

// ---------------------------------------------------------------------------
//  Profile truncation window.
//
//  A Gaussian is dead by ~4 FWHM; a Lorentzian is not. At 8 FWHM a pure
//  Lorentzian still loses ~4% of its area to truncation and another ~4% to the
//  pedestal subtraction in prepareVoigt. That in itself is survivable, because
//  the Le Bail extraction measures the area of the profile as actually drawn --
//  but the loss depends on eta, and eta moves as the widths refine, so a FIXED
//  window leaks truncation error straight into the intensity/width
//  correlation. The window is therefore scaled with the Lorentzian content,
//  and the exact analytic area of the UNtruncated shape is carried on the prep
//  object (prep.analyticArea) so integrated intensities can be reported free
//  of the truncation bias.
// ---------------------------------------------------------------------------
const CALCULATION_WINDOW_MULTIPLIER_G = 8.0;    // pure Gaussian
const CALCULATION_WINDOW_MULTIPLIER_L = 25.0;   // pure Lorentzian
const PROFILE_WINDOW_MAX_DEG = 10.0;            // hard cap, keeps cost bounded
const CALCULATION_WINDOW_MULTIPLIER = CALCULATION_WINDOW_MULTIPLIER_G; // legacy alias

// ---------------------------------------------------------------------------
//  Minimum resolvable profile FWHM, degrees 2-theta. Set from the real data
//  step (setMinProfileFwhmFromAxis) as soon as the pattern is known.
//
//  A peak narrower than a couple of data steps is not measurable, and letting
//  the refiner go there is exactly how a fit with a bad cell "improves" Rwp:
//  a spike in the wrong place costs less than a broad peak in the wrong place,
//  so the profile collapses instead of the cell correcting.
// ---------------------------------------------------------------------------
let MIN_PROFILE_FWHM_DEG = 1e-3;
function setMinProfileFwhmFromAxis(tthAxis) {
    if (!tthAxis || tthAxis.length < 3) return;
    const step = (tthAxis[tthAxis.length - 1] - tthAxis[0]) / (tthAxis.length - 1);
    if (isFinite(step) && step > 0) MIN_PROFILE_FWHM_DEG = Math.max(1e-4, 2 * step);
}

/**
 * Smooth positive branch: >= floor, strictly increasing, derivative never 0.
 *     softPositive(x, f) = (x + sqrt(x^2 + 4 f^2)) / 2
 * f at x = 0, ~x for x >> f, ~f^2/|x| for x << -f, d/dx > 0 everywhere.
 *
 * FIX: replaces the `if (sigma_sq_cdeg2 > 0) { ... }` guards in the width
 * calculators. Those guards were a hard dead zone -- once a width combination
 * went non-positive the width stuck at its default value, the one-sided
 * finite-difference column came out identically zero, and the parameter was
 * frozen there for the rest of the fit. It could not recover even after the
 * cell and zero-point had been corrected.
 */
function softPositive(x, floor) {
    if (!isFinite(x)) return floor;
    return 0.5 * (x + Math.sqrt(x * x + 4 * floor * floor));
}

// Peak heights of the UNIT-AREA Gaussian and Lorentzian, in units of 1/FWHM.
const PV_G0 = 2 * Math.sqrt(Math.LN2 / Math.PI);   // 0.9394372787
const PV_L0 = 2 / Math.PI;                         // 0.6366197724
// Areas of the UNIT-HEIGHT Gaussian and Lorentzian, in units of FWHM.
const PV_GAUSS_AREA   = 1 / PV_G0;                 // 1.0644670194
const PV_LORENTZ_AREA = 1 / PV_L0;                 // 1.5707963268

/**
 * Mixing coefficients for a pseudo-Voigt evaluated on UNIT-HEIGHT components.
 *
 * FIX (eta convention). The pseudo-Voigt is defined as eta*L + (1-eta)*G with
 * L and G normalised to unit AREA -- that is what makes eta a mixing fraction,
 * and it is the convention the Thompson-Cox-Hastings eta polynomial was fitted
 * for. This code carries a peak HEIGHT per reflection, so the components are
 * evaluated unit-height and the area normalisation is folded into these
 * coefficients, which are then renormalised so the mixture is still 1 at the
 * peak centre.
 *
 * Mixing unit-height components directly (the previous behaviour) is a
 * DIFFERENT function: L(0)/G(0) = 0.637/0.939 = 0.678 at equal FWHM, so a TCH
 * eta of 0.578 (Gamma_L/Gamma = 0.5) was rendered as a shape whose true area
 * mixing fraction is 0.48. For simple_pvoigt, where eta is refined, that was
 * only a mislabelling; for tch_aniso, where eta is COMPUTED from the
 * polynomial, the modelled Lorentzian content was systematically too high and
 * X / Y had to absorb the difference.
 */
function pvMixCoeffs(eta, hg, hl) {
    const e = Math.max(0, Math.min(1, eta));
    const wL = e * PV_L0 / Math.max(1e-12, hl);
    const wG = (1 - e) * PV_G0 / Math.max(1e-12, hg);
    const s = wL + wG;
    if (!(s > 0) || !isFinite(s)) return { cL: e, cG: 1 - e };
    return { cL: wL / s, cG: wG / s };
}

/**
 * First index i with arr[i] >= target, on an ASCENDING array. O(log n).
 * PERF FIX: replaces `let s=0; while (s<n && arr[s]<target) s++;`, which made
 * every pattern evaluation O(n_peaks * n_points) instead of
 * O(n_peaks * points_per_peak). This is the single hottest path in the fit.
 */
function lowerBound(arr, target) {
    let lo = 0, hi = arr.length;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (arr[mid] < target) lo = mid + 1; else hi = mid;
    }
    return lo;
}
// ---------------------------------------------------------------------------
//  Le Bail revival floor.
//
//  BUG FIX (absorbing zero). calculatePattern skips reflections with
//  intensity <= 0, and the decomposition weight is
//      f_j(i) = I_j * phi_j(i) / SUM_k I_k * phi_k(i),
//  so a reflection sitting at exactly zero gets f = 0 at every point, is
//  handed no intensity, and stays at zero forever. Since the extracted area is
//  clipped at zero and the net observed profile is (correctly) allowed to go
//  negative, a reflection whose window happens to sit over background -- which
//  is precisely what happens while the cell or the zero-point is still wrong --
//  had a ~50% chance of dying on the FIRST pass, including during the seeding
//  cycles, and could never come back once the cell was corrected.
//
//  The heights used for modelling are therefore floored at a tiny fraction of
//  the mean, which keeps f_j > 0 so the partition can revive, while
//  contributing nothing measurable to the pattern (1e-5 of the mean height).
//  The MEASURED quantities -- I_integrated and I_sigma -- are reported
//  unclipped, so a genuinely absent reflection is still detectable as
//  I/sigma ~ 0 (or negative) rather than being hidden behind a floor.
// ---------------------------------------------------------------------------
const LEBAIL_REVIVAL_FRACTION = 1e-5;

let workerWorkingData = null; // To store the sliced data sent from the main thread
let workerBackgroundAnchors = []; 


//   START: Monotonic Cubic Spline Helper Functions  

/**
 * Creates a monotonic cubic Hermite spline interpolation function.
 * Uses the Fritsch-Carlson method to determine tangents.
 * @param {Array<object>} points - Array of {tth, y} points, must be sorted by tth.
/**
 * Calculates the peak position shift due to sample displacement or transparency.
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

/*  ============================================================
 *  GSAS profile-function conventions (CW types 3 & 4)
 *  ============================================================
 *  Units of the profile parameters follow GSAS exactly:
 *
 *    Gaussian (Caglioti):   σ²[cdeg²]/(8·ln2)
 *                          = GU·tan²θ + GV·tanθ + GW + GP/cos²θ
 *      → H_G[deg] = √(8·ln2 · σ²) / 100
 *
 *    Lorentzian:            γ_L[cdeg]
 *                          = LX/cosθ + LY·tanθ                (TCH)
 *                          = LX/cosθ                          (Simple/Split)
 *      → γ_L[deg] = γ_L[cdeg] / 100
 *
 *  In this UI:
 *    • Simple pVoigt  GU,GV,GW,GP,LX           = GSAS GU,GV,GW,GP,LX
 *    • TCH            U,V,W                    = GSAS GU,GV,GW
 *                     X = GSAS LY (strain),  Y = GSAS LX (size)
 *    • Split pVoigt   GU_*,GV_*,GW_*,LX_*      = GSAS units
 *
 *  Stephens (1999) anisotropic strain (TOPAS / GSAS-II form):
 *    M_HKL = ΣS_HKL · h^H k^K l^L
 *    pp    = d² · √max(M_HKL, 0) / 1000
 *    γ_L_aniso[deg] = (1.8/π) · pp · η_aniso     · tanθ
 *    γ_G_aniso[deg] = (1.8/π) · pp · (1−η_aniso) · tanθ
 *
 *  All width helpers below return widths in DEGREES 2θ.
 *  These constants & helpers MUST stay in lock-step with the
 *  main-thread copy in powder5.html (calculateProfileWidths).
 */
const GSAS_GAUSSIAN_TO_DEG   = Math.sqrt(8 * Math.LN2) / 100.0;
const GSAS_LORENTZIAN_TO_DEG = 1.0 / 100.0;
const STEPHENS_PREFACTOR     = 1.8 / Math.PI;
const STEPHENS_DEFAULT_ETA   = 1.0;

/**
 * Caglioti sigma^2 (GSAS centidegree^2 / 8ln2) -> Gaussian FWHM in degrees.
 * Component floors are a QUARTER of the minimum resolvable FWHM, so neither
 * component alone pins the total; the combined floor in calculateProfileWidths
 * does that, and it preserves the Gaussian/Lorentzian ratio (hence eta).
 */
function gaussianFromSigmaSq(sigma_sq_cdeg2) {
    const floorSq = Math.pow(0.25 * MIN_PROFILE_FWHM_DEG / GSAS_GAUSSIAN_TO_DEG, 2);
    return Math.sqrt(softPositive(sigma_sq_cdeg2, floorSq)) * GSAS_GAUSSIAN_TO_DEG;
}
/** Lorentzian gamma (GSAS centidegrees) -> FWHM in degrees. */
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
    function updateHklPositions(hklList, params, system) {
    const { a, b, c, alpha, beta, gamma, lambda } = params;
    if (!hklList || hklList.length === 0) return;
    if (!params || !lambda || lambda <= 0 || !a || a <= 0) {
         hklList.forEach(peak => { if(peak) { peak.tth = null; peak.d = null; } });
         return;
    }

    const deg2rad = Math.PI / 180;
    const lambda_sq_over_4 = (lambda * lambda) / 4.0;
    const a_sq = a * a;

    let b_sq, c_sq;
    let a_star_sq, b_star_sq, c_star_sq, ab_star, bc_star, ac_star;

    // Monoclinic shares the general reciprocal metric with triclinic. With
    // the two fixed angles at 90 it reduces to the closed form it replaces,
    // exactly (agreement to 1 part in 1e16), and unique axis a, b and c all
    // come out right without a branch each.
    if (system === 'triclinic' || system === 'monoclinic') {
        const al = (alpha || 90) * deg2rad;
        const be = (beta || 90) * deg2rad;
        const ga = (gamma || 90) * deg2rad;
        const b_val = b || a;
        const c_val = c || a;
        
        const cosA = Math.cos(al), cosB = Math.cos(be), cosG = Math.cos(ga);
        const sinA = Math.sin(al), sinB = Math.sin(be), sinG = Math.sin(ga);
        
        const volume_factor = 1 - cosA*cosA - cosB*cosB - cosG*cosG + 2*cosA*cosB*cosG;
        if (volume_factor <= 1e-9 || Math.abs(sinA) < 1e-9 || Math.abs(sinB) < 1e-9 || Math.abs(sinG) < 1e-9) {
             hklList.forEach(peak => { if(peak) { peak.tth = null; peak.d = null; } });
             return;
        }
        const V = a * b_val * c_val * Math.sqrt(volume_factor);
        
        const a_star = b_val * c_val * sinA / V;
        const b_star = a * c_val * sinB / V;
        const c_star = a * b_val * sinG / V;

        const cosA_star = (cosB * cosG - cosA) / (sinB * sinG);
        const cosB_star = (cosA * cosG - cosB) / (sinA * sinG);
        const cosG_star = (cosA * cosB - cosG) / (sinA * sinB);

        a_star_sq = a_star * a_star;
        b_star_sq = b_star * b_star;
        c_star_sq = c_star * c_star;
        ab_star = 2 * a_star * b_star * cosG_star;
        bc_star = 2 * b_star * c_star * cosA_star;
        ac_star = 2 * a_star * c_star * cosB_star;
        
        } else if (system === 'tetragonal' || system === 'hexagonal' || system === 'rhombohedral' || system === 'trigonal') {
         c_sq = (c && c > 0) ? (c * c) : a_sq;
    } else if (system === 'orthorhombic') {
         b_sq = (b && b > 0) ? (b * b) : a_sq;
         c_sq = (c && c > 0) ? (c * c) : a_sq;
    }

    hklList.forEach(peak => {
        if (!peak || peak.h_orig === undefined || peak.k_orig === undefined || peak.l_orig === undefined) {
             if(peak) { peak.tth = null; peak.d = null; }
             return;
        }
        const h = peak.h_orig;
        const k = peak.k_orig;
        const l = peak.l_orig;
        const h2 = h * h;
        const k2 = k * k;
        const l2 = l * l;

        let inv_d_sq = 0;
        try {
            switch(system) {
                case 'cubic':
                    inv_d_sq = (h2 + k2 + l2) / a_sq;
                    break;
                case 'tetragonal':
                    inv_d_sq = (h2 + k2) / a_sq + l2 / c_sq;
                    break;
                case 'orthorhombic':
                    inv_d_sq = h2/a_sq + k2/b_sq + l2/c_sq;
                    break;
                case 'hexagonal':
                case 'rhombohedral':
                case 'trigonal':
                     if (a_sq <= 0 || c_sq <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = 4 * (h2 + h*k + k2) / (3 * a_sq) + l2 / c_sq;
                    break;
                case 'monoclinic':
                case 'triclinic':
                    if (a_star_sq === undefined) throw new Error("Invalid lattice param");
                    inv_d_sq = h2 * a_star_sq + k2 * b_star_sq + l2 * c_star_sq + h * k * ab_star + k * l * bc_star + h * l * ac_star;
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
                     peak.tth = 2 * thetaRad / deg2rad;
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
module.exports = { calculatePeakShift, calculateProfileWidths, getPeakFWHM,
  prepareVoigt, evalVoigt, updateHklPositions, setMinProfileFwhmFromAxis };
