// refinement_worker.js
// version 130, 26 july 2026, fixed some format in read data, changed the zoom
// version 129, 25 july 2026, Fixed TCH error
// version 123 (fixes: init race, tetragonal/triclinic HKL generation, LM step clamp)
// Spline Background in nov 2025, version 115
// FIX: must match the file the main thread loads (was stale "_v2", -> 404).
const initPromise = fetch('cctbx_space_groups_all_settings_v4.json')
    .then(res => {
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        return res.json();
    })
    .then(data => { 
        globalThis.SG_DATABASE = data; 
        // FIX: mathjs is no longer imported here. The worker only ever used it
        // for lusolve/clone/flatten/add; those are now local (solveLinearSystem
        // below, plus the hand-rolled JtJ). That removes a ~700 kB third-party
        // download from the critical path of every refinement, makes the worker
        // work offline, and eliminates a failure mode where a CDN hiccup killed
        // all fitting. importScripts cannot carry a Subresource Integrity hash,
        // so there was no way to pin it safely either.
        importScripts('sg_engine.js');

        if (typeof SG_ENGINE === 'undefined') {
            throw new Error('sg_engine.js did not load correctly.');
        }
        console.log("Worker dependencies loaded successfully.");
    })
    .catch(err => {
        console.error("Worker Error:", err);
        postMessage({ type: 'error', message: `Worker initialization failed: ${err.message}` });
        // FIX: rethrow. A bare .catch() returns a RESOLVED promise, so
        // `await initPromise` in onmessage never threw and execution fell
        // through to `SG_ENGINE.resolve()` -> ReferenceError. Rethrowing
        // makes the guard in onmessage reachable. Mark the promise handled
        // so the rethrow does not surface as an unhandled rejection.
        throw err;
    });
initPromise.catch(() => {});



// The HKL list (indices, multiplicities, allowed reflections) is generated
// on the main thread and shipped in with the refinement request, so the
// worker no longer duplicates that logic. SG_ENGINE is still needed here
// for SG_ENGINE.resolve() when the request arrives.
function enforceSymmetryConstraintsWorker(params, system) {
    if (!params || params.profileType !== "tch_aniso") return;
    switch (system) {
        case 'cubic':
            if (params.S400 !== undefined) { params.S040 = params.S400; params.S004 = params.S400; }
            if (params.S220 !== undefined) { params.S202 = params.S220; params.S022 = params.S220; }
            break;
        case 'hexagonal': case 'tetragonal': case 'rhombohedral': case 'trigonal':
            if (params.S400 !== undefined) params.S040 = params.S400;
            if (params.S202 !== undefined) params.S022 = params.S202;
            break;
    }
}

/**
 * Solves A x = b for dense square A by LU decomposition with partial pivoting.
 * Replaces math.lusolve. A is consumed (modified in place); pass a copy.
 * Returns a plain Array, or null if A is singular to working precision.
 */
function solveLinearSystem(A, b) {
    const n = A.length;
    if (n === 0) return [];
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = b[i];

    for (let col = 0; col < n; col++) {
        let piv = col, maxAbs = Math.abs(A[col][col]);
        for (let r = col + 1; r < n; r++) {
            const v = Math.abs(A[r][col]);
            if (v > maxAbs) { maxAbs = v; piv = r; }
        }
        if (!(maxAbs > 0) || !isFinite(maxAbs)) return null;   // singular column

        if (piv !== col) {
            const tr = A[piv]; A[piv] = A[col]; A[col] = tr;
            const tx = x[piv]; x[piv] = x[col]; x[col] = tx;
        }

        const d = A[col][col];
        for (let r = col + 1; r < n; r++) {
            const f = A[r][col] / d;
            if (f === 0) continue;
            A[r][col] = 0;
            for (let c = col + 1; c < n; c++) A[r][c] -= f * A[col][c];
            x[r] -= f * x[col];
        }
    }

    for (let r = n - 1; r >= 0; r--) {
        let s = x[r];
        for (let c = r + 1; c < n; c++) s -= A[r][c] * x[c];
        const d = A[r][r];
        if (!isFinite(d) || d === 0) return null;
        x[r] = s / d;
    }

    for (let i = 0; i < n; i++) if (!isFinite(x[i])) return null;
    return Array.from(x);
}

//   2. Define Constants & Global Worker State  
const CALCULATION_WINDOW_MULTIPLIER = 8.0;

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
let workerWorkingData = null; // To store the sliced data sent from the main thread
let workerBackgroundAnchors = []; 


//   START: Monotonic Cubic Spline Helper Functions  

/**
 * Creates a monotonic cubic Hermite spline interpolation function.
 * Uses the Fritsch-Carlson method to determine tangents.
 * @param {Array<object>} points - Array of {tth, y} points, must be sorted by tth.
 * @returns {function(number): number | null} - A function that takes a tth value and returns the interpolated y value, or null if spline calculation fails.
 */
function createMonotonicCubicSplineInterpolator(points) {
    const n = points.length;
    if (n < 2) return null;

    // Step 1: Calculate interval widths (h) and slopes (delta)
    const h = new Array(n - 1);
    const delta = new Array(n - 1);
    for (let i = 0; i < n - 1; i++) {
        h[i] = points[i + 1].tth - points[i].tth;
        if (h[i] <= 0) {
            console.error("Monotonic spline failed: Points must have strictly increasing tth values.");
            return null;
        }
        delta[i] = (points[i + 1].y - points[i].y) / h[i];
    }

    // Step 2: Calculate tangents (m) using Fritsch-Carlson method
    const m = new Array(n);
    // Endpoint tangents (can be simple estimates)
    m[0] = delta[0];
    m[n - 1] = delta[n - 2];
    // Internal tangents
    for (let i = 1; i < n - 1; i++) {
        if (delta[i - 1] * delta[i] <= 0) {
            m[i] = 0; // Slope changes sign, force tangent to 0
        } else {
            // Weighted average, biased towards shorter interval
             m[i] = (h[i] * delta[i - 1] + h[i-1] * delta[i]) / (h[i-1] + h[i]);
            // Alternative simple average: m[i] = (delta[i - 1] + delta[i]) / 2;
        }
    }

    // Step 3: Enforce monotonicity constraint on tangents
    for (let i = 0; i < n - 1; i++) {
        if (delta[i] === 0) { // Flat segment
            m[i] = 0;
            m[i + 1] = 0;
        } else {
            const alpha = m[i] / delta[i];
            const beta = m[i + 1] / delta[i];
            const tau = alpha * alpha + beta * beta;
            // If condition violated, scale tangents to preserve monotonicity
            if (tau > 9) { // Fritsch & Carlson condition
                const factor = 3.0 / Math.sqrt(tau);
                m[i] = alpha * delta[i] * factor;
                m[i + 1] = beta * delta[i] * factor;
            }
        }
    }

    // Step 4: Return the interpolation function using Hermite basis functions
    return function(tthValue) {
        // Flat extrapolation outside the defined range. This MUST be checked
        // before evaluating the cubic (it used to run after, wasting the work
        // and reading like a bug).
        if (tthValue <= points[0].tth)     return points[0].y;
        if (tthValue >= points[n - 1].tth) return points[n - 1].y;

        // Find the interval [i, i+1] that contains tthValue
        let i = 0;
        if (tthValue >= points[n - 1].tth) {
            i = n - 2; // Handle edge case: exactly the last point or beyond
        } else {
            while (i < n - 1 && points[i + 1].tth <= tthValue) { // Find segment where tthValue >= start
                i++;
            }
        }
         // Handle edge case: before the first point
         if (i < 0) i = 0;
         if (i >= n - 1) i = n - 2; // Should not happen with above checks, but safety

        const x_i = points[i].tth;
        const x_ip1 = points[i + 1].tth;
        const y_i = points[i].y;
        const y_ip1 = points[i + 1].y;
        const h_i = h[i]; // Use pre-calculated h[i] = x_ip1 - x_i
        const m_i = m[i];
        const m_ip1 = m[i + 1];

        // Normalized position within interval [0, 1]
        const t = (h_i > 1e-9) ? (tthValue - x_i) / h_i : 0; // Avoid division by zero

        // Cubic Hermite spline basis functions
        const h00 = 2 * t * t * t - 3 * t * t + 1;
        const h10 = t * t * t - 2 * t * t + t;
        const h01 = -2 * t * t * t + 3 * t * t;
        const h11 = t * t * t - t * t;

        // Interpolation formula
        return h00 * y_i + h10 * h_i * m_i + h01 * y_ip1 + h11 * h_i * m_ip1;
    };
}


    function calculateTotalBackground(tthAxis, params, splinePoints, outArr = null) {
        const n = tthAxis.length;
        if (n === 0 || !splinePoints || splinePoints.length < 2) {
            if (outArr) { outArr.fill(0); return outArr; }
            return new Float64Array(n);
        }

        const background = outArr || new Float64Array(n);
        if (outArr) background.fill(0);
        
        const sortedPoints = [...splinePoints].sort((a, b) => a.tth - b.tth);
        
        let interpolate = null;
        try {
            // Attempt to create the monotonic cubic spline interpolator
            interpolate = createMonotonicCubicSplineInterpolator(sortedPoints);
            if (!interpolate) throw new Error("Interpolator creation returned null."); // Explicitly check
        } catch (splineError) {
             console.error("Failed to create monotonic cubic spline, falling back to linear:", splineError);
             interpolate = null; // Ensure it's null if creation failed
        }

        // If spline creation failed, use linear interpolation as fallback
        if (!interpolate) {
            console.warn("Using linear interpolation for background (Monotonic spline failed).");
            let p_idx = 0;
            const numSplinePoints = sortedPoints.length; // Use length of sortedPoints
            for (let i = 0; i < n; i++) {
                const tth = tthAxis[i];
                 // Find segment using linear scan
                 while (p_idx < numSplinePoints - 2 && sortedPoints[p_idx + 1].tth < tth) {
                     p_idx++;
                 }
                const p1 = sortedPoints[p_idx];
                const p2 = sortedPoints[Math.min(p_idx + 1, numSplinePoints - 1)]; // Ensure p2 index is valid

                if (tth <= p1.tth) {
                    background[i] = p1.y;
                } else if (tth >= p2.tth) {
                    background[i] = p2.y;
                } else {
                    const tth_diff = p2.tth - p1.tth;
                    background[i] = (tth_diff > 1e-9) ? p1.y + (p2.y - p1.y) * (tth - p1.tth) / tth_diff : p1.y;
                }
                // Ensure non-negative background
                if (background[i] < 0) background[i] = 0;
            }
        } else {
            // Use the monotonic cubic spline interpolator
            for (let i = 0; i < n; i++) {
                background[i] = interpolate(tthAxis[i]);
                 // Ensure non-negative background
                 if (background[i] < 0) background[i] = 0;
            }
        }

        return background;
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

    let b_sq, c_sq, sin_beta_sq, cos_beta;
    let a_star_sq, b_star_sq, c_star_sq, ab_star, bc_star, ac_star;

    if (system === 'triclinic') {
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
        
    } else if (system === 'monoclinic') {
        const beta_rad = (beta || 90) * deg2rad;
        sin_beta_sq = Math.sin(beta_rad);
        sin_beta_sq *= sin_beta_sq;
        cos_beta = Math.cos(beta_rad);
        b_sq = (b && b > 0) ? (b * b) : a_sq;
        c_sq = (c && c > 0) ? (c * c) : a_sq;
         if (Math.abs(sin_beta_sq) < 1e-9) {
             hklList.forEach(peak => { if(peak) { peak.tth = null; peak.d = null; } });
             return;
         }
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
                     if (a_sq <= 0 || b_sq <= 0 || c_sq <= 0 || a <= 0 || c <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = (1/sin_beta_sq) * (h2/a_sq + k2*sin_beta_sq/b_sq + l2/c_sq - (2*h*l*cos_beta)/(a*c));
                    break;
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
 * Calculates widths for simple_pvoigt (GSAS units → degrees 2θ).
 */
function _calculateWidths_Simple(tth, hkl, params, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe) {
    let gamma_G = 1e-4;
    let gamma_L = 1e-4;

    const GU = params.GU || 0;
    const GV = params.GV || 0;
    const GW = params.GW || 0;
    const GP = params.GP || 0;
    const LX = params.LX || 0;

    //   Caglioti in cdeg²/(8 ln2). Convert √σ² to FWHM in degrees.
    const sigma_sq_cdeg2 = GU * tanTheta * tanTheta + GV * tanTheta + GW + GP / cosThetaSq_safe;
    if (sigma_sq_cdeg2 > 0 && isFinite(sigma_sq_cdeg2)) {
        gamma_G = Math.sqrt(sigma_sq_cdeg2) * GSAS_GAUSSIAN_TO_DEG;
    }
    //   Lorentzian (size only) in centidegrees → degrees.
    const gL_cdeg = LX / cosTheta_safe;
    if (gL_cdeg > 0 && isFinite(gL_cdeg)) {
        gamma_L = gL_cdeg * GSAS_LORENTZIAN_TO_DEG;
    }

    return { gamma_G, gamma_L };
}

/**
 * Calculates widths for split_pvoigt (GSAS units → degrees 2θ).
 */
function _calculateWidths_Split(tth, hkl, params, side, safeThetaRad, tanTheta, cosTheta_safe, cosThetaSq_safe) {
    let gamma_G = 1e-4;
    let gamma_L = 1e-4;
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
    if (sigma_sq_cdeg2 > 0 && isFinite(sigma_sq_cdeg2)) {
        gamma_G = Math.sqrt(sigma_sq_cdeg2) * GSAS_GAUSSIAN_TO_DEG;
    }
    const gL_cdeg = LX / cosTheta_safe;
    if (gL_cdeg > 0 && isFinite(gL_cdeg)) {
        gamma_L = gL_cdeg * GSAS_LORENTZIAN_TO_DEG;
    }

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
    let gamma_G = 1e-4;
    let gamma_L = 1e-4;

    const U = params.U || 0;
    const V = params.V || 0;
    const W = params.W || 0;
    const X = params.X || 0;     //  GSAS LY  (strain · tanθ)
    const Y = params.Y || 0;     //  GSAS LX  (size   / cosθ)

    const sigma_sq_cdeg2 = U * tanTheta * tanTheta + V * tanTheta + W;
    if (sigma_sq_cdeg2 > 0 && isFinite(sigma_sq_cdeg2)) {
        gamma_G = Math.sqrt(sigma_sq_cdeg2) * GSAS_GAUSSIAN_TO_DEG;
    }
    const gL_cdeg = X * tanTheta + Y / cosTheta_safe;
    if (gL_cdeg > 0 && isFinite(gL_cdeg)) {
        gamma_L = gL_cdeg * GSAS_LORENTZIAN_TO_DEG;
    }

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
    if (!params || !params.profileType) return { gamma_G: 1e-4, gamma_L: 1e-4 };
    
    const profileType = params.profileType;
    const thetaRad = tth * (Math.PI / 180) / 2;

    const MAX_ANGLE_RAD = Math.PI / 2.0 - 1e-6;
    const safeThetaRad = Math.min(thetaRad, MAX_ANGLE_RAD);
     if (safeThetaRad < 1e-6) {
         return { gamma_G: 1e-4, gamma_L: 1e-4 };
     }

    const tanTheta = Math.tan(safeThetaRad);
    const cosTheta = Math.cos(safeThetaRad);
    const cosTheta_safe = Math.max(cosTheta, 1e-9);
    const cosThetaSq_safe = Math.max(cosTheta * cosTheta, 1e-9);

    let widths = { gamma_G: 1e-4, gamma_L: 1e-4 };

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

    return {
        gamma_G: Math.max(1e-4, isFinite(widths.gamma_G) ? widths.gamma_G : 1e-4),
        gamma_L: Math.max(1e-4, isFinite(widths.gamma_L) ? widths.gamma_L : 1e-4)
    };
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
        prep = { type: 1, x0, asym_param, H_G, H_L, eta, Cg };
    } else if (profileType === "split_pvoigt") {
        const wL = calculateProfileWidths(tth_peak, hkl, params, 'left');
        const wR = calculateProfileWidths(tth_peak, hkl, params, 'right');
        const H_G_L = Math.max(1e-9, wL.gamma_G), H_L_L = Math.max(1e-9, wL.gamma_L);
        const H_G_R = Math.max(1e-9, wR.gamma_G), H_L_R = Math.max(1e-9, wR.gamma_L);
        const eta = Math.max(0, Math.min(1, params.eta_split || 0.5));
        fwhm_total = Math.max(getPeakFWHM(H_G_L, H_L_L), getPeakFWHM(H_G_R, H_L_R));
        prep = { type: 2, x0, asym_param, H_G_L, H_L_L, H_G_R, H_L_R, eta, Cg };
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
        prep = { type: 3, x0, asym_param, fwhm, eta, Cg };
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
    const truncationRadius = CALCULATION_WINDOW_MULTIPLIER * Math.max(0.01, prep.fwhm_total)
                             * (asym_param !== 0 ? 2.0 : 1.0);

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
        // FIX: use the SIGNED delta. The old code took Math.abs(delta) inside
        // the correction, so BOTH flanks were stretched by the same factor and
        // SL/HL produced pure symmetric broadening -- no asymmetry whatsoever.
        // With the signed form one flank is broadened (factor > 1) and the
        // other sharpened (factor < 1), which is the point of the correction.
        const t = Math.max(-0.95, Math.min(0.95, p.asym_param * delta));
        delta /= (1.0 - t);   // 1 - t is in [0.05, 1.95], always positive
    }

    let hg, hl;
    if (p.type === 1) {
        hg = p.H_G; hl = p.H_L;
    } else if (p.type === 2) {
        hg = delta < 0 ? p.H_G_L : p.H_G_R;
        hl = delta < 0 ? p.H_L_L : p.H_L_R;
    } else {
        hg = hl = p.fwhm;
    }

    const dg = (delta / hg) * (delta / hg);
    const dl = (delta / hl) * (delta / hl);
    const g = Math.exp(-p.Cg * dg);
    const l = 1.0 / (1.0 + 4.0 * dl);
    const v = p.eta * l + (1.0 - p.eta) * g - p.pedestal;
    return v > 0 ? v : 0.0;
}

    function calculatePattern(tthAxis, hklList, params, outArr = null) {
    const n_points = tthAxis ? tthAxis.length : 0;
    if (n_points === 0 || !hklList || hklList.length === 0 || !params) {
        if (outArr) { outArr.fill(0); return outArr; }
        return new Float64Array(n_points);
    }

    const pattern = outArr || new Float64Array(n_points);
    if (outArr) pattern.fill(0);
    
    const deg2rad = Math.PI / 180;
    
    const lambda1 = params.lambda || 1.54056;
    const lambda2 = params.lambda2 || 0;
    const ratio21 = params.ratio || 0;
    const zeroShift = params.zeroShift || 0;
    //   K-alpha 1  
    // The window now comes from prep.max_d, so the loop bounds and the
    // profile's own truncation radius can never disagree. This also removes
    // the duplicate calculateProfileWidths/getPeakFWHM call that used to be
    // made purely to size the window, and it evaluates the widths at the same
    // angle prepareVoigt uses (basePos1) rather than at peak.tth.
    hklList.forEach(peak => {
        if (!peak || !peak.intensity || peak.intensity <= 0 || !peak.tth || peak.tth < 0 || peak.tth > 180) return;

        const basePos1 = peak.tth + zeroShift;
        const peakPos1 = basePos1 + calculatePeakShift(basePos1, params);

        const prep1 = prepareVoigt(basePos1, peakPos1, peak, params);
        const min_tth1 = peakPos1 - prep1.max_d;
        const max_tth1 = peakPos1 + prep1.max_d;

        const startIndex = lowerBound(tthAxis, min_tth1);
        if (startIndex === n_points) return;

        for (let i = startIndex; i < n_points; i++) {
            const current_tth = tthAxis[i];
            if (current_tth > max_tth1) break;
            pattern[i] += peak.intensity * evalVoigt(current_tth, prep1);
        }
    });

    //   K-alpha 2  
    const doubletEnabled = ratio21 > 1e-6 && lambda2 > 1e-6 && Math.abs(lambda1 - lambda2) > 1e-6;
    if (doubletEnabled) {
        const lambdaRatio = lambda2 / lambda1;
        hklList.forEach(peak => {
            if (!peak || !peak.intensity || peak.intensity <= 0 || !peak.tth || peak.tth <= 0 || peak.tth >= 180) return;

            const sinTheta1 = Math.sin(peak.tth * deg2rad / 2.0);
            const sinTheta2 = sinTheta1 * lambdaRatio;
            if (Math.abs(sinTheta2) >= 1) return;

            const tth2 = 2 * Math.asin(sinTheta2) / deg2rad;
            const basePos2 = tth2 + zeroShift;
            const peakPos2 = basePos2 + calculatePeakShift(basePos2, params);

            const prep2 = prepareVoigt(basePos2, peakPos2, peak, params);
            const min_tth2 = peakPos2 - prep2.max_d;
            const max_tth2 = peakPos2 + prep2.max_d;

            const startIndex2 = lowerBound(tthAxis, min_tth2);
            if (startIndex2 === n_points) return;

            for (let i = startIndex2; i < n_points; i++) {
                const current_tth = tthAxis[i];
                if (current_tth > max_tth2) break;
                pattern[i] += peak.intensity * ratio21 * evalVoigt(current_tth, prep2);
            }
        });
    }

    for (let i = 0; i < n_points; i++) {
        if (!isFinite(pattern[i])) {
            pattern[i] = 0;
        }
    }
    return pattern;
}


function leBailIntensityExtraction(expData, hklList, params, scratch_calc = null) {
    if (!expData || !expData.tth || !expData.intensity || !expData.background || !hklList) return;
    const n_points = expData.tth.length;
    if (n_points === 0) return;

    // If all peak intensities are 0, seed them to 1000 so total_calc > 0
    let maxInt = 0;
    for (let i = 0; i < hklList.length; i++) {
        if ((hklList[i].intensity || 0) > maxInt) maxInt = hklList[i].intensity;
    }
    if (maxInt <= 1e-6) {
        hklList.forEach(p => { if (p) p.intensity = 1000; });
    }
    // -----------------------------------

    const deg2rad = Math.PI / 180;
    const lambda1 = params.lambda || 1.54056;
    const lambda2 = params.lambda2 || 0;
    const ratio21 = params.ratio || 0;
    const zeroShift = params.zeroShift || 0;
    const doubletEnabled = ratio21 > 1e-6 && lambda2 > 1e-6 && Math.abs(lambda1 - lambda2) > 1e-6;
    const lambdaRatio = doubletEnabled ? lambda2 / lambda1 : 1.0;

    // 1. Compute net observed
    const y_obs_net = new Float64Array(n_points);
    for (let i = 0; i < n_points; i++) {
        // FIX: do NOT clip at zero here. Clipping away the negative half of the
        // noise leaves only positive excursions, which the decomposition then
        // hands out as real intensity -- a systematic positive bias that lands
        // hardest on weak and absent reflections, exactly the ones you most
        // need to be honest about. Measured on a synthetic pattern, a
        // reflection with true I = 4.2 came back as 14.5. The partition is
        // linear in y, so letting it go negative is unbiased; the resulting
        // intensity is clipped once, at the end.
        y_obs_net[i] = expData.intensity[i] - expData.background[i];
    }

    // 2. Compute total calculated pattern with current peak intensities
    const total_calc = calculatePattern(expData.tth, hklList, params, scratch_calc);

    // 3. Extract area per peak without storing M arrays of size N
    const currentCycleIntensities = new Float64Array(hklList.length);
    const currentCycleShapeAreas   = new Float64Array(hklList.length);
    const currentCycleVariances    = new Float64Array(hklList.length);

    hklList.forEach((peak, j) => {
        if (!peak || !peak.tth || peak.tth <= 0 || peak.tth >= 180) return;

        const basePos1 = peak.tth + zeroShift;
        const peakPos1 = basePos1 + calculatePeakShift(basePos1, params);
        const prep1 = prepareVoigt(basePos1, peakPos1, peak, params);
        const min_tth1 = peakPos1 - prep1.max_d;
        const max_tth1 = peakPos1 + prep1.max_d;

        let tth2 = null, min_tth2 = 0, max_tth2 = 0, prep2 = null;
        if (doubletEnabled) {
             const sinTheta1 = Math.sin(peak.tth * deg2rad / 2.0);
             const sinTheta2 = sinTheta1 * lambdaRatio;
             if (Math.abs(sinTheta2) < 1) {
                 tth2 = 2 * Math.asin(sinTheta2) / deg2rad;
                 const basePos2 = tth2 + zeroShift;
                 const peakPos2 = basePos2 + calculatePeakShift(basePos2, params);
                 prep2 = prepareVoigt(basePos2, peakPos2, peak, params);
                 min_tth2 = peakPos2 - prep2.max_d;
                 max_tth2 = peakPos2 + prep2.max_d;
             }
        }

        // Find loop bounds (binary search; see lowerBound)
        const actualMin = tth2 ? Math.min(min_tth1, min_tth2) : min_tth1;
        const actualMax = tth2 ? Math.max(max_tth1, max_tth2) : max_tth1;

        const startIdx = lowerBound(expData.tth, actualMin);
        
        // Le Bail (1988) decomposition, exactly as written:
        //
        //   I_hkl = SUM_i  [ y_obs(i) - y_bkg(i) ] * -------------------------
        //                                            I_hkl*W_hkl(i)
        //                                            SUM_k I_k*W_k(i)
        //
        // The profile shape's own integral is accumulated over the SAME points
        // in the SAME loop. That is what converts the integrated intensity to
        // the height this code stores, and doing it here rather than from the
        // analytic pseudo-Voigt area matters: the rendered profile is truncated
        // at max_d and pedestal-subtracted, so its true area is a few percent
        // below the analytic value. Using the analytic area made the extracted
        // intensity inconsistent with the peak actually drawn -- and it costs
        // nothing to do it right, because we are already walking these points.
        let extracted_area = 0;
        let shape_area = 0;
        let variance = 0;
        let prev_partitioned_I = 0;
        let prev_profile_val = 0;

        // Uncertainty on the extracted intensity, by propagating counting
        // statistics through the decomposition. Writing
        //     I = SUM_i c_i * f_i * (y_i - b_i)
        // with c_i the trapezoid weight of point i and f_i its partition
        // fraction (held fixed, as is standard), Poisson counts give
        //     sigma^2(I) = SUM_i (c_i * f_i)^2 * y_i .
        // c_i spans half of each adjacent interval, so a point's weight is only
        // known once the NEXT point is seen; carry it forward and settle up.
        let pend_c = 0, pend_f = 0, pend_var = 0, havePending = false;

        for (let i = startIdx; i < n_points; i++) {
            const current_tth = expData.tth[i];
            if (current_tth > actualMax) break;

            let profile_val = 0;
            if (current_tth >= min_tth1 && current_tth <= max_tth1) {
                profile_val += evalVoigt(current_tth, prep1);
            }
            if (tth2 && current_tth >= min_tth2 && current_tth <= max_tth2) {
                profile_val += ratio21 * evalVoigt(current_tth, prep2);
            }

            let current_fraction = 0;
            if (profile_val > 0 && total_calc[i] > 1e-9) {
                current_fraction = (peak.intensity * profile_val) / total_calc[i];
            }

            const current_partitioned_I = y_obs_net[i] * current_fraction;

            let half = 0;
            if (i > startIdx) {
                const step_width = expData.tth[i] - expData.tth[i-1];
                if (step_width > 0) {
                    extracted_area += ((prev_partitioned_I + current_partitioned_I) / 2) * step_width;
                    shape_area     += ((prev_profile_val   + profile_val)         / 2) * step_width;
                    half = step_width / 2;
                }
            }

            if (havePending) {
                const c = pend_c + half;                 // now complete
                variance += (c * pend_f) * (c * pend_f) * pend_var;
            }
            pend_c   = half;
            pend_f   = current_fraction;
            pend_var = Math.max(1, expData.intensity[i]);   // Poisson on raw counts
            havePending = true;

            prev_partitioned_I = current_partitioned_I;
            prev_profile_val   = profile_val;
        }
        if (havePending) {
            variance += (pend_c * pend_f) * (pend_c * pend_f) * pend_var;
        }

        currentCycleIntensities[j] = extracted_area;
        currentCycleShapeAreas[j]  = shape_area;
        currentCycleVariances[j]   = variance;
    });

    hklList.forEach((peak, idx) => {
        if (!peak) return;
        const extracted_area = Math.max(0, currentCycleIntensities[idx] || 0);
        const shape_area     = currentCycleShapeAreas[idx] || 0;

        // FIX: no `Math.max(1.0, ...)` floor. That floor stopped any reflection
        // from ever going to zero, so systematically absent reflections kept a
        // permanent phantom intensity. Le Bail intensities are allowed to
        // refine to nothing; only negatives are clipped.
        peak.intensity  = (shape_area > 1e-12) ? (extracted_area / shape_area) : 0;
        peak.shapeArea  = shape_area;        // integral of the rendered Ka1+Ka2 shape
        peak.I_integrated = extracted_area;  // the Le Bail intensity itself
        peak.I_sigma = Math.sqrt(Math.max(0, currentCycleVariances[idx] || 0));
        delete peak.intensity_previous;
    });
}


//   Statistics  
function calculateStatistics(localWorkingData, netCalcPattern, fitFlags, finalBackground, params, hklList, refinementMode) {
    const y_obs = localWorkingData.intensity;
    const y_bkg = finalBackground;
    const weights = localWorkingData.weights;
    const N = y_obs.length;

    if (!y_obs || !netCalcPattern || !y_bkg || !weights ||
        N === 0 || y_obs.length !== netCalcPattern.length || y_obs.length !== y_bkg.length || y_obs.length !== weights.length) {
        console.error("Statistics calculation error: Mismatched or invalid array inputs.");
        return { r_p: -1, rwp: -1, chi2: -1, scaleFactor: 1, sum_w_res_sq: 0 };
    }

    // No separate scale factor: after the Le Bail fix, `netCalcPattern` is
    // already per-peak-correct (intensities were extracted inside the
    // refiner's objective before calculating the pattern).
    const y_calc = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        y_calc[i] = netCalcPattern[i] + y_bkg[i];
    }

    // FIX: Rp and Rwp now use the SAME convention (both on the full observed
    // profile, as GSAS/FullProf report them). Previously Rp was divided by
    // the BACKGROUND-SUBTRACTED sum while Rwp used the raw profile, so on
    // high-background data Rp looked far worse than Rwp for no real reason.
    // The background-subtracted variants are computed alongside and reported
    // as r_p_net / rwp_net for anyone who wants them.
    let sum_w_res_sq = 0, sum_w_obs_sq = 0, sum_w_obs_net_sq = 0;
    let sum_abs_res = 0, sum_abs_obs = 0, sum_abs_obs_net = 0;
    for (let i = 0; i < N; i++) {
        const obs_i = y_obs[i];
        const calc_i = y_calc[i];
         const w_i = (weights[i] !== undefined && isFinite(weights[i])) ? weights[i] : 1.0;


         if (isFinite(obs_i) && isFinite(calc_i)) {
            const res = obs_i - calc_i;
            const obs_net = obs_i - y_bkg[i];

            sum_w_res_sq += w_i * res * res;
            sum_w_obs_sq += w_i * obs_i * obs_i;
            sum_w_obs_net_sq += w_i * obs_net * obs_net;
            sum_abs_res += Math.abs(res);
            sum_abs_obs += Math.abs(obs_i);
            sum_abs_obs_net += Math.abs(obs_net);
         }
    }

    const Rp  = (sum_abs_obs   > 1e-9) ? 100 * (sum_abs_res / sum_abs_obs) : 0;
    const Rwp = (sum_w_obs_sq  > 1e-9) ? 100 * Math.sqrt(Math.max(0, sum_w_res_sq / sum_w_obs_sq)) : 0;
    const RpNet  = (sum_abs_obs_net  > 1e-9) ? 100 * (sum_abs_res / sum_abs_obs_net) : 0;
    const RwpNet = (sum_w_obs_net_sq > 1e-9) ? 100 * Math.sqrt(Math.max(0, sum_w_res_sq / sum_w_obs_net_sq)) : 0;


     const { paramMapping } = getParameterMapping(fitFlags || {}, params || {}, hklList || [], refinementMode || 'le-bail');
     const P_base = paramMapping.filter(m => !m.isIntensity).length;
     let P = P_base;
     if (refinementMode === 'pawley' && hklList && localWorkingData && localWorkingData.tth && localWorkingData.tth.length > 0) {
          const tthMin = localWorkingData.tth[0];
          const tthMax = localWorkingData.tth[localWorkingData.tth.length - 1];
         const refinedIntensitiesCount = hklList.filter(hkl =>
             hkl && hkl.tth && hkl.tth >= tthMin && hkl.tth <= tthMax
         ).length;
         P += refinedIntensitiesCount;
     }

    const degreesOfFreedom = N - P;

    let chi2 = 0;
    if (degreesOfFreedom > 0) {
        chi2 = sum_w_res_sq / degreesOfFreedom;
        if (!isFinite(chi2) || chi2 < 0) chi2 = 0;
    }

    return {
        r_p: isFinite(Rp) ? Rp : -1,
        rwp: isFinite(Rwp) ? Rwp : -1,
        r_p_net: isFinite(RpNet) ? RpNet : -1,
        rwp_net: isFinite(RwpNet) ? RwpNet : -1,
        chi2: isFinite(chi2) ? chi2 : -1,
        scaleFactor: 1.0,
        sum_w_res_sq: isFinite(sum_w_res_sq) ? sum_w_res_sq : 0
    };
}

//   Refinement Algorithms (LM, PT)  
async function refineParametersLM(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
        const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode);
        if (paramMapping.length === 0) {
             const finalProgress = (leBailCycle + 1) / totalLeBailCycles;
             postMessage({ type: 'progress', value: Math.min(1.0, finalProgress) });
            return { params: initialParams, hklList: hklList, ss_res: 0 };
        }

        let params = initialParams;
        let workingHklList = JSON.parse(JSON.stringify(hklList));

        let finalJtJ = null, ss_res = Infinity;
        let lambda = 0.001;

        const y_obs = workerWorkingData.intensity;
        const sqrt_weights = workerWorkingData.weights.map(w => Math.sqrt(w));
        const n_points = y_obs.length;
        const n_params = paramMapping.length;

        const residuals = new Float64Array(n_points);
        const y_calc_total = new Float64Array(n_points);
        const y_calc_baseline = new Float64Array(n_points);
        const jacobian_col = new Float64Array(n_points);


        const scratch_bkg = new Float64Array(n_points);
        const scratch_pattern = new Float64Array(n_points);

        const baseProgress = leBailCycle / totalLeBailCycles;
        const cycleProgressSpan = 1 / totalLeBailCycles;

        // ===================================================================
        // CRITICAL FIX -- Le Bail intensities are extracted ONCE PER CYCLE and
        // then held FIXED for every function evaluation inside this call.
        //
        // The previous version re-extracted them inside the objective, on every
        // single evaluation. That makes profile-WIDTH refinement degenerate:
        // the extraction re-partitions the observed profile against whatever
        // shape it is handed, so it silently absorbs the width error and
        // flattens the very curvature the refiner needs. Measured on a
        // synthetic 9-reflection cubic pattern (true LX = 6):
        //
        //     LX      Rwp re-extracting     Rwp with fixed intensities
        //      6           5.19 %                    5.19 %
        //     12          20.15 %                   32.34 %
        //     25          40.11 %                  103.39 %
        //    200          56.66 %                  658.34 %
        //
        // Re-extracting flattens the landscape and plateaus around 56 %, so LM
        // wanders into the flat region and never returns; with fixed
        // intensities the minimum at the true width is sharp and it converges.
        // Extract-then-refine, repeated over a few outer cycles, is the actual
        // Le Bail method.
        // ===================================================================
        // Re-extract the Le Bail intensities ONCE at the top of each LM
        // iteration. They are then constants for that whole iteration --
        // baseline, every Jacobian column, and the trial step -- which is what
        // makes the derivatives and the accept/reject test meaningful. This is
        // the standard extract-then-least-squares cycle.
        const refreshLeBailIntensities = () => {
            if (refinementMode !== 'le-bail') return;
            enforceSymmetryConstraintsWorker(params, system);
            updateHklPositions(workingHklList, params, system);
            const bkg0 = calculateTotalBackground(workerWorkingData.tth, params,
                                                  workerBackgroundAnchors, scratch_bkg);
            leBailIntensityExtraction(
                { tth: workerWorkingData.tth, intensity: y_obs, background: bkg0 },
                workingHklList, params, scratch_pattern
            );
        };

        const calculateTotalPattern = (targetArray) => {
            enforceSymmetryConstraintsWorker(params, system);
            updateHklPositions(workingHklList, params, system);

            const y_bkg = calculateTotalBackground(workerWorkingData.tth, params, workerBackgroundAnchors, scratch_bkg);
            const netCalcPattern = calculatePattern(workerWorkingData.tth, workingHklList, params, scratch_pattern);
            for (let i = 0; i < n_points; i++) {
                targetArray[i] = netCalcPattern[i] + y_bkg[i];
            }
            return 1.0;   // no separate scale factor
        };


        //   LM main loop                  -
        //   Fixes (April 2026):
        //     B. Convergence test moved to the END of the iteration, and
        //        only triggers on an ACCEPTED step. The old placement
        //        (before the step) could confuse "last step was rejected"
        //        with "we've converged" and exit prematurely.
        //     C. On rejection, we restore oldParams/oldHklList correctly,
        //        AND we do NOT advance ss_res — it still reflects the last
        //        accepted cost.
        //     D. Finite-difference step now uses a larger relative floor and,
        //        when the parameter is at (or very near) its lower bound,
        //        flips sign so the perturbed point stays feasible.
        //        This avoids silent clamping by set() during Jacobian build,
        //        which would have zeroed that parameter's column.
        //
        //   Also: since get/set are now raw-space, the LM solve itself is in
        //   raw space. Damping uses Marquardt's lambda*diag(JtJ) (see below),
        //   which is the local-curvature scaling that replaces the old
        //   initial-value normalization.
        let lastAcceptedCost = Infinity;
        const deadParamWarned = new Set();
        let lastProgressPct = -1;

        for (let iter = 0; iter < maxIter; iter++) {
             let oldParams = null;
             let oldHklList = null;
            try {
                 refreshLeBailIntensities();
                 calculateTotalPattern(y_calc_baseline);

                 let cost = 0;
                 let residualFailed = false;
                 for (let i = 0; i < n_points; i++) {
                     residuals[i] = (y_obs[i] - y_calc_baseline[i]) * sqrt_weights[i];
                     if (isFinite(residuals[i])) {
                          cost += residuals[i] * residuals[i];
                     } else {
                          residualFailed = true;
                          break;
                     }
                 }
                 if (residualFailed) {
                     lambda = Math.min(1e9, lambda * 10);
                     console.warn(`LM iter ${iter}: Residual calculation failed. Increased lambda to ${lambda}.`);
                     continue;
                 }

                 // First-iteration seed for the accepted cost (no step yet).
                 if (iter === 0) lastAcceptedCost = cost;
                 ss_res = cost;



                 //   Build Jacobian column by column (one-sided FD, bound-aware)  
                 const jacobian_T = [];
                 const savedIntensities = workingHklList.map(h => h ? h.intensity : 0);

                 for (let p = 0; p < n_params; p++) {

                     const mapping = paramMapping[p];
                     const originalValue = mapping.get(params, workingHklList);

                     const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
                     const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;

                     // Base step: relative, with an absolute floor and a floor
                     // tied to the parameter's declared typical magnitude (so a
                     // parameter sitting AT zero still gets a usable probe).
                     const baseStep = Math.max(1e-7,
                                               Math.abs(originalValue) * 1e-4,
                                               (mapping.defaultScale || 0) * 1e-4);

                     // ADAPTIVE FD.
                     // A width parameter can fall into a dead zone: gamma is
                     // clamped by `Math.max(1e-4, ...)` inside
                     // calculateProfileWidths, so once e.g. GW drops near its
                     // lower bound the clamp swallows any small perturbation,
                     // the column comes out identically zero, and the parameter
                     // is frozen there for the rest of the fit -- it cannot
                     // recover even after the cell and zero-point are corrected.
                     // Escalate the probe until the model actually responds.
                     let fd_step = baseStep;
                     let fd_sign = 1;
                     let responded = false;

                     for (let attempt = 0; attempt < 5; attempt++) {
                         // Step away from whichever bound we are against, so
                         // set() cannot silently clamp the probe and zero the
                         // column. If hemmed in on both sides, fall back to +.
                         fd_sign = (originalValue + fd_step > maxV && originalValue - fd_step >= minV) ? -1 : 1;

                         mapping.set(params, workingHklList, originalValue + fd_sign * fd_step);
                         const applied = mapping.get(params, workingHklList);
                         const actual_h = applied - originalValue;   // what set() really did

                         if (actual_h !== 0) {
                             calculateTotalPattern(y_calc_total);
                             for (let i = 0; i < n_points; i++) {
                                 if (!isFinite(y_calc_baseline[i]) || !isFinite(y_calc_total[i])) {
                                     jacobian_col[i] = 0;
                                 } else {
                                     // Divide by the step actually taken, not the
                                     // one requested -- they differ whenever set()
                                     // clamps at a bound.
                                     jacobian_col[i] = (y_calc_total[i] - y_calc_baseline[i])
                                                       / actual_h * sqrt_weights[i];
                                 }
                             }
                             for (let i = 0; i < n_points; i++) {
                                 if (jacobian_col[i] !== 0) { responded = true; break; }
                             }
                         }

                         if (responded) break;
                         fd_step *= 100;    // 1e2 .. 1e8 x the base step
                     }

                     if (!responded) jacobian_col.fill(0);

                     mapping.set(params, workingHklList, originalValue);
                     workingHklList.forEach((h, idx) => { if (h) h.intensity = savedIntensities[idx]; });
                     // PERF FIX: keep the column typed. `[...jacobian_col]`
                     // unboxed a Float64Array into a plain Array of n_points
                     // heap numbers, n_params times per iteration.
                     jacobian_T.push(jacobian_col.slice());
                 }

                 // PERF FIX: form JtJ and Jtr directly instead of
                 // math.transpose() + math.multiply(). The transpose allocated
                 // a throwaway n_points x n_params boxed matrix purely to be
                 // consumed once, and mathjs's generic multiply dispatch is
                 // far slower than a typed inner loop. JtJ is symmetric, so
                 // only the upper triangle is computed.
                 const JtJ = new Array(n_params);
                 for (let i = 0; i < n_params; i++) JtJ[i] = new Array(n_params).fill(0);
                 const Jtr = new Array(n_params).fill(0);
                 let productFailed = false;

                 for (let i = 0; i < n_params; i++) {
                     const ci = jacobian_T[i];
                     let s_r = 0;
                     for (let k = 0; k < n_points; k++) s_r += ci[k] * residuals[k];
                     if (!isFinite(s_r)) productFailed = true;
                     Jtr[i] = s_r;

                     for (let j = i; j < n_params; j++) {
                         const cj = jacobian_T[j];
                         let s = 0;
                         for (let k = 0; k < n_points; k++) s += ci[k] * cj[k];
                         if (!isFinite(s)) productFailed = true;
                         JtJ[i][j] = s;
                         JtJ[j][i] = s;
                     }
                 }

                 if (productFailed) {
                      console.warn(`LM iter ${iter}: JtJ/Jtr contained non-finite entries.`);
                      finalJtJ = null;
                      lambda = Math.min(1e9, lambda * 5);
                      continue;
                 }
                 finalJtJ = JtJ;

                 // A parameter with an identically-zero Jacobian column has no
                 // effect on the model -- it is either symmetry-slaved or its
                 // set() is not reaching the model. It would silently consume a
                 // degree of freedom and produce a meaningless ESD, so say so.
                 for (let i = 0; i < n_params; i++) {
                     if (JtJ[i][i] === 0 && !deadParamWarned.has(paramMapping[i].name)) {
                         deadParamWarned.add(paramMapping[i].name);
                         console.warn(`LM: parameter "${paramMapping[i].name}" has no effect on the calculated pattern (zero Jacobian column).`);
                         postMessage({ type: 'warning', message: `Parameter "${paramMapping[i].name}" has no effect on the fit and will not move.` });
                     }
                 }

                 // Marquardt diagonal scaling: A = JtJ + lambda * diag(JtJ).
                 // This replaces the old scheme where parameters were
                 // pre-normalized by their initial magnitudes; diag(JtJ)
                 // carries the correct CURRENT local curvature for each
                 // parameter, so damping is applied consistently whether
                 // a parameter has moved a lot or a little from its start.
                 // solveLinearSystem consumes its matrix, so hand it a copy and
                 // keep finalJtJ intact for the main thread's ESD calculation.
                 const A_lm = finalJtJ.map(row => Float64Array.from(row));
                 for (let i = 0; i < n_params; i++) {
                     const d = finalJtJ[i][i];
                     const diag = (isFinite(d) && d > 1e-30) ? d : 1e-6;
                     A_lm[i][i] += lambda * diag;
                 }

                 const p_step = solveLinearSystem(A_lm, Jtr);
                 if (!p_step) {
                      console.warn(`LM iter ${iter}: normal equations are singular. Increasing lambda.`);
                      lambda = Math.min(1e9, lambda * 10);
                      continue;
                 }

                 if (p_step.some(v => !isFinite(v))) {
                     console.warn(`LM iter ${iter}: Step calculation resulted in NaN/Infinity. Increasing lambda.`);
                     lambda = Math.min(1e9, lambda * 5);
                     continue;
                 }

                 // Snapshot BEFORE applying the trial step, so we can revert on rejection.
                 oldParams = JSON.parse(JSON.stringify(params));
                 oldHklList = JSON.parse(JSON.stringify(workingHklList));

                 // Apply trial step in raw parameter space.
                 const p_current = paramMapping.map(m => m.get(params, workingHklList));
                 
                 // ===============================================================
                 // Trust region on the LM step.
                 //
                 // CRITICAL FIX: the old code clipped each component of p_step
                 // INDEPENDENTLY. The Gauss-Newton step is a joint solution of
                 // coupled normal equations, so truncating one component while
                 // leaving the others intact rotates the vector -- and a rotated
                 // Gauss-Newton step is no longer guaranteed to point downhill.
                 // Observed directly: refining `a` together with `zeroShift`
                 // (strongly anti-correlated, both move peak positions), the
                 // solve asked for da=+1.9e-3 with dz=+4.1e-2; dz was clipped to
                 // +4.0e-3 while da passed through untouched, and the resulting
                 // step went UPHILL. Every iteration was rejected and the fit sat
                 // still while lambda climbed. Refining widths alone, where
                 // nothing clipped, converged normally.
                 //
                 // Scaling the whole vector by one factor keeps the direction and
                 // only shortens the step, which is what a trust region is meant
                 // to do.
                 // ===============================================================
                 let stepScale = 1.0;
                 for (let idx = 0; idx < p_step.length; idx++) {
                     const m = paramMapping[idx];
                     // `scale` is |initial value|, which collapses for parameters
                     // that legitimately start near zero (zeroShift, GV, the
                     // Stephens terms). Fall back to the declared default so such
                     // a parameter is not frozen by an absurdly small budget.
                     const magnitude = Math.max(Math.abs(m.scale || 0), m.defaultScale || 0, 1e-6);
                     const maxAllowedStep = (m.step || 0.2) * 4.0 * magnitude;
                     const want = Math.abs(p_step[idx]);
                     if (want > maxAllowedStep && want > 0) {
                         stepScale = Math.min(stepScale, maxAllowedStep / want);
                     }
                 }
                 const p_step_clamped = p_step.map(v => v * stepScale);

                 const p_new = p_current.map((v, i) => v + p_step_clamped[i]);




                 paramMapping.forEach((m, i) => m.set(params, workingHklList, p_new[i]));

                 calculateTotalPattern(y_calc_total);
                 let new_cost = 0;
                 for (let i = 0; i < n_points; i++) {
                     const res = (y_obs[i] - y_calc_total[i]) * sqrt_weights[i];
                      if (isFinite(res)) {
                          new_cost += res * res;
                      } else {
                          new_cost = Infinity;
                          break;
                      }
                 }

                 let stepAccepted = false;
                 if (new_cost < cost && isFinite(new_cost)) {
                     lambda = Math.max(1e-9, lambda / 3);
                     stepAccepted = true;
                 } else {
                     // Reject: restore previous state. Do NOT update lastAcceptedCost.
                     params = oldParams;
                     workingHklList = oldHklList;
                     lambda = Math.min(1e9, lambda * 2);
                 }

                 // Throttle: only post when the whole-percent value changes.
                 // Posting every iteration floods the main thread's event loop
                 // with messages it can only render once per frame anyway.
                 const overallProgress = Math.min(1.0, baseProgress + ((iter + 1) / maxIter) * cycleProgressSpan);
                 const pct = Math.floor(overallProgress * 100);
                 if (pct !== lastProgressPct) {
                     lastProgressPct = pct;
                     postMessage({ type: 'progress', value: overallProgress });
                 }

                 // Convergence test: only compares ACCEPTED costs, so a
                 // rejected step can never look like "converged" here.
                 if (stepAccepted) {
                     const denom = Math.max(Math.abs(lastAcceptedCost), 1.0);
                     if (Math.abs(lastAcceptedCost - new_cost) < 1e-9 * denom) {
                         lastAcceptedCost = new_cost;
                         ss_res = new_cost;
                         break;
                     }
                     lastAcceptedCost = new_cost;
                     ss_res = new_cost;
                 }


             } catch (error) {
                  console.error("Error during LM iteration:", iter, error);
                  postMessage({ type: 'error', message: `Error in LM iter ${iter}: ${error.message}` });
                  // Build parameterInfo locally — the outer-scope one isn't in scope yet.
                  const errorParamInfo = paramMapping.map(m => ({
                      name: m.name,
                      scale: 1.0,
                      typicalMagnitude: m.scale,
                      isIntensity: m.isIntensity
                  }));
                  return { params: oldParams || initialParams, hklList: oldHklList || JSON.parse(JSON.stringify(hklList)), ss_res: ss_res, error: true, JtJ: finalJtJ, parameterInfo: errorParamInfo, algorithm: 'lm', fitFlags };
             }

        }

         const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
         postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });

          // LM now operates in raw parameter space, so JtJ, Jtr, and the
          // resulting covariance matrix are all in raw units. The main
          // thread multiplies reported sigma by `scale` to un-normalize, so
          // we report scale=1.0 here to keep that code path correct without
          // touching it. We still include `.typicalMagnitude` in case the UI
          // wants to display a "typical value" for the parameter.
          const parameterInfoForMainThread = paramMapping.map(m => ({
               name: m.name,
               scale: 1.0,
               typicalMagnitude: m.scale,
               isIntensity: m.isIntensity
          }));


        return {
             params,
             hklList: workingHklList,
             JtJ: finalJtJ,
             parameterInfo: parameterInfoForMainThread,
             ss_res,
             algorithm: 'lm',
             fitFlags
        };
}

//  
async function refineParametersPT(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode);
    if (paramMapping.length === 0) {
         const finalProgress = (leBailCycle + 1) / totalLeBailCycles;
         postMessage({ type: 'progress', value: Math.min(1.0, finalProgress) });
        return { params: initialParams, hklList: hklList, ss_res: 0 };
    }

    const numReplicas = 8;
    const maxTemp = 1.0;
    const minTemp = 1e-5;
    const swapInterval = 10;


    const n_points = workerWorkingData.tth.length;
    const scratch_bkg = new Float64Array(n_points);
    const scratch_pattern = new Float64Array(n_points);

    // Le Bail intensities are extracted once, below, and held fixed for the
    // whole run -- see the long note in refineParametersLM for why
    // re-extracting inside the objective breaks width refinement.
    if (refinementMode === 'le-bail') {
        enforceSymmetryConstraintsWorker(initialParams, system);
        updateHklPositions(hklList, initialParams, system);
        const bkg0 = calculateTotalBackground(workerWorkingData.tth, initialParams, workerBackgroundAnchors);
        leBailIntensityExtraction(
            { tth: workerWorkingData.tth, intensity: workerWorkingData.intensity, background: bkg0 },
            hklList,
            initialParams
        );
    }

    const objective = (p_obj, hkl_list_obj) => {
         try {
            enforceSymmetryConstraintsWorker(p_obj, system); 
            updateHklPositions(hkl_list_obj, p_obj, system);
            const y_bkg = calculateTotalBackground(workerWorkingData.tth, p_obj, workerBackgroundAnchors, scratch_bkg);
            const netCalcPattern = calculatePattern(workerWorkingData.tth, hkl_list_obj, p_obj, scratch_pattern);
             let sum_w_res_sq = 0;
             for (let i = 0; i < workerWorkingData.tth.length; i++) {
                  const w_i = workerWorkingData.weights[i];
                  const res = workerWorkingData.intensity[i] - (netCalcPattern[i] + y_bkg[i]);
                   if (isFinite(w_i) && isFinite(res)) {
                        sum_w_res_sq += w_i * res * res;
                   } else {
                       return 1e12;
                   }
             }
             return isFinite(sum_w_res_sq) ? sum_w_res_sq : 1e12;
         } catch (err) {
              console.warn("PT objective function error:", err);
              return 1e12;
         }
    };


    const temperatures = Array.from({ length: numReplicas }, (_, i) =>
        maxTemp * Math.pow(minTemp / maxTemp, i / (numReplicas - 1 || 1))
    );

    const initialCost = objective(initialParams, hklList);

    // FIX: one reference energy scale shared by BOTH acceptance tests.
    // Previously the local Metropolis step divided by the replica's own
    // CURRENT cost (relative scale) while the replica-exchange test used the
    // raw cost difference (absolute scale). With chi-square values of 1e5-1e7
    // counts the swap exponent saturated, so exchanges were accepted either
    // always or never and the run degenerated into 8 independent annealers
    // rather than parallel tempering. Using a fixed scale also keeps the
    // local acceptance ratio well defined (a cost-dependent denominator
    // breaks detailed balance).
    const costScale = Math.max(1e-9, Math.abs(initialCost));

    let replicas = temperatures.map(temp => ({
        params: JSON.parse(JSON.stringify(initialParams)),
        hklList: JSON.parse(JSON.stringify(hklList)),
        cost: initialCost,
        temp: temp
    }));

    let bestOverallParams = JSON.parse(JSON.stringify(initialParams));
    let bestOverallHklList = JSON.parse(JSON.stringify(hklList));
    let bestOverallCost = initialCost;

    const baseProgress = leBailCycle / totalLeBailCycles;
    const cycleProgressSpan = 1 / totalLeBailCycles;
    let lastProgressPct = -1;
    let consecutiveErrors = 0;

    for (let iter = 0; iter < maxIter; iter++) {
         try {
             for (let i = 0; i < numReplicas; i++) {
                 let replica = replicas[i];

                 const p_idx = Math.floor(Math.random() * paramMapping.length);
                 const mapping = paramMapping[p_idx];
                 const original_val = mapping.get(replica.params, replica.hklList);

                 // PERF FIX: a move perturbs exactly ONE scalar, so rolling it
                 // back only needs that scalar plus the Le Bail intensities
                 // (everything else objective() touches -- symmetry-slaved
                 // params, peak positions -- is recomputed from scratch on the
                 // next call). The old code deep-cloned params AND the whole
                 // hklList via JSON twice per replica per iteration: 16 full
                 // JSON round-trips of every reflection, every iteration.
                 // In Pawley mode the perturbed scalar may itself be an
                 // intensity, but mapping.set restores it; nothing else in
                 // objective() mutates intensities any more.
                 const step_scale = Math.max(0.01, replica.temp);
                 const step_size = (mapping.step || 0.05) * step_scale * mapping.scale;
                 const random_step = (Math.random() - 0.5) * 2 * step_size;
                 mapping.set(replica.params, replica.hklList, original_val + random_step);

                 const neighbor_cost = objective(replica.params, replica.hklList);
                 const delta_cost = neighbor_cost - replica.cost;

                 const acceptance_prob = (replica.temp > 1e-9)
                     ? Math.exp(Math.max(-700, -delta_cost / (costScale * replica.temp)))
                     : 0;
                 if (delta_cost < 0 || acceptance_prob > Math.random()) {
                     replica.cost = neighbor_cost;
                 } else {
                     mapping.set(replica.params, replica.hklList, original_val);
                 }

                 if (replica.cost < bestOverallCost) {
                     bestOverallCost = replica.cost;
                     bestOverallParams = JSON.parse(JSON.stringify(replica.params));
                     bestOverallHklList = JSON.parse(JSON.stringify(replica.hklList));
                 }
             }

             if (iter > 0 && iter % swapInterval === 0) {
                 for (let i = 0; i < numReplicas - 1; i++) {
                     const rep1 = replicas[i];
                     const rep2 = replicas[i + 1];

                     if (Math.abs(rep1.temp - rep2.temp) < 1e-9) continue;

                     const delta_beta = (1 / rep1.temp) - (1 / rep2.temp);
                     const delta_cost = (rep1.cost - rep2.cost) / costScale;   // same scale as the local step
                     const acceptance_prob_swap = Math.exp(Math.max(-700, Math.min(50, delta_beta * delta_cost)));

                     if (acceptance_prob_swap > Math.random()) {
                         [rep1.params, rep2.params] = [rep2.params, rep1.params];
                         [rep1.hklList, rep2.hklList] = [rep2.hklList, rep1.hklList];
                         [rep1.cost, rep2.cost] = [rep2.cost, rep1.cost];
                     }
                 }
             }

             // Throttled, as in LM: only post on whole-percent changes.
             const overallProgress = Math.min(1.0, baseProgress + ((iter + 1) / maxIter) * cycleProgressSpan);
             const pct = Math.floor(overallProgress * 100);
             if (pct !== lastProgressPct) {
                 lastProgressPct = pct;
                 postMessage({ type: 'progress', value: overallProgress });
             }
             consecutiveErrors = 0;

         } catch (error) {
              // FIX: the old handler posted an error EVERY iteration, so a
              // persistent fault produced maxIter error toasts and kept
              // burning CPU. Give up after a few consecutive failures.
              console.error("Error during PT iteration:", iter, error);
              consecutiveErrors++;
              if (consecutiveErrors >= 5) {
                  postMessage({ type: 'error', message: `PT aborted at iteration ${iter}: ${error.message}` });
                  break;
              }
         }

    }

     updateHklPositions(bestOverallHklList, bestOverallParams, system);
     // After the per-peak Le Bail extraction inside the objective, the
     // intensities in bestOverallHklList are already correct. No final
     // global rescaling is needed or appropriate.

     const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
     postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });

      const parameterInfoForMainThread = paramMapping.map(m => ({
           name: m.name,
           scale: m.scale,
           isIntensity: m.isIntensity
      }));


    return {
        params: bestOverallParams,
        hklList: bestOverallHklList,
        algorithm: 'pt',
        parameterInfo: parameterInfoForMainThread,
        fitFlags,
        ss_res: bestOverallCost
    };
}


//   Parameter Mapping (Helper)  
function getParameterMapping(fitFlags, initialParams, hklList, refinementMode) {
    const mappings = [];

    // NOTE (fix, April 2026):
    //   get/set now operate in RAW parameter space. Previously they
    //   normalized by a `scale` frozen at the INITIAL value, which caused
    //   LM to size its steps wrong once a parameter had moved appreciably
    //   from its starting point. The LM routine now uses Marquardt's own
    //   diagonal scaling (lambda * diag(JtJ)), which is the correct
    //   curvature-aware scaling.
    //
    //   `scale` is still exposed on the mapping because:
    //     - PT uses it to size random proposal steps (typical magnitude).
    //     - The main thread previously multiplied ESDs by `scale` to
    //       un-normalize. Since LM is now raw-space, that multiplication
    //       is a no-op when `scale === 1.0`. We therefore report
    //       `scale: 1.0` in the LM `parameterInfo` so existing main-thread
    //       code keeps working unchanged. PT still reports its typical
    //       magnitude for diagnostic display.
    const createMapping = (flag, name, defaultScale = 1.0, minVal = -Infinity, maxVal = Infinity, step = 0.2) => {
        if (!flag) return null;
        if (!initialParams || !(name in initialParams) || !isFinite(initialParams[name])) {
            // Refining a parameter the model never received (or received as
            // NaN) can only produce a dead column and a bogus ESD. Drop it and
            // tell the user rather than failing silently.
            console.warn(`Refinement flag set for "${name}" but no finite initial value was supplied; skipping.`);
            postMessage({ type: 'warning', message: `Parameter "${name}" has no valid value and was excluded from the refinement.` });
            return null;
        }
        const initialValue = initialParams[name] ?? 0;
        const scale = Math.abs(initialValue) > 1e-9 ? Math.abs(initialValue) : defaultScale;

        return {
            name: name,
            scale: scale,          // typical magnitude, for PT step-sizing & display
            defaultScale: defaultScale,  // floor for the trust-region budget when
                                         // the parameter starts at (or near) zero
            step: step,
            isIntensity: false,
            minVal: minVal,
            maxVal: maxVal,
            get: (p_obj, hkl_list_obj) => (p_obj[name] ?? 0),
            // FIX: the old body only assigned `if (p_obj.hasOwnProperty(name))`.
            // When the key was absent, get() still returned 0 via `?? 0`, so the
            // parameter got a mapping slot, an identically-zero Jacobian column
            // and a consumed degree of freedom, while being physically unable
            // to move. Assign unconditionally.
            set: (p_obj, hkl_list_obj, rawValue) => {
                if (!p_obj) return;
                if (!(rawValue > minVal)) rawValue = Math.max(rawValue, minVal);
                if (!(rawValue < maxVal)) rawValue = Math.min(rawValue, maxVal);
                if (!isFinite(rawValue)) return;
                p_obj[name] = rawValue;
            }
        };
    };


    //   Pawley Intensity Parameters  
    // Raw-space (see NOTE above createMapping).
     if (refinementMode === 'pawley' && hklList && workerWorkingData && workerWorkingData.tth && workerWorkingData.tth.length > 0) {
          const tthMin = workerWorkingData.tth[0];
          const tthMax = workerWorkingData.tth[workerWorkingData.tth.length - 1];

          hklList.forEach((hkl, index) => {
               if (hkl && hkl.tth && hkl.tth >= tthMin && hkl.tth <= tthMax) {
                    const hkl_name = `I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`;
                    const initialIntensity = (hkl.intensity !== undefined && hkl.intensity > 1e-6) ? hkl.intensity : 1000.0;
                    const scale = initialIntensity;  // typical magnitude, used by PT

                    mappings.push({
                        name: hkl_name,
                        scale: scale,
                        defaultScale: 1000.0,
                        step: 0.3,
                        isIntensity: true,
                        index: index,
                        minVal: 0,
                        maxVal: Infinity,
                        get: (p_obj, hkl_list_obj) => {
                            return hkl_list_obj?.[index]?.intensity ?? 0;
                        },
                        set: (p_obj, hkl_list_obj, rawValue) => {
                             if (hkl_list_obj?.[index]) {
                                hkl_list_obj[index].intensity = Math.max(0, rawValue);
                             }
                        }
                    });
               }
          });
     }


    //   Other Parameters  
    const profileType = String(initialParams.profileType || "simple_pvoigt");

    // ---------------------------------------------------------------------
    // Data-derived lower bounds for the width parameters.
    //
    // A peak narrower than a couple of data steps is not resolvable and is not
    // physically meaningful, but nothing stopped the refiner from driving the
    // widths there. When the starting cell is poor, shrinking the peaks to
    // spikes genuinely lowers the weighted residual (a spike in the wrong place
    // costs less than a broad peak in the wrong place), so the fit would
    // collapse the profile instead of correcting the cell -- and once the
    // widths hit the internal 1e-4 deg clamp inside calculateProfileWidths the
    // Jacobian column went dead and they could never recover.
    //
    // Floor them at a FWHM of two data steps, expressed back in GSAS units.
    let GW_MIN = 1e-6, LX_MIN = 1e-6;
    if (workerWorkingData && workerWorkingData.tth && workerWorkingData.tth.length > 2) {
        const t = workerWorkingData.tth;
        const stepDeg = (t[t.length - 1] - t[0]) / (t.length - 1);
        const minFwhmDeg = Math.max(1e-3, 2 * stepDeg);
        // Gaussian: FWHM[deg] = sqrt(sigma_sq) * GSAS_GAUSSIAN_TO_DEG
        GW_MIN = Math.pow(minFwhmDeg / GSAS_GAUSSIAN_TO_DEG, 2);
        // Lorentzian: FWHM[deg] = (LX / cos(theta)) * GSAS_LORENTZIAN_TO_DEG.
        // Evaluate at the middle of the scanned range.
        const midTheta = ((t[0] + t[t.length - 1]) / 2) * (Math.PI / 180) / 2;
        const cosMid = Math.max(0.1, Math.cos(midTheta));
        LX_MIN = minFwhmDeg * cosMid / GSAS_LORENTZIAN_TO_DEG;
    }

    mappings.push(createMapping(fitFlags.a, 'a', 4.0, 0.1, Infinity, 0.01));
    mappings.push(createMapping(fitFlags.b, 'b', 4.0, 0.1, Infinity, 0.01));
    mappings.push(createMapping(fitFlags.c, 'c', 6.0, 0.1, Infinity, 0.01));
    mappings.push(createMapping(fitFlags.alpha, 'alpha', 90.0, 0, 180, 0.05));
    mappings.push(createMapping(fitFlags.beta, 'beta', 90.0, 0, 180, 0.05));
    mappings.push(createMapping(fitFlags.gamma, 'gamma', 120.0, 0, 180, 0.05));
    mappings.push(createMapping(fitFlags.zeroShift, 'zeroShift', 0.01, -Infinity, Infinity, 0.1));

    if (profileType === "simple_pvoigt") {
        mappings.push(createMapping(fitFlags.GU, 'GU', 0.01, 0, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GV, 'GV', 0.01, -Infinity, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GW, 'GW', 0.01, GW_MIN, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GP, 'GP', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.LX, 'LX', 0.01, LX_MIN, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.eta, 'eta', 0.5, 0, 1, 0.1));
        mappings.push(createMapping(fitFlags.shft, 'shft', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns, 'trns', 0.01, -Infinity, Infinity, 0.1));
    } 
    else if (profileType === "tch_aniso") {
        mappings.push(createMapping(fitFlags.U, 'U', 0.01, 0, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.V, 'V', 0.01, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.W, 'W', 0.01, GW_MIN, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.X, 'X', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.Y, 'Y', 0.01, LX_MIN, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.SL, 'SL', 0.001, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.HL, 'HL', 0.001, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.S400, 'S400', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S040, 'S040', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S004, 'S004', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S220, 'S220', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S202, 'S202', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S022, 'S022', 0.1, -Infinity, Infinity, 0.2));
    }
    else if (profileType === "split_pvoigt") {
        mappings.push(createMapping(fitFlags.GU_L, 'GU_L', 0.01, 0, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GV_L, 'GV_L', 0.01, -Infinity, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GW_L, 'GW_L', 0.01, GW_MIN, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.LX_L, 'LX_L', 0.01, LX_MIN, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.GU_R, 'GU_R', 0.01, 0, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GV_R, 'GV_R', 0.01, -Infinity, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GW_R, 'GW_R', 0.01, GW_MIN, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.LX_R, 'LX_R', 0.01, LX_MIN, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.eta_split, 'eta_split', 0.5, 0, 1, 0.1));
        mappings.push(createMapping(fitFlags.shft_split, 'shft_split', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns_split, 'trns_split', 0.01, -Infinity, Infinity, 0.1));
    }
   
    const paramMapping = mappings.filter(Boolean);
    return { paramMapping };
}


//   4. Worker Message Handler  
self.onmessage = async function(e) {
    // FIX: wait for async init (SG database + mathjs + sg_engine). Without
    // this, a fit posted before init completed crashed with
    // "SG_ENGINE is not defined".
    try {
        await initPromise;
    } catch (initErr) {
        postMessage({ type: 'error', message: `Worker initialization failed: ${initErr && initErr.message ? initErr.message : initErr}` });
        return;
    }

    const {
        initialParams,
        fitFlags,
        workingData, // contains tth, intensity, weights, startIndex
        masterHklList,
        selectedSgNumber,
        selectedSgQuery, // original user input string, for setting preservation
        system,
        maxIterations,
        algorithm,
        refinementMode,
        backgroundAnchors 
    } = e.data;

     // Store working data globally in the worker
     workerWorkingData = workingData;
     workerBackgroundAnchors = backgroundAnchors; 
     // Resolve the space group via the shared engine. Prefer the original
     // user query (which preserves any non-standard setting, e.g. "P1211");
     // fall back to the numeric lookup.
     const selectedSg = SG_ENGINE.resolve(selectedSgQuery) || SG_ENGINE.resolve(selectedSgNumber);
     if (!selectedSg) {
         postMessage({ type: 'error', message: `Worker Error: Space group ${selectedSgNumber} could not be resolved.` });
         return;
     }


    let finalResults;
    let currentHklList = JSON.parse(JSON.stringify(masterHklList || [])); // Start with master list
    let currentParams = initialParams; // Use initial params for both modes initially
     currentParams.system = system;
    try {
        //   Le Bail Mode  
        // Extract-then-refine, repeated. Each cycle re-extracts the peak
        // intensities from the observed profile using the parameters the
        // previous cycle arrived at, then refines with those intensities held
        // fixed. The outer loop is what lets the intensities catch up; holding
        // them fixed *within* a cycle is what keeps the width parameters
        // determinate (see the note in refineParametersLM).
        if (refinementMode === 'le-bail') {
            // Textbook Le Bail start: every reflection gets the SAME intensity.
            // The method's whole point is that no structural model is assumed,
            // so the first decomposition partitions purely on profile overlap.
            // Carrying intensities over from a previous run (or from imported
            // structure factors) biases the result and makes it depend on
            // history; starting flat makes a given dataset give a given answer.
            currentHklList.forEach(peak => { if (peak) peak.intensity = 100.0; });

            // Then iterate the decomposition alone a few times so the
            // intensities settle before any parameter is allowed to move. The
            // partitioning is a fixed-point iteration and converges in a
            // handful of passes.
            try {
                const bkgSeed = calculateTotalBackground(workerWorkingData.tth, currentParams, workerBackgroundAnchors);
                enforceSymmetryConstraintsWorker(currentParams, system);
                updateHklPositions(currentHklList, currentParams, system);
                for (let seedCycle = 0; seedCycle < 5; seedCycle++) {
                    leBailIntensityExtraction(
                        { tth: workerWorkingData.tth,
                          intensity: workerWorkingData.intensity,
                          background: bkgSeed },
                        currentHklList, currentParams);
                }
            } catch (seedErr) {
                console.warn("Le Bail seeding pass failed; continuing from flat intensities.", seedErr);
            }

            // refineParametersLM then re-extracts at the top of every iteration
            // (refreshLeBailIntensities), so the extract/least-squares
            // alternation happens inside the iteration loop.
            if (algorithm === 'lm') {
                finalResults = await refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { // pt
                finalResults = await refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            }
            if (!finalResults || finalResults.error) {
                throw new Error(`Refinement algorithm (${algorithm}) failed during Le Bail fit. ${finalResults?.error ? 'See previous error.' : ''}`);
            }
        }
        //   Pawley Mode  
        else { // refinementMode === 'pawley'
             // Perform a SINGLE Le Bail extraction first to initialize intensities
             postMessage({ type: 'progress', value: 0.05, message: 'Initializing Pawley intensities...' }); // Small progress update
             try {
                // Ensure HKL positions are correct with initial params before extraction
                updateHklPositions(currentHklList, currentParams, system);

                const backgroundForInit = calculateTotalBackground(workerWorkingData.tth, currentParams, workerBackgroundAnchors);
                const expDataForInit = {
                    tth: workerWorkingData.tth,
                    intensity: workerWorkingData.intensity,
                    background: backgroundForInit
                };
                 leBailIntensityExtraction(expDataForInit, currentHklList, currentParams);
                 // currentHklList now contains initial intensity estimates (as heights)
                 console.log("Pawley intensities initialized via Le Bail extraction.");
             } catch (initError) {
                  console.warn("Could not initialize Pawley intensities using Le Bail extraction:", initError);
                  // Fallback: Initialize with 1000 if extraction fails
                  currentHklList.forEach(peak => {
                       if (peak) peak.intensity = 1000.0;
                  });
             }
            // Now run the chosen algorithm ONCE for Pawley, refining initialized intensities simultaneously
            if (algorithm === 'lm') {
                finalResults = await refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { // pt
                finalResults = await refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            }

             if (!finalResults || finalResults.error) {
                  throw new Error(`Refinement algorithm (${algorithm}) failed during Pawley fit. ${finalResults?.error ? 'See previous error.' : ''}`);
             }

        } // End Pawley Mode

        //   Final Calculations & Posting Result  
        if (!finalResults || !finalResults.params || !finalResults.hklList) {
            throw new Error("Refinement finished but produced invalid results.");
        }

        //   (Rest of the onmessage function remains the same: calculating stats, sending results)  
        const finalParams = finalResults.params;
        const finalHklList = finalResults.hklList;

        const finalNetPatternForStats = calculatePattern(workerWorkingData.tth, finalHklList, finalParams);
        const finalBackgroundForStats = calculateTotalBackground(workerWorkingData.tth, finalParams, workerBackgroundAnchors);

        const finalStats = calculateStatistics(
            workerWorkingData,
            finalNetPatternForStats,
            finalResults.fitFlags || fitFlags,
            finalBackgroundForStats,
            finalParams,
            finalHklList,
            refinementMode
        );

        const resultPayload = {
            params: finalParams,
            hklList: finalHklList,
            stats: finalStats,
            algorithm: algorithm,
            refinementMode: refinementMode,
            fitFlags: finalResults.fitFlags || fitFlags,
            parameterInfo: finalResults.parameterInfo || [],
            JtJ: finalResults.JtJ || null,
            ss_res: finalResults.ss_res
        };

        postMessage({ type: 'result', results: resultPayload });


    } catch (error) {
        console.error("Worker Error during refinement:", error);
        postMessage({ type: 'error', message: `Refinement failed: ${error.message}` });
    } finally {

         workerWorkingData = null;
         workerBackgroundAnchors = [];
    }
};

//   error handler  
self.onerror = function(event) {
     console.error("Unhandled Worker Error:", event.message, event);
     postMessage({ type: 'error', message: `Unhandled Worker Error: ${event.message}` });
};