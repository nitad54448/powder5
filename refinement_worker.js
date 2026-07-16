// refinement_worker.js
// version 122, 22 april 202, 
// Spline Background in nov 2025, version 115
const initPromise = fetch('cctbx_space_groups_all_settings_v2.json')
    .then(res => {
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        return res.json();
    })
    .then(data => { 
        globalThis.SG_DATABASE = data; 
        importScripts('https://cdnjs.cloudflare.com/ajax/libs/mathjs/12.4.3/math.min.js', 'sg_engine.js');
        
        if (typeof math === 'undefined') {
            throw new Error('math.js did not load correctly.');
        }
        if (typeof SG_ENGINE === 'undefined') {
            throw new Error('sg_engine.js did not load correctly.');
        }
        console.log("Worker dependencies loaded successfully.");
    })
    .catch(err => {
        console.error("Worker Error:", err);
        postMessage({ type: 'error', message: `Worker initialization failed: ${err.message}` });
        self.close();
    });



// Engine shims so legacy call sites inside this file still work. The
// authoritative implementations live in sg_engine.js.
function isReflectionAllowed(h, k, l, spaceGroup) {
    return SG_ENGINE.isReflectionAllowed(h, k, l, spaceGroup);
}

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

//   2. Define Constants & Global Worker State  
const CALCULATION_WINDOW_MULTIPLIER = 8.0;
const PEAK_HEIGHT_CUTOFF = 0.002;
const HIGH_WEIGHT_MULTIPLIER = 50.0;
let workerWorkingData = null; // To store the sliced data sent from the main thread
let hklIndexCache = {}; 
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
        const interpolatedY = h00 * y_i + h10 * h_i * m_i + h01 * y_ip1 + h11 * h_i * m_ip1;

         // Extrapolation: use flat extrapolation outside the defined range
         if (tthValue < points[0].tth) return points[0].y;
         if (tthValue > points[n-1].tth) return points[n-1].y;

        return interpolatedY;
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


// Engine shim: delegate multiplicity to sg_engine.js so main thread
// and worker share one authoritative implementation.
function getMultiplicityAndCanonicalHKL(h, k, l, laue_class) {
    return SG_ENGINE.getMultiplicity(h, k, l, laue_class);
}


/**
 * Generates the list of raw HKL indices {h, k, l, multiplicity} for PREVIEW/REFINEMENT.
 * Checks cache first in the worker context.
 */
function generateAndCacheHklIndices(spaceGroup, maxTth, params) {
    let sgNumber = null;
    if (spaceGroup && typeof spaceGroup.number === 'number') {
        sgNumber = spaceGroup.number;
        if (typeof self !== 'undefined' && typeof importScripts === 'function' && hklIndexCache[sgNumber]) {
             return hklIndexCache[sgNumber];
        }
    }

    const { a, b, c, lambda } = params;
    const sgSystem = spaceGroup ? spaceGroup.system : null;
    const sgLaueClass = spaceGroup ? spaceGroup.laue_class : null;

    if (!lambda || typeof lambda !== 'number' || lambda <= 0 ||
        !sgLaueClass || typeof sgLaueClass !== 'string' || sgLaueClass.length === 0 ||
        !a || typeof a !== 'number' || a <= 0 ||
        !sgSystem || typeof sgSystem !== 'string' || sgSystem.length === 0 ) {
         console.error(`HKL Gen Error: Invalid parameters provided. Lambda: ${lambda}, Laue Class: ${sgLaueClass}, System: ${sgSystem}, a: ${a}`);
         return [];
     }

    const maxDim = Math.max(a || 0, b || a || 0, c || a || 0);
    if (maxDim <= 0) return [];

    const sinThetaMax = Math.sin(maxTth * Math.PI / 360);
    if (sinThetaMax <= 0) return [];
    const dMin = lambda / (2 * sinThetaMax);
    if (dMin <= 0) return [];
    const inv_d_sq_max = 1.0 / (dMin * dMin);
    const baseMaxIndex = Math.ceil(maxDim / dMin) + 5;
    const maxIndex = Math.min(baseMaxIndex, 50);

    const a_sq = a * a;
    const b_sq = (sgSystem === 'orthorhombic' || sgSystem === 'monoclinic' || sgSystem === 'triclinic') ? ((b && b > 0) ? b * b : a_sq) : a_sq;
    const c_sq = (sgSystem !== 'cubic') ? ((c && c > 0) ? c * c : a_sq) : a_sq;

    let rawReflections = [];
    const addedHKLs = new Set();
    const getKey = (h, k, l) => `${h},${k},${l}`;

    const loopAndAdd = (h, k, l) => {
        if (h === 0 && k === 0 && l === 0) return;
        const key = getKey(h, k, l);
        if (addedHKLs.has(key)) return;

        if (isReflectionAllowed(h, k, l, spaceGroup)) {
            const { multiplicity } = getMultiplicityAndCanonicalHKL(h, k, l, sgLaueClass);
             if (multiplicity > 0) {
                 rawReflections.push({
                    h_orig: h, k_orig: k, l_orig: l,
                    hkl_list: [`(${h},${k},${l})`],
                    multiplicity: multiplicity
                });
                addedHKLs.add(key);
             }
        }
    };

    const maxI = maxIndex;

    if (sgSystem === 'monoclinic' || sgSystem === 'triclinic') {
        for (let h = -maxI; h <= maxI; h++) {
            for (let k = 0; k <= maxI; k++) {
                for (let l = 0; l <= maxI; l++) {
                     if (k === 0 && l === 0 && h <= 0) continue;
                     if (k === 0 && h === 0 && l === 0) continue;
                    loopAndAdd(h, k, l);
                }
            }
        }
    } else if (sgSystem === 'orthorhombic') {
        if (a_sq <= 0 || b_sq <= 0 || c_sq <= 0) return [];
        for (let h = 0; h <= maxI; h++) {
            const h_term = (h*h) / a_sq;
            if (h_term > inv_d_sq_max && h > 0) break;

            for (let k = 0; k <= maxI; k++) {
                const hk_term = h_term + (k*k) / b_sq;
                if (hk_term > inv_d_sq_max && k > 0) break;

                for (let l = 0; l <= maxI; l++) {
                    if (h === 0 && k === 0 && l === 0) continue;
                    const hkl_term = hk_term + (l*l) / c_sq;
                    if (hkl_term > inv_d_sq_max) {
                         break;
                    }
                    loopAndAdd(h, k, l);
                }
            }
        }
    } else if (sgSystem === 'hexagonal' || sgSystem === 'trigonal' || sgSystem === 'rhombohedral') {
        if (a_sq <= 0 || c_sq <= 0) return [];
        const a_term_prefactor = 4.0 / (3.0 * a_sq);
        const l_limit = (sgLaueClass === '6/m' || sgLaueClass === '6/mmm') ? 0 : -maxI;

        for (let h = 0; h <= maxI; h++) {
            const h_term_only = a_term_prefactor * (h*h);
             if (h_term_only > inv_d_sq_max && h > 0) break;

            for (let k = 0; k <= h; k++) {
                const hk_term = a_term_prefactor * (h*h + h*k + k*k);
                 if (hk_term > inv_d_sq_max && k > 0) break;

                for (let l = l_limit; l <= maxI; l++) {
                    if (h === 0 && k === 0 && l === 0) continue;
                    const hkl_term = hk_term + (l*l) / c_sq;
                    if (hkl_term > inv_d_sq_max) {
                        if (l >= 0) break;
                    }
                     if (h === 0 && k === 0 && l === 0) continue;
                    loopAndAdd(h, k, l);
                }
            }
        }
    } else { // Cubic, Tetragonal
        const isTetragonal = (sgSystem === 'tetragonal');
        if (a_sq <= 0 || (isTetragonal && c_sq <= 0)) return [];

        for (let h = 0; h <= maxI; h++) {
            const h_term_factor = (h*h);

            for (let k = 0; k <= h; k++) {
                const hk_term_factor = h_term_factor + (k*k);

                for (let l = 0; l <= k; l++) {
                    if (h === 0 && k === 0 && l === 0) continue;

                    let inv_d_sq_term;
                    if (isTetragonal) {
                        inv_d_sq_term = hk_term_factor / a_sq + (l*l) / c_sq;
                    } else { // Cubic
                        inv_d_sq_term = (hk_term_factor + l*l) / a_sq;
                    }

                    if (inv_d_sq_term > inv_d_sq_max) {
                         break;
                    }
                    loopAndAdd(h, k, l);
                }
                 let l0_term;
                 if(isTetragonal) l0_term = hk_term_factor / a_sq;
                 else l0_term = hk_term_factor / a_sq;
                 if (l0_term > inv_d_sq_max && k > 0) {
                     break;
                 }

            }
             let k0l0_term;
             if(isTetragonal) k0l0_term = h_term_factor / a_sq;
             else k0l0_term = h_term_factor / a_sq;
             if (k0l0_term > inv_d_sq_max && h > 0) {
                 break;
             }
        }
    }

    if (sgNumber !== null && typeof self !== 'undefined' && typeof importScripts === 'function') {
         hklIndexCache[sgNumber] = rawReflections;
    }
    return rawReflections;
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

/**
 * Calculates area for simple_pvoigt
 */
function _getArea_Simple(tth_peak, hkl, params) {
    const GAUSS_AREA_CONST = 1.0644677;
    const LORENTZ_AREA_CONST = 1.5707963;

    const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
    const gG = Math.max(1e-9, gamma_G);
    const gL = Math.max(1e-9, gamma_L);
    const area_G = gG * GAUSS_AREA_CONST;
    const area_L = gL * LORENTZ_AREA_CONST;
    const currentEta = Math.max(0, Math.min(1, params.eta || 0.5));
    return currentEta * area_L + (1 - currentEta) * area_G;
}

/**
 * Calculates area for split_pvoigt
 */
function _getArea_Split(tth_peak, hkl, params) {
    const GAUSS_AREA_CONST = 1.0644677;
    const LORENTZ_AREA_CONST = 1.5707963;

    const { gamma_G: gG_L, gamma_L: gL_L } = calculateProfileWidths(tth_peak, hkl, params, 'left');
    const { gamma_G: gG_R, gamma_L: gL_R } = calculateProfileWidths(tth_peak, hkl, params, 'right');
    
    const currentEta = Math.max(0, Math.min(1, params.eta_split || 0.5));
    
    const area_G_L = Math.max(1e-9, gG_L) * GAUSS_AREA_CONST;
    const area_L_L = Math.max(1e-9, gL_L) * LORENTZ_AREA_CONST;
    const totalArea_L = currentEta * area_L_L + (1 - currentEta) * area_G_L;

    const area_G_R = Math.max(1e-9, gG_R) * GAUSS_AREA_CONST;
    const area_L_R = Math.max(1e-9, gL_R) * LORENTZ_AREA_CONST;
    const totalArea_R = currentEta * area_L_R + (1 - currentEta) * area_G_R;
    
    // Average area of the two halves
    return (totalArea_L + totalArea_R) / 2.0; 
}

/**
 * [NEW HELPER] Calculates area for tch_aniso
 */
function _getArea_TCH(tth_peak, hkl, params) {
    const GAUSS_AREA_CONST = 1.0644677;
    const LORENTZ_AREA_CONST = 1.5707963;

    const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
    const gG = Math.max(1e-9, gamma_G);
    const gL = Math.max(1e-9, gamma_L);
    
    const fwhm = getPeakFWHM(gG, gL);
    const ratio = (fwhm > 1e-9) ? gL / fwhm : 0;
    const eta_calc = 1.36603 * ratio - 0.47719 * (ratio * ratio) + 0.11116 * Math.pow(ratio, 3);
    const currentEta = Math.max(0, Math.min(1, eta_calc));
    
    // TCH area uses the convoluted FWHM
    const area_G_combined = fwhm * GAUSS_AREA_CONST;
    const area_L_combined = fwhm * LORENTZ_AREA_CONST;
    
    return currentEta * area_L_combined + (1 - currentEta) * area_G_combined;
}



/**
 * Calculates the integrated area under a pseudo-Voigt peak shape.
 * This is now a DISPATCHER function.
 */
function getPseudoVoigtArea(tth_peak, hkl, params) {
    if (!params || !params.profileType) return 1.0;

    let totalArea = 1.0;
    switch (params.profileType) {
        case "simple_pvoigt":
            totalArea = _getArea_Simple(tth_peak, hkl, params);
            break;
        case "split_pvoigt":
            totalArea = _getArea_Split(tth_peak, hkl, params);
            break;
        case "tch_aniso":
            totalArea = _getArea_TCH(tth_peak, hkl, params);
            break;
    }
    return (isFinite(totalArea) && totalArea > 0) ? totalArea : 1.0;
}



/**
 * Applies asymmetry correction (FCJ model for TCH).
 */
function applyAsymmetry(x, x0, tth_peak, params) {
    // Only apply for TCH profile, and only if params are set
    if (!params || params.profileType !== "tch_aniso" || (!params.SL && !params.HL)) {
        return x - x0;
    }
            
    const delta_2theta = x - x0;
    if (Math.abs(delta_2theta) < 1e-9) return 0;
    if (tth_peak < 0.1 || tth_peak >= 180) return delta_2theta;

    const theta_rad = tth_peak * (Math.PI / 180) / 2.0;
    const safe_theta_rad = Math.max(1e-6, Math.min(Math.PI / 2.0 - 1e-6, theta_rad));
    const tan_theta = Math.tan(safe_theta_rad);
    if (Math.abs(tan_theta) < 1e-9) return delta_2theta;
    const cot_theta = 1.0 / tan_theta;

    const SL = params.SL || 0;
    const HL = params.HL || 0;
    const asymmetry_param = SL * cot_theta + HL;
    if (!isFinite(asymmetry_param)) return delta_2theta;

    const correction_term = asymmetry_param * Math.abs(delta_2theta);
    const MAX_CORRECTION_EFFECT = 0.95;
    const asymmetry_factor = Math.max(1e-6, 1.0 - Math.min(Math.abs(correction_term), MAX_CORRECTION_EFFECT));
    const corrected = delta_2theta / asymmetry_factor;
    return isFinite(corrected) ? corrected : delta_2theta;
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

    if (profileType === "simple_pvoigt") {
        const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
        const H_G = Math.max(1e-9, gamma_G), H_L = Math.max(1e-9, gamma_L);
        const eta = Math.max(0, Math.min(1, params.eta || 0.5));
        return { type: 1, x0, asym_param, H_G, H_L, eta, max_d: 10 * (H_G + H_L), Cg };
    } else if (profileType === "split_pvoigt") {
        const wL = calculateProfileWidths(tth_peak, hkl, params, 'left');
        const wR = calculateProfileWidths(tth_peak, hkl, params, 'right');
        const H_G_L = Math.max(1e-9, wL.gamma_G), H_L_L = Math.max(1e-9, wL.gamma_L);
        const H_G_R = Math.max(1e-9, wR.gamma_G), H_L_R = Math.max(1e-9, wR.gamma_L);
        const eta = Math.max(0, Math.min(1, params.eta_split || 0.5));
        return { type: 2, x0, asym_param, H_G_L, H_L_L, H_G_R, H_L_R, eta, max_d: 10 * Math.max(H_G_L + H_L_L, H_G_R + H_L_R), Cg };
    } else { // tch_aniso
        const { gamma_G, gamma_L } = calculateProfileWidths(tth_peak, hkl, params, 'center');
        const H_G = Math.max(1e-9, gamma_G), H_L = Math.max(1e-9, gamma_L);
        const fwhm = getPeakFWHM(H_G, H_L);
        let eta = 0;
        if (fwhm > 1e-9) {
            const ratio = H_L / fwhm;
            eta = Math.max(0, Math.min(1, 1.36603 * ratio - 0.47719 * (ratio * ratio) + 0.11116 * Math.pow(ratio, 3)));
        }
        return { type: 3, x0, asym_param, fwhm: Math.max(1e-9, fwhm), eta, max_d: 10 * (H_G + H_L), Cg };
    }
}

function evalVoigt(x, p) {
    if (!p) return 0.0;
    let delta = x - p.x0;
    if (p.asym_param !== 0 && Math.abs(delta) >= 1e-9) {
        const factor = Math.max(1e-6, 1.0 - Math.min(Math.abs(p.asym_param * Math.abs(delta)), 0.95));
        delta /= factor;
    }
    if (Math.abs(delta) > p.max_d) return 0.0;

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
    return p.eta * l + (1.0 - p.eta) * g;
}

// Keep pseudoVoigt as a backward-compatible wrapper for any legacy call sites
function pseudoVoigt(x, x0, tth_peak, hkl, params) {
    return evalVoigt(x, prepareVoigt(tth_peak, x0, hkl, params));
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
    const WINDOW_MULT = CALCULATION_WINDOW_MULTIPLIER;
    const HEIGHT_CUTOFF = PEAK_HEIGHT_CUTOFF;

    //   K-alpha 1  
    hklList.forEach(peak => {
        if (!peak || !peak.intensity || peak.intensity <= 0 || !peak.tth || peak.tth < 0 || peak.tth > 180) return;

        const basePos1 = peak.tth + zeroShift;
        const shift1 = calculatePeakShift(basePos1, params);
        const peakPos1 = basePos1 + shift1;
        
        const { gamma_G, gamma_L } = calculateProfileWidths(peak.tth, peak, params, 'center');
        const fwhm_approx1 = getPeakFWHM(gamma_G, gamma_L);
        
        const window1 = WINDOW_MULT * Math.max(0.01, fwhm_approx1);
        const min_tth1 = peakPos1 - window1;
        const max_tth1 = peakPos1 + window1;

        let startIndex = 0;
        while (startIndex < n_points && tthAxis[startIndex] < min_tth1) startIndex++;
        if (startIndex === n_points) return;

        const prep1 = prepareVoigt(basePos1, peakPos1, peak, params);
        for (let i = startIndex; i < n_points; i++) {
            const current_tth = tthAxis[i];
            if (current_tth > max_tth1) break;
            const intensityAtPoint = evalVoigt(current_tth, prep1);

            if (intensityAtPoint > HEIGHT_CUTOFF) {
                pattern[i] += peak.intensity * intensityAtPoint;
            }
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
            const shift2 = calculatePeakShift(basePos2, params);
            const peakPos2 = basePos2 + shift2;
            
            const { gamma_G: gG2, gamma_L: gL2 } = calculateProfileWidths(tth2, peak, params, 'center');
            const fwhm_approx2 = getPeakFWHM(gG2, gL2);

            const window2 = WINDOW_MULT * Math.max(0.01, fwhm_approx2);
            const min_tth2 = peakPos2 - window2;
            const max_tth2 = peakPos2 + window2;

            let startIndex2 = 0;
            while (startIndex2 < n_points && tthAxis[startIndex2] < min_tth2) startIndex2++;
            if (startIndex2 === n_points) return;

          const prep2 = prepareVoigt(basePos2, peakPos2, peak, params);
            for (let i = startIndex2; i < n_points; i++) {
                const current_tth = tthAxis[i];
                if (current_tth > max_tth2) break;
                const intensityAtPoint = evalVoigt(current_tth, prep2);

                if (intensityAtPoint > HEIGHT_CUTOFF) {
                    pattern[i] += peak.intensity * ratio21 * intensityAtPoint;
                }
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
        y_obs_net[i] = Math.max(0, expData.intensity[i] - expData.background[i]);
    }

    // 2. Compute total calculated pattern with current peak intensities
    const total_calc = calculatePattern(expData.tth, hklList, params, scratch_calc);

    // 3. Extract area per peak without storing M arrays of size N
    const currentCycleIntensities = new Float64Array(hklList.length);
    const WINDOW_MULT = CALCULATION_WINDOW_MULTIPLIER + 2; 
    
    hklList.forEach((peak, j) => {
        if (!peak || !peak.tth || peak.tth <= 0 || peak.tth >= 180) return;
        
        const basePos1 = peak.tth + zeroShift;
        const peakPos1 = basePos1 + calculatePeakShift(basePos1, params);
        const { gamma_G, gamma_L } = calculateProfileWidths(basePos1, peak, params);
        const window1 = WINDOW_MULT * Math.max(0.01, getPeakFWHM(gamma_G, gamma_L));
        
        const min_tth1 = peakPos1 - window1;
        const max_tth1 = peakPos1 + window1;
        const prep1 = prepareVoigt(basePos1, peakPos1, peak, params);
        
        let tth2 = null, basePos2 = null, peakPos2 = null, min_tth2 = 0, max_tth2 = 0, prep2 = null;
        if (doubletEnabled) {
             const sinTheta1 = Math.sin(peak.tth * deg2rad / 2.0);
             const sinTheta2 = sinTheta1 * lambdaRatio;
             if (Math.abs(sinTheta2) < 1) {
                 tth2 = 2 * Math.asin(sinTheta2) / deg2rad;
                 basePos2 = tth2 + zeroShift;
                 peakPos2 = basePos2 + calculatePeakShift(basePos2, params);
                 const widths2 = calculateProfileWidths(basePos2, peak, params);
                 const window2 = WINDOW_MULT * Math.max(0.01, getPeakFWHM(widths2.gamma_G, widths2.gamma_L));
                 min_tth2 = peakPos2 - window2;
                 max_tth2 = peakPos2 + window2;
                 prep2 = prepareVoigt(basePos2, peakPos2, peak, params);
             }
        }
        
        // Find loop bounds
        let startIdx = 0;
        const actualMin = tth2 ? Math.min(min_tth1, min_tth2) : min_tth1;
        const actualMax = tth2 ? Math.max(max_tth1, max_tth2) : max_tth1;
        
        while (startIdx < n_points && expData.tth[startIdx] < actualMin) startIdx++;
        
        let extracted_area = 0;
        let prev_partitioned_I = 0;
        
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
            if (profile_val > (PEAK_HEIGHT_CUTOFF / 10)) { 
                 const calc_k = peak.intensity * profile_val;
                 if (total_calc[i] > 1e-9) {
                     current_fraction = calc_k / total_calc[i];
                 }
            }
            
            const current_partitioned_I = y_obs_net[i] * current_fraction;
            
            if (i > startIdx) {
                const step_width = expData.tth[i] - expData.tth[i-1];
                if (step_width > 0) {
                    extracted_area += ((prev_partitioned_I + current_partitioned_I) / 2) * step_width;
                }
            }
            prev_partitioned_I = current_partitioned_I;
        }
        currentCycleIntensities[j] = extracted_area;
    });

    hklList.forEach((peak, idx) => {
        if (!peak) return;
        const extracted_area = Math.max(0, currentCycleIntensities[idx] || 0);
        let peak_height = 0;
        if (extracted_area > 0 && peak.tth) {
            const shapeArea_Ka1 = getPseudoVoigtArea(peak.tth, peak, params);
            let total_area_factor = shapeArea_Ka1;
            if (doubletEnabled) {
                const sinTheta1 = Math.sin(peak.tth * deg2rad / 2.0);
                const sinTheta2 = sinTheta1 * (lambda1 > 1e-6 ? (lambda2 / lambda1) : 1.0);
                if (Math.abs(sinTheta2) < 1) {
                    const tth2 = 2 * Math.asin(sinTheta2) / deg2rad;
                    const shapeArea_Ka2 = getPseudoVoigtArea(tth2, peak, params);
                    total_area_factor += ratio21 * shapeArea_Ka2;
                }
            }
            if (total_area_factor > 1e-9) {
                peak_height = extracted_area / total_area_factor;
            }
        }
        peak.intensity = Math.max(1.0, peak_height);
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

    let sum_w_res_sq = 0, sum_w_obs_sq = 0, sum_abs_res = 0, sum_abs_obs_net = 0;
    for (let i = 0; i < N; i++) {
        const obs_i = y_obs[i];
        const calc_i = y_calc[i];
         const w_i = (weights[i] !== undefined && isFinite(weights[i])) ? weights[i] : 1.0;


         if (isFinite(obs_i) && isFinite(calc_i)) {
            const res = obs_i - calc_i;
            const obs_net = obs_i - y_bkg[i];

            sum_w_res_sq += w_i * res * res;
            sum_w_obs_sq += w_i * obs_i * obs_i;
            sum_abs_res += Math.abs(res);
            sum_abs_obs_net += Math.abs(obs_net);
         }
    }

    const Rp = (sum_abs_obs_net > 1e-9) ? 100 * (sum_abs_res / sum_abs_obs_net) : 0;
    const Rwp = (sum_w_obs_sq > 1e-9) ? 100 * Math.sqrt(Math.max(0, sum_w_res_sq / sum_w_obs_sq)) : 0;


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

        const calculateTotalPattern = (targetArray) => {
            enforceSymmetryConstraintsWorker(params, system);
            updateHklPositions(workingHklList, params, system);

            const y_bkg = calculateTotalBackground(workerWorkingData.tth, params, workerBackgroundAnchors, scratch_bkg);
            // In Le Bail mode, extract per-peak intensities from the observed
            // pattern (minus background) using the current profile shapes BEFORE
            // computing the calculated pattern. This replaces the old
            // single-global-scale-factor approximation, which was mathematically
            // wrong: it tried to match per-peak intensities by scaling the whole
            // pattern by one number, so the refiner had to move lattice/profile
            // parameters to compensate for intensities it was forbidden to adjust.
            // Per-peak extraction is the actual Le Bail method.
            if (refinementMode === 'le-bail') {
                leBailIntensityExtraction(
                    { tth: workerWorkingData.tth, intensity: y_obs, background: y_bkg },
                    workingHklList,
                    params,
                    scratch_pattern
                );
            }

const netCalcPattern = calculatePattern(workerWorkingData.tth, workingHklList, params, scratch_pattern);
            for (let i = 0; i < n_points; i++) {
                targetArray[i] = netCalcPattern[i] + y_bkg[i];
            }
            return 1.0;   // no separate scale factor anymore
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

        for (let iter = 0; iter < maxIter; iter++) {
             let oldParams = null;
             let oldHklList = null;
            try {
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

                     // Larger floor than before (1e-7 absolute, 1e-4 relative)
                     // to keep the FD signal above numerical noise for small
                     // profile params like GW, GP, SL, HL.
                     let fd_step = Math.max(1e-7, Math.abs(originalValue) * 1e-4);

                     // If we're at (or below) the lower bound, flipping sign
                     // keeps the perturbed point inside the feasible set so
                     // set() won't silently clamp it and zero the column.
                     const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
                     const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;
                     let fd_sign = 1;
                     if (isFinite(minV) && (originalValue + fd_step > maxV || originalValue - fd_step < minV)) {
                         // near upper bound -> step down; near lower bound -> step up.
                         // If hemmed in on both sides (unusual), fall back to +.
                         if (originalValue + fd_step > maxV && originalValue - fd_step >= minV) {
                             fd_sign = -1;
                         } else if (originalValue - fd_step < minV && originalValue + fd_step <= maxV) {
                             fd_sign = 1;
                         } else {
                             fd_sign = 1;
                         }
                     }

                     mapping.set(params, workingHklList, originalValue + fd_sign * fd_step);
                     calculateTotalPattern(y_calc_total);

                     for (let i = 0; i < n_points; i++) {
                          if (!isFinite(y_calc_baseline[i]) || !isFinite(y_calc_total[i])) {
                              jacobian_col[i] = 0;
                          } else {
                              // Note the fd_sign: df/dx ≈ (f(x+h) - f(x)) / h
                              // with h = fd_sign * fd_step, i.e. signed step.
                              jacobian_col[i] = (y_calc_total[i] - y_calc_baseline[i])
                                                / (fd_sign * fd_step) * sqrt_weights[i];
                          }
                     }

                     mapping.set(params, workingHklList, originalValue);
                     workingHklList.forEach((h, idx) => { if (h) h.intensity = savedIntensities[idx]; });
                     jacobian_T.push([...jacobian_col]);
                 }

                 const jacobian = math.transpose(jacobian_T);
                 const JtJ = math.multiply(jacobian_T, jacobian);
                 const Jtr = math.multiply(jacobian_T, Array.from(residuals));

                  if (n_params === 1 && typeof JtJ === 'number' && isFinite(JtJ)) {
                       finalJtJ = [[JtJ]];
                  } else if (n_params > 0 && Array.isArray(JtJ) && JtJ.length === n_params && Array.isArray(JtJ[0]) && JtJ[0].length === n_params) {
                       finalJtJ = JtJ;
                  } else {
                       console.warn(`LM iter ${iter}: JtJ calculation failed or format invalid.`);
                       finalJtJ = null;
                       lambda = Math.min(1e9, lambda * 5);
                       continue;
                  }

                 // Marquardt diagonal scaling: A = JtJ + lambda * diag(JtJ).
                 // This replaces the old scheme where parameters were
                 // pre-normalized by their initial magnitudes; diag(JtJ)
                 // carries the correct CURRENT local curvature for each
                 // parameter, so damping is applied consistently whether
                 // a parameter has moved a lot or a little from its start.
                 const A_lm = math.clone(finalJtJ);
                 for (let i = 0; i < n_params; i++) {
                     const d = finalJtJ[i][i];
                     const diag = (isFinite(d) && d > 1e-30) ? d : 1e-6;
                     A_lm[i][i] += lambda * diag;
                 }

                 let p_step;
                 try {
                    const p_step_matrix = math.lusolve(A_lm, Jtr);
                    p_step = math.flatten(p_step_matrix);
                 } catch (solveError) {
                      console.warn(`LM iter ${iter}: Solve failed: ${solveError.message}. Increasing lambda.`);
                      lambda = Math.min(1e9, lambda * 10);
                      continue;
                 }

                 if (!p_step || p_step.some(v => !isFinite(v))) {
                     console.warn(`LM iter ${iter}: Step calculation resulted in NaN/Infinity. Increasing lambda.`);
                     lambda = Math.min(1e9, lambda * 5);
                     continue;
                 }

                 // Snapshot BEFORE applying the trial step, so we can revert on rejection.
                 oldParams = JSON.parse(JSON.stringify(params));
                 oldHklList = JSON.parse(JSON.stringify(workingHklList));

                 // Apply trial step in raw parameter space.
                 const p_current = paramMapping.map(m => m.get(params, workingHklList));
                 
                 // CRITICAL IMPROVEMENT: Prevent wild steps from locking params at boundaries
const p_step_clamped = p_step.map((step, idx) => {
    const maxAllowedStep = (paramMapping[idx].step || 0.2) * 4.0; 
    return Math.max(-maxAllowedStep, Math.min(maxAllowedStep, step));
});
                 
                 const p_new = math.add(p_current, p_step_clamped);




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

                 const progressWithinCycle = (iter + 1) / maxIter;
                 const overallProgress = baseProgress + progressWithinCycle * cycleProgressSpan;
                 postMessage({ type: 'progress', value: Math.min(1.0, overallProgress) });

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

    const objective = (p_obj, hkl_list_obj) => {
         try {
            enforceSymmetryConstraintsWorker(p_obj, system); 
            updateHklPositions(hkl_list_obj, p_obj, system);
const y_bkg = calculateTotalBackground(workerWorkingData.tth, p_obj, workerBackgroundAnchors, scratch_bkg);
             // Le Bail: extract per-peak intensities against the current shapes
             // before building the calculated pattern (see LM fix notes above).
             if (refinementMode === 'le-bail') {
                  leBailIntensityExtraction(
                      { tth: workerWorkingData.tth,
                        intensity: workerWorkingData.intensity,
                        background: y_bkg },
                      hkl_list_obj,
                      p_obj,
                      scratch_pattern
                  );
             }

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


    for (let iter = 0; iter < maxIter; iter++) {
         try {
             for (let i = 0; i < numReplicas; i++) {
                 let replica = replicas[i];
                 const originalParams = JSON.parse(JSON.stringify(replica.params));
                 const originalHklList = JSON.parse(JSON.stringify(replica.hklList));

                 const p_idx = Math.floor(Math.random() * paramMapping.length);
                 const mapping = paramMapping[p_idx];
                 const original_val = mapping.get(replica.params, replica.hklList);
                 const step_scale = Math.max(0.01, replica.temp);
                 const base_step_suggestion = (refinementMode === 'pawley' && mapping.isIntensity) ? 0.05 : 0.05;
                 const step_size = (mapping.step || base_step_suggestion) * step_scale * mapping.scale;
                 const random_step = (Math.random() - 0.5) * 2 * step_size;
                 const new_val = original_val + random_step;
                 mapping.set(replica.params, replica.hklList, new_val);

                 const neighbor_cost = objective(replica.params, replica.hklList);
                 const delta_cost = neighbor_cost - replica.cost;

                  const acceptance_prob = (replica.temp > 1e-9 && replica.cost > 0) ? Math.exp(-delta_cost / (replica.cost * replica.temp)) : 0;
                 if (delta_cost < 0 || acceptance_prob > Math.random()) {
                     replica.cost = neighbor_cost;
                 } else {
                     replica.params = originalParams;
                     replica.hklList = originalHklList;
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
                     const delta_cost = rep1.cost - rep2.cost;
                     const acceptance_prob_swap = Math.exp(Math.min(50, delta_beta * delta_cost));

                     if (acceptance_prob_swap > Math.random()) {
                         [rep1.params, rep2.params] = [rep2.params, rep1.params];
                         [rep1.hklList, rep2.hklList] = [rep2.hklList, rep1.hklList];
                         [rep1.cost, rep2.cost] = [rep2.cost, rep1.cost];
                     }
                 }
             }

             const progressWithinCycle = (iter + 1) / maxIter;
             const overallProgress = baseProgress + progressWithinCycle * cycleProgressSpan;
             postMessage({ type: 'progress', value: Math.min(1.0, overallProgress) });


         } catch (error) {
              console.error("Error during PT iteration:", iter, error);
              postMessage({ type: 'error', message: `Error in PT iter ${iter}: ${error.message}` });
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
        const initialValue = initialParams[name] ?? 0;
        const scale = Math.abs(initialValue) > 1e-9 ? Math.abs(initialValue) : defaultScale;

        return {
            name: name,
            scale: scale,        // typical magnitude, for PT step-sizing & display
            step: step,
            isIntensity: false,
            minVal: minVal,
            maxVal: maxVal,
            get: (p_obj, hkl_list_obj) => (p_obj[name] ?? 0),
            set: (p_obj, hkl_list_obj, rawValue) => {
                if (rawValue < minVal) rawValue = minVal;
                if (rawValue > maxVal) rawValue = maxVal;
                if (p_obj && p_obj.hasOwnProperty(name)) {
                    p_obj[name] = rawValue;
                }
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
        mappings.push(createMapping(fitFlags.GW, 'GW', 0.01, 1e-6, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GP, 'GP', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.LX, 'LX', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.eta, 'eta', 0.5, 0, 1, 0.1));
        mappings.push(createMapping(fitFlags.shft, 'shft', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns, 'trns', 0.01, -Infinity, Infinity, 0.1));
    } 
    else if (profileType === "tch_aniso") {
        mappings.push(createMapping(fitFlags.U, 'U', 0.01, 0, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.V, 'V', 0.01, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.W, 'W', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.X, 'X', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.Y, 'Y', 0.01, 1e-6, Infinity, 0.2));
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
        mappings.push(createMapping(fitFlags.GW_L, 'GW_L', 0.01, 1e-6, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.LX_L, 'LX_L', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.GU_R, 'GU_R', 0.01, 0, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GV_R, 'GV_R', 0.01, -Infinity, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GW_R, 'GW_R', 0.01, 1e-6, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.LX_R, 'LX_R', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.eta_split, 'eta_split', 0.5, 0, 1, 0.1));
        mappings.push(createMapping(fitFlags.shft_split, 'shft_split', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns_split, 'trns_split', 0.01, -Infinity, Infinity, 0.1));
    }
   
    const paramMapping = mappings.filter(Boolean);
    return { paramMapping };
}


//   4. Worker Message Handler  
self.onmessage = async function(e) {
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
        // With per-peak intensity extraction done inside the refiner's
        // objective (see refineParametersLM/PT), Le Bail is just a single
        // refinement call. The old 4-cycle outer loop was a workaround
        // for a buggy objective that used a global scale factor instead
        // of per-peak extraction; it is no longer needed.
        if (refinementMode === 'le-bail') {
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
         hklIndexCache = {};
    }
};

//   error handler  
self.onerror = function(event) {
     console.error("Unhandled Worker Error:", event.message, event);
     postMessage({ type: 'error', message: `Unhandled Worker Error: ${event.message}` });
};