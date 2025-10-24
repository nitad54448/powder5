// refinement_worker.js
// version 111, 24 oct 2025
try {
    importScripts('numeric.min.js', 'rules_spaceGroups.js');
} catch (e) {
    console.error("Worker Error: Failed to import scripts.", e);
    // Post an error back immediately if scripts fail to load
    postMessage({ type: 'error', message: `Worker failed to load scripts: ${e.message}` });
    self.close(); // Terminate the worker if essential scripts are missing
}

// Check if numeric library loaded
if (typeof numeric === 'undefined') {
    postMessage({ type: 'error', message: 'Worker Error: numeric.min.js did not load correctly.' });
    self.close();
}
// Check if spaceGroups data loaded
if (typeof spaceGroups === 'undefined') {
    postMessage({ type: 'error', message: 'Worker Error: rules_spaceGroups.js did not load correctly.' });
    self.close();
}


// --- 2. Define Constants & Global Worker State ---
const NUM_BACKGROUND_PARAMS = 9; // Make sure this matches the main script
const CALCULATION_WINDOW_MULTIPLIER = 6.0;
const PEAK_HEIGHT_CUTOFF = 0.002;
const HIGH_WEIGHT_MULTIPLIER = 50.0;

let workerWorkingData = null; // To store the sliced data sent from the main thread
let hklIndexCache = {}; // Cache for HKL indices within the worker



/**
 * Calculates the total background contribution from Chebyshev polynomials and an amorphous hump.
 * @param {Float64Array} tthAxis - The array of 2-theta values.
 * @param {object} params - The object containing all refinement parameters.
 * @returns {Float64Array} A new array containing the calculated background intensity at each point.
 */
function calculateTotalBackground(tthAxis, params) {
    // --- 1. Pre-computation and setup ---
    const n = tthAxis.length;
    if (n === 0) {
        return new Float64Array();
    }

    const background = new Float64Array(n); // Allocate the typed array upfront

    // Extract and check Chebyshev parameters
    const chebyshevCoefficients = [];
    let hasChebyshev = false;
    for (let i = 0; i < NUM_BACKGROUND_PARAMS; i++) {
        const coeff = params[`B${i}`] || 0;
        chebyshevCoefficients[i] = coeff;
        if (Math.abs(coeff) > 1e-9) hasChebyshev = true;
    }

    // Extract and pre-calculate for Amorphous Hump
    const humpHeight = params.hump_H || 0;
    let hasHump = humpHeight > 1e-9;
    let humpPosition, hwhm_sq;
    if (hasHump) {
        humpPosition = params.hump_P || 0;
        const fwhm = params.hump_W || 1;
        hwhm_sq = (fwhm / 2) * (fwhm / 2);
        if (hwhm_sq < 1e-9) hasHump = false; // Avoid division by zero
    }

    // Pre-calculate scaling factors for Chebyshev
    let tthMin, tthRange;
    if (hasChebyshev && n > 0) { // Check length
        tthMin = tthAxis[0];
        tthRange = tthAxis[tthAxis.length - 1] - tthMin;
        if (tthRange <= 0) hasChebyshev = false; // Avoid division by zero
    } else {
        hasChebyshev = false;
    }


    // Early exit if there's nothing to calculate
    if (!hasChebyshev && !hasHump) {
        return background; // Return the zero-filled Float64Array
    }

    // --- 2. Single loop calculation ---
    for (let i = 0; i < n; i++) {
        const tth = tthAxis[i];
        let totalBackgroundValue = 0;

        // Calculate Chebyshev part if needed
        if (hasChebyshev) {
            const x_prime = (2 * (tth - tthMin) / tthRange) - 1;

            let t_n_minus_1 = 1; // T_0(x)
            let t_n = x_prime; // T_1(x)

            let chebyshevValue = chebyshevCoefficients[0] * t_n_minus_1;
            if (chebyshevCoefficients.length > 1) {
                chebyshevValue += (chebyshevCoefficients[1] || 0) * t_n;
            }

            for (let n = 2; n < chebyshevCoefficients.length; n++) {
                const t_n_plus_1 = 2 * x_prime * t_n - t_n_minus_1;
                chebyshevValue += (chebyshevCoefficients[n] || 0) * t_n_plus_1;
                t_n_minus_1 = t_n;
                t_n = t_n_plus_1;
            }
            totalBackgroundValue += chebyshevValue;
        }

        // Calculate Amorphous Hump part if needed
        if (hasHump) {
            const diff = tth - humpPosition;
            const humpValue = humpHeight / (1 + (diff * diff) / hwhm_sq);
            totalBackgroundValue += humpValue;
        }

        background[i] = totalBackgroundValue;
    }

    return background;
}

// --- Space Group / HKL Functions (from rules_spaceGroups.js, now assumed global) ---
// These functions like isReflectionAllowed, evaluateRuleTree, etc.,
// are expected to be available globally because rules_spaceGroups.js was imported... see script in powder5.html

// --- HKL List Generation & Position Update ---
function updateHklPositions(hklList, params, system) {
    const { a, b, c, alpha, beta, lambda } = params;
    if (!a || !lambda || a <= 0 || lambda <= 0 || !hklList) return; // Added !hklList check

    const deg2rad = Math.PI / 180;
    const lambda_sq_over_4 = (lambda * lambda) / 4.0;
    const a_sq = a * a;

    let b_sq, c_sq, sin_beta_sq, cos_beta;
    if (system === 'monoclinic') {
        const beta_rad = (beta || 90) * deg2rad;
        sin_beta_sq = Math.sin(beta_rad);
        sin_beta_sq *= sin_beta_sq;
        cos_beta = Math.cos(beta_rad);
        b_sq = (b || a) * (b || a);
        c_sq = (c || a) * (c || a);
         // Add check for sin_beta_sq being too small (beta near 0 or 180)
         if (Math.abs(sin_beta_sq) < 1e-9) {
             hklList.forEach(peak => { if(peak) { peak.tth = null; peak.d = null; } });
             return;
         }
    } else if (system === 'tetragonal' || system === 'hexagonal' || system === 'rhombohedral' || system === 'trigonal') {
         c_sq = (c && c > 0) ? (c * c) : a_sq; // Use 'a' if c is invalid
    } else if (system === 'orthorhombic') {
         b_sq = (b && b > 0) ? (b * b) : a_sq;
         c_sq = (c && c > 0) ? (c * c) : a_sq;
    }


    hklList.forEach(peak => {
        if (!peak || peak.h_orig === undefined || peak.k_orig === undefined || peak.l_orig === undefined) {
             if(peak) { peak.tth = null; peak.d = null; } // Invalidate if indices missing
             return;
        }
        const h = peak.h_orig;
        const k = peak.k_orig;
        const l = peak.l_orig;
        const h2 = h * h;
        const k2 = k * k;
        const l2 = l * l;

        let inv_d_sq = 0;
        try { // Add try-catch for safety
            switch(system) {
                case 'cubic':
                     if (a_sq <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = (h2 + k2 + l2) / a_sq;
                    break;
                case 'tetragonal':
                     if (a_sq <= 0 || c_sq <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = (h2 + k2) / a_sq + l2 / c_sq;
                    break;
                case 'orthorhombic':
                     if (a_sq <= 0 || b_sq <= 0 || c_sq <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = h2/a_sq + k2/b_sq + l2/c_sq;
                    break;
                case 'hexagonal':
                case 'rhombohedral': // Using hexagonal axes
                case 'trigonal':     // Using hexagonal axes
                     if (a_sq <= 0 || c_sq <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = 4 * (h2 + h*k + k2) / (3 * a_sq) + l2 / c_sq;
                    break;
                case 'monoclinic':
                     if (a_sq <= 0 || b_sq <= 0 || c_sq <= 0 || a <= 0 || c <= 0) throw new Error("Invalid lattice param");
                    inv_d_sq = (1/sin_beta_sq) * (h2/a_sq + k2*sin_beta_sq/b_sq + l2/c_sq - (2*h*l*cos_beta)/(a*c));
                    break;
                 // Add triclinic if needed, complex formula
                 // case 'triclinic': ...
                default:
                    throw new Error(`Unknown system: ${system}`);
            }

            // Validate inv_d_sq
            if (!isFinite(inv_d_sq) || inv_d_sq <= 1e-12) {
                peak.tth = null;
                peak.d = null;
            } else {
                const sinThetaSq = lambda_sq_over_4 * inv_d_sq;
                if (sinThetaSq <= 1 && sinThetaSq > 0) {
                     const thetaRad = Math.asin(Math.sqrt(sinThetaSq));
                     peak.tth = 2 * thetaRad / deg2rad; // Convert back to degrees
                    peak.d = 1 / Math.sqrt(inv_d_sq);
                } else {
                    peak.tth = null; // Cannot calculate angle (sin^2 > 1 or <= 0)
                    peak.d = null;
                }
            }
        } catch (error) {
             // console.error(`Error calculating d-spacing for HKL (${h},${k},${l}) in ${system}: ${error.message}`);
             peak.tth = null;
             peak.d = null;
        }

    }); // end forEach peak
}


function getMultiplicityAndCanonicalHKL(h, k, l, laue_class) {
    if (h === 0 && k === 0 && l === 0) {
        return { multiplicity: 1, canonical_hkl_obj: [0, 0, 0] };
    }

    let m = 0;
     // Use absolute values for comparisons to handle negative indices correctly
     const abs_h = Math.abs(h);
     const abs_k = Math.abs(k);
     const abs_l = Math.abs(l);

     // Sort absolute values: h' >= k' >= l'
     let [h_p, k_p, l_p] = [abs_h, abs_k, abs_l].sort((a, b) => b - a);


    switch (laue_class) {
         case 'm-3m':
             if (h_p > k_p && k_p > l_p && l_p >= 0) m = 48;        // h > k > l > 0
             else if (h_p === k_p && k_p > l_p && l_p >= 0) m = 24;  // h = k > l >= 0
             else if (h_p > k_p && k_p === l_p && l_p >= 0) m = 24;  // h > k = l >= 0
             // else if (h_p > l_p && k_p > l_p && h_p === k_p) m = 24; // Redundant with h=k>l
             else if (h_p === k_p && k_p === l_p && l_p > 0) m = 8;   // h = k = l > 0
             else if (h_p > 0 && k_p === 0 && l_p === 0) m = 6;       // h00
             else if (h_p === k_p && l_p === 0 && h_p > 0) m = 12;      // hh0
             else if (h_p > k_p && k_p > 0 && l_p === 0) m = 24;      // hk0 (h>k>0)
             else m = 1; // Should not happen for non-zero hkl
             break;
         case 'm-3':
              if (h_p > k_p && k_p > l_p && l_p >= 0) m = 24;        // h > k > l >= 0
              else if ((h_p === k_p && k_p > l_p && l_p >= 0) || (h_p > k_p && k_p === l_p && l_p >= 0)) m = 12; // h=k>l>=0 or h>k=l>=0
              else if (h_p === k_p && k_p === l_p && l_p > 0) m = 8;   // h=k=l>0
              else if (h_p > 0 && k_p === 0 && l_p === 0) m = 6;       // h00
              else if (h_p === k_p && l_p === 0 && h_p > 0) m = 12;      // hh0 - Correction: m-3 has 12 for hh0
              else if (h_p > k_p && k_p > 0 && l_p === 0) m = 12;      // hk0 (h>k>0) - Correction: m-3 has 12 for hk0
              else m = 1;
              break;
         case '6/mmm':
             if (l_p > 0) { // l != 0
                  if (abs_h === 0 && abs_k === 0) m = 2; // 00l
                  else if (abs_h > 0 && abs_k === 0) m = 12; // h0l
                  else if (abs_h === abs_k && abs_k > 0) m = 12; // hhl
                  else if (abs_h > abs_k && abs_k >= 0) m = 24; // hkl, h>k>=0
                   // Need to consider 2hk.l, h h 2h .l etc for hexagonal? Usually handled by checking h,k,i=-h-k
                   // Assuming standard 3-index notation for now.
                  else m = 24; // Fallback general
              } else { // l=0 plane
                  if (abs_h === 0 && abs_k === 0) m = 1; // Origin
                  else if (abs_h > 0 && abs_k === 0) m = 6; // h00
                  else if (abs_h === abs_k && abs_k > 0) m = 6; // hh0
                  else if (abs_h > abs_k && abs_k >= 0) m = 12; // hk0, h>k>=0
                  else m = 12; // Fallback general hk0
              }
              break;
         case '6/m':
              if (l_p > 0) m = (abs_h > 0 || abs_k > 0) ? 12 : 2; // 00l vs hkl/h0l/hhl etc.
              else m = (abs_h > 0 || abs_k > 0) ? 6 : 1; // Origin vs hk0/h00/hh0
              break;
        case '-3m': // Needs careful check for rhombohedral vs hexagonal indexing if applicable
             // Using hexagonal axes convention
             if (l_p !== 0) { // l != 0
                 if (abs_h === 0 && abs_k === 0) { m = 2; } // 00l
                 // Special conditions like h-k+l = 3n or -h+k+l = 3n might matter for R centering, handled by isReflectionAllowed
                 // Multiplicity based purely on Laue group symmetry:
                 else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) { m = 12; } // h0l, hhl, 0kl forms
                 else { m = 24; } // General hkl
             } else { // l=0 plane
                 if (abs_h === 0 && abs_k === 0) { m = 1; } // Origin
                 else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) { m = 6; } // h00, hh0 forms
                 else { m = 12; } // General hk0
             }
             break;
         case '-3':
             if (abs_h === 0 && abs_k === 0 && l_p === 0) m = 1; // Origin
             else if (abs_h === 0 && abs_k === 0) m = 2; // 00l
             else m = 6; // General hkl and hk0
             break;
         case '4/mmm':
             if (l_p > 0) { // l != 0
                 if (abs_h === 0 && abs_k === 0) m = 2; // 00l
                 else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) m = 8; // h0l, hhl forms
                 else m = 16; // General hkl
             } else { // l=0 plane
                 if (abs_h === 0 && abs_k === 0) m = 1; // Origin
                 else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) m = 4; // h00, hh0 forms
                 else m = 8; // General hk0
             }
             break;
         case '4/m':
              if (l_p > 0) m = (abs_h > 0 || abs_k > 0) ? 8 : 2; // 00l vs hkl/h0l/hhl
              else m = (abs_h > 0 || abs_k > 0) ? 4 : 1; // Origin vs hk0/h00/hh0
              break;
         case 'mmm':
              if (abs_h > 0 && abs_k > 0 && l_p > 0) m = 8; // hkl
              else if ((abs_h > 0 && abs_k > 0 && l_p === 0) || (abs_h > 0 && abs_k === 0 && l_p > 0) || (abs_h === 0 && abs_k > 0 && l_p > 0)) m = 4; // hk0, h0l, 0kl
              else if (abs_h > 0 || abs_k > 0 || l_p > 0) m = 2; // h00, 0k0, 00l
              else m = 1; // Origin
              break;
         case '2/m': // Assumes unique axis b
              if (abs_k > 0) m = 4; // hkl, 0kl
              else if (abs_k === 0 && (abs_h !== 0 || l_p !== 0)) m = 2; // h0l (including h00, 00l)
              else m = 1; // Origin
              break;
         case '-1':
              if (abs_h === 0 && abs_k === 0 && l_p === 0) m = 1; // Origin
              else m = 2; // General hkl
              break;
        default:
            console.warn("Unknown Laue class:", laue_class, "- assuming multiplicity 1");
            m = 1;
            break;
    }
    // Simple canonical: return original h,k,l for now
    return { multiplicity: m, canonical_hkl_obj: [h, k, l] };
}


function generateAndCacheHklIndices(spaceGroup, maxTth, params) {
    const sgNumber = spaceGroup.number;
    // 1. Check the cache first! If found, return the cached list immediately.
    if (hklIndexCache[sgNumber]) {
        return hklIndexCache[sgNumber];
    }

    // 2. If not in cache, perform the expensive generation.
    const { a, b, c, lambda } = params; // Include b, c
    const { system, laue_class } = spaceGroup;

    // Basic validation
     if (!lambda || lambda <= 0 || !laue_class) {
          console.error("Cannot generate HKL: Missing lambda or laue_class.");
          return [];
     }
     // Ensure at least 'a' is valid
     if (!a || a <=0) {
         console.error("Cannot generate HKL: Invalid lattice parameter 'a'.");
         return [];
     }


    // Use a robust max dimension, considering potentially missing b or c
    const maxDim = Math.max(a || 0, b || a || 0, c || a || 0);
    if (maxDim <= 0) {
        console.error("Cannot generate HKL: All lattice dimensions are invalid.");
        return [];
    }

     // Estimate max index needed based on maxTth and smallest possible d-spacing (related to max dimension)
     const sinThetaMax = Math.sin(maxTth * Math.PI / 360);
      if (sinThetaMax <= 0) return []; // Avoid division by zero/negative dMin
     const dMinEstimate = lambda / (2 * sinThetaMax);
     // Calculate max index somewhat generously. This is an approximation.
     const maxIndex = Math.ceil(maxDim / dMinEstimate) + 5; // Add buffer

    let rawReflections = [];

    // Use a Set to avoid adding duplicate HKLs if the loop logic generates equivalent ones
    const addedHKLs = new Set();
    const getKey = (h, k, l) => `${h},${k},${l}`;


    const loopAndAdd = (h, k, l) => {
        if (h === 0 && k === 0 && l === 0) return; // Skip origin

        // Normalize indices based on symmetry for caching/checking duplicates if needed
        const key = getKey(h, k, l);
        if (addedHKLs.has(key)) return; // Already processed equivalent reflection

        if (isReflectionAllowed(h, k, l, spaceGroup)) {
            // Get multiplicity based on the specific Laue class
            const { multiplicity } = getMultiplicityAndCanonicalHKL(h, k, l, laue_class);
            if (multiplicity > 0) {
                 rawReflections.push({
                    h_orig: h,
                    k_orig: k,
                    l_orig: l,
                    hkl_list: [`(${h},${k},${l})`], // Store original for reference
                    multiplicity: multiplicity
                });
                addedHKLs.add(key); // Mark as added
            }
        }
    };

    // Loops need to cover the necessary range based on system symmetry.
    // These ranges might need adjustment for specific symmetries (e.g., negative indices).
    // The current loops might miss reflections if negative indices are required and
    // not generated by symmetry operations implicitly handled by laue class logic.
    const maxI = maxIndex; // Use calculated maxIndex

    // Adjust loop ranges based on expected index limits per system
    // These are simplified; more complex symmetries might need +/- ranges.
    if (system === 'monoclinic' || system === 'triclinic') {
        // These often require negative indices
        for (let h = -maxI; h <= maxI; h++) {
            for (let k = 0; k <= maxI; k++) { // Often k>=0 is sufficient by convention
                for (let l = -maxI; l <= maxI; l++) {
                     loopAndAdd(h, k, l);
                }
            }
        }
    } else if (system === 'orthorhombic') {
         for (let h = 0; h <= maxI; h++) {
             for (let k = 0; k <= maxI; k++) {
                 for (let l = 0; l <= maxI; l++) {
                      loopAndAdd(h, k, l);
                 }
             }
         }
    }
     else if (system === 'hexagonal' || system === 'trigonal' || system === 'rhombohedral') {
         // Hexagonal indices often have specific ranges, h,k >=0, l can be +/-
          for (let h = 0; h <= maxI; h++) {
              for (let k = 0; k <= maxI; k++) {
                  // Ensure h, k are not both zero unless l is non-zero
                  if (h === 0 && k === 0) {
                       for (let l = -maxI; l <= maxI; l++) loopAndAdd(h, k, l);
                  } else {
                       // Constraint for hexagonal indices like h >= k >= 0 might apply depending on convention
                       // This loop is broader:
                       for (let l = -maxI; l <= maxI; l++) {
                           loopAndAdd(h, k, l);
                           loopAndAdd(k, h, l); // Cover permutations if loop is restricted
                           // Consider i = -(h+k) constraints/equivalences if needed
                       }
                  }
              }
          }
     }
     else { // Cubic, Tetragonal (often h>=k>=l>=0 or similar conventions)
          for (let h = 0; h <= maxI; h++) {
              for (let k = 0; k <= h; k++) {
                  for (let l = 0; l <= k; l++) {
                       loopAndAdd(h, k, l);
                       // Add permutations if loop uses restricted ranges like h>=k>=l
                       if (h !== k || k !== l) {
                           loopAndAdd(h, l, k);
                           loopAndAdd(k, h, l);
                           loopAndAdd(k, l, h);
                           loopAndAdd(l, h, k);
                           loopAndAdd(l, k, h);
                       }
                  }
              }
          }
     }


    // 3. Store the newly generated list in the cache before returning it.
    hklIndexCache[sgNumber] = rawReflections;
    return rawReflections;
}


// --- Profile Shape & Width ---
function calculatePeakShift(tth, params) {
    // Ensure params and profileType exist
    if (!params || !params.profileType) return 0;
    const profileType = String(params.profileType || "4"); // Ensure string

    if (profileType === "4") {
        const thetaRad = tth * (Math.PI / 180) / 2;
        // Handle edge case near 90 degrees 2theta (theta near PI/2)
        if (Math.abs(thetaRad - Math.PI / 2.0) < 1e-6) {
             return 0; // Avoid potential division by zero or large shifts
        }
        const cosTheta = Math.cos(thetaRad);
        const sin2Theta = Math.sin(2 * thetaRad);
        // Ensure parameters exist, default to 0
        const shftParam = params.shft || 0;
        const trnsParam = params.trns || 0;

        const displacementShift = -(shftParam / 1000) * cosTheta * (180 / Math.PI);
        const transparencyShift = trnsParam * sin2Theta * (180 / Math.PI);
         const totalShift = displacementShift + transparencyShift;
         return isFinite(totalShift) ? totalShift : 0; // Ensure finite result
    } else {
        return 0; // No shift for TCH profile in this implementation
    }
}


function calculateProfileWidths(tth, hkl, params) {
    // Ensure params and profileType exist
    if (!params || !params.profileType) return { gamma_G: 1e-4, gamma_L: 1e-4 };
    const profileType = String(params.profileType || "4"); // Ensure string type
    const thetaRad = tth * (Math.PI / 180) / 2;

    // Prevent division by zero or tan(Infinity) at 90 degrees
    const MAX_ANGLE_RAD = Math.PI / 2.0 - 1e-6;
    const safeThetaRad = Math.min(thetaRad, MAX_ANGLE_RAD);
     if (safeThetaRad < 1e-6) { // Also handle angles near zero if needed
         return { gamma_G: 1e-4, gamma_L: 1e-4 };
     }

    const tanTheta = Math.tan(safeThetaRad);
    const cosTheta = Math.cos(safeThetaRad);

    let gamma_G = 1e-4; // Default tiny width, strictly positive
    let gamma_L = 1e-4; // Default tiny width, strictly positive

    if (profileType === "4") {
        // --- Robust Simple pVoigt Calculation ---
        const GU = params.GU || 0;
        const GV = params.GV || 0;
        const GW = params.GW || 0;
        const GP = params.GP || 0;
        const LX = params.LX || 0;

        // Ensure cosTheta is not zero before dividing
        const cosThetaSq_safe = Math.max(cosTheta * cosTheta, 1e-9);
        const gamma_G_sq = GU * tanTheta * tanTheta + GV * tanTheta + GW + GP / cosThetaSq_safe;

        if (gamma_G_sq > 0 && isFinite(gamma_G_sq)) {
            gamma_G = Math.sqrt(gamma_G_sq);
        } else {
            gamma_G = 1e-4; // Fallback width
        }

        const cosTheta_safe = Math.max(cosTheta, 1e-9);
        const calculated_L = LX / cosTheta_safe;

        if (calculated_L > 0 && isFinite(calculated_L)) {
            gamma_L = calculated_L;
        } else {
             gamma_L = 1e-4; // Fallback width
        }

    } else { // Profile 3 (TCH)
        // --- TCH Profile Calculation ---
         const U = params.U || 0;
         const V = params.V || 0;
         const W = params.W || 0;
         const X = params.X || 0;
         const Y = params.Y || 0;

        const gamma_G_sq = U * tanTheta * tanTheta + V * tanTheta + W;

        if (gamma_G_sq > 0 && isFinite(gamma_G_sq)) {
            gamma_G = Math.sqrt(gamma_G_sq);
        } else {
            gamma_G = 1e-4; // Fallback
        }

        const cosTheta_safe = Math.max(cosTheta, 1e-9);
        const calculated_L = X * tanTheta + Y / cosTheta_safe;

        if (calculated_L > 0 && isFinite(calculated_L)) {
            gamma_L = calculated_L;
        } else {
             gamma_L = 1e-4; // Fallback
        }

        // Anisotropic broadening (Stephens) - Profile 3 only
         if (hkl && hkl.d && hkl.h_orig !== undefined) {
             const d_sq = hkl.d * hkl.d;
             if (d_sq > 1e-9) { // Avoid division by zero
                const d_inv_sq = 1 / d_sq;
                 // Use original hkl indices directly
                 const h_val = hkl.h_orig;
                 const k_val = hkl.k_orig;
                 const l_val = hkl.l_orig;

                 const h2 = h_val * h_val, k2 = k_val * k_val, l2 = l_val * l_val;
                 const h4 = h2*h2, k4 = k2*k2, l4 = l2*l2; // Need 4th powers for Stephens

                 // Ensure anisotropy parameters exist, default to 0
                 const S400 = params.S400 || 0;
                 const S040 = params.S040 || 0;
                 const S004 = params.S004 || 0;
                 const S220 = params.S220 || 0;
                 const S202 = params.S202 || 0;
                 const S022 = params.S022 || 0;

                 // Stephens formula - check exact form if needed, often involves h^4, k^4, l^4 and cross terms
                 let H_aniso = S400 * h4 + S040 * k4 + S004 * l4 +
                               S220 * h2 * k2 + S202 * h2 * l2 + S022 * k2 * l2;
                 // The formula might need normalization or scaling factors depending on definition
                 // Often involves division by (h^2+k^2+l^2)^2 or similar, or uses d^-4 scaling
                 H_aniso *= d_inv_sq * d_inv_sq; // Example using d^-4 scaling

                 // Only add anisotropic broadening if it's a valid, finite, positive number
                 if(isFinite(H_aniso) && H_aniso > 0) {
                      // The division by 1000 might be specific to implementation/units
                      gamma_L += H_aniso / 1000.0;
                 }
            }
        }
    }

    // Final safeguard: ensure widths are always positive and finite
    return {
        gamma_G: Math.max(1e-4, isFinite(gamma_G) ? gamma_G : 1e-4),
        gamma_L: Math.max(1e-4, isFinite(gamma_L) ? gamma_L : 1e-4)
    };
}


function getPeakFWHM(gamma_G, gamma_L) {
    // Ensure inputs are positive numbers
    const gG = Math.max(1e-9, gamma_G || 1e-9);
    const gL = Math.max(1e-9, gamma_L || 1e-9);

    // Thompson-Cox-Hastings approximation
    const fwhm_g_5 = Math.pow(gG, 5);
    const fwhm_l_5 = Math.pow(gL, 5);
    const fwhm_g_4_l = Math.pow(gG, 4) * gL;
    const fwhm_g_3_l_2 = Math.pow(gG, 3) * Math.pow(gL, 2);
    const fwhm_g_2_l_3 = Math.pow(gG, 2) * Math.pow(gL, 3);
    const fwhm_g_l_4 = gG * Math.pow(gL, 4);

    const fwhm_pow5 = fwhm_g_5 +
                      2.69269 * fwhm_g_4_l +
                      2.42843 * fwhm_g_3_l_2 +
                      4.47163 * fwhm_g_2_l_3 +
                      0.07842 * fwhm_g_l_4 +
                      fwhm_l_5;

     if (fwhm_pow5 < 0 || !isFinite(fwhm_pow5)) return Math.max(gG, gL, 1e-6); // Fallback

     const fwhm = Math.pow(fwhm_pow5, 0.2);
     return Math.max(1e-6, fwhm); // Ensure minimum positive FWHM
}


function getPseudoVoigtArea(gamma_G, gamma_L, params) {
    const GAUSS_AREA_CONST = 1.0644677; // sqrt(pi / (4*ln(2)))
    const LORENTZ_AREA_CONST = 1.5707963; // pi / 2

    // Ensure params and profileType exist
     if (!params || !params.profileType) return 1.0; // Default area? Or throw error?
    const profileType = String(params.profileType || "4");
    const gG = Math.max(1e-9, gamma_G || 1e-9); // Ensure positive widths
    const gL = Math.max(1e-9, gamma_L || 1e-9);


    if (profileType === "4") {
        const area_G = gG * GAUSS_AREA_CONST;
        const area_L = gL * LORENTZ_AREA_CONST;
        const currentEta = Math.max(0, Math.min(1, params.eta || 0.5));
        const totalArea = currentEta * area_L + (1 - currentEta) * area_G;
         return isFinite(totalArea) && totalArea > 0 ? totalArea : 1.0; // Ensure positive finite area

    } else { // TCH Profile (Type 3)
        const fwhm = getPeakFWHM(gG, gL);

        const ratio = (fwhm > 1e-9) ? gL / fwhm : 0;
        const eta_calc = 1.36603 * ratio - 0.47719 * (ratio * ratio) + 0.11116 * Math.pow(ratio, 3);
        const currentEta = Math.max(0, Math.min(1, eta_calc));

        const area_G_combined = fwhm * GAUSS_AREA_CONST;
        const area_L_combined = fwhm * LORENTZ_AREA_CONST;
        const totalArea = currentEta * area_L_combined + (1 - currentEta) * area_G_combined;
         return isFinite(totalArea) && totalArea > 0 ? totalArea : 1.0; // Ensure positive finite area
    }
}


function applyAsymmetry(x, x0, tth_peak, params) {
    // Ensure params exist
    if (!params) return x - x0;
    const profileType = String(params.profileType || "4");

    // No asymmetry for profile 4 or if SL/HL are zero/missing
    if (profileType === "4" || (!params.SL && !params.HL)) {
        return x - x0;
    }

    const delta_2theta = x - x0;
    if (Math.abs(delta_2theta) < 1e-9) return 0;
    if (tth_peak < 0.1 || tth_peak >= 180) return delta_2theta; // Invalid peak position

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
    // Ensure factor is positive and not too close to zero
    const asymmetry_factor = Math.max(1e-6, 1.0 - Math.min(Math.abs(correction_term), MAX_CORRECTION_EFFECT));

    const corrected = delta_2theta / asymmetry_factor;
     return isFinite(corrected) ? corrected : delta_2theta; // Fallback if calculation fails
}

function pseudoVoigt(x, x0, tth_peak, gamma_G, gamma_L, params) {
     // Ensure params exist
     if (!params) return 0.0;
     const H_G = Math.max(1e-9, gamma_G || 1e-9);
     const H_L = Math.max(1e-9, gamma_L || 1e-9);

    const corrected_delta = applyAsymmetry(x, x0, tth_peak, params);
    if (Math.abs(corrected_delta) > 10 * (H_G + H_L)) return 0.0;

    const profileType = String(params.profileType || "4");
    const Cg = 2.772588722239781; // 4 * ln(2)

    let result = 0.0;
    try { // Add try-catch for safety during calculation
         if (profileType === "4") {
             const currentEta = Math.max(0, Math.min(1, params.eta || 0.5));
             const delta_over_Hg_sq = Math.pow(corrected_delta / H_G, 2);
             const delta_over_Hl_sq = Math.pow(corrected_delta / H_L, 2);
             const gaussianShape = Math.exp(-Cg * delta_over_Hg_sq);
             const lorentzianShape = 1 / (1 + 4 * delta_over_Hl_sq);
             result = currentEta * lorentzianShape + (1 - currentEta) * gaussianShape;
         } else { // TCH (Profile 3)
             const fwhm = getPeakFWHM(H_G, H_L);
             if (fwhm <= 1e-9) return Math.abs(corrected_delta) < 1e-6 ? 1.0 : 0.0; // Approx delta

             const ratio = H_L / fwhm;
             const eta_calc = 1.36603 * ratio - 0.47719 * (ratio * ratio) + 0.11116 * Math.pow(ratio, 3);
             const currentEta = Math.max(0, Math.min(1, eta_calc));

             const delta_over_fwhm_sq = Math.pow(corrected_delta / fwhm, 2);
             const gaussianShape = Math.exp(-Cg * delta_over_fwhm_sq);
             const lorentzianShape = 1 / (1 + 4 * delta_over_fwhm_sq);
             result = currentEta * lorentzianShape + (1 - currentEta) * gaussianShape;
         }
    } catch (calcError) {
         console.warn("Error in pseudoVoigt calculation:", calcError);
         return 0.0; // Return 0 if any math error occurs
    }

     return (isFinite(result) && result >= 0) ? result : 0.0; // Ensure non-negative finite result
}


function calculatePattern(tthAxis, hklList, params) {
     // Ensure params exist
     if (!params) return new Float64Array(tthAxis ? tthAxis.length : 0);

    const n_points = tthAxis ? tthAxis.length : 0;
    if (n_points === 0 || !hklList || hklList.length === 0) {
        return new Float64Array(n_points);
    }

    const pattern = new Float64Array(n_points);
    const deg2rad = Math.PI / 180;
    const lambda1 = params.lambda || 1.54056;
    const lambda2 = params.lambda2 || 0;
    const ratio21 = params.ratio || 0;
    const zeroShift = params.zeroShift || 0;
    // Use global constants for calculation window and cutoff
    const WINDOW_MULT = CALCULATION_WINDOW_MULTIPLIER;
    const HEIGHT_CUTOFF = PEAK_HEIGHT_CUTOFF;


    // --- K-alpha 1 ---
    hklList.forEach(peak => {
        if (!peak || !peak.intensity || peak.intensity <= 0 || !peak.tth || peak.tth < 0 || peak.tth > 180) return;

        const basePos1 = peak.tth + zeroShift;
        const shift1 = calculatePeakShift(basePos1, params);
        const peakPos1 = basePos1 + shift1;
        const { gamma_G, gamma_L } = calculateProfileWidths(basePos1, peak, params);
        const shapeArea1 = getPseudoVoigtArea(gamma_G, gamma_L, params);
        if (shapeArea1 <= 1e-9 || !isFinite(shapeArea1)) return;

        const fwhm_approx1 = getPeakFWHM(gamma_G, gamma_L);
        const window1 = WINDOW_MULT * Math.max(0.01, fwhm_approx1);
        const min_tth1 = peakPos1 - window1;
        const max_tth1 = peakPos1 + window1;

        let startIndex = 0;
        while (startIndex < n_points && tthAxis[startIndex] < min_tth1) startIndex++;
        if (startIndex === n_points) return;

        for (let i = startIndex; i < n_points; i++) {
            const current_tth = tthAxis[i];
            if (current_tth > max_tth1) break;
            const intensityAtPoint = pseudoVoigt(current_tth, peakPos1, basePos1, gamma_G, gamma_L, params);
            if (intensityAtPoint > HEIGHT_CUTOFF) {
                pattern[i] += peak.intensity * (intensityAtPoint / shapeArea1);
            }
        }
    });

    // --- K-alpha 2 ---
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
            const { gamma_G: gG2, gamma_L: gL2 } = calculateProfileWidths(basePos2, peak, params);
            const shapeArea2 = getPseudoVoigtArea(gG2, gL2, params);
            if (shapeArea2 <= 1e-9 || !isFinite(shapeArea2)) return;

            const fwhm_approx2 = getPeakFWHM(gG2, gL2);
            const window2 = WINDOW_MULT * Math.max(0.01, fwhm_approx2);
            const min_tth2 = peakPos2 - window2;
            const max_tth2 = peakPos2 + window2;

            let startIndex2 = 0;
            while (startIndex2 < n_points && tthAxis[startIndex2] < min_tth2) startIndex2++;
            if (startIndex2 === n_points) return;

            for (let i = startIndex2; i < n_points; i++) {
                const current_tth = tthAxis[i];
                if (current_tth > max_tth2) break;
                const intensityAtPoint = pseudoVoigt(current_tth, peakPos2, basePos2, gG2, gL2, params);
                if (intensityAtPoint > HEIGHT_CUTOFF) {
                    pattern[i] += peak.intensity * ratio21 * (intensityAtPoint / shapeArea2);
                }
            }
        });
    }

    // Final check for NaNs/Infinities in the pattern
    for (let i = 0; i < n_points; i++) {
        if (!isFinite(pattern[i])) {
             console.warn(`calculatePattern: Found non-finite value at index ${i}, setting to 0.`);
            pattern[i] = 0; // Replace invalid numbers with zero
        }
    }

    return pattern;
}



// --- Le Bail Intensity Extraction 
function leBailIntensityExtraction(expData, hklList, params) {
    // Basic validation
    if (!expData || !expData.tth || !expData.intensity || !expData.background || !hklList ||
        expData.tth.length !== expData.intensity.length || expData.tth.length !== expData.background.length) {
        console.error("leBailIntensityExtraction: Invalid input data.");
        if (hklList) hklList.forEach(p => { if(p) p.intensity = 0; }); // Reset intensities on error
        return;
    }
     const n_points = expData.tth.length;
     if (n_points === 0) return; // Nothing to do

    const deg2rad = Math.PI / 180;
    const lambda1 = params.lambda || 1.54056;
    const lambda2 = params.lambda2 || 0;
    const ratio21 = params.ratio || 0;
    const zeroShift = params.zeroShift || 0;
    const doubletEnabled = ratio21 > 1e-6 && lambda2 > 1e-6 && Math.abs(lambda1 - lambda2) > 1e-6;
    const lambdaRatio = doubletEnabled ? lambda2 / lambda1 : 1.0;

    // Use a temporary array to store calculated intensities for this cycle
    const currentCycleIntensities = new Array(hklList.length).fill(0.0);

    for (let i = 0; i < n_points; i++) {
        let totalProfileAtPoint = 0;
        const contributingPeaks = []; // Stores { peakIndex: index, profileValue: value }

        const current_tth = expData.tth[i];

        // --- Calculate contributions from ALL peaks (Ka1 and Ka2) at this point ---
        hklList.forEach((peak, peakIdx) => {
            // --- Ka1 contribution ---
            if (peak && peak.tth && peak.tth > 0 && peak.tth < 180) {
                 const basePos1 = peak.tth + zeroShift;
                 const shift1 = calculatePeakShift(basePos1, params);
                 const peakPos1 = basePos1 + shift1;
                 const { gamma_G, gamma_L } = calculateProfileWidths(basePos1, peak, params);
                 const fwhm_approx1 = getPeakFWHM(gamma_G, gamma_L);
                 // Use a slightly larger window maybe for extraction robustness?
                 const window1 = (CALCULATION_WINDOW_MULTIPLIER + 2) * Math.max(0.01, fwhm_approx1);

                if (Math.abs(current_tth - peakPos1) < window1) {
                    const profileVal1 = pseudoVoigt(current_tth, peakPos1, basePos1, gamma_G, gamma_L, params);
                    // Use a lower cutoff for contributions?
                    if (profileVal1 > PEAK_HEIGHT_CUTOFF / 10) {
                        totalProfileAtPoint += profileVal1;
                        contributingPeaks.push({ idx: peakIdx, val: profileVal1 });
                    }
                }
            }

            // --- Ka2 contribution ---
            if (doubletEnabled && peak && peak.tth && peak.tth > 0 && peak.tth < 180) {
                 const sinTheta1 = Math.sin(peak.tth * deg2rad / 2.0);
                 const sinTheta2 = sinTheta1 * lambdaRatio;

                 if (Math.abs(sinTheta2) < 1) {
                      const tth2 = 2 * Math.asin(sinTheta2) / deg2rad;
                      const basePos2 = tth2 + zeroShift;
                      const shift2 = calculatePeakShift(basePos2, params);
                      const peakPos2 = basePos2 + shift2;
                      const { gamma_G: gamma_G2, gamma_L: gamma_L2 } = calculateProfileWidths(basePos2, peak, params);
                      const fwhm_approx2 = getPeakFWHM(gamma_G2, gamma_L2);
                      const window2 = (CALCULATION_WINDOW_MULTIPLIER + 2) * Math.max(0.01, fwhm_approx2);

                     if (Math.abs(current_tth - peakPos2) < window2) {
                         const profileVal2 = pseudoVoigt(current_tth, peakPos2, basePos2, gamma_G2, gamma_L2, params);
                         if (profileVal2 > PEAK_HEIGHT_CUTOFF / 10) { // Lower cutoff?
                             const weightedProfileVal2 = profileVal2 * ratio21;
                             totalProfileAtPoint += weightedProfileVal2;
                             contributingPeaks.push({ idx: peakIdx, val: weightedProfileVal2 });
                          }
                     }
                 }
            }
        }); // End hklList.forEach

        // --- Distribute observed intensity ---
        if (totalProfileAtPoint > 1e-9 && contributingPeaks.length > 0) {
            const obsIntensityNet = expData.intensity[i] - (expData.background[i] || 0);

            if (obsIntensityNet > 0) {
                // Distribute intensity based on relative profile contribution at this point
                contributingPeaks.forEach(contribution => {
                    if (contribution.idx >= 0 && contribution.idx < currentCycleIntensities.length) {
                        currentCycleIntensities[contribution.idx] += obsIntensityNet * (contribution.val / totalProfileAtPoint);
                    }
                });
            }
        }
    } // End loop over data points

    hklList.forEach((peak, idx) => {
        if (!peak) return; // Safety check
        // Directly assign the calculated intensity for this cycle
        const newIntensity = Math.max(0, currentCycleIntensities[idx] || 0); // Ensure non-negative and handle undefined
        peak.intensity = newIntensity;

        // Ensure intensity_previous is removed if it exists
        delete peak.intensity_previous;
    });

} 


// --- Statistics ---
function calculateStatistics(localWorkingData, netCalcPattern, fitFlags, finalBackground, params, hklList, refinementMode) {
    const y_obs = localWorkingData.intensity;
    const y_bkg = finalBackground;
    const weights = localWorkingData.weights;
    const N = y_obs.length;

    // Validate inputs
    if (!y_obs || !netCalcPattern || !y_bkg || !weights ||
        N === 0 || y_obs.length !== netCalcPattern.length || y_obs.length !== y_bkg.length || y_obs.length !== weights.length) {
        console.error("Statistics calculation error: Mismatched or invalid array inputs.");
        return { r_p: -1, rwp: -1, chi2: -1, scaleFactor: 1, sum_w_res_sq: 0 };
    }

    let scaleFactor = 1.0;
    if (refinementMode === 'le-bail') {
        let sum_w_y_net_y_calc = 0;
        let sum_w_y_calc_sq = 0;
        for (let i = 0; i < N; i++) {
             // Ensure weights are valid numbers
             const w_i = (weights[i] !== undefined && isFinite(weights[i])) ? weights[i] : 1.0;
            const y_obs_net = y_obs[i] - y_bkg[i];
            const y_calc_i = netCalcPattern[i] || 0; // Ensure y_calc is numeric
             if (isFinite(y_obs_net) && isFinite(y_calc_i)) {
                sum_w_y_net_y_calc += w_i * y_obs_net * y_calc_i;
                sum_w_y_calc_sq += w_i * y_calc_i * y_calc_i;
             }
        }
        // Avoid division by zero, ensure scale factor is positive and finite
        scaleFactor = (sum_w_y_calc_sq > 1e-12) ? Math.max(0, sum_w_y_net_y_calc / sum_w_y_calc_sq) : 1.0;
         if (!isFinite(scaleFactor)) scaleFactor = 1.0; // Fallback
    }

    // Calculate total calculated pattern y_calc = scale * net + background
    const y_calc = numeric.add(numeric.mul(netCalcPattern, scaleFactor), y_bkg);

    let sum_w_res_sq = 0, sum_w_obs_sq = 0, sum_abs_res = 0, sum_abs_obs_net = 0; // Use net observed for Rp
    for (let i = 0; i < N; i++) {
        const obs_i = y_obs[i];
        const calc_i = y_calc[i];
         const w_i = (weights[i] !== undefined && isFinite(weights[i])) ? weights[i] : 1.0;


         if (isFinite(obs_i) && isFinite(calc_i)) {
            const res = obs_i - calc_i;
            const obs_net = obs_i - y_bkg[i]; // Net observed intensity for Rp

            sum_w_res_sq += w_i * res * res;
            sum_w_obs_sq += w_i * obs_i * obs_i; // Use gross observed for Rwp denominator
            sum_abs_res += Math.abs(res);
            sum_abs_obs_net += Math.abs(obs_net); // Sum of |Yobs_net| for Rp
         }
    }

    // Calculate Rp using net intensities
    const Rp = (sum_abs_obs_net > 1e-9) ? 100 * (sum_abs_res / sum_abs_obs_net) : 0;

    // Calculate Rwp using gross intensities
    const Rwp = (sum_w_obs_sq > 1e-9) ? 100 * Math.sqrt(Math.max(0, sum_w_res_sq / sum_w_obs_sq)) : 0; // Ensure sqrt argument is non-negative


    // --- Calculate Chi-squared (Goodness of Fit) , attention style GSAS---
    // Count number of refined parameters (P)
     const { paramMapping } = getParameterMapping(fitFlags || {}, params || {}, hklList || [], refinementMode || 'le-bail');
     const P_base = paramMapping.filter(m => !m.isIntensity).length; // Count non-intensity parameters
     let P = P_base;
     if (refinementMode === 'pawley' && hklList && localWorkingData && localWorkingData.tth && localWorkingData.tth.length > 0) { // Check workingData
         // Count intensities actually refined (within the working data range)
          const tthMin = localWorkingData.tth[0];
          const tthMax = localWorkingData.tth[localWorkingData.tth.length - 1];
         const refinedIntensitiesCount = hklList.filter(hkl =>
             hkl && hkl.tth && hkl.tth >= tthMin && hkl.tth <= tthMax
         ).length;
         P += refinedIntensitiesCount;
     }

    // Degrees of freedom
    const degreesOfFreedom = N - P;

    let chi2 = 0;
    if (degreesOfFreedom > 0) {
        // Chi^2 = Sum(w * (y_obs - y_calc)^2) / (N - P)
        chi2 = sum_w_res_sq / degreesOfFreedom;
        if (!isFinite(chi2) || chi2 < 0) chi2 = 0; // Ensure valid value
    }


    return {
        r_p: isFinite(Rp) ? Rp : -1,
        rwp: isFinite(Rwp) ? Rwp : -1,
        chi2: isFinite(chi2) ? chi2 : -1,
        scaleFactor: scaleFactor,
        sum_w_res_sq: isFinite(sum_w_res_sq) ? sum_w_res_sq : 0
    };
}


// --- Refinement Algorithms (LM, SA, PT) ---
// These need the worker's postMessage for progress updates
async function refineParametersLM(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
        const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode);
        if (paramMapping.length === 0) {
             const finalProgress = (leBailCycle + 1) / totalLeBailCycles;
             postMessage({ type: 'progress', value: Math.min(1.0, finalProgress) }); // Ensure 100% for this cycle
            return { params: initialParams, hklList: hklList, ss_res: 0 };
        }

        let params = initialParams;
        let workingHklList = JSON.parse(JSON.stringify(hklList)); // Deep copy

        let finalJtJ = null, ss_res = Infinity;
        let lambda = 0.001; // LM damping parameter

        // Use the worker's workingData
        const y_obs = workerWorkingData.intensity;
        const sqrt_weights = workerWorkingData.weights.map(w => Math.sqrt(w)); // Assumes weights is Float64Array or similar
        const n_points = y_obs.length;
        const n_params = paramMapping.length;

        // Pre-allocate arrays
        const residuals = new Float64Array(n_points);
        const y_calc_total = new Float64Array(n_points);
        const y_calc_baseline = new Float64Array(n_points);
        const jacobian_col = new Float64Array(n_points);

        // --- Progress Scaling Setup ---
        const baseProgress = leBailCycle / totalLeBailCycles;
        const cycleProgressSpan = 1 / totalLeBailCycles;
        // --- End Progress Scaling Setup ---

        const calculateTotalPattern = (targetArray) => {
            updateHklPositions(workingHklList, params, system);
            const y_bkg = calculateTotalBackground(workerWorkingData.tth, params);
            const netCalcPattern = calculatePattern(workerWorkingData.tth, workingHklList, params);

            let scaleFactor = 1.0;
            if (refinementMode === 'le-bail') {
                let num = 0, den = 0;
                for (let i = 0; i < n_points; i++) {
                     const w_i = workerWorkingData.weights[i];
                     const y_obs_net_i = y_obs[i] - y_bkg[i];
                     const net_calc_i = netCalcPattern[i] || 0;
                     if (isFinite(w_i) && isFinite(y_obs_net_i) && isFinite(net_calc_i)) {
                        num += w_i * y_obs_net_i * net_calc_i;
                        den += w_i * net_calc_i * net_calc_i;
                     }
                }
                 scaleFactor = (den > 1e-12) ? Math.max(0, num / den) : 1.0;
                 if (!isFinite(scaleFactor)) scaleFactor = 1.0;
            }

             const totalPattern = numeric.add(numeric.mul(netCalcPattern, scaleFactor), y_bkg);
             for (let i = 0; i < n_points; i++) {
                 targetArray[i] = totalPattern[i]; // Copy result to targetArray
             }
            return scaleFactor; // Return scale factor used
        };


        for (let iter = 0; iter < maxIter; iter++) {
             let oldParams = null;       // Store previous state only if needed for rollback
             let oldHklList = null;
            try {
                 calculateTotalPattern(y_calc_baseline);

                 // Calculate residuals and cost
                 let cost = 0;
                 for (let i = 0; i < n_points; i++) {
                     residuals[i] = (y_obs[i] - y_calc_baseline[i]) * sqrt_weights[i];
                     if (isFinite(residuals[i])) { // Check for NaN/Infinity
                          cost += residuals[i] * residuals[i];
                     } else {
                         // If cost calculation fails early, try increasing lambda? Or break?
                         // Let's try increasing lambda and retrying this iter
                          lambda = Math.min(1e9, lambda * 10);
                          console.warn(`LM iter ${iter}: Residual calculation failed. Increased lambda to ${lambda}.`);
                          if(oldParams) { // Rollback if we have a previous state
                               params = oldParams;
                               workingHklList = oldHklList;
                          }
                         continue; // Skip rest of iter, retry with higher damping
                         // throw new Error(`Residual calculation failed at index ${i}`);
                     }
                 }

                 // Check convergence
                 if (iter > 0 && Math.abs(ss_res - cost) < 1e-9 * ss_res) {
                     break; // Converged
                 }
                 ss_res = cost;

                 // Calculate Jacobian using finite differences
                 const jacobian_T = []; // Store as columns (transposed)

                 for (let p = 0; p < n_params; p++) {
                     const mapping = paramMapping[p];
                     const originalValue = mapping.get(params, workingHklList);
                     const fd_step = Math.max(1e-6, Math.abs(originalValue) * 1e-5);

                     mapping.set(params, workingHklList, originalValue + fd_step);
                     calculateTotalPattern(y_calc_total);

                     for (let i = 0; i < n_points; i++) {
                          if (!isFinite(y_calc_baseline[i]) || !isFinite(y_calc_total[i])) {
                              jacobian_col[i] = 0;
                          } else {
                              jacobian_col[i] = (y_calc_total[i] - y_calc_baseline[i]) / fd_step * sqrt_weights[i];
                          }
                     }
                     mapping.set(params, workingHklList, originalValue); // Restore
                     jacobian_T.push([...jacobian_col]);
                 }

                 // Calculate JtJ and Jtr
                 const jacobian = numeric.transpose(jacobian_T);
                 const JtJ = numeric.dot(jacobian_T, jacobian);
                 const Jtr = numeric.dot(jacobian_T, residuals);

                 // --- Ensure JtJ is Correct Format ---
                  if (n_params === 1 && typeof JtJ === 'number' && isFinite(JtJ)) {
                       finalJtJ = [[JtJ]];
                  } else if (n_params > 0 && Array.isArray(JtJ) && JtJ.length === n_params && Array.isArray(JtJ[0]) && JtJ[0].length === n_params) {
                       finalJtJ = JtJ;
                  } else {
                       console.warn("LM iter ${iter}: JtJ calculation failed or format invalid.");
                       finalJtJ = null; // Mark as invalid
                       // Try increasing lambda and continuing?
                       lambda = Math.min(1e9, lambda * 5);
                       continue; // Skip solve step for this iteration
                       // throw new Error("JtJ calculation failed or format invalid.");
                  }
                 // --- End Format Check ---

                 // Apply Levenberg-Marquardt damping
                 const A_lm = numeric.clone(finalJtJ);
                 for (let i = 0; i < n_params; i++) {
                     A_lm[i][i] += lambda * (finalJtJ[i][i] || 1e-6);
                 }

                 // Solve for parameter step
                 let p_step;
                 try {
                     p_step = numeric.solve(A_lm, Jtr);
                 } catch (solveError) {
                      console.warn(`LM iter ${iter}: Solve failed: ${solveError.message}. Increasing lambda.`);
                      lambda = Math.min(1e9, lambda * 10);
                       if(oldParams) { // Rollback if possible
                            params = oldParams;
                            workingHklList = oldHklList;
                       }
                      continue; // Retry iter with higher lambda
                 }

                 if (!p_step || p_step.some(v => !isFinite(v))) {
                     console.warn(`LM iter ${iter}: Step calculation resulted in NaN/Infinity. Increasing lambda.`);
                     lambda = Math.min(1e9, lambda * 5);
                     if(oldParams) { // Rollback if possible
                          params = oldParams;
                          workingHklList = oldHklList;
                     }
                     continue; // Retry iter with higher lambda
                 }

                 // Store current state for potential rollback *before* applying step
                 oldParams = JSON.parse(JSON.stringify(params));
                 oldHklList = JSON.parse(JSON.stringify(workingHklList));

                 // Apply step to get new parameters
                 const p_current = paramMapping.map(m => m.get(params, workingHklList));
                 const p_new = numeric.add(p_current, p_step);
                 paramMapping.forEach((m, i) => m.set(params, workingHklList, p_new[i]));

                 // Calculate new cost
                 calculateTotalPattern(y_calc_total);
                 let new_cost = 0;
                 for (let i = 0; i < n_points; i++) {
                     const res = (y_obs[i] - y_calc_total[i]) * sqrt_weights[i];
                      if (isFinite(res)) {
                          new_cost += res * res;
                      } else {
                          new_cost = Infinity; // Invalid cost
                          break;
                      }
                 }

                 // Accept or reject step
                 if (new_cost < cost && isFinite(new_cost)) {
                     lambda = Math.max(1e-9, lambda / 3); // Decrease damping
                     // Keep new params
                 } else {
                     // Restore old state
                     params = oldParams;
                     workingHklList = oldHklList;
                     lambda = Math.min(1e9, lambda * 2); // Increase damping
                     // ss_res remains the same as the previous iteration
                 }

                 // --- MODIFIED PROGRESS UPDATE ---
                 const progressWithinCycle = (iter + 1) / maxIter;
                 const overallProgress = baseProgress + progressWithinCycle * cycleProgressSpan;
                 postMessage({ type: 'progress', value: Math.min(1.0, overallProgress) });
                 // --- END MODIFICATION ---


             } catch (error) {
                  console.error("Error during LM iteration:", iter, error);
                  postMessage({ type: 'error', message: `Error in LM iter ${iter}: ${error.message}` });
                   // Return the state before the error if possible
                  return { params: oldParams || initialParams, hklList: oldHklList || JSON.parse(JSON.stringify(hklList)), ss_res: ss_res, error: true, JtJ: finalJtJ, parameterInfo: parameterInfoForMainThread, algorithm: 'lm', fitFlags };
             }

        } // End of iteration loop

         // --- Ensure final progress update for the cycle ---
         const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
         postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });
         // --- End final progress update ---

         // --- Create serializable parameter info ---
          const parameterInfoForMainThread = paramMapping.map(m => ({
               name: m.name,
               scale: m.scale,
               isIntensity: m.isIntensity
          }));
          // --- End serializable info ---


        return {
             params,
             hklList: workingHklList,
             JtJ: finalJtJ,
             parameterInfo: parameterInfoForMainThread, // Send serializable info
             ss_res,
             algorithm: 'lm',
             fitFlags
        };
} // End refineParametersLM


async function refineParametersSA(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
     const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode);
     if (paramMapping.length === 0) {
          const finalProgress = (leBailCycle + 1) / totalLeBailCycles;
          postMessage({ type: 'progress', value: Math.min(1.0, finalProgress) });
         return { params: initialParams, hklList: hklList, ss_res: 0 };
     }

     const objective = (p_obj, hkl_list_obj) => {
         try {
             updateHklPositions(hkl_list_obj, p_obj, system);
             const netCalcPattern = calculatePattern(workerWorkingData.tth, hkl_list_obj, p_obj);
             const y_bkg = calculateTotalBackground(workerWorkingData.tth, p_obj);

             let scaleFactor = 1.0;
             if (refinementMode === 'le-bail') {
                  let num = 0, den = 0;
                  for (let i = 0; i < workerWorkingData.tth.length; i++) {
                       const w_i = workerWorkingData.weights[i];
                       const y_obs_net_i = workerWorkingData.intensity[i] - y_bkg[i];
                       const net_calc_i = netCalcPattern[i] || 0;
                       if (isFinite(w_i) && isFinite(y_obs_net_i) && isFinite(net_calc_i)) {
                            num += w_i * y_obs_net_i * net_calc_i;
                            den += w_i * net_calc_i * net_calc_i;
                       }
                  }
                  scaleFactor = (den > 1e-12) ? Math.max(0, num / den) : 1.0;
                  if (!isFinite(scaleFactor)) scaleFactor = 1.0;
             }

             const y_calc_total = numeric.add(numeric.mul(netCalcPattern, scaleFactor), y_bkg);

             let sum_w_res_sq = 0;
             for (let i = 0; i < workerWorkingData.tth.length; i++) {
                  const w_i = workerWorkingData.weights[i];
                  const res = workerWorkingData.intensity[i] - y_calc_total[i];
                  if (isFinite(w_i) && isFinite(res)) {
                      sum_w_res_sq += w_i * res * res;
                  } else {
                      return 1e12; // Return large cost if calculation fails
                  }
             }
             return isFinite(sum_w_res_sq) ? sum_w_res_sq : 1e12; // Ensure finite cost
         } catch (err) {
              console.warn("SA objective function error:", err);
              return 1e12; // Large cost on error
         }
     };

     // SA parameters
     let currentParams = JSON.parse(JSON.stringify(initialParams));
     let currentHklList = JSON.parse(JSON.stringify(hklList));
     let current_cost = objective(currentParams, currentHklList);

     let bestParams = JSON.parse(JSON.stringify(currentParams));
     let bestHklList = JSON.parse(JSON.stringify(currentHklList));
     let best_cost = current_cost;

     let T = (refinementMode === 'pawley') ? 0.1 : 1.0; // Start cooler for Pawley
     const T_min = 1e-7;
     const coolingRate = (maxIter > 1 && T > T_min) ? Math.pow(T_min / T, 1.0 / (maxIter -1)) : 0.99;

     const baseProgress = leBailCycle / totalLeBailCycles;
     const cycleProgressSpan = 1 / totalLeBailCycles;


     for (let step = 0; step < maxIter; step++) {
          try {
              const originalParams = JSON.parse(JSON.stringify(currentParams));
              const originalHklList = JSON.parse(JSON.stringify(currentHklList));

              // Perturb a random subset of parameters
              const paramsToChange = Math.max(1, Math.floor(paramMapping.length * 0.1));
              const changedIndices = new Set();
               while (changedIndices.size < paramsToChange && changedIndices.size < paramMapping.length) {
                    changedIndices.add(Math.floor(Math.random() * paramMapping.length));
               }

               changedIndices.forEach(p_idx => {
                    const mapping = paramMapping[p_idx];
                    const original_val = mapping.get(currentParams, currentHklList);
                    const step_scale = Math.max(0.01, T);
                    // Use smaller steps for Pawley intensities
                    const base_step_suggestion = (refinementMode === 'pawley' && mapping.isIntensity) ? 0.005 : 0.05;
                    const step_size = (mapping.step || base_step_suggestion) * step_scale * mapping.scale;
                    const random_step = (Math.random() - 0.5) * 2 * step_size;
                    const new_val = original_val + random_step;
                    mapping.set(currentParams, currentHklList, new_val);
               });


              const neighbor_cost = objective(currentParams, currentHklList);
              const delta_cost = neighbor_cost - current_cost;

              // Metropolis acceptance criterion
              const acceptance_prob = (T > 1e-9 && current_cost > 0) ? Math.exp(-delta_cost / (current_cost * T)) : 0;
              const should_accept = (delta_cost < 0) || (acceptance_prob > Math.random());

              if (should_accept) {
                  current_cost = neighbor_cost;
              } else {
                  currentParams = originalParams;
                  currentHklList = originalHklList;
              }

              if (current_cost < best_cost) {
                  best_cost = current_cost;
                  bestParams = JSON.parse(JSON.stringify(currentParams));
                  bestHklList = JSON.parse(JSON.stringify(currentHklList));
              }

              T = Math.max(T_min, T * coolingRate);

              const progressWithinCycle = (step + 1) / maxIter;
              const overallProgress = baseProgress + progressWithinCycle * cycleProgressSpan;
              postMessage({ type: 'progress', value: Math.min(1.0, overallProgress) });


          } catch (error) {
               console.error("Error during SA iteration:", step, error);
               postMessage({ type: 'error', message: `Error in SA iter ${step}: ${error.message}` });
               currentParams = JSON.parse(JSON.stringify(bestParams)); // Revert to best known
               currentHklList = JSON.parse(JSON.stringify(bestHklList));
               current_cost = best_cost;
               T = Math.max(T_min, T * coolingRate); // Still cool down
          }

     } // End of loop

     // --- Finalization ---
     updateHklPositions(bestHklList, bestParams, system);
     const finalNetCalcPattern = calculatePattern(workerWorkingData.tth, bestHklList, bestParams);
     const finalY_bkg = calculateTotalBackground(workerWorkingData.tth, bestParams);
     let finalScaleFactor = 1.0;
     if (refinementMode === 'le-bail') {
          let num = 0, den = 0;
          for (let i = 0; i < workerWorkingData.tth.length; i++) {
               const w_i = workerWorkingData.weights[i];
               const y_obs_net_i = workerWorkingData.intensity[i] - finalY_bkg[i];
               const net_calc_i = finalNetCalcPattern[i] || 0;
               if (isFinite(w_i) && isFinite(y_obs_net_i) && isFinite(net_calc_i)) {
                    num += w_i * y_obs_net_i * net_calc_i;
                    den += w_i * net_calc_i * net_calc_i;
               }
          }
          finalScaleFactor = (den > 1e-12) ? Math.max(0, num / den) : 1.0;
          if (!isFinite(finalScaleFactor)) finalScaleFactor = 1.0;
           bestHklList.forEach(hkl => {
               if (hkl) hkl.intensity *= finalScaleFactor;
           });
     }

     // --- Ensure final progress update for the cycle ---
     const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
     postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });
     // --- End final progress update ---

     // --- Create serializable parameter info ---
      const parameterInfoForMainThread = paramMapping.map(m => ({
           name: m.name,
           scale: m.scale,
           isIntensity: m.isIntensity
      }));
      // --- End serializable info ---


     return {
         params: bestParams,
         hklList: bestHklList,
         algorithm: 'sa',
         parameterInfo: parameterInfoForMainThread, // Send serializable info
         fitFlags,
         ss_res: best_cost
     };
} // End refineParametersSA


// ---
async function refineParametersPT(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode);
    if (paramMapping.length === 0) {
         const finalProgress = (leBailCycle + 1) / totalLeBailCycles;
         postMessage({ type: 'progress', value: Math.min(1.0, finalProgress) });
        return { params: initialParams, hklList: hklList, ss_res: 0 };
    }

    // --- PT Configuration ---
    const numReplicas = 8;
    const maxTemp = 1.0;
    const minTemp = 1e-5;
    const swapInterval = 10;

    // --- Objective function ---
    const objective = (p_obj, hkl_list_obj) => {
         try {
             updateHklPositions(hkl_list_obj, p_obj, system);
             const netCalcPattern = calculatePattern(workerWorkingData.tth, hkl_list_obj, p_obj);
             const y_bkg = calculateTotalBackground(workerWorkingData.tth, p_obj);

             let scaleFactor = 1.0;
             if (refinementMode === 'le-bail') {
                  let num = 0, den = 0;
                  for (let i = 0; i < workerWorkingData.tth.length; i++) {
                       const w_i = workerWorkingData.weights[i];
                       const y_obs_net_i = workerWorkingData.intensity[i] - y_bkg[i];
                       const net_calc_i = netCalcPattern[i] || 0;
                       if (isFinite(w_i) && isFinite(y_obs_net_i) && isFinite(net_calc_i)) {
                            num += w_i * y_obs_net_i * net_calc_i;
                            den += w_i * net_calc_i * net_calc_i;
                       }
                  }
                  scaleFactor = (den > 1e-12) ? Math.max(0, num / den) : 1.0;
                  if (!isFinite(scaleFactor)) scaleFactor = 1.0;
             }

             const y_calc_total = numeric.add(numeric.mul(netCalcPattern, scaleFactor), y_bkg);

             let sum_w_res_sq = 0;
             for (let i = 0; i < workerWorkingData.tth.length; i++) {
                  const w_i = workerWorkingData.weights[i];
                  const res = workerWorkingData.intensity[i] - y_calc_total[i];
                   if (isFinite(w_i) && isFinite(res)) {
                        sum_w_res_sq += w_i * res * res;
                   } else {
                       return 1e12; // High cost if calculation fails
                   }
             }
             return isFinite(sum_w_res_sq) ? sum_w_res_sq : 1e12;
         } catch (err) {
              console.warn("PT objective function error:", err);
              return 1e12;
         }
    };


    // --- Initialization ---
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

    // --- Progress Scaling Setup ---
    const baseProgress = leBailCycle / totalLeBailCycles;
    const cycleProgressSpan = 1 / totalLeBailCycles;
    // --- End Progress Scaling Setup ---


    // --- Main PT Loop ---
    for (let iter = 0; iter < maxIter; iter++) {
         try {
             // --- Part 1: Standard Monte Carlo step for each replica ---
             for (let i = 0; i < numReplicas; i++) {
                 let replica = replicas[i];
                 const originalParams = JSON.parse(JSON.stringify(replica.params));
                 const originalHklList = JSON.parse(JSON.stringify(replica.hklList));

                 // Perturb a random parameter
                 const p_idx = Math.floor(Math.random() * paramMapping.length);
                 const mapping = paramMapping[p_idx];
                 const original_val = mapping.get(replica.params, replica.hklList);
                 const step_scale = Math.max(0.01, replica.temp);
                  // Use smaller steps for Pawley intensities
                 const base_step_suggestion = (refinementMode === 'pawley' && mapping.isIntensity) ? 0.005 : 0.05;
                 const step_size = (mapping.step || base_step_suggestion) * step_scale * mapping.scale;
                 const random_step = (Math.random() - 0.5) * 2 * step_size;
                 const new_val = original_val + random_step;
                 mapping.set(replica.params, replica.hklList, new_val);

                 // Calculate new cost and decide whether to accept
                 const neighbor_cost = objective(replica.params, replica.hklList);
                 const delta_cost = neighbor_cost - replica.cost;

                 // Metropolis criterion, scaled by temperature
                  const acceptance_prob = (replica.temp > 1e-9 && replica.cost > 0) ? Math.exp(-delta_cost / (replica.cost * replica.temp)) : 0;
                 if (delta_cost < 0 || acceptance_prob > Math.random()) {
                     replica.cost = neighbor_cost; // Accept
                 } else {
                     replica.params = originalParams; // Reject
                     replica.hklList = originalHklList;
                 }

                 // Update the global best solution
                 if (replica.cost < bestOverallCost) {
                     bestOverallCost = replica.cost;
                     bestOverallParams = JSON.parse(JSON.stringify(replica.params));
                     bestOverallHklList = JSON.parse(JSON.stringify(replica.hklList));
                 }
             }

             // --- Part 2: Attempt swaps between adjacent replicas ---
             if (iter > 0 && iter % swapInterval === 0) {
                 for (let i = 0; i < numReplicas - 1; i++) {
                     const rep1 = replicas[i];
                     const rep2 = replicas[i + 1];

                     if (Math.abs(rep1.temp - rep2.temp) < 1e-9) continue;

                     const delta_beta = (1 / rep1.temp) - (1 / rep2.temp);
                     const delta_cost = rep1.cost - rep2.cost;
                     const acceptance_prob_swap = Math.exp(Math.min(50, delta_beta * delta_cost)); // Cap exponent

                     if (acceptance_prob_swap > Math.random()) {
                         [rep1.params, rep2.params] = [rep2.params, rep1.params];
                         [rep1.hklList, rep2.hklList] = [rep2.hklList, rep1.hklList];
                         [rep1.cost, rep2.cost] = [rep2.cost, rep1.cost];
                     }
                 }
             }

             // --- MODIFIED PROGRESS UPDATE ---
             const progressWithinCycle = (iter + 1) / maxIter;
             const overallProgress = baseProgress + progressWithinCycle * cycleProgressSpan;
             postMessage({ type: 'progress', value: Math.min(1.0, overallProgress) });
             // --- END MODIFICATION ---


         } catch (error) {
              console.error("Error during PT iteration:", iter, error);
              postMessage({ type: 'error', message: `Error in PT iter ${iter}: ${error.message}` });
              // Continue loop, relying on bestOverall state
         }

    } // End of loop

    // --- Finalization ---
     updateHklPositions(bestOverallHklList, bestOverallParams, system);
     const finalNetCalcPattern = calculatePattern(workerWorkingData.tth, bestOverallHklList, bestOverallParams);
     const finalY_bkg = calculateTotalBackground(workerWorkingData.tth, bestOverallParams);
     let finalScaleFactor = 1.0;
     if (refinementMode === 'le-bail') {
          let num = 0, den = 0;
          for (let i = 0; i < workerWorkingData.tth.length; i++) {
               const w_i = workerWorkingData.weights[i];
               const y_obs_net_i = workerWorkingData.intensity[i] - finalY_bkg[i];
               const net_calc_i = finalNetCalcPattern[i] || 0;
               if (isFinite(w_i) && isFinite(y_obs_net_i) && isFinite(net_calc_i)) {
                    num += w_i * y_obs_net_i * net_calc_i;
                    den += w_i * net_calc_i * net_calc_i;
               }
          }
          finalScaleFactor = (den > 1e-12) ? Math.max(0, num / den) : 1.0;
          if (!isFinite(finalScaleFactor)) finalScaleFactor = 1.0;

           bestOverallHklList.forEach(hkl => {
               if(hkl) hkl.intensity *= finalScaleFactor;
           });
     }

     // --- Ensure final progress update for the cycle ---
     const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
     postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });
     // --- End final progress update ---

     // --- Create serializable parameter info ---
      const parameterInfoForMainThread = paramMapping.map(m => ({
           name: m.name,
           scale: m.scale,
           isIntensity: m.isIntensity
      }));
      // --- End serializable info ---


    return {
        params: bestOverallParams,
        hklList: bestOverallHklList,
        algorithm: 'pt',
        parameterInfo: parameterInfoForMainThread, // Send serializable info
        fitFlags,
        ss_res: bestOverallCost
    };
} // End refineParametersPT


// --- Parameter Mapping (Helper) ---

function getParameterMapping(fitFlags, initialParams, hklList, refinementMode) {
    const mappings = [];

    const createMapping = (flag, name, defaultScale = 1.0, minVal = -Infinity, maxVal = Infinity, step = 0.2) => {
        if (!flag) return null;
        const initialValue = initialParams[name] ?? 0;
         const scale = Math.abs(initialValue) > 1e-9 ? Math.abs(initialValue) : defaultScale;

        return {
            name: name,
            scale: scale,
             step: step,
            isIntensity: false,
            get: (p_obj, hkl_list_obj) => (p_obj[name] ?? 0) / scale,
            set: (p_obj, hkl_list_obj, normalizedValue) => {
                let rawValue = normalizedValue * scale;
                if (rawValue < minVal) rawValue = minVal;
                if (rawValue > maxVal) rawValue = maxVal;
                 // Ensure the parameter exists before setting
                 if (p_obj && p_obj.hasOwnProperty(name)) {
                     p_obj[name] = rawValue;
                 }
            }
        };
    };


    // --- Pawley Intensity Parameters ---
     if (refinementMode === 'pawley' && hklList && workerWorkingData && workerWorkingData.tth && workerWorkingData.tth.length > 0) {
          const tthMin = workerWorkingData.tth[0];
          const tthMax = workerWorkingData.tth[workerWorkingData.tth.length - 1];

          hklList.forEach((hkl, index) => {
               if (hkl && hkl.tth && hkl.tth >= tthMin && hkl.tth <= tthMax) {
                    const hkl_name = `I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`;
                    const initialIntensity = (hkl.intensity !== undefined && hkl.intensity > 1e-6) ? hkl.intensity : 1000.0;
                    const scale = initialIntensity;

                    mappings.push({
                        name: hkl_name,
                        scale: scale,
                        step: 0.3, // Step suggestion for intensity
                        isIntensity: true, // Mark as intensity
                        index: index,
                        get: (p_obj, hkl_list_obj) => {
                            const intensity = hkl_list_obj?.[index]?.intensity ?? 0;
                            return intensity / scale;
                        },
                        set: (p_obj, hkl_list_obj, normalizedValue) => {
                             if (hkl_list_obj?.[index]) {
                                let rawValue = normalizedValue * scale;
                                hkl_list_obj[index].intensity = Math.max(0, rawValue); // Ensure >= 0
                             }
                        }
                    });
               }
          });
     }


    // --- Other Parameters ---
    const profileType = String(initialParams.profileType || "4");

    // Lattice Parameters
    mappings.push(createMapping(fitFlags.a, 'a', 4.0, 0.1, Infinity, 0.01));
    mappings.push(createMapping(fitFlags.b, 'b', 4.0, 0.1, Infinity, 0.01));
    mappings.push(createMapping(fitFlags.c, 'c', 6.0, 0.1, Infinity, 0.01));
    mappings.push(createMapping(fitFlags.alpha, 'alpha', 90.0, 0, 180, 0.05));
    mappings.push(createMapping(fitFlags.beta, 'beta', 90.0, 0, 180, 0.05));
    mappings.push(createMapping(fitFlags.gamma, 'gamma', 120.0, 0, 180, 0.05));

    // Instrumental
    mappings.push(createMapping(fitFlags.zeroShift, 'zeroShift', 0.01, -Infinity, Infinity, 0.1));

    // Profile Parameters
    if (profileType === "4") { // Simple pVoigt
        mappings.push(createMapping(fitFlags.GU, 'GU', 0.01, 0, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GV, 'GV', 0.01, -Infinity, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GW, 'GW', 0.01, 1e-6, Infinity, 0.05));
        mappings.push(createMapping(fitFlags.GP, 'GP', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.LX, 'LX', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.eta, 'eta', 0.5, 0, 1, 0.1));
        mappings.push(createMapping(fitFlags.shft, 'shft', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns, 'trns', 0.01, -Infinity, Infinity, 0.1));
    } else { // TCH (Profile 3)
        mappings.push(createMapping(fitFlags.U, 'U', 0.01, 0, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.V, 'V', 0.01, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.W, 'W', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.X, 'X', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.Y, 'Y', 0.01, 1e-6, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.SL, 'SL', 0.001, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.HL, 'HL', 0.001, -Infinity, Infinity, 0.1));
         if (profileType === "3") {
            mappings.push(createMapping(fitFlags.S400, 'S400', 0.1, -Infinity, Infinity, 0.2));
            mappings.push(createMapping(fitFlags.S040, 'S040', 0.1, -Infinity, Infinity, 0.2));
            mappings.push(createMapping(fitFlags.S004, 'S004', 0.1, -Infinity, Infinity, 0.2));
            mappings.push(createMapping(fitFlags.S220, 'S220', 0.1, -Infinity, Infinity, 0.2));
            mappings.push(createMapping(fitFlags.S202, 'S202', 0.1, -Infinity, Infinity, 0.2));
            mappings.push(createMapping(fitFlags.S022, 'S022', 0.1, -Infinity, Infinity, 0.2));
        }
    }

    // Background Parameters
    for (let i = 0; i < NUM_BACKGROUND_PARAMS; i++) {
        mappings.push(createMapping(fitFlags[`B${i}`], `B${i}`, i === 0 ? 100 : 1.0, -Infinity, Infinity, 0.3));
    }
    mappings.push(createMapping(fitFlags.hump_H, 'hump_H', 100.0, 0, Infinity, 0.2));
    mappings.push(createMapping(fitFlags.hump_P, 'hump_P', 30.0, -Infinity, Infinity, 0.1));
    mappings.push(createMapping(fitFlags.hump_W, 'hump_W', 5.0, 0.01, Infinity, 0.1));

    // Filter out null mappings (where flag was false)
    const paramMapping = mappings.filter(Boolean);
    return { paramMapping };
}



// --- 4. Worker Message Handler ---
self.onmessage = async function(e) {
    const {
        initialParams,
        fitFlags,
        workingData, // contains tth, intensity, weights, startIndex
        masterHklList,
        spaceGroupsData, // Pass the whole array if needed by HKL generation logic
        selectedSgNumber,
        system,
        maxIterations,
        algorithm,
        refinementMode
    } = e.data;

     // Store working data globally in the worker
     workerWorkingData = workingData;
     // If spaceGroups data is passed, potentially assign it globally if needed
     // self.spaceGroups = spaceGroupsData; // Uncomment if needed globally

    // Find the selected space group within the worker context
     const selectedSg = spaceGroups.find(sg => sg.number === selectedSgNumber);
     if (!selectedSg) {
         postMessage({ type: 'error', message: `Worker Error: Space group ${selectedSgNumber} not found.` });
         return;
     }


    let finalResults;
    let currentHklList = JSON.parse(JSON.stringify(masterHklList || [])); // Start with master list

    try {
        // --- Le Bail Mode ---
        if (refinementMode === 'le-bail') {
            const LE_BAIL_CYCLES = 4; // Or make configurable
            let currentParams = initialParams; // Start with initial params

            for (let cycle = 0; cycle < LE_BAIL_CYCLES; cycle++) {
                // --- REMOVED intermediate progress post ---
                // postMessage({ type: 'progress', value: (cycle + 0.5) / LE_BAIL_CYCLES, message: `Le Bail Cycle ${cycle + 1}/${LE_BAIL_CYCLES}` });

                let refinementResults;
                // --- Pass cycle info to refinement functions ---
                if (algorithm === 'lm') {
                    refinementResults = await refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, cycle, LE_BAIL_CYCLES);
                } else if (algorithm === 'sa') {
                    refinementResults = await refineParametersSA(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, cycle, LE_BAIL_CYCLES);
                } else { // pt
                    refinementResults = await refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, cycle, LE_BAIL_CYCLES);
                }
                // --- End Pass cycle info ---


                if (refinementResults && refinementResults.params && refinementResults.hklList && !refinementResults.error) {
                    currentParams = refinementResults.params; // Update params for next cycle/extraction
                    currentHklList = refinementResults.hklList; // Update HKL intensities
                    finalResults = refinementResults; // Store latest results
                } else {
                     throw new Error(`Refinement algorithm (${algorithm}) failed during Le Bail cycle ${cycle + 1}. ${refinementResults?.error ? 'See previous error.' : ''}`);
                }

                // --- Intensity Extraction Step ---
                const backgroundForExtraction = calculateTotalBackground(workerWorkingData.tth, currentParams);
                const expDataForExtraction = {
                    tth: workerWorkingData.tth,
                    intensity: workerWorkingData.intensity,
                    background: backgroundForExtraction
                };
                 leBailIntensityExtraction(expDataForExtraction, currentHklList, currentParams);


            } // End Le Bail cycle loop

        }
        // --- Pawley Mode ---
        else {
             let currentParams = initialParams; // Start with initial params

             // Initialize intensities to a constant value instead
             currentHklList.forEach(peak => {
                 if (peak) {
                     peak.intensity = 1000.0; // Or another reasonable starting value
                 }
             });

            // Run chosen algorithm ONCE for Pawley, refining intensities simultaneously
            // --- Pass default cycle info (0, 1) for Pawley ---
            if (algorithm === 'lm') {
                finalResults = await refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else if (algorithm === 'sa') {
                finalResults = await refineParametersSA(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { // pt
                finalResults = await refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            }
             // --- End Pass cycle info ---

             if (!finalResults || finalResults.error) {
                  throw new Error(`Refinement algorithm (${algorithm}) failed during Pawley fit. ${finalResults?.error ? 'See previous error.' : ''}`);
             }

        } // End Pawley Mode

        // --- Final Calculations & Posting Result ---
        if (!finalResults || !finalResults.params || !finalResults.hklList) {
            throw new Error("Refinement finished but produced invalid results.");
        }

        const finalParams = finalResults.params;
        const finalHklList = finalResults.hklList;

        // Calculate final stats using the results from the refinement
        const finalNetPatternForStats = calculatePattern(workerWorkingData.tth, finalHklList, finalParams);
        const finalBackgroundForStats = calculateTotalBackground(workerWorkingData.tth, finalParams);

        // Pass workerWorkingData to calculateStatistics
        const finalStats = calculateStatistics(
            workerWorkingData, // Use the data stored in the worker
            finalNetPatternForStats,
            finalResults.fitFlags || fitFlags, // Use flags from result if available
            finalBackgroundForStats,
            finalParams,
            finalHklList,
            refinementMode
        );


        // Prepare result object to send back
        const resultPayload = {
            params: finalParams,
            hklList: finalHklList,
            stats: finalStats,
            algorithm: algorithm,
            refinementMode: refinementMode,
            fitFlags: finalResults.fitFlags || fitFlags,
            parameterInfo: finalResults.parameterInfo || [], // Send serializable info
            // paramMapping removed
            JtJ: finalResults.JtJ || null,
            ss_res: finalResults.ss_res // Include final cost
        };

        // Post the final result
        postMessage({ type: 'result', results: resultPayload });

    } catch (error) {
        console.error("Worker Error during refinement:", error);
        postMessage({ type: 'error', message: `Refinement failed: ${error.message}` });
    } finally {

         workerWorkingData = null;
         hklIndexCache = {}; // Clear HKL cache after run? Maybe keep it?
    }
}; 

// --- error handler -
self.onerror = function(event) {
     console.error("Unhandled Worker Error:", event.message, event);
     postMessage({ type: 'error', message: `Unhandled Worker Error: ${event.message}` });
};

