// refinement_worker.js
// version 131, 29 july 2026, LM: Marquardt-scaled solve + lambda no longer
//   relaxed on a trust-region-clamped step (fixes stalling / uphill drift)
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
    .then(async data => { 
        globalThis.SG_DATABASE = data; 
        importScripts('sg_engine.js');

        if (typeof SG_ENGINE === 'undefined') {
            throw new Error('sg_engine.js did not load correctly.');
        }
        await initWebGPU();
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


let gpuDevice = null;
let computePipeline = null;
let gpuBuffers = { tth: null, hkl: null, out: null, read: null, tthSize: 0, hklSize: 0, outSize: 0 };
let reusableHklData = new Float32Array(1024);

async function initWebGPU() {
    if (!navigator.gpu) {
        console.warn("WebGPU not supported. Falling back to CPU.");
        return false;
    }
    try {
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) return false;
        gpuDevice = await adapter.requestDevice();
        
        gpuDevice.lost.then((info) => {
            console.warn("WebGPU device lost:", info);
            gpuDevice = null;
        });

        const wgslCode = await fetch('xrd_compute.wgsl').then(res => res.text());
        const shaderModule = gpuDevice.createShaderModule({ code: wgslCode });

        computePipeline = gpuDevice.createComputePipeline({
            layout: 'auto',
            compute: { module: shaderModule, entryPoint: 'main' },
        });
        return true;
    } catch (e) {
        console.warn("Failed to initialize WebGPU:", e);
        gpuDevice = null;
        return false;
    }
}

async function calculatePattern(tthAxis, hklList, params, outArr = null) {
    if (!gpuDevice || !computePipeline || tthAxis.length === 0) {
        return calculatePatternCPU(tthAxis, hklList, params, outArr);
    }

    const n_points = tthAxis.length;
    const pattern = outArr || new Float64Array(n_points);
    const deg2rad = Math.PI / 180;
    
    const lambda1 = params.lambda || 1.54056;
    const lambda2 = params.lambda2 || 0;
    const ratio21 = params.ratio || 0;
    const zeroShift = params.zeroShift || 0;
    const doubletEnabled = ratio21 > 1e-6 && lambda2 > 1e-6 && Math.abs(lambda1 - lambda2) > 1e-6;
    const lambdaRatio = doubletEnabled ? (lambda2 / lambda1) : 1.0;

    let totalPeaks = 0;
    hklList.forEach(peak => {
        if (!peak || !peak.intensity || peak.intensity <= 0 || !peak.tth || peak.tth < 0 || peak.tth > 180) return;
        totalPeaks++;
        if (doubletEnabled) {
            const sinTheta2 = Math.sin(peak.tth * deg2rad / 2.0) * lambdaRatio;
            if (Math.abs(sinTheta2) < 1) totalPeaks++;
        }
    });


    if (totalPeaks === 0) {
        pattern.fill(0);
        return pattern;
    }

    const requiredSize = totalPeaks * 16;
    if (reusableHklData.length < requiredSize) {
        reusableHklData = new Float32Array(requiredSize * 2); // Grow with headroom
    }
    const hklData = reusableHklData.subarray(0, requiredSize);

    let offset = 0;

    const writePrep = (intensity, prep) => {
        hklData[offset++] = intensity; hklData[offset++] = prep.x0;
        hklData[offset++] = prep.max_d; hklData[offset++] = prep.pedestal;
        hklData[offset++] = prep.type; hklData[offset++] = prep.asym_param || 0.0;
        hklData[offset++] = prep.eta || 0.0; hklData[offset++] = prep.fwhm_total || 0.0;
        
        if (prep.type === 2) {
            hklData[offset++] = prep.H_G_L; hklData[offset++] = prep.H_L_L;
            hklData[offset++] = prep.cL_L;  hklData[offset++] = prep.cG_L;
            hklData[offset++] = prep.H_G_R; hklData[offset++] = prep.H_L_R;
            hklData[offset++] = prep.cL_R;  hklData[offset++] = prep.cG_R;
        } else if (prep.type === 1) {
            hklData[offset++] = prep.H_G; hklData[offset++] = prep.H_L;
            hklData[offset++] = prep.cL; hklData[offset++] = prep.cG;
            hklData[offset++] = 0; hklData[offset++] = 0; hklData[offset++] = 0; hklData[offset++] = 0;
        } else {
            hklData[offset++] = prep.fwhm; hklData[offset++] = prep.fwhm;
            hklData[offset++] = prep.cL; hklData[offset++] = prep.cG;
            hklData[offset++] = 0; hklData[offset++] = 0; hklData[offset++] = 0; hklData[offset++] = 0;
        }
    };

    hklList.forEach(peak => {
        if (!peak || !peak.intensity || peak.intensity <= 0 || !peak.tth || peak.tth < 0 || peak.tth > 180) return;
        const basePos1 = peak.tth + zeroShift;
        const peakPos1 = basePos1 + calculatePeakShift(basePos1, params);
        writePrep(peak.intensity, prepareVoigt(basePos1, peakPos1, peak, params));

        if (doubletEnabled) {
            const sinTheta2 = Math.sin(peak.tth * deg2rad / 2.0) * lambdaRatio;
            if (Math.abs(sinTheta2) < 1) {
                const tth2 = 2 * Math.asin(sinTheta2) / deg2rad;
                const basePos2 = tth2 + zeroShift;
                const peakPos2 = basePos2 + calculatePeakShift(basePos2, params);
                writePrep(peak.intensity * ratio21, prepareVoigt(basePos2, peakPos2, peak, params));
            }
        }
    });

    try {
        const tthData = new Float32Array(tthAxis);
        if (!gpuBuffers.tth || gpuBuffers.tthSize < tthData.byteLength) {
            if (gpuBuffers.tth) gpuBuffers.tth.destroy();
            gpuBuffers.tth = gpuDevice.createBuffer({ size: tthData.byteLength, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
            gpuBuffers.tthSize = tthData.byteLength;
        }
        gpuDevice.queue.writeBuffer(gpuBuffers.tth, 0, tthData);

        const hklByteLength = Math.max(hklData.byteLength, 64); // Min 64 bytes (1 PeakPrep struct)
        if (!gpuBuffers.hkl || gpuBuffers.hklSize < hklByteLength) {
            if (gpuBuffers.hkl) gpuBuffers.hkl.destroy();
            gpuBuffers.hkl = gpuDevice.createBuffer({ size: hklByteLength, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
            gpuBuffers.hklSize = hklByteLength;
        }
        if (hklData.length > 0) gpuDevice.queue.writeBuffer(gpuBuffers.hkl, 0, hklData);

        if (!gpuBuffers.out || gpuBuffers.outSize < tthData.byteLength) {
            if (gpuBuffers.out) gpuBuffers.out.destroy();
            if (gpuBuffers.read) gpuBuffers.read.destroy();
            gpuBuffers.out = gpuDevice.createBuffer({ size: tthData.byteLength, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
            gpuBuffers.read = gpuDevice.createBuffer({ size: tthData.byteLength, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
            gpuBuffers.outSize = tthData.byteLength;
        }

        const bindGroup = gpuDevice.createBindGroup({
            layout: computePipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: { buffer: gpuBuffers.tth } },
                { binding: 1, resource: { buffer: gpuBuffers.hkl, size: hklByteLength } },
                { binding: 2, resource: { buffer: gpuBuffers.out } },
            ],
        });

        const commandEncoder = gpuDevice.createCommandEncoder();
        const passEncoder = commandEncoder.beginComputePass();
        passEncoder.setPipeline(computePipeline);
        passEncoder.setBindGroup(0, bindGroup);
        passEncoder.dispatchWorkgroups(Math.ceil(n_points / 256));
        passEncoder.end();

        commandEncoder.copyBufferToBuffer(gpuBuffers.out, 0, gpuBuffers.read, 0, tthData.byteLength);
        gpuDevice.queue.submit([commandEncoder.finish()]);
        
        await gpuBuffers.read.mapAsync(GPUMapMode.READ);
        const resultArray = new Float32Array(gpuBuffers.read.getMappedRange(0, tthData.byteLength));
        pattern.set(resultArray); 
        gpuBuffers.read.unmap();

        return pattern;
    } catch (e) {
        console.warn("GPU computation failed, falling back to CPU.", e);
        gpuDevice = null; 
        return calculatePatternCPU(tthAxis, hklList, params, outArr);
    }
}


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
        prep = { type: 1, x0, asym_param, H_G, H_L, eta,
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
        prep = { type: 2, x0, asym_param, H_G_L, H_L_L, H_G_R, H_L_R, eta,
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
        prep = { type: 3, x0, asym_param, fwhm, eta,
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
        // FIX: use the SIGNED delta. The old code took Math.abs(delta) inside
        // the correction, so BOTH flanks were stretched by the same factor and
        // SL/HL produced pure symmetric broadening -- no asymmetry whatsoever.
        // With the signed form one flank is broadened (factor > 1) and the
        // other sharpened (factor < 1), which is the point of the correction.
        const t = Math.max(-0.95, Math.min(0.95, p.asym_param * delta));
        delta /= (1.0 - t);   // 1 - t is in [0.05, 1.95], always positive
    }

    let hg, hl, cL, cG;
    if (p.type === 1) {
        hg = p.H_G; hl = p.H_L; cL = p.cL; cG = p.cG;
    } else if (p.type === 2) {
        if (delta < 0) { hg = p.H_G_L; hl = p.H_L_L; cL = p.cL_L; cG = p.cG_L; }
        else           { hg = p.H_G_R; hl = p.H_L_R; cL = p.cL_R; cG = p.cG_R; }
    } else {
        hg = hl = p.fwhm; cL = p.cL; cG = p.cG;
    }

    const dg = (delta / hg) * (delta / hg);
    const dl = (delta / hl) * (delta / hl);
    const g = Math.exp(-p.Cg * dg);
    const l = 1.0 / (1.0 + 4.0 * dl);
    // cL / cG are the area-correct mixing weights (see pvMixCoeffs),
    // renormalised so this expression is 1 at delta = 0.
    const v = cL * l + cG * g - p.pedestal;
    return v > 0 ? v : 0.0;
}

    function calculatePatternCPU(tthAxis, hklList, params, outArr = null) {
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


async function leBailIntensityExtraction(expData, hklList, params, scratch_calc = null) {
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
    // Any INDIVIDUAL reflection at zero is also lifted onto the revival floor
    // relative to the current maximum, for the same reason (see
    // LEBAIL_REVIVAL_FRACTION): a zero weight can never grow back.
    else {
        const floorNow = LEBAIL_REVIVAL_FRACTION * maxInt;
        if (floorNow > 0) {
            for (let i = 0; i < hklList.length; i++) {
                const p = hklList[i];
                if (p && !(p.intensity > floorNow)) p.intensity = floorNow;
            }
        }
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
    const total_calc = await calculatePattern(expData.tth, hklList, params, scratch_calc);

    // 3. Extract area per peak without storing M arrays of size N
    const currentCycleIntensities = new Float64Array(hklList.length);
    const currentCycleShapeAreas   = new Float64Array(hklList.length);
    const currentCycleVariances    = new Float64Array(hklList.length);
    // Area of the UNtruncated, un-pedestalled shape, accumulated alongside the
    // rendered one. The rendered area is what keeps the stored height
    // consistent with the pattern actually drawn; the analytic area is what
    // gives an integrated intensity free of the few-percent truncation bias.
    const currentCycleAnalyticAreas = new Float64Array(hklList.length);
    // Point variances. Prefer the weights the main thread computed (which carry
    // its own floor convention) over a bare Poisson count, so the extraction's
    // sigma is on the same footing as the weighted residual it came from.
    const pointWeights = expData.weights || null;

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
                if (current_fraction > 1.0) current_fraction = 1.0; // Prevent float drift explosion
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
            pend_var = (pointWeights && isFinite(pointWeights[i]) && pointWeights[i] > 0)
                       ? (1 / pointWeights[i])                  // sigma^2 as weighted
                       : Math.max(1, expData.intensity[i]);     // Poisson fallback
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
        currentCycleAnalyticAreas[j] = (prep1.analyticArea || 0)
                                     + (prep2 ? ratio21 * (prep2.analyticArea || 0) : 0);
    });

    // Mean of the positive heights, for the revival floor. Computed before
    // anything is written back so the floor does not chase its own tail.
    let heightSum = 0, heightCount = 0;
    for (let idx = 0; idx < hklList.length; idx++) {
        const sa = currentCycleShapeAreas[idx] || 0;
        if (sa <= 1e-12) continue;
        const h = (currentCycleIntensities[idx] || 0) / sa;
        if (h > 0) { heightSum += h; heightCount++; }
    }
    const meanHeight = heightCount > 0 ? heightSum / heightCount : 0;
    const revivalFloor = meanHeight > 0 ? LEBAIL_REVIVAL_FRACTION * meanHeight : 0;

    hklList.forEach((peak, idx) => {
        if (!peak) return;
        const raw_area       = currentCycleIntensities[idx] || 0;   // may be negative
        const shape_area     = currentCycleShapeAreas[idx] || 0;
        const analytic_area  = currentCycleAnalyticAreas[idx] || 0;

        // The modelling height. Floored, not clipped -- see the note on
        // LEBAIL_REVIVAL_FRACTION. A hard zero here is an absorbing state.
        const height = (shape_area > 1e-12) ? (raw_area / shape_area) : 0;
        peak.intensity = Math.max(revivalFloor, height);

        peak.shapeArea = shape_area;            // integral of the rendered Ka1+Ka2 shape
        peak.shapeAreaAnalytic = analytic_area; // same shape, untruncated

        // Reported, MEASURED quantities: deliberately not clipped, so that a
        // systematically absent reflection shows up as I/sigma near zero (or
        // negative) instead of being masked by the floor above.
        peak.I_integrated = raw_area;
        peak.I_sigma = Math.sqrt(Math.max(0, currentCycleVariances[idx] || 0));
        // Truncation-corrected integrated intensity: the same decomposition
        // result rescaled from the drawn area to the full analytic area. This
        // is the number to quote; peak.I_integrated is the number consistent
        // with the drawn pattern.
        peak.I_integrated_full = (shape_area > 1e-12)
            ? raw_area * (analytic_area / shape_area)
            : raw_area;
        peak.I_sigma_full = (shape_area > 1e-12)
            ? peak.I_sigma * (analytic_area / shape_area)
            : peak.I_sigma;
        delete peak.intensity_previous;
    });
}


//   Statistics  
function calculateStatistics(localWorkingData, netCalcPattern, fitFlags, finalBackground, params, hklList, refinementMode, system) {
    const y_obs = localWorkingData.intensity;
    const y_bkg = finalBackground;
    const weights = localWorkingData.weights;
    const N = y_obs.length;

    if (!y_obs || !netCalcPattern || !y_bkg || !weights ||
        N === 0 || y_obs.length !== netCalcPattern.length || y_obs.length !== y_bkg.length || y_obs.length !== weights.length) {
        console.error("Statistics calculation error: Mismatched or invalid array inputs.");
        return { r_p: -1, rwp: -1, chi2: -1, scaleFactor: 1, sum_w_res_sq: 0,
                 nParams: 0, nIntensities: 0, dof: 0 };
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


     const { paramMapping } = getParameterMapping(fitFlags || {}, params || {}, hklList || [], refinementMode || 'le-bail', system);
     const P_base = paramMapping.filter(m => !m.isIntensity).length;

    // ---------------------------------------------------------------------
    //  FIX: the extracted Le Bail intensities are FREE PARAMETERS and are now
    //  counted.
    //
    //  Only Pawley mode used to add them, because only Pawley puts them in the
    //  parameter mapping. But a Le Bail fit also determines one intensity per
    //  reflection from the same data -- they simply come from the decomposition
    //  instead of from the normal equations. Leaving them out made N - P far too
    //  large, so chi-square came out too small and every ESD derived from it was
    //  correspondingly optimistic.
    //
    //  Reflections are counted only if they actually overlap the fitted range
    //  (shapeArea > 0 once an extraction has run), since a reflection outside
    //  the window is determined by nothing and costs no degree of freedom.
    // ---------------------------------------------------------------------
    let nIntensities = 0;
    if ((refinementMode === 'pawley' || refinementMode === 'le-bail') &&
        hklList && localWorkingData && localWorkingData.tth && localWorkingData.tth.length > 0) {
        const tthMin = localWorkingData.tth[0];
        const tthMax = localWorkingData.tth[localWorkingData.tth.length - 1];
        nIntensities = hklList.filter(hkl => {
            if (!hkl || !hkl.tth) return false;
            if (typeof hkl.shapeArea === 'number') return hkl.shapeArea > 0;
            return hkl.tth >= tthMin && hkl.tth <= tthMax;
        }).length;
    }

    const P = P_base + nIntensities;
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
        sum_w_res_sq: isFinite(sum_w_res_sq) ? sum_w_res_sq : 0,
        nParams: P_base,
        nIntensities: nIntensities,
        dof: degreesOfFreedom
    };
}

//   Refinement Algorithms (LM, PT)  
// ---------------------------------------------------------------------------
//  LM tuning constants.
// ---------------------------------------------------------------------------
//  Finite differences. The escalation ladder is BOUNDED now: the old loop
//  multiplied the probe by 100 up to five times (1e8 x the base step). A width
//  parameter sitting in a numerical dead zone could be driven to a nonsensical
//  value, the pattern then changed wildly, `responded` flipped true, and the
//  resulting secant -- across eight orders of magnitude -- was accepted as a
//  derivative. Cap the probe at a few percent of the parameter's own magnitude:
//  large enough to escape a dead zone, small enough to still be a derivative.
const FD_BASE_FRACTION     = 1e-4;
const FD_MAX_PROBE_FRACTION = 0.02;
const FD_MAX_ATTEMPTS      = 4;
//  Convergence. The old threshold was 1e-9 relative, which in Le Bail mode can
//  never fire -- the intensity re-extraction at the top of the next iteration
//  moves the cost by far more than that -- so LM always burned maxIter. Require
//  a realistic relative drop, twice in a row, so a single flat step does not
//  end the fit.
//  Relative tolerance below which a Jacobian column counts as empty, i.e. the
//  parameter is not determined by the data at this point. Compared against the
//  LARGEST diagonal of JtJ, so it is scale-free.
const JTJ_RANK_TOL        = 1e-12;
const LM_COST_TOL         = 1e-7;
const LM_CONVERGED_REPEATS = 2;
//  Outer (extract-then-refine) monotonicity. Relative slack before an iteration
//  is judged to have made the TRUE, re-extracted cost worse, and how many such
//  rejections to tolerate before giving up. See the long note in the main loop.
const LM_OUTER_TOL         = 1e-6;
const LM_MAX_OUTER_REJECTS = 6;
//  Trust-region scaling on the step caps: where it starts, and how fast it grows
//  after an iteration that improved the re-extracted cost.
const LM_TRUST_INITIAL     = 0.05;
const LM_TRUST_GROWTH      = 1.6;
//  ---------------------------------------------------------------------------
//  How many decomposition passes to run when re-extracting, and whether to start
//  each re-extraction from flat intensities.
//
//  BUG FIX. A single pass does not settle, so the intensities carry HISTORY: the
//  cost at a given set of parameters depends on how the fit arrived there. That
//  makes the Le Bail objective path-dependent, and a path-dependent objective
//  cannot be minimised meaningfully -- the decomposition keeps adapting to cover
//  mispositioned peaks, so the cost can fall while the cell drifts away. Measured
//  on a synthetic cubic pattern, the SETTLED objective is emphatic (Rwp 9.4% at
//  the true cell against 91.3% at the aliased cell+zero-point minimum), but with
//  one pass per iteration the refiner slid into the aliased minimum anyway and
//  reported 81% -- lower than a settled extraction would ever give there.
//
//  Restarting from flat makes the extraction a deterministic function of the
//  parameters, which is what the outer monotonicity guard needs in order to mean
//  anything. The decomposition is a fixed-point iteration and converges in a
//  handful of passes, and the cost is small next to building n_params Jacobian
//  columns.
//  ---------------------------------------------------------------------------
const LEBAIL_PASSES_PER_ITER = 6;
const LEBAIL_RESET_EACH_ITER = true;

async function refineParametersLM(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
    const cloneParams = (p) => ({ ...p });
    const cloneHkl = (list) => list.map(h => (h ? { ...h } : null));
    
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode, system);
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
    calculateTotalBackground(workerWorkingData.tth, initialParams, workerBackgroundAnchors, scratch_bkg);
    const scratch_pattern = new Float64Array(n_points);

    const baseProgress = leBailCycle / totalLeBailCycles;
    const cycleProgressSpan = 1 / totalLeBailCycles;

    const refreshLeBailIntensities = async () => {
        if (refinementMode !== 'le-bail') return;
        enforceSymmetryConstraintsWorker(params, system);

        updateHklPositions(workingHklList, params, system);
        const bkg0 = scratch_bkg;
        
        if (LEBAIL_RESET_EACH_ITER) {



            for (let i = 0; i < workingHklList.length; i++) {
                if (workingHklList[i]) workingHklList[i].intensity = 100.0;
            }
        }
        const expData = { tth: workerWorkingData.tth, intensity: y_obs, background: bkg0, weights: workerWorkingData.weights };
        for (let pass = 0; pass < LEBAIL_PASSES_PER_ITER; pass++) {
            await leBailIntensityExtraction(expData, workingHklList, params, scratch_pattern);
        }
    };

    const calculateTotalPattern = async (targetArray) => {
        enforceSymmetryConstraintsWorker(params, system);




            updateHklPositions(workingHklList, params, system);

        const netCalcPattern = await calculatePattern(workerWorkingData.tth, workingHklList, params, scratch_pattern);
        for (let i = 0; i < n_points; i++) {
            targetArray[i] = netCalcPattern[i] + scratch_bkg[i];
        }


        return 1.0; 
    };

    const evaluateResiduals = async (target) => {
        await calculateTotalPattern(target);
        let cost = 0;
        for (let i = 0; i < n_points; i++) {
            residuals[i] = (y_obs[i] - target[i]) * sqrt_weights[i];
            if (!isFinite(residuals[i])) return null;
            cost += residuals[i] * residuals[i];
        }
        return isFinite(cost) ? cost : null;
    };

    const deadParamWarned = new Set();

    const buildNormalEquations = async () => {
        const jacobian_T = [];
        const savedIntensities = workingHklList.map(h => h ? h.intensity : 0);

        for (let p = 0; p < n_params; p++) {
            const mapping = paramMapping[p];
            const originalValue = mapping.get(params, workingHklList);

            const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
            const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;

            const magnitude = Math.max(Math.abs(originalValue), mapping.defaultScale || 0, 1e-9);
            const baseStep  = Math.max(1e-9, magnitude * FD_BASE_FRACTION);
            const maxProbe  = Math.max(baseStep, FD_MAX_PROBE_FRACTION * magnitude);

            let fd_step = baseStep;
            let fd_sign = 1;
            let responded = false;

            for (let attempt = 0; attempt < FD_MAX_ATTEMPTS; attempt++) {
                fd_sign = (originalValue + fd_step > maxV && originalValue - fd_step >= minV) ? -1 : 1;

                mapping.set(params, workingHklList, originalValue + fd_sign * fd_step);
                const applied  = mapping.get(params, workingHklList);
                const actual_h = applied - originalValue;   

                if (actual_h !== 0) {
                    await calculateTotalPattern(y_calc_total);
                    for (let i = 0; i < n_points; i++) {
                        if (!isFinite(y_calc_baseline[i]) || !isFinite(y_calc_total[i])) {
                            jacobian_col[i] = 0;
                        } else {
                            jacobian_col[i] = (y_calc_total[i] - y_calc_baseline[i])
                                              / actual_h * sqrt_weights[i];
                        }
                    }
                    for (let i = 0; i < n_points; i++) {
                        if (jacobian_col[i] !== 0) { responded = true; break; }
                    }
                }

                if (responded || fd_step >= maxProbe) break;
                fd_step = Math.min(maxProbe, fd_step * 30);
            }

            if (!responded) jacobian_col.fill(0);

            mapping.set(params, workingHklList, originalValue);
            workingHklList.forEach((h, idx) => { if (h) h.intensity = savedIntensities[idx]; });
            jacobian_T.push(jacobian_col.slice());
        }

        const JtJ = new Array(n_params);
        for (let i = 0; i < n_params; i++) JtJ[i] = new Array(n_params).fill(0);
        const Jtr = new Array(n_params).fill(0);
        let ok = true;

        for (let i = 0; i < n_params; i++) {
            const ci = jacobian_T[i];
            let s_r = 0;
            for (let k = 0; k < n_points; k++) s_r += ci[k] * residuals[k];
            if (!isFinite(s_r)) ok = false;
            Jtr[i] = s_r;

            for (let j = i; j < n_params; j++) {
                const cj = jacobian_T[j];
                let s = 0;
                for (let k = 0; k < n_points; k++) s += ci[k] * cj[k];
                if (!isFinite(s)) ok = false;
                JtJ[i][j] = s;
                JtJ[j][i] = s;
            }
        }
        if (!ok) return null;

        let maxDiagWarn = 0;
        for (let i = 0; i < n_params; i++) {
            const d = JtJ[i][i];
            if (isFinite(d) && d > maxDiagWarn) maxDiagWarn = d;
        }
        for (let i = 0; i < n_params; i++) {
            if (!(JtJ[i][i] > JTJ_RANK_TOL * maxDiagWarn) && !deadParamWarned.has(paramMapping[i].name)) {
                deadParamWarned.add(paramMapping[i].name);
                console.warn(`LM: parameter "${paramMapping[i].name}" has no effect on the calculated pattern (zero Jacobian column).`);
                postMessage({ type: 'warning', message: `Parameter "${paramMapping[i].name}" has no effect on the fit and will not move.` });
            }
        }
        return { JtJ, Jtr };
    };

    let lastAcceptedCost = Infinity;
    let convergedRepeats = 0;
    let lastProgressPct = -1;
    let bestTrueCost = Infinity, bestParams = null, bestHkl = null;
    let trustFactor = LM_TRUST_INITIAL;
    let outerRejects = 0;

    for (let iter = 0; iter < maxIter; iter++) {
         let oldParams = null;
         let oldHklList = null;
        try {
             await refreshLeBailIntensities();
             const cost = await evaluateResiduals(y_calc_baseline);
             if (cost === null) {
                 lambda = Math.min(1e9, lambda * 10);
                 console.warn(`LM iter ${iter}: Residual calculation failed. Increased lambda to ${lambda}.`);
                 continue;
             }


if (cost < bestTrueCost) {
                 bestTrueCost = cost;
                 bestParams = cloneParams(params);
                 bestHkl = cloneHkl(workingHklList);
                 trustFactor = Math.min(1.0, trustFactor * LM_TRUST_GROWTH);
             } else if (isFinite(bestTrueCost) && bestParams &&
                        cost > bestTrueCost * (1 + LM_OUTER_TOL)) {
                 if (++outerRejects > LM_MAX_OUTER_REJECTS) {
                     console.warn(`LM iter ${iter}: no further progress in the re-extracted cost; stopping.`);
                     break;
                 }
                 params = cloneParams(bestParams);
                 workingHklList = cloneHkl(bestHkl);

                 lambda = Math.min(1e9, lambda * 10);
                 trustFactor = Math.max(1e-3, trustFactor * 0.25);
                 lastAcceptedCost = bestTrueCost;
                 convergedRepeats = 0;
                 continue;
             }

             if (!isFinite(lastAcceptedCost)) lastAcceptedCost = cost;
             ss_res = cost;

             const normals = await buildNormalEquations();
             if (!normals) {
                  console.warn(`LM iter ${iter}: JtJ/Jtr contained non-finite entries.`);
                  finalJtJ = null;
                  lambda = Math.min(1e9, lambda * 5);
                  continue;
             }
             finalJtJ = normals.JtJ;
             const Jtr = normals.Jtr;

             let maxDiag = 0;
             for (let i = 0; i < n_params; i++) {
                 const d = finalJtJ[i][i];
                 if (isFinite(d) && d > maxDiag) maxDiag = d;
             }
             const active = [];
             for (let i = 0; i < n_params; i++) {
                 if (finalJtJ[i][i] > JTJ_RANK_TOL * maxDiag) active.push(i);
             }
             if (active.length === 0) {
                  console.warn(`LM iter ${iter}: no parameter has any effect on the pattern.`);
                  break;
             }

             const nA = active.length;
             const Dscale = new Float64Array(nA);
             for (let ii = 0; ii < nA; ii++) {
                 Dscale[ii] = Math.sqrt(finalJtJ[active[ii]][active[ii]]);
             }
             const A_lm = new Array(nA);
             for (let ii = 0; ii < nA; ii++) {
                 const row = new Float64Array(nA);
                 const di = Dscale[ii];
                 for (let jj = 0; jj < nA; jj++) {
                     row[jj] = finalJtJ[active[ii]][active[jj]] / (di * Dscale[jj]);
                 }
                 row[ii] += lambda;
                 A_lm[ii] = row;
             }
             const rhs_scaled = new Array(nA);
             for (let ii = 0; ii < nA; ii++) rhs_scaled[ii] = Jtr[active[ii]] / Dscale[ii];

             const step_active = solveLinearSystem(A_lm, rhs_scaled);
             if (!step_active) {
                  console.warn(`LM iter ${iter}: normal equations are singular. Increasing lambda.`);
                  lambda = Math.min(1e9, lambda * 10);
                  continue;
             }
             const p_step = new Array(n_params).fill(0);
             for (let ii = 0; ii < nA; ii++) p_step[active[ii]] = step_active[ii] / Dscale[ii];

             if (p_step.some(v => !isFinite(v))) {
                 console.warn(`LM iter ${iter}: Step calculation resulted in NaN/Infinity. Increasing lambda.`);
                 lambda = Math.min(1e9, lambda * 5);
                 continue;
             }

             const p_current = paramMapping.map(m => m.get(params, workingHklList));

             let stepScale = 1.0;
             for (let idx = 0; idx < p_step.length; idx++) {
                 const m = paramMapping[idx];
                 const cap = trustFactor * Math.min(
                     (m.maxStepAbs !== undefined) ? m.maxStepAbs : Infinity,
                     Math.max(2 * Math.abs(p_current[idx]), 4 * (m.defaultScale || 1e-9))
                 );
                 const want = Math.abs(p_step[idx]);
                 if (want > cap && want > 0) stepScale = Math.min(stepScale, cap / want);
             }
             const p_step_clamped = p_step.map(v => v * stepScale);

             let predicted = 0;
             for (let i = 0; i < n_params; i++) {
                 predicted += 2 * p_step_clamped[i] * Jtr[i];
                 let q = 0;
                 for (let j = 0; j < n_params; j++) q += finalJtJ[i][j] * p_step_clamped[j];
                 predicted -= p_step_clamped[i] * q;
             }

            oldParams = cloneParams(params);
             oldHklList = cloneHkl(workingHklList);

             const p_new = p_current.map((v, i) => v + p_step_clamped[i]);
             paramMapping.forEach((m, i) => m.set(params, workingHklList, p_new[i]));

             await calculateTotalPattern(y_calc_total);
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

             const actual = cost - new_cost;
             const rho = (predicted > 0 && isFinite(predicted))
                         ? actual / predicted
                         : (actual > 0 ? 1 : -1);

             let stepAccepted = false;
             if (new_cost < cost && isFinite(new_cost)) {
                 stepAccepted = true;
                 if (rho > 0.75 && stepScale > 0.99) lambda = Math.max(1e-9, lambda / 3);
                 else if (rho < 0.25)                lambda = Math.min(1e9, lambda * 2);
             } else {
                 params = oldParams;
                 workingHklList = oldHklList;
                 lambda = Math.min(1e9, lambda * 3);
             }

             const overallProgress = Math.min(1.0, baseProgress + ((iter + 1) / maxIter) * cycleProgressSpan);
             const pct = Math.floor(overallProgress * 100);
             if (pct !== lastProgressPct) {
                 lastProgressPct = pct;
                 postMessage({ type: 'progress', value: overallProgress });
             }

             if (stepAccepted) {
                 const denom = Math.max(Math.abs(lastAcceptedCost), 1.0);
                 const relDrop = Math.abs(lastAcceptedCost - new_cost) / denom;
                 lastAcceptedCost = new_cost;
                 ss_res = new_cost;
                 if (relDrop < LM_COST_TOL) {
                     if (++convergedRepeats >= LM_CONVERGED_REPEATS) break;
                 } else {
                     convergedRepeats = 0;
                 }
             }

         } catch (error) {
              console.error("Error during LM iteration:", iter, error);
              postMessage({ type: 'error', message: `Error in LM iter ${iter}: ${error.message}` });
              const errorParamInfo = paramMapping.map(m => ({
                  name: m.name,
                  scale: 1.0,
                  typicalMagnitude: m.scale,
                  isIntensity: m.isIntensity
              }));
              return { params: oldParams || params || initialParams,
                       hklList: oldHklList || workingHklList || JSON.parse(JSON.stringify(hklList)),
                       ss_res: ss_res, error: true, JtJ: finalJtJ,
                       parameterInfo: errorParamInfo, algorithm: 'lm', fitFlags };
         }
    }

    try {
        if (bestParams && isFinite(bestTrueCost)) {
            const currentCost = await evaluateResiduals(y_calc_baseline);
            if (currentCost === null || currentCost > bestTrueCost) {
                params = bestParams;
                workingHklList = bestHkl;
            }
        }
    } catch (err) {
        console.warn("LM: could not compare the final point with the best one.", err);
    }

    try {
        await refreshLeBailIntensities();
        const finalCost = await evaluateResiduals(y_calc_baseline);
        if (finalCost !== null) ss_res = finalCost;
        const finalNormals = await buildNormalEquations();
        if (finalNormals) finalJtJ = finalNormals.JtJ;
    } catch (err) {
        console.warn("LM: could not rebuild the normal equations at the final point; ESDs come from the last iteration.", err);
    }

     const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
     postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });

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
// ---------------------------------------------------------------------------
//  How often PT re-extracts the Le Bail intensities, in iterations.
//
//  BUG FIX. Holding the intensities fixed for a whole LM iteration is correct
//  (see the long note in refineParametersLM), but that argument does NOT carry
//  over to a Monte Carlo run of thousands of moves. Intensities frozen at the
//  INITIAL parameters are the decomposition of the observed profile against the
//  starting cell, so the objective's minimum sits exactly where the search
//  started -- which makes a global search over the cell unable to find anything
//  else, the one job PT exists to do. Re-extract periodically, per replica, at
//  that replica's own parameters.
//
//  Every refresh changes the objective function, so every recorded cost has to
//  be requoted on the new footing at the same time.
// ---------------------------------------------------------------------------
const LEBAIL_PT_REFRESH_INTERVAL = 40;

async function refineParametersPT(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode, system);
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
    calculateTotalBackground(workerWorkingData.tth, initialParams, workerBackgroundAnchors, scratch_bkg);
    const scratch_pattern = new Float64Array(n_points);

    const objective = async(p_obj, hkl_list_obj) => {
         try {
        
 enforceSymmetryConstraintsWorker(p_obj, system); 
            updateHklPositions(hkl_list_obj, p_obj, system);
            const netCalcPattern = await calculatePattern(workerWorkingData.tth, hkl_list_obj, p_obj, scratch_pattern);
             let sum_w_res_sq = 0;
             for (let i = 0; i < workerWorkingData.tth.length; i++) {
                  const w_i = workerWorkingData.weights[i];
                  const res = workerWorkingData.intensity[i] - (netCalcPattern[i] + scratch_bkg[i]);             



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

    // Re-extract this replica's Le Bail intensities at its own parameters, then
    // requote its cost against the new intensities. In Pawley mode the
    // intensities are refined parameters, so there is nothing to extract and
    // this just refreshes the cost.
    const refreshReplica = async(rep) => {


if (refinementMode === 'le-bail') {
            enforceSymmetryConstraintsWorker(rep.params, system);
            updateHklPositions(rep.hklList, rep.params, system);
            await leBailIntensityExtraction(
                { tth: workerWorkingData.tth, intensity: workerWorkingData.intensity,
                  background: scratch_bkg, weights: workerWorkingData.weights },
                rep.hklList, rep.params, scratch_pattern);
        }
       
       
       
        rep.cost = await objective(rep.params, rep.hklList);
        return rep.cost;
    };

    // Templates, captured before anything mutates them.
    const paramTemplate = JSON.parse(JSON.stringify(initialParams));
    const hklTemplate   = JSON.parse(JSON.stringify(hklList));

    // Seed the intensities once from the starting parameters.
    const seed = { params: initialParams, hklList: hklList, cost: Infinity, temp: 0 };
    const initialCost = await refreshReplica(seed);

    const temperatures = Array.from({ length: numReplicas }, (_, i) =>
        maxTemp * Math.pow(minTemp / maxTemp, i / (numReplicas - 1 || 1))
    );

    // One reference energy scale shared by BOTH acceptance tests.
    // Previously the local Metropolis step divided by the replica's own CURRENT
    // cost (relative scale) while the replica-exchange test used the raw cost
    // difference (absolute scale). With chi-square values of 1e5-1e7 counts the
    // swap exponent saturated, so exchanges were accepted either always or never
    // and the run degenerated into 8 independent annealers rather than parallel
    // tempering. A fixed scale also keeps the local acceptance ratio well
    // defined -- a cost-dependent denominator breaks detailed balance.
    const costScale = Math.max(1e-9, Math.abs(initialCost));

    let replicas = temperatures.map(temp => ({
        params: JSON.parse(JSON.stringify(initialParams)),
        hklList: JSON.parse(JSON.stringify(hklList)),
        cost: initialCost,
        temp: temp
    }));

    // PERF FIX: the best state is stored as the mapped PARAMETER VECTOR, not as
    // a deep clone of params + hklList. Only mapped scalars are ever perturbed,
    // so the vector plus the templates reconstructs the state exactly (slaved
    // parameters and peak positions are recomputed from it). The old code ran
    // `JSON.parse(JSON.stringify(replica.hklList))` on every improvement, which
    // early in a run is most moves of most replicas -- a full JSON round-trip of
    // every reflection, dominating the whole search.
    let bestVector = paramMapping.map(m => m.get(initialParams, hklList));
    let bestOverallCost = initialCost;

    // Scratch replica used to requote the best cost after a refresh.
    const probe = {
        params: JSON.parse(JSON.stringify(paramTemplate)),
        hklList: JSON.parse(JSON.stringify(hklTemplate)),
        cost: Infinity, temp: 0
    };
    const requoteBest = async() => {
        paramMapping.forEach((m, i) => m.set(probe.params, probe.hklList, bestVector[i]));
        return await refreshReplica(probe);
    };

    const baseProgress = leBailCycle / totalLeBailCycles;
    const cycleProgressSpan = 1 / totalLeBailCycles;
    let lastProgressPct = -1;
    let consecutiveErrors = 0;

    for (let iter = 0; iter < maxIter; iter++) {
         try {
             // Periodic Le Bail re-extraction (see LEBAIL_PT_REFRESH_INTERVAL).
             if (refinementMode === 'le-bail' && iter > 0 &&
                 iter % LEBAIL_PT_REFRESH_INTERVAL === 0) {
                 for (let i = 0; i < numReplicas; i++) await refreshReplica(replicas[i]);
                 // The objective just changed, so the stored best cost is no
                 // longer comparable with the replicas'. Requote it.
                 bestOverallCost = await requoteBest();
             }

             for (let i = 0; i < numReplicas; i++) {
                 let replica = replicas[i];

                 const p_idx = Math.floor(Math.random() * paramMapping.length);
                 const mapping = paramMapping[p_idx];
                 const original_val = mapping.get(replica.params, replica.hklList);

                 // A move perturbs exactly ONE scalar, so rolling it back only
                 // needs that scalar: everything else objective() touches
                 // (symmetry-slaved params, peak positions) is recomputed from
                 // scratch on the next call. In Pawley mode the perturbed scalar
                 // may itself be an intensity, but mapping.set restores it.
                 const step_scale = Math.max(0.01, replica.temp);
                 const step_size = (mapping.step || 0.05) * step_scale * mapping.scale;
                 const random_step = (Math.random() - 0.5) * 2 * step_size;
                 const proposed = original_val + random_step;

                 // FIX: an out-of-bounds proposal is REJECTED, not clamped.
                 // mapping.set clamps, which makes the proposal distribution
                 // asymmetric, breaks detailed balance and piles samples up
                 // against the bound.
                 const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
                 const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;
                 if (!(proposed >= minV && proposed <= maxV)) continue;

                 mapping.set(replica.params, replica.hklList, proposed);

                 const neighbor_cost = await objective(replica.params, replica.hklList);
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
                     for (let k = 0; k < paramMapping.length; k++) {
                         bestVector[k] = paramMapping[k].get(replica.params, replica.hklList);
                     }
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
              // The old handler posted an error EVERY iteration, so a persistent
              // fault produced maxIter error toasts and kept burning CPU. Give up
              // after a few consecutive failures.
              console.error("Error during PT iteration:", iter, error);
              consecutiveErrors++;
              if (consecutiveErrors >= 5) {
                  postMessage({ type: 'error', message: `PT aborted at iteration ${iter}: ${error.message}` });
                  break;
              }
         }

    }

     // Reconstruct the best state from the stored vector and re-extract its
     // intensities, so the returned hklList belongs to the returned parameters.
     const bestOverallParams  = JSON.parse(JSON.stringify(paramTemplate));
     const bestOverallHklList = JSON.parse(JSON.stringify(hklTemplate));
     paramMapping.forEach((m, i) => m.set(bestOverallParams, bestOverallHklList, bestVector[i]));
     const finalCost = await refreshReplica({ params: bestOverallParams, hklList: bestOverallHklList,
                                        cost: Infinity, temp: 0 });

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
        ss_res: isFinite(finalCost) ? finalCost : bestOverallCost
    };
}

/**
 * PT is a global search: it stops at the best SAMPLED point, which is not a
 * stationary point, and it builds no normal-equations matrix -- so a PT-only run
 * reported no ESDs at all. Follow it with a short Levenberg-Marquardt run from
 * the PT optimum to converge properly and to produce JtJ.
 */
async function polishWithLM(ptResults, fitFlags, hklList, system, refinementMode, maxIterations) {
    if (!ptResults || !ptResults.params || !ptResults.hklList) return ptResults;
    const polishIters = Math.max(10, Math.min(40, Math.round(maxIterations / 5)));
    try {
        postMessage({ type: 'progress', value: 0.98, message: 'Converging with Levenberg-Marquardt...' });
        const lm = await refineParametersLM(
            JSON.parse(JSON.stringify(ptResults.params)), fitFlags, polishIters,
            ptResults.hklList, system, refinementMode, 0, 1);
        if (lm && !lm.error && lm.params && lm.hklList) {
            return Object.assign({}, lm, { algorithm: 'pt+lm' });
        }
        // LM failed outright; keep PT's answer, but pass its normal-equations
        // matrix through if it managed to build one, so ESDs are still available.
        if (lm && lm.JtJ) {
            return Object.assign({}, ptResults, {
                JtJ: lm.JtJ, parameterInfo: lm.parameterInfo, algorithm: 'pt'
            });
        }
    } catch (err) {
        console.warn("LM polish after PT failed; returning the PT result.", err);
    }
    return ptResults;
}


//   Parameter Mapping (Helper)  
function getParameterMapping(fitFlags, initialParams, hklList, refinementMode, system) {
    const mappings = [];

    // -----------------------------------------------------------------------
    //  Symmetry-slaved Stephens terms are DROPPED up front.
    //
    //  enforceSymmetryConstraintsWorker overwrites them from their parent term
    //  on every pattern evaluation, so a slaved parameter's finite-difference
    //  probe is undone before the pattern is calculated: the column comes out
    //  identically zero, the escalation ladder runs its full length for nothing
    //  (several full pattern evaluations per slaved term per iteration -- four
    //  of them in a cubic tch_aniso fit), and the parameter still consumes a
    //  degree of freedom and reports a meaningless ESD.
    //
    //  M_HKL in _calculateWidths_TCH also simply does not reference some terms
    //  in the higher-symmetry branches, so those are dead too.
    // -----------------------------------------------------------------------
    const slaved = new Set();
    switch (system) {
        case 'cubic':
            ['S040', 'S004', 'S202', 'S022'].forEach(n => slaved.add(n));
            break;
        case 'tetragonal':
            ['S040', 'S022'].forEach(n => slaved.add(n));
            break;
        case 'hexagonal': case 'trigonal': case 'rhombohedral':
            // M_HKL uses only S400, S004, S202 in this branch.
            ['S040', 'S220', 'S022'].forEach(n => slaved.add(n));
            break;
    }

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
    const createMapping = (flag, name, defaultScale = 1.0, minVal = -Infinity, maxVal = Infinity, step = 0.2, maxStepAbs = Infinity) => {
        if (!flag) return null;
        if (slaved.has(name)) {
            console.warn(`Refinement flag set for "${name}", which is fixed by ${system} symmetry; skipping.`);
            postMessage({ type: 'warning', message: `Parameter "${name}" is determined by ${system} symmetry and was excluded from the refinement.` });
            return null;
        }
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
            // Absolute, PHYSICAL cap on a single LM step for this parameter.
            // Used only as a safeguard against absurd steps; lambda does the
            // actual step-length control. See the note in refineParametersLM.
            maxStepAbs: maxStepAbs,
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
                        maxStepAbs: Infinity,
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
    //  A LOOSE lower bound on the Caglioti constant term.
    //
    //  Not a resolution floor -- that lives in applyFwhmFloor now -- but a
    //  containment bound. Once the widths are pinned at the resolvable-FWHM
    //  floor, GW has no effect on the pattern at all, and an inert parameter
    //  with no lower bound drifts to nonsense (measured: -78732). The LM solve
    //  now holds undetermined directions, which is the real fix; this bound is
    //  the second line of defence, and it is set ten times looser than the
    //  sigma^2 that corresponds to the floor, so a genuinely negative GW -- which
    //  is physical in GSAS when GU is large -- is never obstructed.
    // ---------------------------------------------------------------------
    const sigmaSqAtFloor = Math.pow(MIN_PROFILE_FWHM_DEG / GSAS_GAUSSIAN_TO_DEG, 2);
    const GW_FLOOR = -10 * sigmaSqAtFloor;
    //  Absolute step caps for the profile coefficients. Without one, the cap is
    //  purely relative (2 x |value|), which permits exponential growth.
    const WIDTH_MAX_STEP = Math.max(25, 5 * sigmaSqAtFloor);

    // ---------------------------------------------------------------------
    // FIX: the per-parameter width floors (GW_MIN, LX_MIN) that used to live
    // here are gone.
    //
    // They bounded the WRONG QUANTITY. The Gaussian width comes from
    //     sigma^2 = GU tan^2(th) + GV tan(th) + GW + GP/cos^2(th),
    // so a bound on GW alone is simultaneously too weak -- GV < 0 can still
    // collapse the peak at low angle -- and too strong, because it ignores what
    // the other terms already contribute. Worse, a HARD bound plus the internal
    // `Math.max(1e-4, ...)` clamp created a dead zone: once a width parameter
    // reached its floor the finite-difference column came out identically zero
    // and the parameter was frozen there for the rest of the fit, unable to
    // recover even after the cell and zero-point were corrected.
    //
    // The floor now lives where the physics is: applyFwhmFloor bounds the
    // COMBINED FWHM at two data steps, smoothly (sqrt(t^2 + f^2)), scaling both
    // components equally so the Gaussian/Lorentzian ratio -- and hence the TCH
    // eta -- is untouched, and softPositive keeps every derivative alive. The
    // bounds below are therefore only physical sanity limits.
    // ---------------------------------------------------------------------

    //  createMapping(flag, name, defaultScale, minVal, maxVal, ptStep, maxStepAbs)
    //
    //  NOTE on defaultScale: it must be the parameter's REAL typical magnitude.
    //  The UI ships GU = 2, GW = 5, LX = 1, so the old blanket 0.01 was two to
    //  three orders of magnitude low. It matters for any parameter that starts
    //  at exactly zero -- GV, GP, SL, HL, the Stephens terms -- because `scale`
    //  then falls back to it and both the finite-difference probe and the LM
    //  step cap are derived from it.
    //
    //  defaultScale is the parameter's TYPICAL magnitude; it sizes the PT
    //  proposal and the finite-difference probe for a parameter that starts at
    //  zero, and it sets the relative floor on the LM step cap. zeroShift's used
    //  to be 0.01, which combined with step=0.1 gave it an LM budget of
    //  0.004 deg per iteration -- the tightest in the whole model, on the one
    //  parameter most likely to need a large first move.
    //
    //  maxStepAbs is a physical limit on one LM step: half an angstrom on a cell
    //  edge, a couple of degrees on a cell angle, 0.2 deg on the zero-point.
    mappings.push(createMapping(fitFlags.a, 'a', 4.0, 0.1, Infinity, 0.01, 0.5));
    mappings.push(createMapping(fitFlags.b, 'b', 4.0, 0.1, Infinity, 0.01, 0.5));
    mappings.push(createMapping(fitFlags.c, 'c', 6.0, 0.1, Infinity, 0.01, 0.5));
    mappings.push(createMapping(fitFlags.alpha, 'alpha', 90.0, 1, 179, 0.05, 2.0));
    mappings.push(createMapping(fitFlags.beta, 'beta', 90.0, 1, 179, 0.05, 2.0));
    mappings.push(createMapping(fitFlags.gamma, 'gamma', 120.0, 1, 179, 0.05, 2.0));
    mappings.push(createMapping(fitFlags.zeroShift, 'zeroShift', 0.05, -Infinity, Infinity, 0.1, 0.2));

    //  Width lower bounds are now only sanity limits: the resolvable-FWHM floor
    //  is applied smoothly to the COMBINED width in applyFwhmFloor, and
    //  softPositive keeps the derivative alive even at the bound, so sitting on
    //  a bound no longer kills the Jacobian column.
    if (profileType === "simple_pvoigt") {
        mappings.push(createMapping(fitFlags.GU, 'GU', 1.0, 0, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GV, 'GV', 1.0, -Infinity, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GW, 'GW', 1.0, GW_FLOOR, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GP, 'GP', 1.0, -Infinity, Infinity, 0.1, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.LX, 'LX', 1.0, 0, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.eta, 'eta', 0.5, 0, 1, 0.1, 0.25));
        mappings.push(createMapping(fitFlags.shft, 'shft', 0.1, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns, 'trns', 0.1, -Infinity, Infinity, 0.1));
    } 
    else if (profileType === "tch_aniso") {
        mappings.push(createMapping(fitFlags.U, 'U', 1.0, 0, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.V, 'V', 1.0, -Infinity, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.W, 'W', 1.0, GW_FLOOR, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.X, 'X', 1.0, 0, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.Y, 'Y', 1.0, 0, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.SL, 'SL', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.HL, 'HL', 0.01, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.S400, 'S400', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S040, 'S040', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S004, 'S004', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S220, 'S220', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S202, 'S202', 0.1, -Infinity, Infinity, 0.2));
        mappings.push(createMapping(fitFlags.S022, 'S022', 0.1, -Infinity, Infinity, 0.2));
    }
    else if (profileType === "split_pvoigt") {
        mappings.push(createMapping(fitFlags.GU_L, 'GU_L', 1.0, 0, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GV_L, 'GV_L', 1.0, -Infinity, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GW_L, 'GW_L', 1.0, GW_FLOOR, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.LX_L, 'LX_L', 1.0, 0, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GU_R, 'GU_R', 1.0, 0, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GV_R, 'GV_R', 1.0, -Infinity, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.GW_R, 'GW_R', 1.0, GW_FLOOR, Infinity, 0.05, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.LX_R, 'LX_R', 1.0, 0, Infinity, 0.2, WIDTH_MAX_STEP));
        mappings.push(createMapping(fitFlags.eta_split, 'eta_split', 0.5, 0, 1, 0.1, 0.25));
        mappings.push(createMapping(fitFlags.shft_split, 'shft_split', 0.1, -Infinity, Infinity, 0.1));
        mappings.push(createMapping(fitFlags.trns_split, 'trns_split', 0.1, -Infinity, Infinity, 0.1));
    }
   
    const paramMapping = mappings.filter(Boolean);
    return { paramMapping };
}


// 4. Worker Message Handler  
self.onmessage = async function(e) {
    try {
        await initPromise;
    } catch (initErr) {
        postMessage({ type: 'error', message: `Worker initialization failed: ${initErr && initErr.message ? initErr.message : initErr}` });
        return;
    }

    const {
        initialParams, fitFlags, workingData, masterHklList,
        selectedSgNumber, selectedSgQuery, system, maxIterations,
        algorithm, refinementMode, backgroundAnchors 
    } = e.data;

    workerWorkingData = workingData;
    workerBackgroundAnchors = backgroundAnchors; 
    setMinProfileFwhmFromAxis(workerWorkingData && workerWorkingData.tth);
    const selectedSg = SG_ENGINE.resolve(selectedSgQuery) || SG_ENGINE.resolve(selectedSgNumber);
    if (!selectedSg) {
        postMessage({ type: 'error', message: `Worker Error: Space group ${selectedSgNumber} could not be resolved.` });
        return;
    }

    let finalResults;
    let currentHklList = JSON.parse(JSON.stringify(masterHklList || []));
    let currentParams = initialParams;
    currentParams.system = system;
    
    try {
        //   Le Bail Mode  
        if (refinementMode === 'le-bail') {
            currentHklList.forEach(peak => { if (peak) peak.intensity = 100.0; });
            try {
                const bkgSeed = calculateTotalBackground(workerWorkingData.tth, currentParams, workerBackgroundAnchors);
                enforceSymmetryConstraintsWorker(currentParams, system);
                updateHklPositions(currentHklList, currentParams, system);
                for (let seedCycle = 0; seedCycle < 5; seedCycle++) {
                    await leBailIntensityExtraction(
                        { tth: workerWorkingData.tth,
                          intensity: workerWorkingData.intensity,
                          background: bkgSeed,
                          weights: workerWorkingData.weights },
                        currentHklList, currentParams);
                }
            } catch (seedErr) {
                console.warn("Le Bail seeding pass failed; continuing from flat intensities.", seedErr);
            }

            if (algorithm === 'lm') {
                finalResults = await refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { 
                finalResults = await refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
                finalResults = await polishWithLM(finalResults, fitFlags, currentHklList, system, refinementMode, maxIterations);
            }
            if (!finalResults || finalResults.error) {
                throw new Error(`Refinement algorithm (${algorithm}) failed during Le Bail fit. ${finalResults?.error ? 'See previous error.' : ''}`);
            }
        }
        //   Pawley Mode  
        else { 
             postMessage({ type: 'progress', value: 0.05, message: 'Initializing Pawley intensities...' });
             try {
                updateHklPositions(currentHklList, currentParams, system);
                const backgroundForInit = calculateTotalBackground(workerWorkingData.tth, currentParams, workerBackgroundAnchors);
                const expDataForInit = {
                    tth: workerWorkingData.tth,
                    intensity: workerWorkingData.intensity,
                    background: backgroundForInit,
                    weights: workerWorkingData.weights
                };
                 await leBailIntensityExtraction(expDataForInit, currentHklList, currentParams);
             } catch (initError) {
                  console.warn("Could not initialize Pawley intensities using Le Bail extraction:", initError);
                  currentHklList.forEach(peak => { if (peak) peak.intensity = 1000.0; });
             }
            
            if (algorithm === 'lm') {
                finalResults = await refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { 
                finalResults = await refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
                finalResults = await polishWithLM(finalResults, fitFlags, currentHklList, system, refinementMode, maxIterations);
            }

             if (!finalResults || finalResults.error) {
                  throw new Error(`Refinement algorithm (${algorithm}) failed during Pawley fit. ${finalResults?.error ? 'See previous error.' : ''}`);
             }
        }

        if (!finalResults || !finalResults.params || !finalResults.hklList) {
            throw new Error("Refinement finished but produced invalid results.");
        }

        const finalParams = finalResults.params;
        const finalHklList = finalResults.hklList;

        const finalNetPatternForStats = await calculatePattern(workerWorkingData.tth, finalHklList, finalParams);
        const finalBackgroundForStats = calculateTotalBackground(workerWorkingData.tth, finalParams, workerBackgroundAnchors);

        const finalStats = calculateStatistics(
            workerWorkingData,
            finalNetPatternForStats,
            finalResults.fitFlags || fitFlags,
            finalBackgroundForStats,
            finalParams,
            finalHklList,
            refinementMode,
            system
        );

        const resultPayload = {
            params: finalParams,
            hklList: finalHklList,
            stats: finalStats,
            algorithm: algorithm,
            effectiveAlgorithm: finalResults.algorithm || algorithm,
            refinementMode: refinementMode,
            fitFlags: finalResults.fitFlags || fitFlags,
            parameterInfo: finalResults.parameterInfo || [],
            JtJ: finalResults.JtJ || null,
            ss_res: finalResults.ss_res,
            dof: finalStats.dof,
            nIntensities: finalStats.nIntensities
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