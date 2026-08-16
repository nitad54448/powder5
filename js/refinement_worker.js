// refinement_worker.js
//version 143, 3 aug 2026... changes in the monoclinic settings, all accepted. added marching cubes 
//version 145, 15 august 2026: WebGPU removed from the fit path (CPU/f64 only)
//   kernel (arrayLength no longer reports a stale, larger buffer) and the
//   readback staging buffer is now serialised against overlapping calls.
//version 133, 31 july 2026, changed to cctbx v5, the json files holds the Harker data for further analysis 
// Harker planes saved in the text report
// version 131, 29 july 2026, LM: Marquardt-scaled solve + lambda no longer
//   relaxed on a trust-region-clamped step (fixes stalling / uphill drift)
// version 130, 26 july 2026, fixed some format in read data, changed the zoom
// version 129, 25 july 2026, Fixed TCH error
// version 123 (fixes: init race, tetragonal/triclinic HKL generation, LM step clamp)
// Spline Background in nov 2025, version 115
// ---------------------------------------------------------------------------
//  Shared modules. Loaded SYNCHRONOUSLY and FIRST, before anything else in
//  this file runs, because every profile and geometry routine below now lives
//  in them. Load order matters: profile.js reads constants.js's declarations
//  as bare globals.
//
//  Anything that used to be declared twice -- once here and once in
//  powder5.html -- now has exactly one definition. See the header of
//  constants.js for why that is a correctness issue and not a tidiness one.
// ---------------------------------------------------------------------------
importScripts('constants.js', 'crystal.js', 'profile.js');


const initPromise = fetch('../sg/all.json')
    .then(res => {
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        return res.json();
    })
    .then(async all => { 
        // Format uniformly: ensure SG_DATABASE has a 'space_groups' property
        if (all.space_groups) {
            globalThis.SG_DATABASE = all;
        } else if (Array.isArray(all)) {
            const space_groups = {};
            all.forEach(sg => { space_groups[sg.number] = sg; });
            globalThis.SG_DATABASE = { space_groups };
        } else {
            globalThis.SG_DATABASE = { space_groups: all };
        }
        
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


// ---------------------------------------------------------------------------
//  NO WEBGPU IN THE FIT PATH (v145).
//
//  The pattern kernel, its device/adapter setup, the peak-prep packing, the
//  per-workgroup range table, the buffer and bind-group caches, the mapAsync
//  serialisation queue and xrd_compute.wgsl are all gone. The refinement runs
//  on the CPU, in f64, always.
//
//  WHY, given the shader was measurably faster per evaluation:
//
//    * It was a SECOND implementation of the model, and the whole point of
//      the v144 consolidation was that a second implementation is a place for
//      the two to drift apart silently -- which had already happened once.
//      The shader could only be kept honest by a test that cannot run in this
//      environment.
//    * It returned f32 into an f64 pipeline. Every derivative, every ESD and
//      therefore every refined parameter depended on whether WebGPU happened
//      to initialise, and a transient device loss switched precision MID-RUN
//      because the catch block set gpuDevice = null and fell through to the
//      CPU. Two users with the same data and the same version could get
//      different numbers and nothing in the report would say why.
//    * Its reason for existing was the finite-difference Jacobian, which
//      called it once per parameter per iteration. The Pawley intensities are
//      now analytic, so the call count collapsed and with it most of the
//      benefit.
//
//  WebGPU is still used, where the arithmetic is genuinely heavy and the
//  result is not a refined number with an uncertainty attached to it: charge
//  flipping (charge_flipping_worker.js) and the Wyckoff swarm
//  (swarm_wyckoff.js).
// ---------------------------------------------------------------------------


// The HKL list (indices, multiplicities, allowed reflections) is generated
// on the main thread and shipped in with the refinement request, so the
// worker no longer duplicates that logic. SG_ENGINE is still needed here
// for SG_ENGINE.resolve() when the request arrives.
function enforceSymmetryConstraintsWorker(params, system) {
    if (!params) return;

    // --- cell metric -------------------------------------------------------
    // updateHklPositions DERIVES the slaved edges: for cubic it reads only
    // `a`, and for the uniaxial systems only `a` and `c`. So b (and c) can sit
    // at a stale value for a whole refinement without the pattern ever
    // noticing -- and they did. Refining `a` on a cubic cell from 8.14626
    // returned a = 8.13000 with b = c = 8.14626 still on the starting value.
    //
    // The FIT was right, because the calculation never looked at them. But the
    // REPORT is wrong, and so is anything downstream that does look: the
    // charge-flipping job is handed {a, b, c, alpha, beta, gamma} and builds
    // its metric from all six, so it would have solved on a cell that is not
    // the one that was refined, silently.
    //
    // Mirroring them here -- on every pattern evaluation, which is where the
    // symmetry constraints already live -- keeps the reported cell identical
    // to the cell that was actually used. The rules below match exactly what
    // updateHklPositions assumes; nothing is invented, and monoclinic and
    // triclinic are deliberately left alone because their unique-axis
    // convention is decided elsewhere (see monoclinicUniqueAxis).
    if (params.a > 0) {
        switch (system) {
            case 'cubic':
                params.b = params.a; params.c = params.a;
                break;
            case 'tetragonal': case 'hexagonal': case 'trigonal': case 'rhombohedral':
                params.b = params.a;
                break;
        }
    }

    // --- Stephens anisotropic broadening terms -----------------------------
    if (params.profileType !== "tch_aniso") return;
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
 * Solves A x = b for a dense SYMMETRIC POSITIVE DEFINITE A by Cholesky
 * factorisation, falling back to LU with partial pivoting if the factorisation
 * runs into a non-positive pivot.
 *
 * WHY CHOLESKY. The only system this function is ever asked to solve is the
 * Marquardt-scaled normal-equation matrix
 *
 *     A = D^-1 (J^T J) D^-1 + lambda I,
 *
 * which is symmetric and, for lambda > 0, positive definite by construction.
 * LU with partial pivoting ignores both facts: it does twice the arithmetic,
 * it destroys the symmetry while pivoting, and -- the part that actually
 * mattered -- it happily returns a finite "solution" for an indefinite matrix.
 * An indefinite A means the normal equations have been corrupted (a NaN that
 * slipped through, a Jacobian column that is not what it claims to be), and
 * the step computed from it points in an arbitrary direction. Cholesky fails
 * loudly on exactly that case, which is the behaviour the LM loop wants: it
 * raises lambda and tries again.
 *
 * A is READ ONLY here (unlike the previous version, which consumed it), so the
 * caller can retry with a larger lambda without rebuilding the matrix.
 *
 * @param {Array<Float64Array|number[]>} A  n x n, symmetric.
 * @param {Float64Array|number[]} b
 * @returns {number[]|null} Solution, or null if A is not usable.
 */
function solveLinearSystem(A, b) {
    const n = A.length;
    if (n === 0) return [];

    const x = solveCholesky(A, b, n);
    if (x) return x;
    return solveLU(A, b, n);
}

/**
 * In-place-free Cholesky: A = L L^T, then two triangular solves.
 * @returns {number[]|null} null if A is not positive definite.
 */
function solveCholesky(A, b, n) {
    // Packed lower triangle, row-major: L[i*(i+1)/2 + j] for j <= i.
    const L = new Float64Array((n * (n + 1)) / 2);

    for (let i = 0; i < n; i++) {
        const rowA = A[i];
        const bi = (i * (i + 1)) / 2;
        for (let j = 0; j <= i; j++) {
            const bj = (j * (j + 1)) / 2;
            let s = rowA[j];
            if (!isFinite(s)) return null;
            for (let k = 0; k < j; k++) s -= L[bi + k] * L[bj + k];
            if (i === j) {
                // A tolerance relative to the diagonal entry, not an absolute
                // one: the scaled matrix has a unit diagonal plus lambda, so
                // anything at or below rounding of the original entry means
                // the direction is not determined by the data.
                const floor = 1e-14 * Math.max(1, Math.abs(rowA[j]));
                if (!(s > floor)) return null;
                L[bi + j] = Math.sqrt(s);
            } else {
                L[bi + j] = s / L[bj + j];
            }
        }
    }

    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        const bi = (i * (i + 1)) / 2;
        let s = b[i];
        for (let k = 0; k < i; k++) s -= L[bi + k] * y[k];
        y[i] = s / L[bi + i];
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let s = y[i];
        for (let k = i + 1; k < n; k++) s -= L[(k * (k + 1)) / 2 + i] * x[k];
        x[i] = s / L[(i * (i + 1)) / 2 + i];
    }
    for (let i = 0; i < n; i++) if (!isFinite(x[i])) return null;
    return Array.from(x);
}

/**
 * LU with partial pivoting. Retained as the fallback for the case Cholesky
 * rejects, and for any caller that hands over a matrix that is not SPD.
 * Works on a copy, so the caller's matrix survives.
 */
function solveLU(A, b, n) {
    const M = new Array(n);
    for (let i = 0; i < n; i++) M[i] = Float64Array.from(A[i]);
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = b[i];

    for (let col = 0; col < n; col++) {
        let piv = col, maxAbs = Math.abs(M[col][col]);
        for (let r = col + 1; r < n; r++) {
            const v = Math.abs(M[r][col]);
            if (v > maxAbs) { maxAbs = v; piv = r; }
        }
        if (!(maxAbs > 0) || !isFinite(maxAbs)) return null;   // singular column

        if (piv !== col) {
            const tr = M[piv]; M[piv] = M[col]; M[col] = tr;
            const tx = x[piv]; x[piv] = x[col]; x[col] = tx;
        }

        const d = M[col][col];
        for (let r = col + 1; r < n; r++) {
            const f = M[r][col] / d;
            if (f === 0) continue;
            M[r][col] = 0;
            for (let c = col + 1; c < n; c++) M[r][c] -= f * M[col][c];
            x[r] -= f * x[col];
        }
    }

    for (let r = n - 1; r >= 0; r--) {
        let s = x[r];
        for (let c = r + 1; c < n; c++) s -= M[r][c] * x[c];
        const d = M[r][r];
        if (!isFinite(d) || d === 0) return null;
        x[r] = s / d;
    }

    for (let i = 0; i < n; i++) if (!isFinite(x[i])) return null;
    return Array.from(x);
}

//   2. Define Constants & Global Worker State  
// ---------------------------------------------------------------------------
//  MOVED. The profile-window multipliers, MIN_PROFILE_FWHM_DEG,
//  setMinProfileFwhmFromAxis, softPositive, the PV_* pseudo-Voigt
//  normalisation constants, pvMixCoeffs, lowerBound and the Le Bail revival
//  floor now live in constants.js, which is loaded by importScripts at the top
//  of this file AND by <script src> in powder5.html.
//
//  They used to be declared here AND in powder5.html. Two independent copies
//  of a physical constant is not a style problem: if one is edited and the
//  other is not, the preview and the fit compute different profile functions,
//  the refinement converges, the report is self-consistent, and the numbers
//  are wrong -- with nothing in the program able to detect it.
// ---------------------------------------------------------------------------


let workerWorkingData = null; // To store the sliced data sent from the main thread
let workerBackgroundAnchors = []; 


// ---------------------------------------------------------------------------
//  MOVED to profile.js: createMonotonicCubicSplineInterpolator,
//  calculateTotalBackground and the CPU pattern calculation
//  (calculatePatternCPU, plus buildPeakContributions / accumulatePattern that
//  it and the GPU upload now share).
//
//  These were the second copies of code that also lived in powder5.html; the
//  spline had a third copy in spline_background.js. profile.js is imported at
//  the top of this file, so every name below is already in scope.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//  MOVED. updateHklPositions now lives in crystal.js, together with the metric
//  tensors, cell volume and the monoclinic unique-axis convention. The copy
//  that used to be here was byte-identical to the one in powder5.html.
//
//  Note for anyone tempted to "fix" the monoclinic branch: it is correct. The
//  general triclinic reciprocal metric reduces to the closed monoclinic form
//  EXACTLY for unique axis a, b and c alike (agreement to 1 part in 1e15,
//  asserted in test/test_refactor_equivalence.js). The unique axis is carried
//  by whichever angle is not 90, which is the caller's contract; use
//  normaliseCell() from crystal.js to guarantee it.
// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------
//  MOVED. calculatePeakShift, the width calculators, getPeakFWHM,
//  prepareVoigt, evalVoigt and getPseudoVoigtArea now live in profile.js.
//
//  There were two copies of this, ~200 lines each, here and in powder5.html,
//  and they HAD ALREADY DIVERGED: this one carried a reciprocal-precompute
//  optimisation in evalVoigt the main thread never received, and the main
//  thread carried getPseudoVoigtArea this one never had. profile.js is the
//  merge of the two, and test/test_refactor_equivalence.js asserts it is
//  bit-for-bit identical to what this file used to compute.
// ---------------------------------------------------------------------------

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

    // 2. Compute total calculated pattern with current peak intensities.
    //    The contributions are built once and reused for the per-peak walk
    //    below, so prepareVoigt runs once per reflection per extraction pass
    //    instead of twice.
    const contributions = buildPeakContributions(expData.tth, hklList, params);
    // scratch_calc is optional and IS null on the seeding pass, which runs
    // before the refinement has allocated its scratch buffers.
    const total_calc = accumulatePattern(expData.tth, hklList, contributions,
                                         scratch_calc || new Float64Array(n_points));

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
        // Number.isFinite, not truthiness: tth === 0 is falsy but is a position.
        if (!hasRenderablePosition(peak)) return;

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


     // quiet = true: this call COUNTS parameters for the degrees of freedom. The
     // refinement has already run and already warned about every excluded
     // parameter; repeating the whole set now is noise the user cannot act on.
     const { paramMapping } = getParameterMapping(fitFlags || {}, params || {}, hklList || [],
                                                  refinementMode || 'le-bail', system, true);
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
            if (!hkl || !Number.isFinite(hkl.tth)) return false;
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

    // -----------------------------------------------------------------------
    //  THE BACKGROUND IS FIXED, AND THE REPORT SAYS SO.
    //
    //  It is interpolated through anchor points the user places, evaluated
    //  once before the refinement starts and held constant for every residual
    //  and every Jacobian column. That is a deliberate choice -- an anchored
    //  spline has no meaningful curvature model to refine against and letting
    //  it float would absorb the weak reflections -- but it is not a free one,
    //  and reporting it is the difference between a limitation and a hidden
    //  assumption:
    //
    //    * No ESD anywhere in the output includes background uncertainty. The
    //      covariance comes from J^T J, and the background contributes no
    //      column to J, so every sigma is conditional on the anchors being
    //      exactly right.
    //    * The anchors are not counted in P, so the degrees of freedom, and
    //      hence chi-square and every ESD scaled by it, are optimistic by
    //      roughly the number of anchors.
    //    * Any error in the anchors is absorbed wholesale into the extracted
    //      intensities, because the decomposition works on y_obs - y_bkg.
    // -----------------------------------------------------------------------
    const nBackgroundAnchors = Array.isArray(workerBackgroundAnchors) ? workerBackgroundAnchors.length : 0;

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
        dof: degreesOfFreedom,
        backgroundRefined: false,
        nBackgroundAnchors,
        backgroundNote: nBackgroundAnchors > 0
            ? `Background: fixed spline through ${nBackgroundAnchors} user anchor point`
              + `${nBackgroundAnchors === 1 ? '' : 's'}, not refined. Its ${nBackgroundAnchors} `
              + `degree${nBackgroundAnchors === 1 ? '' : 's'} of freedom are not counted in P, and no `
              + `reported ESD includes background uncertainty.`
            : 'Background: none defined (zero baseline), not refined.'
    };
}

//   Refinement Algorithms (LM, PT)
// ---------------------------------------------------------------------------
//  MOVED. Every LM, finite-difference and Le Bail tuning constant
//  (FD_BASE_FRACTION, FD_MAX_PROBE_FRACTION, FD_MAX_ATTEMPTS, JTJ_RANK_TOL,
//  LM_COST_TOL, LM_CONVERGED_REPEATS, LM_OUTER_TOL, LM_MAX_OUTER_REJECTS,
//  LM_TRUST_INITIAL, LM_TRUST_GROWTH, LEBAIL_PASSES_PER_ITER,
//  LEBAIL_RESET_EACH_ITER, LEBAIL_FLAT_SEED) now lives in constants.js, where
//  each one carries its numerical or physical justification rather than a bare
//  literal.
// ---------------------------------------------------------------------------

function refineParametersLM(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
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

    const refreshLeBailIntensities = () => {
        if (refinementMode !== 'le-bail') return;
        enforceSymmetryConstraintsWorker(params, system);

        updateHklPositions(workingHklList, params, system);
        const bkg0 = scratch_bkg;
        
        if (LEBAIL_RESET_EACH_ITER) {



            for (let i = 0; i < workingHklList.length; i++) {
                if (workingHklList[i]) workingHklList[i].intensity = LEBAIL_FLAT_SEED;
            }
        }
        const expData = { tth: workerWorkingData.tth, intensity: y_obs, background: bkg0, weights: workerWorkingData.weights };

        // Run the decomposition TO CONVERGENCE rather than for a fixed six
        // passes. Six was chosen as the count that settles "an ordinary
        // pattern", which means it is wasteful on an easy one and short on a
        // heavily overlapped one -- and the overlapped case is the one where
        // being short matters, because the LM then fits its profile parameters
        // against intensities that are still moving. Starting from flat keeps
        // the result a deterministic function of the parameters, which is what
        // LEBAIL_RESET_EACH_ITER exists to guarantee; stopping on a measured
        // change rather than a counter just makes it a settled one.
        let prev = null;
        for (let pass = 0; pass < LEBAIL_MAX_PASSES; pass++) {
            leBailIntensityExtraction(expData, workingHklList, params, scratch_pattern);

            const cur = new Float64Array(workingHklList.length);
            for (let i = 0; i < workingHklList.length; i++) {
                cur[i] = workingHklList[i] ? (workingHklList[i].intensity || 0) : 0;
            }
            if (prev && pass + 1 >= LEBAIL_PASSES_PER_ITER) {
                let num = 0, den = 0;
                for (let i = 0; i < cur.length; i++) {
                    num += Math.abs(cur[i] - prev[i]);
                    den += Math.abs(cur[i]);
                }
                if (den > 0 && num / den < LEBAIL_PASS_TOL) break;
            }
            prev = cur;
        }
    };



    const calculateTotalPattern = (targetArray, requiresHklUpdate = true) => {
        enforceSymmetryConstraintsWorker(params, system);

        if (requiresHklUpdate) {
            updateHklPositions(workingHklList, params, system);
        }

        const netCalcPattern = calculatePatternCPU(workerWorkingData.tth, workingHklList, params, scratch_pattern);
      
        for (let i = 0; i < n_points; i++) {
            targetArray[i] = netCalcPattern[i] + scratch_bkg[i];
        }


        return 1.0; 
    };

    const evaluateResiduals = (target) => {
        calculateTotalPattern(target);
        let cost = 0;
        for (let i = 0; i < n_points; i++) {
            residuals[i] = (y_obs[i] - target[i]) * sqrt_weights[i];
            if (!isFinite(residuals[i])) return null;
            cost += residuals[i] * residuals[i];
        }
        return isFinite(cost) ? cost : null;
    };

    const deadParamWarned = new Set();

    // ===================================================================
    //  THE JACOBIAN
    //
    //  Every column is stored as a SLICE, { start, values }, covering only
    //  the range of data points the parameter can affect. Profile, cell and
    //  zero-point parameters move every peak, so their slice is the whole
    //  axis; a Pawley intensity moves exactly one reflection, so its slice is
    //  that reflection's profile window and is typically a few dozen points
    //  out of ten thousand.
    //
    //  THIS IS THE WHOLE POINT. The old code computed every column by finite
    //  difference: perturb the parameter, recompute the ENTIRE pattern (all
    //  peaks, all points), and
    //  subtract. For a Pawley fit with n reflections that is n full pattern
    //  evaluations per iteration to obtain n columns whose derivative is
    //  d(y_i)/d(I_j) = profile_j(x_i) -- a quantity that is analytic, exact,
    //  local, and which buildPeakContributions has already prepared. It also
    //  made J^T J a dense n x n accumulation over all N points, when in truth
    //  two reflections that do not overlap have an inner product of exactly
    //  zero.
    //
    //  Both facts are now used: intensity columns are analytic, and the inner
    //  product runs over the INTERSECTION of two slices, which for
    //  non-overlapping reflections is empty. The finite-difference machinery
    //  below is still there, and still correct, for the parameters that
    //  genuinely need it.
    // ===================================================================

    /** Inner product of two column slices, over their overlap only. */
    const colDot = (a, b) => {
        const lo = Math.max(a.start, b.start);
        const hi = Math.min(a.start + a.values.length, b.start + b.values.length);
        let acc = 0;
        for (let k = lo; k < hi; k++) acc += a.values[k - a.start] * b.values[k - b.start];
        return acc;
    };

    /** Inner product of a column slice with a full-length vector. */
    const colDotFull = (a, v) => {
        let acc = 0;
        const end = a.start + a.values.length;
        for (let k = a.start; k < end; k++) acc += a.values[k - a.start] * v[k];
        return acc;
    };

    const buildNormalEquations = () => {
        const columns = new Array(n_params).fill(null);
        const savedIntensities = workingHklList.map(h => h ? h.intensity : 0);

        // Contributions at the CURRENT parameters, shared by every analytic
        // column. Built once; prepareVoigt is the expensive part of a pattern
        // evaluation and this is the only place it needs to run for them.
        const haveIntensityParams = paramMapping.some(m => m.isIntensity);
        const contributions = haveIntensityParams
            ? buildPeakContributions(workerWorkingData.tth, workingHklList, params)
            : null;

        for (let p = 0; p < n_params; p++) {
            const mapping = paramMapping[p];

            // ---- analytic column: a Pawley intensity ----------------------
            if (mapping.isIntensity && contributions) {
                const col = intensityDerivativeColumn(
                    workerWorkingData.tth, contributions, mapping.index);
                if (col) {
                    // The residual carries sqrt(w), so the Jacobian must too.
                    for (let k = 0; k < col.values.length; k++) {
                        col.values[k] *= sqrt_weights[col.start + k];
                    }
                    columns[p] = col;
                } else {
                    // The reflection falls outside the fitted range: no data
                    // point constrains it. A zero-length slice is the honest
                    // representation and it is what makes the parameter show
                    // up as undetermined rather than as an ESD of zero.
                    columns[p] = { start: 0, values: new Float64Array(0) };
                }
                continue;
            }

            // ---- finite-difference column --------------------------------
            const originalValue = mapping.get(params, workingHklList);
            const isStructural = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].includes(mapping.name);

            const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
            const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;

            const magnitude = Math.max(Math.abs(originalValue), mapping.defaultScale || 0, 1e-9);
            const baseStep  = Math.max(1e-9, magnitude * FD_BASE_FRACTION);
            const maxProbe  = Math.max(baseStep, FD_MAX_PROBE_FRACTION * magnitude);

            let fd_step = baseStep;
            let fd_sign = 1;
            let responded = false;
            jacobian_col.fill(0);

            for (let attempt = 0; attempt < FD_MAX_ATTEMPTS; attempt++) {
                fd_sign = (originalValue + fd_step > maxV && originalValue - fd_step >= minV) ? -1 : 1;

                mapping.set(params, workingHklList, originalValue + fd_sign * fd_step);
                const applied  = mapping.get(params, workingHklList);
                const actual_h = applied - originalValue;

                if (actual_h !== 0) {
                    calculateTotalPattern(y_calc_total, isStructural);

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

            // Restore the peak positions if a structural parameter shifted them.
            if (isStructural) {
                updateHklPositions(workingHklList, params, system);
            }

            columns[p] = { start: 0, values: jacobian_col.slice() };
        }

        const JtJ = new Array(n_params);
        for (let i = 0; i < n_params; i++) JtJ[i] = new Array(n_params).fill(0);
        const Jtr = new Array(n_params).fill(0);
        let ok = true;

        for (let i = 0; i < n_params; i++) {
            const ci = columns[i];
            const s_r = colDotFull(ci, residuals);
            if (!isFinite(s_r)) ok = false;
            Jtr[i] = s_r;

            for (let j = i; j < n_params; j++) {
                const s = colDot(ci, columns[j]);
                if (!isFinite(s)) ok = false;
                JtJ[i][j] = s;
                JtJ[j][i] = s;
            }
        }
        if (!ok) return null;

        // A parameter has no effect iff ITS OWN column is zero. Full stop.
        //
        // This used to be `JtJ[i][i] > JTJ_RANK_TOL * maxDiag`, i.e. each
        // diagonal was compared against the LARGEST diagonal in the matrix --
        // and that is meaningless, because the diagonals are in different
        // physical units. d(pattern)/d(a) is counts per angstrom on a cell
        // edge; d(pattern)/d(I) is dimensionless. Measured on a real cubic
        // problem the cell diagonal came out at 1.3e10 against 1e-2 for the
        // intensities: a ratio of 1.3e12, past the 1e-12 tolerance, so the
        // moment a lattice parameter was refined EVERY Pawley intensity was
        // reported as "has no effect on the fit and will not move".
        //
        // They moved perfectly well -- the `active` filter below selects on
        // JtJ[i][i] > 0, which is scale-free and was always right -- but the
        // user was told a hundred times over that the fit was doing nothing.
        // A warning that fires on a healthy refinement is worse than no
        // warning, because it teaches people to ignore the real ones.
        //
        // The test now matches the `active` filter exactly, so a parameter is
        // warned about precisely when it is dropped from the solve. Genuine
        // rank deficiency -- parameters that are individually alive but
        // collectively degenerate -- is a different condition and is reported
        // where it can be reported properly: the ESDs, and the overlap-cluster
        // section of the report.
        for (let i = 0; i < n_params; i++) {
            const d = JtJ[i][i];
            if (!(isFinite(d) && d > 0) && !deadParamWarned.has(paramMapping[i].name)) {
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
             refreshLeBailIntensities();
             const cost = evaluateResiduals(y_calc_baseline);
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

             const normals = buildNormalEquations();
             if (!normals) {
                  console.warn(`LM iter ${iter}: JtJ/Jtr contained non-finite entries.`);
                  finalJtJ = null;
                  lambda = Math.min(1e9, lambda * 5);
                  continue;
             }
             finalJtJ = normals.JtJ;
             const Jtr = normals.Jtr;



             const active = [];
             for (let i = 0; i < n_params; i++) {
                 // Only drop genuinely dead parameters. 
                 // Marquardt diagonal scaling handles the scale differences safely.
                 if (finalJtJ[i][i] > 0) active.push(i);
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

            oldParams = cloneParams(params);
             oldHklList = cloneHkl(workingHklList);

             const p_new = p_current.map((v, i) => v + p_step_clamped[i]);
             paramMapping.forEach((m, i) => m.set(params, workingHklList, p_new[i]));

             // FIX: set() silently clamps to [minVal, maxVal] (and returns
             // without assigning on a non-finite value), so the step that is
             // actually taken can be shorter than p_step_clamped -- zero, for a
             // parameter already sitting on a bound with the step pointing
             // outwards. `predicted` used to be computed from the PROPOSED
             // step while `actual` measured the APPLIED one, so every bound
             // contact depressed rho through pure bookkeeping.
             //
             // That never rejected a step -- acceptance below is decided by
             // cost alone -- but it did drive lambda. With rho stuck under
             // 0.75 the decrease branch is unreachable while the increase
             // branch fires every iteration, so lambda ratchets towards 1e9
             // and the whole refinement grinds to a halt, not just the
             // parameter on the bound. Pawley mode is the bad case: intensity
             // parameters are bounded below by zero and weak reflections park
             // there permanently, throttling the cell and profile terms too.
             //
             // Read back what the bounds allowed and measure the same step
             // the cost function sees. Same convention buildNormalEquations()
             // already uses for the finite-difference denominator. This has to
             // happen BEFORE calculateTotalPattern(), which runs
             // enforceSymmetryConstraintsWorker() and, in Le Bail mode,
             // re-extracts the intensities.
             const p_applied  = paramMapping.map(m => m.get(params, workingHklList));
             const p_step_real = p_applied.map((v, i) => v - p_current[i]);

             let predicted = 0;
             for (let i = 0; i < n_params; i++) {
                 predicted += 2 * p_step_real[i] * Jtr[i];
                 let q = 0;
                 for (let j = 0; j < n_params; j++) q += finalJtJ[i][j] * p_step_real[j];
                 predicted -= p_step_real[i] * q;
             }

             // Every component blocked: the model is sitting in a corner of
             // the bound box with the step pointing out of it. Nothing will
             // change, so stop rather than spin up lambda for the remaining
             // iterations.
             if (p_step_real.every(v => v === 0)) {
                 console.warn(`LM iter ${iter}: the proposed step is entirely blocked by parameter bounds; stopping.`);
                 params = oldParams;
                 workingHklList = oldHklList;
                 break;
             }

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

             const actual = cost - new_cost;
             const rho = (predicted > 0 && isFinite(predicted))
                         ? actual / predicted
                         : (actual > 0 ? 1 : -1);

             let stepAccepted = false;
             if (new_cost < cost && isFinite(new_cost)) {
                 stepAccepted = true;
                 // The old gate was `stepScale > 0.99`, which meant any trust-
                 // region clipping at all -- however slight -- also made the
                 // decrease branch unreachable, so lambda could only ever go
                 // up. Now that `predicted` describes the step that was really
                 // taken, a mildly clipped step is still evidence the model is
                 // good; only a heavily clipped one is not.
                 if (rho > 0.75 && stepScale > 0.5) lambda = Math.max(1e-9, lambda / 3);
                 else if (rho < 0.25)               lambda = Math.min(1e9, lambda * 2);
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
            const currentCost = evaluateResiduals(y_calc_baseline);
            if (currentCost === null || currentCost > bestTrueCost) {
                params = bestParams;
                workingHklList = bestHkl;
            }
        }
    } catch (err) {
        console.warn("LM: could not compare the final point with the best one.", err);
    }

    try {
        refreshLeBailIntensities();
        const finalCost = evaluateResiduals(y_calc_baseline);
        if (finalCost !== null) ss_res = finalCost;
        const finalNormals = buildNormalEquations();
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

    // One last constraint pass before the result leaves the worker. `params`
    // may have been rolled back to a saved best after the final evaluation, and
    // a rollback restores the MASTER edge without re-deriving the slaved ones
    // -- which is how a cubic fit came back reporting a = 8.13000 alongside
    // b = c = 8.13081. Harmless to the pattern, which reads only `a`, and
    // exactly the kind of inconsistency someone downstream would trust.
    enforceSymmetryConstraintsWorker(params, system);

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
// ---------------------------------------------------------------------------
//  MOVED. LEBAIL_PT_REFRESH_INTERVAL and the PT ladder constants
//  (PT_NUM_REPLICAS, PT_MAX_TEMP, PT_MIN_TEMP, PT_SWAP_INTERVAL,
//  PT_COST_SCALE_FLOOR, PT_RESCALE_THRESHOLD) now live in constants.js.
// ---------------------------------------------------------------------------


function refineParametersPT(initialParams, fitFlags, maxIter, hklList, system, refinementMode, leBailCycle = 0, totalLeBailCycles = 1) {
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList, refinementMode, system);
    if (paramMapping.length === 0) {
         const finalProgress = (leBailCycle + 1) / totalLeBailCycles;
         postMessage({ type: 'progress', value: Math.min(1.0, finalProgress) });
        return { params: initialParams, hklList: hklList, ss_res: 0 };
    }

    // The ladder itself lives in constants.js, with the reasoning for each
    // value (replica count vs neighbour swap acceptance, what T = 1 and
    // T = 1e-5 mean in units of the relative cost scale below).
    const numReplicas  = PT_NUM_REPLICAS;
    const maxTemp      = PT_MAX_TEMP;
    const minTemp      = PT_MIN_TEMP;
    const swapInterval = PT_SWAP_INTERVAL;


const n_points = workerWorkingData.tth.length;
    const scratch_bkg = new Float64Array(n_points);
    calculateTotalBackground(workerWorkingData.tth, initialParams, workerBackgroundAnchors, scratch_bkg);
    const scratch_pattern = new Float64Array(n_points);

    const objective = (p_obj, hkl_list_obj) => {
         try {
        
 enforceSymmetryConstraintsWorker(p_obj, system); 
            updateHklPositions(hkl_list_obj, p_obj, system);
            const netCalcPattern = calculatePatternCPU(workerWorkingData.tth, hkl_list_obj, p_obj, scratch_pattern);
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
    const refreshReplica = (rep) => {


if (refinementMode === 'le-bail') {
            enforceSymmetryConstraintsWorker(rep.params, system);
            updateHklPositions(rep.hklList, rep.params, system);
            leBailIntensityExtraction(
                { tth: workerWorkingData.tth, intensity: workerWorkingData.intensity,
                  background: scratch_bkg, weights: workerWorkingData.weights },
                rep.hklList, rep.params, scratch_pattern);
        }
       
       
       
        rep.cost = objective(rep.params, rep.hklList);
        return rep.cost;
    };

    // Templates, captured before anything mutates them.
    const paramTemplate = JSON.parse(JSON.stringify(initialParams));
    const hklTemplate   = JSON.parse(JSON.stringify(hklList));

    // Seed the intensities once from the starting parameters.
    const seed = { params: initialParams, hklList: hklList, cost: Infinity, temp: 0 };
    const initialCost = refreshReplica(seed);

    const temperatures = Array.from({ length: numReplicas }, (_, i) =>
        maxTemp * Math.pow(minTemp / maxTemp, i / (numReplicas - 1 || 1))
    );

    // ---------------------------------------------------------------------
    //  One reference energy scale shared by BOTH acceptance tests.
    //
    //  Why a scale at all. Previously the local Metropolis step divided by the
    //  replica's own CURRENT cost (relative scale) while the replica-exchange
    //  test used the raw cost difference (absolute scale). With chi-square
    //  values of 1e5-1e7 counts the swap exponent saturated, so exchanges were
    //  accepted either always or never and the run degenerated into 8
    //  independent annealers rather than parallel tempering. Dividing both
    //  tests by one shared scale makes the temperature ladder RELATIVE:
    //  exp(-(dChi2/scale)/T). At T = PT_MAX_TEMP = 1 a move that doubles
    //  chi-square is accepted with probability e^-1, which is the definition of
    //  a chain that explores freely; at T = PT_MIN_TEMP = 1e-5 a 0.01% increase
    //  is already rejected almost surely.
    //
    //  NOTE, because this looks like a bug and is not: costScale is of order
    //  chi-square itself (1e5-1e7). That is deliberate. Setting it to 1 and
    //  using the raw chi-square difference would put the exponent at
    //  -1e6/T for every uphill move, freezing every chain in the ladder
    //  solid -- including the hot end whose entire job is to move.
    //
    //  BUG FIX (stale scale). The scale used to be frozen at initialCost for
    //  the whole run. A successful global search drops chi-square by one to two
    //  orders of magnitude, and a frozen scale then means the ladder is
    //  silently 10-100x hotter in relative terms than it was designed to be:
    //  the cold replicas stop being cold, and PT loses the exploit end of its
    //  explore/exploit balance exactly when it is closing in on a solution.
    //
    //  It cannot simply be re-quoted continuously: a cost-dependent denominator
    //  that moves under the chain's feet on every step breaks detailed balance.
    //  So it is re-quoted ONLY at the periodic Le Bail refresh boundaries --
    //  points at which the objective function is ALREADY being redefined and
    //  every replica cost re-quoted anyway, so no additional damage is done --
    //  and only when the cost has moved by more than PT_RESCALE_THRESHOLD.
    //  Between two refreshes the acceptance criterion is exactly constant,
    //  which is the condition the Metropolis argument actually needs.
    // ---------------------------------------------------------------------
    let costScale = Math.max(PT_COST_SCALE_FLOOR, Math.abs(initialCost));

    /**
     * Re-quotes costScale from the current replica costs, if they have moved
     * far enough to matter. Call ONLY where every replica cost is being
     * re-quoted anyway.
     * @param {Array<{cost: number}>} reps
     * @returns {boolean} True if the scale was changed.
     */
    const maybeRescaleCost = (reps) => {
        let sum = 0, n = 0;
        for (const r of reps) {
            if (r && isFinite(r.cost)) { sum += Math.abs(r.cost); n++; }
        }
        if (n === 0) return false;
        const proposed = Math.max(PT_COST_SCALE_FLOOR, sum / n);
        // Only act on a change of regime, not on Monte Carlo jitter.
        const ratio = proposed / costScale;
        if (ratio > (1 - PT_RESCALE_THRESHOLD) && ratio < 1 / (1 - PT_RESCALE_THRESHOLD)) {
            return false;
        }
        costScale = proposed;
        return true;
    };

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
    const requoteBest = () => {
        paramMapping.forEach((m, i) => m.set(probe.params, probe.hklList, bestVector[i]));
        return refreshReplica(probe);
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
                 for (let i = 0; i < numReplicas; i++) refreshReplica(replicas[i]);
                 // The objective just changed, so the stored best cost is no
                 // longer comparable with the replicas'. Requote it.
                 bestOverallCost = requoteBest();
                 // This is the ONE point in the run where the acceptance
                 // criterion may change without breaking detailed balance --
                 // the objective has just been redefined regardless. See the
                 // long note on costScale.
                 maybeRescaleCost(replicas);
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
                 
                 if (proposed >= minV && proposed <= maxV) {
                     mapping.set(replica.params, replica.hklList, proposed);

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
                         for (let k = 0; k < paramMapping.length; k++) {
                             bestVector[k] = paramMapping[k].get(replica.params, replica.hklList);
                         }
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
     const finalCost = refreshReplica({ params: bestOverallParams, hklList: bestOverallHklList,
                                        cost: Infinity, temp: 0 });

     const finalCycleProgress = (leBailCycle + 1) / totalLeBailCycles;
     postMessage({ type: 'progress', value: Math.min(1.0, finalCycleProgress) });

      const parameterInfoForMainThread = paramMapping.map(m => ({
           name: m.name,
           scale: m.scale,
           isIntensity: m.isIntensity
      }));


    // Same reason as in refineParametersLM: the returned parameter set is a
    // saved snapshot, so re-derive the slaved cell edges from it.
    enforceSymmetryConstraintsWorker(bestOverallParams, system);

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
function polishWithLM(ptResults, fitFlags, hklList, system, refinementMode, maxIterations) {
    if (!ptResults || !ptResults.params || !ptResults.hklList) return ptResults;
    const polishIters = Math.max(10, Math.min(40, Math.round(maxIterations / 5)));
    try {
        postMessage({ type: 'progress', value: 0.98, message: 'Converging with Levenberg-Marquardt...' });
        const lm = refineParametersLM(
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
/**
 * Build the list of refinable parameters.
 *
 * @param {object} fitFlags
 * @param {ProfileParams} initialParams
 * @param {Reflection[]} hklList
 * @param {string} refinementMode  'pawley' | 'le-bail'
 * @param {string} system
 * @param {boolean} [quiet=false]  Suppress the postMessage warnings.
 *        calculateStatistics calls this purely to COUNT parameters, at the end
 *        of a fit that has already run. Without this flag the user was shown
 *        the whole "parameter X is fixed by symmetry / has no valid value"
 *        set a second time, after the refinement had finished, with nothing
 *        they could do about it -- the same warnings they had already been
 *        given when the refinement started.
 * @returns {{paramMapping: object[]}}
 */
function getParameterMapping(fitFlags, initialParams, hklList, refinementMode, system, quiet = false) {
    const mappings = [];
    const warn = (msg) => {
        if (quiet) return;
        console.warn(msg);
        postMessage({ type: 'warning', message: msg });
    };

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
    // The cell edges that enforceSymmetryConstraintsWorker derives from `a`.
    // They were previously refinable: ticking "b" on a cubic cell created a
    // parameter that moved, was reported, and had no effect whatsoever on the
    // calculated pattern, because updateHklPositions reads only `a`. Now it is
    // excluded with the same message the slaved Stephens terms get.
    switch (system) {
        case 'cubic':
            ['b', 'c'].forEach(n => slaved.add(n));
            break;
        case 'tetragonal': case 'hexagonal': case 'trigonal': case 'rhombohedral':
            slaved.add('b');
            break;
    }
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
            warn(`Parameter "${name}" is determined by ${system} symmetry and was excluded from the refinement.`);
            return null;
        }
        if (!initialParams || !(name in initialParams) || !isFinite(initialParams[name])) {
            // Refining a parameter the model never received (or received as
            // NaN) can only produce a dead column and a bogus ESD. Drop it and
            // tell the user rather than failing silently.
            warn(`Parameter "${name}" has no valid value and was excluded from the refinement.`);
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

          // -----------------------------------------------------------------
          //  A PAWLEY INTENSITY IS NOT BOUNDED BELOW BY ZERO.
          //
          //  It used to be (minVal: 0, plus a Math.max(0, .) in set()), and
          //  that quietly reintroduced exactly the bias the Le Bail extraction
          //  goes out of its way to avoid -- see the long note in
          //  leBailIntensityExtraction about clipping the negative half of the
          //  noise. An intensity is a measurement of a quantity whose TRUE
          //  value is non-negative, but the measurement itself carries noise,
          //  and forcing every realisation to be >= 0 turns a symmetric error
          //  into a one-sided one. It lands hardest on weak and systematically
          //  absent reflections, which are exactly the ones a space-group test
          //  needs to be honest about.
          //
          //  It was also an absorbing state in practice: a reflection clamped
          //  to zero disappeared from the calculated pattern entirely (the old
          //  `intensity <= 0` filter), so nothing downstream could tell a
          //  measured zero from a reflection that was never there.
          //
          //  The intensities are therefore free. calculatePattern renders a
          //  negative one as a negative contribution, and the reported
          //  I / sigma is left unclipped so the reader can see it.
          // -----------------------------------------------------------------
          hklList.forEach((hkl, index) => {
               if (hkl && Number.isFinite(hkl.tth) && hkl.tth >= tthMin && hkl.tth <= tthMax) {
                    const hkl_name = `I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`;
                    const initialIntensity = (Number.isFinite(hkl.intensity) && hkl.intensity > 1e-6) ? hkl.intensity : 1000.0;
                    const scale = initialIntensity;  // typical magnitude, used by PT

                    mappings.push({
                        name: hkl_name,
                        scale: scale,
                        defaultScale: 1000.0,
                        step: 0.3,
                        maxStepAbs: Infinity,
                        isIntensity: true,
                        index: index,
                        minVal: -Infinity,
                        maxVal: Infinity,
                        get: (p_obj, hkl_list_obj) => {
                            return hkl_list_obj?.[index]?.intensity ?? 0;
                        },
                        set: (p_obj, hkl_list_obj, rawValue) => {
                             if (hkl_list_obj?.[index] && Number.isFinite(rawValue)) {
                                hkl_list_obj[index].intensity = rawValue;
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

    mappings.push(createMapping(fitFlags.zeroShift, 'zeroShift', 0.5, -Infinity, Infinity, 0.1, 1.0));
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
            currentHklList.forEach(peak => { if (peak) peak.intensity = LEBAIL_FLAT_SEED; });
            try {
                const bkgSeed = calculateTotalBackground(workerWorkingData.tth, currentParams, workerBackgroundAnchors);
                enforceSymmetryConstraintsWorker(currentParams, system);
                updateHklPositions(currentHklList, currentParams, system);
                for (let seedCycle = 0; seedCycle < 5; seedCycle++) {
                    leBailIntensityExtraction(
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
                finalResults = refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { 
                finalResults = refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
                finalResults = polishWithLM(finalResults, fitFlags, currentHklList, system, refinementMode, maxIterations);
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
                 leBailIntensityExtraction(expDataForInit, currentHklList, currentParams);
             } catch (initError) {
                  console.warn("Could not initialize Pawley intensities using Le Bail extraction:", initError);
                  currentHklList.forEach(peak => { if (peak) peak.intensity = 1000.0; });
             }
            
            if (algorithm === 'lm') {
                finalResults = refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
            } else { 
                finalResults = refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList, system, refinementMode, 0, 1);
                finalResults = polishWithLM(finalResults, fitFlags, currentHklList, system, refinementMode, maxIterations);
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

        const finalNetPatternForStats = calculatePatternCPU(workerWorkingData.tth, finalHklList, finalParams);
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