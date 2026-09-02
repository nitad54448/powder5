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
importScripts('constants.js', 'crystal.js', 'profile.js', 'pawley_linear.js');


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

    // Removed the solveLU fallback. If Cholesky fails, the matrix is ill-conditioned 
    // and LM must cleanly increase lambda. LU returns non-descent directions and takes O(N^3).
    return solveCholesky(A, b, n);
}

/**
 * In-place-free Cholesky: A = L L^T, then two triangular solves.
 * @returns {number[]|null} null if A is not positive definite.
 */
function solveCholesky(A, b, n) {
    // Packed lower triangle, row-major: L[i*(i+1)/2 + j] for j <= i.
    const L = new Float64Array((n * (n + 1)) / 2);
    const firstNonZero = new Int32Array(n);

    // Compute the skyline envelope (first non-zero index of each row)
    for (let i = 0; i < n; i++) {
        let fnz = i;
        for (let j = 0; j <= i; j++) {
            if (A[i][j] !== 0) {
                fnz = j;
                break;
            }
        }
        firstNonZero[i] = fnz;
    }

    for (let i = 0; i < n; i++) {
        const rowA = A[i];
        const bi = (i * (i + 1)) / 2;
        const fi = firstNonZero[i];
        
        for (let j = fi; j <= i; j++) {
            const bj = (j * (j + 1)) / 2;
            const fj = firstNonZero[j];
            let s = rowA[j];
            if (!isFinite(s)) return null;
            
            // Only loop over the overlapping envelope
            const startK = Math.max(fi, fj);
            for (let k = startK; k < j; k++) {
                s -= L[bi + k] * L[bj + k];
            }
            
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
        for (let k = firstNonZero[i]; k < i; k++) {
            s -= L[bi + k] * y[k];
        }
        y[i] = s / L[bi + i];
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let s = y[i];
        for (let k = i + 1; k < n; k++) {
            if (firstNonZero[k] <= i) {
                s -= L[(k * (k + 1)) / 2 + i] * x[k];
            }
        }
        x[i] = s / L[(i * (i + 1)) / 2 + i];
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

    // Only ratio21 is still needed here: the Ka2 position, its window and its
    // prepared profile now come from buildPeakContributions below, which is
    // the single place that resolves them.
    const ratio21 = params.ratio || 0;

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

    // ONE pass over the contributions instead of a second prepareVoigt sweep.
    // The comment above already claimed this ("built once and reused for the
    // per-peak walk below, so prepareVoigt runs once per reflection per
    // extraction pass instead of twice") -- the code did not do it, and rebuilt
    // prep1/prep2 from scratch here. Measured at 2500 reflections that was
    // 8.2 ms of duplicated work per pass, and the pass runs up to
    // LEBAIL_MAX_PASSES = 24 times per LM iteration.
    const buckets = bucketContributionsByPeak(contributions, hklList.length);

    hklList.forEach((peak, j) => {
        // Number.isFinite, not truthiness: tth === 0 is falsy but is a position.
        if (!hasRenderablePosition(peak)) return;
        const bucket = buckets[j];
        if (!bucket || bucket.length === 0) return;

        // The Ka1 contribution carries weight 1; anything else is the Ka2
        // satellite. Read the weight rather than the emission order, so a
        // future change in buildPeakContributions cannot silently swap them.
        let prep1 = null, prep2 = null;
        let startIdx = Infinity, stopIdx = -Infinity;
        for (let c = 0; c < bucket.length; c++) {
            const con = bucket[c];
            if (con.weight === 1.0 && !prep1) prep1 = con.prep; else prep2 = con.prep;
            if (con.start < startIdx) startIdx = con.start;
            if (con.stop  > stopIdx)  stopIdx  = con.stop;
        }
        if (!prep1) prep1 = bucket[0].prep;
        if (!(stopIdx > startIdx)) return;

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
        // below the analytic value.
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

        for (let i = startIdx; i < stopIdx; i++) {
            const current_tth = expData.tth[i];

            let profile_val = 0;
            for (let c = 0; c < bucket.length; c++) {
                const con = bucket[c];
                if (i >= con.start && i < con.stop) {
                    profile_val += con.weight * evalVoigt(current_tth, con.prep);
                }
            }

            let current_fraction = 0;
            if (profile_val > 0 && total_calc[i] > 1e-9) {
                current_fraction = (peak.intensity * profile_val) / total_calc[i];
                if (current_fraction > 1.0) current_fraction = 1.0; // Prevent float drift explosion
            }

            const current_partitioned_I = y_obs_net[i] * current_fraction;

            let half = 0;
            if (i > startIdx) {
                const step_width = expData.tth[i] - expData.tth[i - 1];
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
    const isBgRefined = fitFlags && fitFlags.fitBackground;

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
        backgroundRefined: isBgRefined,
        nBackgroundAnchors,
        backgroundNote: nBackgroundAnchors > 0
            ? `Background: spline through ${nBackgroundAnchors} user anchor point`
              + `${nBackgroundAnchors === 1 ? '' : 's'}, `
              + (isBgRefined ? `Y positions refined (2-theta positions fixed). They are `
              + `counted in P, and reported ESDs include background uncertainty.`
              : `not refined. Its ${nBackgroundAnchors} `
              + `degree${nBackgroundAnchors === 1 ? '' : 's'} of freedom are not counted in P, and no `
              + `reported ESD includes background uncertainty.`)
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

function refineParametersLM(initialParams, fitFlags, maxIter, hklList, system, refinementMode,
                            leBailCycle = 0, totalLeBailCycles = 1, progressWindow = null) {
    const cloneParams = (p) => ({ ...p });
    const cloneHkl = (list) => list.map(h => (h ? { ...h } : null));

    const separable = (refinementMode === 'pawley');

    // The LM search space. In Pawley mode the intensities are excluded here
    // and solved exactly instead; the FULL mapping (including them) is rebuilt
    // once at the end, because the report's ESDs and the charge-flipping
    // hand-off both read the intensity block of J^T J.
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList,
                                                 refinementMode, system, false,
                                                 { excludeIntensities: separable });

    const baseProgress = progressWindow ? progressWindow.base : (leBailCycle / totalLeBailCycles);
    const cycleProgressSpan = progressWindow ? progressWindow.span : (1 / totalLeBailCycles);
    const postProgress = (frac) => {
        const v = Math.min(1.0, baseProgress + Math.max(0, Math.min(1, frac)) * cycleProgressSpan);
        return v;
    };

    if (paramMapping.length === 0 && !separable) {
        postMessage({ type: 'progress', value: postProgress(1) });
        return { params: initialParams, hklList: hklList, ss_res: 0 };
    }

    let params = initialParams;
    let workingHklList = JSON.parse(JSON.stringify(hklList));

    let ss_res = Infinity;
    let lambda = 0.001;

    const tthAxis = workerWorkingData.tth;
    const y_obs = workerWorkingData.intensity;
    const n_points = y_obs.length;
    const n_params = paramMapping.length;
    const n_peaks = workingHklList.length;

    const sqrt_weights = new Float64Array(n_points);
    for (let i = 0; i < n_points; i++) sqrt_weights[i] = Math.sqrt(workerWorkingData.weights[i]);

    const residuals      = new Float64Array(n_points);
    const y_calc_total   = new Float64Array(n_points);
    const y_calc_baseline= new Float64Array(n_points);
    const jacobian_col   = new Float64Array(n_points);
    const scratch_bkg    = new Float64Array(n_points);
    const scratch_pattern= new Float64Array(n_points);
    const scratch_net    = new Float64Array(n_points);   // y_obs - background
    const scratch_cost   = new Float64Array(n_points);

    // ALLOCATION HOISTING. These were `new Array(n)` / `new Float64Array(n*n)`
    // inside the iteration. At the old n_params they cost ~300 MB of churn per
    // iteration; at the new one they are trivial, but there is no reason to
    // reallocate them at all.
    const JtJ = new Array(n_params);
    for (let i = 0; i < n_params; i++) JtJ[i] = new Float64Array(n_params);
    const Jtr = new Float64Array(n_params);
    const A_lm = new Array(n_params);
    for (let i = 0; i < n_params; i++) A_lm[i] = new Float64Array(n_params);
    const Dscale = new Float64Array(n_params);
    const rhs_scaled = new Float64Array(n_params);
    const p_step = new Float64Array(n_params);

    calculateTotalBackground(tthAxis, initialParams, workerBackgroundAnchors, scratch_bkg);

    // -----------------------------------------------------------------------
    //  Cached state for the separable solve.
    //
    //  `lastBuckets` / `lastWindows` describe the profile at `params`;
    //  `lastNetPattern` is the BACKGROUND-FREE calculated pattern that goes
    //  with them. They are invalidated together, by touchProfile(), and
    //  nothing else may write them.
    // -----------------------------------------------------------------------
    const linCache = {};
    let lastBuckets = null, lastWindows = null, lastCols = null;
    let lastDegenerate = [], lastDegenerateGroups = [];
    const lastNetPattern = new Float64Array(n_points);
    let lastNetValid = false;

    /** Mark the cached profile stale. Call after ANY change to params or positions. */
    const touchProfile = () => { lastBuckets = null; lastNetValid = false; };

    const refreshBackground = () => {
        calculateTotalBackground(tthAxis, params, workerBackgroundAnchors, scratch_bkg);
        for (let i = 0; i < n_points; i++) scratch_net[i] = y_obs[i] - scratch_bkg[i];
    };

    /**
     * Solve the Pawley intensities exactly at the current parameters.
     *
     * FIX 1 -- crystal.js signals "this reflection has no position" by writing
     * peak.tth = null, NOT NaN, and it does so for three distinct situations:
     * indices missing, the cell degenerate, and sin(theta) out of range. All
     * three reach hasRenderablePosition() as false, so buildPeakContributions
     * emits nothing and the reflection arrives here with an empty window.
     *
     * The earlier version then wrote sol.I[j] = 0 back into the list for every
     * such reflection. That is wrong in a way that only shows up in the
     * report: a reflection outside the fitted 2-theta range is not measured to
     * be zero, it is not measured at all, and printing 0.0 for it -- with an
     * ESD -- states something the data does not support. The old LM gave those
     * reflections a zero-length Jacobian slice and left their value alone;
     * this now does the same.
     *
     * FIX 2 -- the pattern is accumulated from the columns the solve already
     * built, instead of calling calculatePatternCPU again. That second call
     * re-ran buildPeakContributions, i.e. a full prepareVoigt sweep over every
     * reflection, for a pattern we had all the pieces of: 12 ms per iteration
     * at 2500 reflections, for nothing.
     *
     * @param {boolean} [rebuildProfile=true]
     * @returns {number|null} Weighted cost, or null if the solve failed.
     */
    const solveIntensitiesAndCost = (rebuildProfile = true) => {
        if (rebuildProfile || !lastBuckets) {
            enforceSymmetryConstraintsWorker(params, system);
            updateHklPositions(workingHklList, params, system);
            const contributions = buildPeakContributions(tthAxis, workingHklList, params);
            lastBuckets = bucketContributionsByPeak(contributions, n_peaks);
            lastWindows = peakWindows(lastBuckets, n_peaks);
        }
        const sol = solvePawleyIntensities(tthAxis, lastBuckets, lastWindows,
                                           sqrt_weights, scratch_net, linCache);
        if (!sol) { lastNetValid = false; return null; }
        lastCols = sol.cols;
        lastDegenerate = sol.degenerate;
        lastDegenerateGroups = sol.degenerateGroups || [];

        // Determined reflections take the solved value. Undetermined ones keep
        // whatever they had -- see FIX 1.
        const undet = sol.undetermined;
        let u = 0;
        for (let j = 0; j < n_peaks; j++) {
            if (u < undet.length && undet[u] === j) { u++; continue; }
            if (workingHklList[j]) workingHklList[j].intensity = sol.I[j];
        }

        // Net pattern from the columns, un-weighted, cached for the background
        // Jacobian and for the baseline render.
        lastNetPattern.fill(0);
        for (let j = 0; j < n_peaks; j++) {
            const cj = sol.cols[j];
            if (!cj) continue;
            const amp = sol.I[j];
            if (amp === 0) continue;
            const s = lastWindows.start[j];
            for (let k = 0; k < cj.length; k++) {
                // The column carries sqrt(w); divide it back out.
                const sw = sqrt_weights[s + k];
                if (sw > 0) lastNetPattern[s + k] += amp * cj[k] / sw;
            }
        }
        lastNetValid = true;

        let cost = 0;
        for (let i = 0; i < n_points; i++) {
            const r = (scratch_net[i] - lastNetPattern[i]) * sqrt_weights[i];
            if (!isFinite(r)) return null;
            cost += r * r;
        }
        return cost;
    };

    const refreshLeBailIntensities = () => {
        if (refinementMode !== 'le-bail') return;
        enforceSymmetryConstraintsWorker(params, system);
        updateHklPositions(workingHklList, params, system);
        if (fitFlags && fitFlags.fitBackground) refreshBackground();

        if (LEBAIL_RESET_EACH_ITER) {
            for (let i = 0; i < n_peaks; i++) {
                if (workingHklList[i]) workingHklList[i].intensity = LEBAIL_FLAT_SEED;
            }
        }
        const expData = { tth: tthAxis, intensity: y_obs, background: scratch_bkg,
                          weights: workerWorkingData.weights };

        // Run the decomposition TO CONVERGENCE rather than for a fixed six
        // passes. Six is wasteful on an easy pattern and short on a heavily
        // overlapped one -- and the overlapped case is the one where being
        // short matters, because the LM then fits its profile parameters
        // against intensities that are still moving. Starting from flat keeps
        // the result a deterministic function of the parameters, which is what
        // LEBAIL_RESET_EACH_ITER exists to guarantee.
        let prev = null;
        for (let pass = 0; pass < LEBAIL_MAX_PASSES; pass++) {
            leBailIntensityExtraction(expData, workingHklList, params, scratch_pattern);
            const cur = new Float64Array(n_peaks);
            for (let i = 0; i < n_peaks; i++) {
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

    /**
     * The full calculated pattern into `targetArray`.
     *
     * FIX 3 -- `scratch_pattern` is NOT a safe cache for the peak part.
     * refreshLeBailIntensities passes it straight into leBailIntensityExtraction
     * as its scratch_calc, which overwrites it with the partition's total_calc.
     * The earlier version's `requiresPeakUpdate = false` path therefore read a
     * buffer that, in Le Bail mode, held the pattern at the intensities of the
     * PREVIOUS decomposition pass. It happened to be correct only because of
     * the order the calls came in, which is not a property to rely on. The peak
     * pattern now has its own buffer, written only here.
     */
    const peakCache = new Float64Array(n_points);
    let peakCacheValid = false;

    const calculateTotalPattern = (targetArray, requiresHklUpdate = true,
                                   requiresPeakUpdate = true, requiresBkgUpdate = true) => {
        enforceSymmetryConstraintsWorker(params, system);
        if (requiresHklUpdate) updateHklPositions(workingHklList, params, system);
        if (requiresBkgUpdate) {
            calculateTotalBackground(tthAxis, params, workerBackgroundAnchors, scratch_bkg);
        }
        if (requiresPeakUpdate || !peakCacheValid) {
            calculatePatternCPU(tthAxis, workingHklList, params, peakCache);
            peakCacheValid = true;
        }
        for (let i = 0; i < n_points; i++) targetArray[i] = peakCache[i] + scratch_bkg[i];
    };

    const evaluateResiduals = (target) => {
        calculateTotalPattern(target, true, true, true);
        let cost = 0;
        for (let i = 0; i < n_points; i++) {
            residuals[i] = (y_obs[i] - target[i]) * sqrt_weights[i];
            if (!isFinite(residuals[i])) return null;
            cost += residuals[i] * residuals[i];
        }
        return isFinite(cost) ? cost : null;
    };

    /** One iteration's baseline: intensities settled, residuals and cost current. */
    const settleAndEvaluate = () => {
        if (separable) {
            refreshBackground();
            const c = solveIntensitiesAndCost(true);
            if (c === null) return null;
            // FIX 2 again: the baseline is the cached net pattern plus the
            // background. No third pattern evaluation.
            peakCache.set(lastNetPattern);
            peakCacheValid = true;
            for (let i = 0; i < n_points; i++) {
                y_calc_baseline[i] = lastNetPattern[i] + scratch_bkg[i];
                residuals[i] = (y_obs[i] - y_calc_baseline[i]) * sqrt_weights[i];
                if (!isFinite(residuals[i])) return null;
            }
            return c;
        }
        refreshLeBailIntensities();
        peakCacheValid = false;      // the extraction moved every intensity
        return evaluateResiduals(y_calc_baseline);
    };


    const deadParamWarned = new Set();

    const colDot = (a, b) => {
        const lo = Math.max(a.start, b.start);
        const hi = Math.min(a.start + a.values.length, b.start + b.values.length);
        if (lo >= hi) return 0;
        let acc = 0;
        for (let k = lo; k < hi; k++) acc += a.values[k - a.start] * b.values[k - b.start];
        return acc;
    };
    const colDotFull = (a, v) => {
        let acc = 0;
        const end = a.start + a.values.length;
        for (let k = a.start; k < end; k++) acc += a.values[k - a.start] * v[k];
        return acc;
    };

    /**
     * Analytic derivative of the pattern with respect to one background anchor
     * height: the spline basis function, obtained by differencing the SPLINE
     * only. No peak is touched.
     *
     * The monotonic Hermite spline is piecewise-cubic in the anchor heights
     * and not strictly linear in them (the Fritsch-Carlson tangent limiter is
     * a min), so this is still a difference -- but of a 10 us spline
     * evaluation rather than a 24 ms pattern render.
     */
    const backgroundColumn = (mapping, out) => {
        const original = mapping.get(params, workingHklList);
        const h = Math.max(1e-6, Math.abs(original) * FD_BASE_FRACTION, 1e-3);
        mapping.set(params, workingHklList, original + h);
        const applied = mapping.get(params, workingHklList);
        const actual_h = applied - original;
        if (actual_h === 0) { mapping.set(params, workingHklList, original); out.fill(0); return false; }
        const probe = new Float64Array(n_points);
        calculateTotalBackground(tthAxis, params, workerBackgroundAnchors, probe);
        mapping.set(params, workingHklList, original);
        calculateTotalBackground(tthAxis, params, workerBackgroundAnchors, scratch_bkg);
        let responded = false;
        for (let i = 0; i < n_points; i++) {
            const d = (probe[i] - scratch_bkg[i]) / actual_h;
            out[i] = d * sqrt_weights[i];
            if (d !== 0) responded = true;
        }
        return responded;
    };

    // -----------------------------------------------------------------------
    //  Normal equations over the NON-LINEAR parameters only
    // -----------------------------------------------------------------------
    const buildNormalEquations = () => {
        const columns = new Array(n_params).fill(null);
        const savedIntensities = new Float64Array(n_peaks);
        for (let i = 0; i < n_peaks; i++) savedIntensities[i] = workingHklList[i] ? workingHklList[i].intensity : 0;

        for (let p = 0; p < n_params; p++) {
            const mapping = paramMapping[p];

            // ---- background: analytic, no peak render --------------------
            if (mapping.isBackground) {
                jacobian_col.fill(0);
                const responded = backgroundColumn(mapping, jacobian_col);
                if (!responded && !deadParamWarned.has(mapping.name)) {
                    deadParamWarned.add(mapping.name);
                    postMessage({ type: 'warning', message: `Parameter "${mapping.name}" has no effect on the fit and will not move.` });
                }
                columns[p] = { start: 0, values: jacobian_col.slice() };
                continue;
            }

            // ---- finite difference on a non-linear parameter -------------
            const originalValue = mapping.get(params, workingHklList);
            const isStructural = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].includes(mapping.name);

            const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
            const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;
            const magnitude = Math.max(Math.abs(originalValue), mapping.defaultScale || 0, 1e-9);
            const baseStep = Math.max(1e-9, magnitude * FD_BASE_FRACTION);
            const maxProbe = Math.max(baseStep, FD_MAX_PROBE_FRACTION * magnitude);

            let fd_step = baseStep, responded = false;
            jacobian_col.fill(0);

            for (let attempt = 0; attempt < FD_MAX_ATTEMPTS; attempt++) {
                const fd_sign = (originalValue + fd_step > maxV && originalValue - fd_step >= minV) ? -1 : 1;
                mapping.set(params, workingHklList, originalValue + fd_sign * fd_step);
                const actual_h = mapping.get(params, workingHklList) - originalValue;

                if (actual_h !== 0) {
                    if (separable && PAWLEY_VARPRO_EXACT_JACOBIAN) {
                        // Exact variable-projection column: the intensities
                        // re-optimise with the probe, so this measures the
                        // derivative of the REDUCED cost surface LM is on.
                        solveIntensitiesAndCost(true);
                        calculateTotalPattern(y_calc_total, false, true, false);
                    } else {
                        calculateTotalPattern(y_calc_total, isStructural, true, false);
                    }
                    for (let i = 0; i < n_points; i++) {
                        jacobian_col[i] = (!isFinite(y_calc_baseline[i]) || !isFinite(y_calc_total[i]))
                            ? 0
                            : (y_calc_total[i] - y_calc_baseline[i]) / actual_h * sqrt_weights[i];
                    }
                    for (let i = 0; i < n_points; i++) if (jacobian_col[i] !== 0) { responded = true; break; }
                }
                if (responded || fd_step >= maxProbe) break;
                fd_step = Math.min(maxProbe, fd_step * 30);
            }
            if (!responded) jacobian_col.fill(0);

            mapping.set(params, workingHklList, originalValue);
            for (let i = 0; i < n_peaks; i++) if (workingHklList[i]) workingHklList[i].intensity = savedIntensities[i];
            // FIX 4: re-derive the slaved cell edges BEFORE recomputing
            // positions. A probe on `a` leaves b and c holding a + h once
            // mapping.set() has put `a` back, because only enforce() mirrors
            // them. updateHklPositions reads neither on a cubic cell, so the
            // pattern is unaffected -- but `params` is the object that gets
            // cloned into bestParams and returned, and the worker's own note
            // at line 1412 records this exact inconsistency escaping once
            // already (a = 8.13000 reported next to b = c = 8.13081).
            if (isStructural || (separable && PAWLEY_VARPRO_EXACT_JACOBIAN)) {
                enforceSymmetryConstraintsWorker(params, system);
                updateHklPositions(workingHklList, params, system);
                touchProfile();
            }
            columns[p] = { start: 0, values: jacobian_col.slice() };
        }

        // Restore the profile cache the probes invalidated.
        if (separable && !lastBuckets) solveIntensitiesAndCost(true);

        let ok = true;
        for (let i = 0; i < n_params; i++) {
            const ci = columns[i];
            const s_r = colDotFull(ci, residuals);
            if (!isFinite(s_r)) ok = false;
            Jtr[i] = s_r;
            for (let j = i; j < n_params; j++) {
                const s = colDot(ci, columns[j]);
                if (!isFinite(s)) ok = false;
                JtJ[i][j] = s; JtJ[j][i] = s;
            }
        }
        if (!ok) return false;

        for (let i = 0; i < n_params; i++) {
            const d = JtJ[i][i];
            if (!(isFinite(d) && d > 0) && !deadParamWarned.has(paramMapping[i].name)) {
                deadParamWarned.add(paramMapping[i].name);
                postMessage({ type: 'warning', message: `Parameter "${paramMapping[i].name}" has no effect on the fit and will not move.` });
            }
        }
        return true;
    };

    // -----------------------------------------------------------------------
    //  Iteration
    // -----------------------------------------------------------------------
    let lastAcceptedCost = Infinity;
    let convergedRepeats = 0, lastProgressPct = -1;
    let bestTrueCost = Infinity, bestParams = null, bestHkl = null;
    let trustFactor = LM_TRUST_INITIAL, outerRejects = 0;
    let haveNormals = false;

    for (let iter = 0; iter < maxIter; iter++) {
        let oldParams = null, oldHklList = null;
        try {
            const cost = settleAndEvaluate();
            if (cost === null) {
                lambda = Math.min(1e9, lambda * 10);
                continue;
            }

            if (cost < bestTrueCost) {
                bestTrueCost = cost;
                bestParams = cloneParams(params);
                bestHkl = cloneHkl(workingHklList);
                trustFactor = Math.min(1.0, trustFactor * LM_TRUST_GROWTH);
            } else if (isFinite(bestTrueCost) && bestParams && cost > bestTrueCost * (1 + LM_OUTER_TOL)) {
                if (++outerRejects > LM_MAX_OUTER_REJECTS) break;
                params = cloneParams(bestParams);
                workingHklList = cloneHkl(bestHkl);
                if (separable) touchProfile();
                lambda = Math.min(1e9, lambda * 10);
                trustFactor = Math.max(1e-3, trustFactor * 0.25);
                lastAcceptedCost = bestTrueCost;
                convergedRepeats = 0;
                continue;
            }

            if (!isFinite(lastAcceptedCost)) lastAcceptedCost = cost;
            ss_res = cost;

            if (n_params === 0) break;   // separable with nothing non-linear left

            if (!buildNormalEquations()) {
                haveNormals = false;
                lambda = Math.min(1e9, lambda * 5);
                continue;
            }
            haveNormals = true;

            const active = [];
            for (let i = 0; i < n_params; i++) if (JtJ[i][i] > 0) active.push(i);
            if (active.length === 0) break;

            const nA = active.length;
            for (let ii = 0; ii < nA; ii++) Dscale[ii] = Math.sqrt(JtJ[active[ii]][active[ii]]);
            for (let ii = 0; ii < nA; ii++) {
                const row = A_lm[ii], di = Dscale[ii], ri = JtJ[active[ii]];
                for (let jj = 0; jj < nA; jj++) row[jj] = ri[active[jj]] / (di * Dscale[jj]);
                row[ii] += lambda;
                rhs_scaled[ii] = Jtr[active[ii]] / di;
            }

            const step_active = solveLinearSystem(A_lm.slice(0, nA), rhs_scaled.subarray(0, nA));
            if (!step_active) { lambda = Math.min(1e9, lambda * 10); continue; }

            p_step.fill(0);
            for (let ii = 0; ii < nA; ii++) p_step[active[ii]] = step_active[ii] / Dscale[ii];
            let bad = false;
            for (let i = 0; i < n_params; i++) if (!isFinite(p_step[i])) bad = true;
            if (bad) { lambda = Math.min(1e9, lambda * 5); continue; }

            const p_current = paramMapping.map(m => m.get(params, workingHklList));

            let stepScale = 1.0;
            for (let idx = 0; idx < n_params; idx++) {
                const m = paramMapping[idx];
                const cap = trustFactor * Math.min(
                    (m.maxStepAbs !== undefined) ? m.maxStepAbs : Infinity,
                    Math.max(2 * Math.abs(p_current[idx]), 4 * (m.defaultScale || 1e-9)));
                const want = Math.abs(p_step[idx]);
                if (want > cap && want > 0) stepScale = Math.min(stepScale, cap / want);
            }

            oldParams = cloneParams(params);
            oldHklList = cloneHkl(workingHklList);

            for (let i = 0; i < n_params; i++) {
                paramMapping[i].set(params, workingHklList, p_current[i] + p_step[i] * stepScale);
            }
            const p_applied = paramMapping.map(m => m.get(params, workingHklList));
            const p_step_real = p_applied.map((v, i) => v - p_current[i]);

            let predicted = 0;
            for (let i = 0; i < n_params; i++) {
                const si = p_step_real[i];
                if (si === 0) continue;
                predicted += 2 * si * Jtr[i];
                const ri = JtJ[i];
                let q = 0;
                for (let j = 0; j < n_params; j++) q += ri[j] * p_step_real[j];
                predicted -= si * q;
            }

            if (p_step_real.every(v => v === 0)) {
                params = oldParams; workingHklList = oldHklList;
                if (separable) touchProfile();
                break;
            }

            // TRIAL COST. In separable mode the intensities are re-solved at
            // the trial point, so this is the reduced objective -- the same
            // function `cost` above measured. That consistency is what makes
            // the accept/reject test and rho mean anything.
            if (separable) touchProfile();
            let new_cost;
            if (separable) {
                refreshBackground();
                const c = solveIntensitiesAndCost(true);
                new_cost = (c === null) ? Infinity : c;
            } else {
                calculateTotalPattern(y_calc_total, true, true, true);
                new_cost = 0;
                for (let i = 0; i < n_points; i++) {
                    const res = (y_obs[i] - y_calc_total[i]) * sqrt_weights[i];
                    if (!isFinite(res)) { new_cost = Infinity; break; }
                    new_cost += res * res;
                }
            }

            const actual = cost - new_cost;
            const rho = (predicted > 0 && isFinite(predicted)) ? actual / predicted : (actual > 0 ? 1 : -1);

            let stepAccepted = false;
            if (new_cost < cost && isFinite(new_cost)) {
                stepAccepted = true;
                if (rho > 0.75 && stepScale > 0.5) lambda = Math.max(1e-9, lambda / 3);
                else if (rho < 0.25) lambda = Math.min(1e9, lambda * 2);
            } else {
                params = oldParams; workingHklList = oldHklList;
                if (separable) touchProfile();
                lambda = Math.min(1e9, lambda * 3);
            }

            const pct = Math.floor(postProgress((iter + 1) / maxIter) * 100);
            if (pct !== lastProgressPct) {
                lastProgressPct = pct;
                postMessage({ type: 'progress', value: postProgress((iter + 1) / maxIter) });
            }

            if (stepAccepted) {
                const denom = Math.max(Math.abs(lastAcceptedCost), 1.0);
                const relDrop = Math.abs(lastAcceptedCost - new_cost) / denom;
                lastAcceptedCost = new_cost;
                ss_res = new_cost;
                if (relDrop < LM_COST_TOL) { if (++convergedRepeats >= LM_CONVERGED_REPEATS) break; }
                else convergedRepeats = 0;
            }
        } catch (error) {
            console.error("Error during LM iteration:", iter, error);
            postMessage({ type: 'error', message: `Error in LM iter ${iter}: ${error.message}` });
            return { params: oldParams || params || initialParams,
                     hklList: oldHklList || workingHklList || JSON.parse(JSON.stringify(hklList)),
                     ss_res, error: true, JtJ: null,
                     parameterInfo: paramMapping.map(m => ({ name: m.name, scale: 1.0,
                         typicalMagnitude: m.scale, isIntensity: m.isIntensity })),
                     algorithm: 'lm', fitFlags };
        }
    }

    // Land on the best point seen.
    try {
        if (bestParams && isFinite(bestTrueCost)) {
            if (separable) touchProfile();
            const currentCost = settleAndEvaluate();
            if (currentCost === null || currentCost > bestTrueCost) {
                params = bestParams; workingHklList = bestHkl;
                if (separable) touchProfile();
            }
        }
        const finalCost = settleAndEvaluate();
        if (finalCost !== null) ss_res = finalCost;
    } catch (err) {
        console.warn("LM: could not settle at the final point.", err);
    }

    // -----------------------------------------------------------------------
    //  FULL normal equations, ONCE, at the converged point.
    //
    //  The report's per-reflection ESDs, intensityOverlapClusters and the
    //  charge-flipping hand-off all read the intensity block of J^T J, so it
    //  still has to exist -- it just does not have to be rebuilt on every one
    //  of forty iterations to get there.
    // -----------------------------------------------------------------------
    const fullMapping = getParameterMapping(fitFlags, params, workingHklList,
                                            refinementMode, system, true).paramMapping;
    const finalJtJ = buildFullNormalEquations(fullMapping, params, workingHklList, system,
                                              refinementMode, sqrt_weights, residuals,
                                              tthAxis, y_obs, scratch_bkg, scratch_pattern,
                                              n_points, fitFlags);

    postMessage({ type: 'progress', value: postProgress(1) });
    enforceSymmetryConstraintsWorker(params, system);

    return {
        params,
        hklList: workingHklList,
        JtJ: finalJtJ,
        parameterInfo: fullMapping.map(m => ({ name: m.name, scale: 1.0,
                                               typicalMagnitude: m.scale, isIntensity: m.isIntensity })),
        ss_res,
        degenerateReflections: lastDegenerate,
        degenerateGroups: lastDegenerateGroups,
        algorithm: 'lm',
        fitFlags
    };
}

/**
 * J^T J over the FULL parameter set, built once at a converged point.
 *
 * Intensity columns are analytic and local; everything else is a finite
 * difference. Identical arithmetic to the old buildNormalEquations, minus the
 * quadratic contribution scan (bucketed) and the dense boxed allocation.
 *
 * @returns {Array<Float64Array>|null}
 */
function buildFullNormalEquations(mapping, params, hklList, system, refinementMode,
                                  sqrtW, residuals, tthAxis, y_obs, bkg, patternScratch,
                                  n_points, fitFlags) {
    const n = mapping.length;
    if (n === 0) return null;
    const n_peaks = hklList.length;

    enforceSymmetryConstraintsWorker(params, system);
    updateHklPositions(hklList, params, system);
    const contributions = buildPeakContributions(tthAxis, hklList, params);
    const buckets = bucketContributionsByPeak(contributions, n_peaks);
    const win = peakWindows(buckets, n_peaks);

    const baseline = new Float64Array(n_points);
    calculatePatternCPU(tthAxis, hklList, params, patternScratch);
    for (let i = 0; i < n_points; i++) baseline[i] = patternScratch[i] + bkg[i];
    for (let i = 0; i < n_points; i++) residuals[i] = (y_obs[i] - baseline[i]) * sqrtW[i];

    const columns = new Array(n);
    const probe = new Float64Array(n_points);
    const trial = new Float64Array(n_points);
    const saved = new Float64Array(n_peaks);
    for (let i = 0; i < n_peaks; i++) saved[i] = hklList[i] ? hklList[i].intensity : 0;

    for (let p = 0; p < n; p++) {
        const m = mapping[p];
        if (m.isIntensity) {
            const col = peakProfileColumn(tthAxis, buckets[m.index],
                                          win.start[m.index], win.stop[m.index], sqrtW, null);
            columns[p] = col ? { start: win.start[m.index], values: col }
                             : { start: 0, values: new Float64Array(0) };
            continue;
        }
        const v0 = m.get(params, hklList);
        const mag = Math.max(Math.abs(v0), m.defaultScale || 0, 1e-9);
        const h = Math.max(1e-9, mag * FD_BASE_FRACTION);
        m.set(params, hklList, v0 + h);
        const actual = m.get(params, hklList) - v0;
        probe.fill(0);
        if (actual !== 0) {
            enforceSymmetryConstraintsWorker(params, system);
            updateHklPositions(hklList, params, system);
            if (m.isBackground) {
                calculateTotalBackground(tthAxis, params, workerBackgroundAnchors, trial);
                for (let i = 0; i < n_points; i++) trial[i] += patternScratch[i];
            } else {
                calculatePatternCPU(tthAxis, hklList, params, trial);
                for (let i = 0; i < n_points; i++) trial[i] += bkg[i];
            }
            for (let i = 0; i < n_points; i++) probe[i] = (trial[i] - baseline[i]) / actual * sqrtW[i];
        }
        m.set(params, hklList, v0);
        for (let i = 0; i < n_peaks; i++) if (hklList[i]) hklList[i].intensity = saved[i];
        enforceSymmetryConstraintsWorker(params, system);
        updateHklPositions(hklList, params, system);
        columns[p] = { start: 0, values: probe.slice() };
    }

    // PLAIN Array rows, not Float64Array.
    //
    // This matrix crosses into the main thread and is consumed by four places
    // that all test it with `Array.isArray(JtJ[0])` -- which is FALSE for a
    // typed array. Returning Float64Array rows made generateReportContent
    // throw "JtJ matrix has incorrect dimensions or format", and made
    // intensityEsdsFromCovariance and intensityOverlapClusters return null
    // WITHOUT SAYING ANYTHING. It also matters for serialisation:
    // JSON.stringify turns a Float64Array into {"0":..,"1":..}, not an array,
    // so an exported matrix would come back as objects.
    //
    // There is no cost. A `new Array(n)` holding only doubles is
    // PACKED_DOUBLE_ELEMENTS in V8 -- the same 8 bytes an element -- and this
    // is built once at the converged point, not per iteration.
    const JtJ = new Array(n);
    for (let i = 0; i < n; i++) JtJ[i] = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
        const ci = columns[i];
        for (let j = i; j < n; j++) {
            const cj = columns[j];
            const lo = Math.max(ci.start, cj.start);
            const hi = Math.min(ci.start + ci.values.length, cj.start + cj.values.length);
            let s = 0;
            for (let k = lo; k < hi; k++) s += ci.values[k - ci.start] * cj.values[k - cj.start];
            if (!isFinite(s)) s = 0;
            JtJ[i][j] = s; JtJ[j][i] = s;
        }
    }
    return JtJ;
}

// ---------------------------------------------------------------------------
//  MOVED. LEBAIL_PT_REFRESH_INTERVAL and the PT ladder constants
//  (PT_NUM_REPLICAS, PT_MAX_TEMP, PT_MIN_TEMP, PT_SWAP_INTERVAL,
//  PT_COST_SCALE_FLOOR, PT_RESCALE_THRESHOLD) now live in constants.js.
// ---------------------------------------------------------------------------


function refineParametersPT(initialParams, fitFlags, maxIter, hklList, system, refinementMode,
                            leBailCycle = 0, totalLeBailCycles = 1, progressWindow = null) {
    const separable = (refinementMode === 'pawley');
    const { paramMapping } = getParameterMapping(fitFlags, initialParams, hklList,
                                                 refinementMode, system, false,
                                                 { excludeIntensities: separable });

    const baseProgress = progressWindow ? progressWindow.base : (leBailCycle / totalLeBailCycles);
    const cycleProgressSpan = progressWindow ? progressWindow.span : (1 / totalLeBailCycles);
    const postProgress = (f) => Math.min(1.0, baseProgress + Math.max(0, Math.min(1, f)) * cycleProgressSpan);

    if (paramMapping.length === 0) {
        postMessage({ type: 'progress', value: postProgress(1) });
        return { params: initialParams, hklList, ss_res: 0 };
    }

    // The iteration count is whatever was asked for. PT proposes one scalar
    // per replica per iteration, so each parameter gets
    // maxIter * PT_NUM_REPLICAS / n_params proposals; the report states that
    // number and leaves the judgement to the person who set it.
    const ptIterations = {
        used: maxIter,
        nParams: paramMapping.length,
        proposalsPerParam: (maxIter * PT_NUM_REPLICAS) / paramMapping.length
    };

    const tthAxis = workerWorkingData.tth;
    const y_obs = workerWorkingData.intensity;
    const n_points = tthAxis.length;
    const n_peaks = hklList.length;

    const sqrtW = new Float64Array(n_points);
    for (let i = 0; i < n_points; i++) sqrtW[i] = Math.sqrt(workerWorkingData.weights[i]);

    const scratch_bkg = new Float64Array(n_points);
    const scratch_pattern = new Float64Array(n_points);
    const scratch_net = new Float64Array(n_points);
    const scratch_cost = new Float64Array(n_points);
    calculateTotalBackground(tthAxis, initialParams, workerBackgroundAnchors, scratch_bkg);

    // One linear-solve working set per replica: the skyline envelope depends
    // on that replica's own profile, so they cannot share a cache.
    const caches = Array.from({ length: PT_NUM_REPLICAS + 1 }, () => ({}));

    /**
     * Cost at a replica's parameters. In Pawley mode the intensities are the
     * exact linear solution at those parameters -- i.e. this is the REDUCED
     * objective, the one whose minimiser in theta is the true minimiser.
     */
    const objective = (p_obj, hkl_list_obj, cache) => {
        try {
            enforceSymmetryConstraintsWorker(p_obj, system);
            updateHklPositions(hkl_list_obj, p_obj, system);
            if (fitFlags && fitFlags.fitBackground) {
                calculateTotalBackground(tthAxis, p_obj, workerBackgroundAnchors, scratch_bkg);
            }
            if (separable) {
                for (let i = 0; i < n_points; i++) scratch_net[i] = y_obs[i] - scratch_bkg[i];
                const contributions = buildPeakContributions(tthAxis, hkl_list_obj, p_obj);
                const buckets = bucketContributionsByPeak(contributions, n_peaks);
                const win = peakWindows(buckets, n_peaks);
                cache.first = null; cache.A = null;
                const sol = solvePawleyIntensities(tthAxis, buckets, win, sqrtW, scratch_net, cache);
                if (!sol) return 1e12;
                for (let j = 0; j < n_peaks; j++) if (hkl_list_obj[j]) hkl_list_obj[j].intensity = sol.I[j];
                const c = weightedCostFromColumns(sol.I, sol.cols, win, sqrtW, scratch_net, scratch_cost);
                return isFinite(c) ? c : 1e12;
            }
            const net = calculatePatternCPU(tthAxis, hkl_list_obj, p_obj, scratch_pattern);
            let s = 0;
            for (let i = 0; i < n_points; i++) {
                const w = workerWorkingData.weights[i];
                const r = y_obs[i] - (net[i] + scratch_bkg[i]);
                if (!isFinite(w) || !isFinite(r)) return 1e12;
                s += w * r * r;
            }
            return isFinite(s) ? s : 1e12;
        } catch (err) {
            console.warn("PT objective error:", err);
            return 1e12;
        }
    };

    const refreshReplica = (rep, cache) => {
        if (refinementMode === 'le-bail') {
            enforceSymmetryConstraintsWorker(rep.params, system);
            updateHklPositions(rep.hklList, rep.params, system);
            if (fitFlags && fitFlags.fitBackground) {
                calculateTotalBackground(tthAxis, rep.params, workerBackgroundAnchors, scratch_bkg);
            }
            leBailIntensityExtraction(
                { tth: tthAxis, intensity: y_obs, background: scratch_bkg,
                  weights: workerWorkingData.weights },
                rep.hklList, rep.params, scratch_pattern);
        }
        rep.cost = objective(rep.params, rep.hklList, cache);
        return rep.cost;
    };

    const paramTemplate = JSON.parse(JSON.stringify(initialParams));
    const hklTemplate = JSON.parse(JSON.stringify(hklList));

    const seed = { params: initialParams, hklList, cost: Infinity, temp: 0 };
    const initialCost = refreshReplica(seed, caches[PT_NUM_REPLICAS]);

    const temperatures = Array.from({ length: PT_NUM_REPLICAS }, (_, i) =>
        PT_MAX_TEMP * Math.pow(PT_MIN_TEMP / PT_MAX_TEMP, i / (PT_NUM_REPLICAS - 1 || 1)));

    let costScale = Math.max(PT_COST_SCALE_FLOOR, Math.abs(initialCost));
    const maybeRescaleCost = (reps) => {
        let sum = 0, n = 0;
        for (const r of reps) if (r && isFinite(r.cost)) { sum += Math.abs(r.cost); n++; }
        if (n === 0) return false;
        const proposed = Math.max(PT_COST_SCALE_FLOOR, sum / n);
        const ratio = proposed / costScale;
        if (ratio > (1 - PT_RESCALE_THRESHOLD) && ratio < 1 / (1 - PT_RESCALE_THRESHOLD)) return false;
        costScale = proposed;
        return true;
    };

    const replicas = temperatures.map(temp => ({
        params: JSON.parse(JSON.stringify(initialParams)),
        hklList: JSON.parse(JSON.stringify(hklList)),
        cost: initialCost, temp
    }));

    let bestVector = paramMapping.map(m => m.get(initialParams, hklList));
    let bestOverallCost = initialCost;

    const probe = { params: JSON.parse(JSON.stringify(paramTemplate)),
                    hklList: JSON.parse(JSON.stringify(hklTemplate)), cost: Infinity, temp: 0 };
    const requoteBest = () => {
        paramMapping.forEach((m, i) => m.set(probe.params, probe.hklList, bestVector[i]));
        return refreshReplica(probe, caches[PT_NUM_REPLICAS]);
    };

    let lastProgressPct = -1, consecutiveErrors = 0;

    for (let iter = 0; iter < maxIter; iter++) {
        try {
            if (refinementMode === 'le-bail' && iter > 0 && iter % LEBAIL_PT_REFRESH_INTERVAL === 0) {
                for (let i = 0; i < PT_NUM_REPLICAS; i++) refreshReplica(replicas[i], caches[i]);
                bestOverallCost = requoteBest();
                maybeRescaleCost(replicas);
            }

            for (let i = 0; i < PT_NUM_REPLICAS; i++) {
                const replica = replicas[i];
                const mapping = paramMapping[Math.floor(Math.random() * paramMapping.length)];
                const original_val = mapping.get(replica.params, replica.hklList);

                const step_scale = Math.max(0.01, replica.temp);
                const step_size = (mapping.step || 0.05) * step_scale * mapping.scale;
                const proposed = original_val + (Math.random() - 0.5) * 2 * step_size;

                const minV = (mapping.minVal !== undefined) ? mapping.minVal : -Infinity;
                const maxV = (mapping.maxVal !== undefined) ? mapping.maxVal : Infinity;
                if (!(proposed >= minV && proposed <= maxV)) continue;   // reject, do not clamp

                mapping.set(replica.params, replica.hklList, proposed);
                const neighbor_cost = objective(replica.params, replica.hklList, caches[i]);
                const delta_cost = neighbor_cost - replica.cost;
                const acceptance_prob = (replica.temp > 1e-9)
                    ? Math.exp(Math.max(-700, -delta_cost / (costScale * replica.temp))) : 0;

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

            if (iter > 0 && iter % PT_SWAP_INTERVAL === 0) {
                for (let i = 0; i < PT_NUM_REPLICAS - 1; i++) {
                    const r1 = replicas[i], r2 = replicas[i + 1];
                    if (Math.abs(r1.temp - r2.temp) < 1e-9) continue;
                    const dBeta = (1 / r1.temp) - (1 / r2.temp);
                    const dCost = (r1.cost - r2.cost) / costScale;
                    if (Math.exp(Math.max(-700, Math.min(50, dBeta * dCost))) > Math.random()) {
                        [r1.params, r2.params] = [r2.params, r1.params];
                        [r1.hklList, r2.hklList] = [r2.hklList, r1.hklList];
                        [r1.cost, r2.cost] = [r2.cost, r1.cost];
                        // The caches follow their replicas: each holds a
                        // skyline envelope derived from that state's profile.
                        [caches[i], caches[i + 1]] = [caches[i + 1], caches[i]];
                    }
                }
            }

            const pct = Math.floor(postProgress((iter + 1) / maxIter) * 100);
            if (pct !== lastProgressPct) {
                lastProgressPct = pct;
                postMessage({ type: 'progress', value: postProgress((iter + 1) / maxIter) });
            }
            consecutiveErrors = 0;
        } catch (error) {
            console.error("Error during PT iteration:", iter, error);
            if (++consecutiveErrors >= 5) {
                postMessage({ type: 'error', message: `PT aborted at iteration ${iter}: ${error.message}` });
                break;
            }
        }
    }

    const bestOverallParams = JSON.parse(JSON.stringify(paramTemplate));
    const bestOverallHklList = JSON.parse(JSON.stringify(hklTemplate));
    paramMapping.forEach((m, i) => m.set(bestOverallParams, bestOverallHklList, bestVector[i]));
    const finalCost = refreshReplica({ params: bestOverallParams, hklList: bestOverallHklList,
                                       cost: Infinity, temp: 0 }, caches[PT_NUM_REPLICAS]);

    postMessage({ type: 'progress', value: postProgress(1) });
    enforceSymmetryConstraintsWorker(bestOverallParams, system);

    return {
        params: bestOverallParams,
        hklList: bestOverallHklList,
        algorithm: 'pt',
        parameterInfo: paramMapping.map(m => ({ name: m.name, scale: m.scale, isIntensity: m.isIntensity })),
        fitFlags,
        ptIterations,
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
    const polishIters = Math.max(20, Math.min(120, Math.round(maxIterations / 2)));
    try {
        postMessage({ type: 'progress', value: 0.90, message: 'Converging with Levenberg-Marquardt...' });
        const lm = refineParametersLM(
            JSON.parse(JSON.stringify(ptResults.params)), fitFlags, polishIters,
            ptResults.hklList, system, refinementMode, 0, 1,
            { base: 0.90, span: 0.10 });
        if (lm && !lm.error && lm.params && lm.hklList) {
            // ptIterations belongs to the PT stage and is not regenerated by
            // the polish, so carry it across rather than letting the merge
            // drop it.
            return Object.assign({}, lm, { algorithm: 'pt+lm',
                                           ptIterations: ptResults.ptIterations });
        }
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
 * @param {{excludeIntensities?: boolean}} [opts={}]  excludeIntensities omits
 *        the Pawley intensity parameters; they are solved exactly instead.
 *        calculateStatistics calls this purely to COUNT parameters, at the end
 *        of a fit that has already run. Without this flag the user was shown
 *        the whole "parameter X is fixed by symmetry / has no valid value"
 *        set a second time, after the refinement had finished, with nothing
 *        they could do about it -- the same warnings they had already been
 *        given when the refinement started.
 * @returns {{paramMapping: object[]}}
 */
function getParameterMapping(fitFlags, initialParams, hklList, refinementMode, system,
                             quiet = false, opts = {}) {
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
     // opts.excludeIntensities: the separable (variable-projection) path in
     // refineParametersLM / refineParametersPT solves the intensities exactly
     // instead of iterating them, so they must not appear in the LM search
     // space. Every other caller leaves opts empty and is unaffected.
     if (refinementMode === 'pawley' && !opts.excludeIntensities &&
         hklList && workerWorkingData && workerWorkingData.tth && workerWorkingData.tth.length > 0) {
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
   
    if (fitFlags && fitFlags.fitBackground && workerBackgroundAnchors && workerBackgroundAnchors.length > 0) {
        workerBackgroundAnchors.forEach((anchor, i) => {
            const bgName = `bg_y_${i}`;
            const initialValue = initialParams[bgName] !== undefined ? initialParams[bgName] : anchor.y;
            const scale = Math.max(1.0, Math.abs(initialValue));

            mappings.push({
                name: bgName,
                scale: scale,
                defaultScale: 100.0,
                step: 1.0,
                maxStepAbs: Infinity,
                isIntensity: false,
                isBackground: true,
                index: i,
                minVal: 0,
                maxVal: Infinity,
                get: (p_obj) => p_obj[bgName] ?? 0,
                set: (p_obj, hkl_list_obj, rawValue) => {
                    if (p_obj) p_obj[bgName] = Math.max(0, rawValue);
                }
            });
        });
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
                finalResults = refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList,
                                                  system, refinementMode, 0, 1, { base: 0.05, span: 0.95 });
            } else {
                // PT owns 5%-90%, the LM polish owns 90%-100%. Previously both
                // ran on (iter+1)/maxIter, so the bar reached 98% at the end of
                // PT and then RESET to 2% and crawled up again through the
                // polish. Nothing was wrong; the run just looked wedged, which
                // is what "the polishing LM seems stuck" actually was.
                finalResults = refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList,
                                                  system, refinementMode, 0, 1, { base: 0.05, span: 0.85 });
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
                finalResults = refineParametersLM(currentParams, fitFlags, maxIterations, currentHklList,
                                                  system, refinementMode, 0, 1, { base: 0.05, span: 0.95 });
            } else {
                // PT owns 5%-90%, the LM polish owns 90%-100%. Previously both
                // ran on (iter+1)/maxIter, so the bar reached 98% at the end of
                // PT and then RESET to 2% and crawled up again through the
                // polish. Nothing was wrong; the run just looked wedged, which
                // is what "the polishing LM seems stuck" actually was.
                finalResults = refineParametersPT(currentParams, fitFlags, maxIterations, currentHklList,
                                                  system, refinementMode, 0, 1, { base: 0.05, span: 0.85 });
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
            // Reflections whose pivot in the intensity normal matrix had to be
            // propped up: no data separates them from a coincident neighbour,
            // so only their group SUM is determined. The report names them.
            degenerateReflections: finalResults.degenerateReflections || [],
            ptIterations: finalResults.ptIterations || null,
            degenerateGroups: finalResults.degenerateGroups || [],
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