// refinement_core.js
// ---------------------------------------------------------------------------
//  DROP-IN REPLACEMENTS for refinement_worker.js.
//
//  Apply as follows (line numbers refer to the version reviewed):
//
//    1. Delete the FIRST solveCholesky (lines 190-255) and solveLU (337-376).
//       The first is a byte-identical duplicate of the second definition, which
//       silently wins; solveLU has been unreachable since solveLinearSystem
//       stopped falling back to it.
//
//    2. importScripts('constants.js', 'crystal.js', 'profile.js',
//                     'pawley_linear.js');
//
//    3. Replace refineParametersLM      (826-1429)  with the version below.
//       Replace refineParametersPT      (1437-1765) with the version below.
//       Replace polishWithLM            (1773-1795) with the version below.
//       Replace the hklList.forEach body inside leBailIntensityExtraction
//                                       (509-622)   with the version below.
//       Replace getParameterMapping's signature line (1816) as noted below.
//
//    4. In powder5.html and the worker, load pawley_linear.js after profile.js.
// ---------------------------------------------------------------------------

'use strict';

// ===========================================================================
//  0. New tuning constants -- move these into constants.js
// ===========================================================================

/**
 * Use the EXACT variable-projection Jacobian for the non-linear parameters,
 * i.e. re-solve the intensities inside every finite-difference probe.
 *
 * WHY false (Kaufman's approximation). The exact form measures
 * d/d(theta) of the REDUCED cost, which is what LM is really minimising, and
 * it is the better Jacobian. It also costs one extra banded solve per column
 * per iteration -- 56 ms at 2500 reflections, times ~12 columns, so about a
 * second an iteration. Kaufman (1975) showed that holding the linear
 * parameters fixed inside the probe costs typically one or two extra
 * iterations and never changes the minimum, because the ACCEPTANCE test below
 * still evaluates the trial point with fully re-solved intensities. Cheaper
 * per iteration and the same answer.
 *
 * Set true if a fit is converging slowly on a strongly correlated cell.
 * @type {boolean}
 */
const PAWLEY_VARPRO_EXACT_JACOBIAN = false;

/**
 * PT proposals per non-linear parameter needed for the search to mean
 * anything. Used only to warn.
 *
 * WHY 50: a Metropolis chain needs O(10) accepted moves per coordinate to
 * decorrelate, and acceptance on the cold rungs is well under half.
 * @type {number}
 */
const PT_MIN_PROPOSALS_PER_PARAM = 50;

// ===========================================================================
//  1. getParameterMapping -- one new option
// ===========================================================================
//  Replace the signature at line 1816 with:
//
//    function getParameterMapping(fitFlags, initialParams, hklList,
//                                 refinementMode, system, quiet = false,
//                                 opts = {}) {
//
//  and the Pawley intensity guard at line 1929 with:
//
//    if (refinementMode === 'pawley' && !opts.excludeIntensities &&
//        hklList && workerWorkingData && workerWorkingData.tth &&
//        workerWorkingData.tth.length > 0) {
//
//  Nothing else in that function changes. Every existing call site keeps its
//  behaviour, because opts defaults to {}.
// ---------------------------------------------------------------------------


// ===========================================================================
//  2. leBailIntensityExtraction -- reuse the preps that were already built
// ===========================================================================
/*
 *  The comment above line 486 claims "the contributions are built once and
 *  reused for the per-peak walk below, so prepareVoigt runs once per
 *  reflection per extraction pass instead of twice". The code did not do
 *  that: the forEach at 509 rebuilt prep1 and prep2 from scratch. Measured at
 *  2500 reflections that is 8.2 ms of duplicated work per pass, and the pass
 *  runs up to LEBAIL_MAX_PASSES = 24 times per LM iteration.
 *
 *  Replace lines 509-622 (the whole `hklList.forEach((peak, j) => {...});`)
 *  with the block below. Everything before and after it is unchanged, and the
 *  arithmetic is identical -- the preps are the same objects
 *  buildPeakContributions already made.
 */

// --- BEGIN replacement for leBailIntensityExtraction's per-peak loop -------
/* eslint-disable */
const _leBailPerPeakLoop = function (hklList, contributions, expData, params,
                                     y_obs_net, total_calc, ratio21,
                                     currentCycleIntensities, currentCycleShapeAreas,
                                     currentCycleVariances, currentCycleAnalyticAreas,
                                     pointWeights, n_points) {
    // One pass over the contributions instead of two prepareVoigt sweeps.
    const buckets = bucketContributionsByPeak(contributions, hklList.length);

    hklList.forEach((peak, j) => {
        if (!hasRenderablePosition(peak)) return;
        const bucket = buckets[j];
        if (!bucket || bucket.length === 0) return;

        // The Ka1 contribution is the one with weight 1; any other is Ka2.
        // buildPeakContributions emits them in that order, but read the flag
        // rather than the position so a future change there cannot break this.
        let prep1 = null, prep2 = null;
        let startIdx = Infinity, stopIdx = -Infinity;
        for (let c = 0; c < bucket.length; c++) {
            const con = bucket[c];
            if (con.weight === 1.0 && !prep1) prep1 = con.prep; else prep2 = con.prep;
            if (con.start < startIdx) startIdx = con.start;
            if (con.stop > stopIdx) stopIdx = con.stop;
        }
        if (!prep1) prep1 = bucket[0].prep;
        if (!(stopIdx > startIdx)) return;

        let extracted_area = 0, shape_area = 0, variance = 0;
        let prev_partitioned_I = 0, prev_profile_val = 0;
        let pend_c = 0, pend_f = 0, pend_var = 0, havePending = false;

        for (let i = startIdx; i < stopIdx; i++) {
            const current_tth = expData.tth[i];

            // Sum this reflection's own contributions at this point. Same
            // quantity the two hand-rolled range tests computed, now driven
            // by the windows buildPeakContributions already resolved.
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
                if (current_fraction > 1.0) current_fraction = 1.0;
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
                const cw = pend_c + half;
                variance += (cw * pend_f) * (cw * pend_f) * pend_var;
            }
            pend_c = half;
            pend_f = current_fraction;
            pend_var = (pointWeights && isFinite(pointWeights[i]) && pointWeights[i] > 0)
                       ? (1 / pointWeights[i]) : Math.max(1, expData.intensity[i]);
            havePending = true;
            prev_partitioned_I = current_partitioned_I;
            prev_profile_val = profile_val;
        }
        if (havePending) variance += (pend_c * pend_f) * (pend_c * pend_f) * pend_var;

        currentCycleIntensities[j] = extracted_area;
        currentCycleShapeAreas[j]  = shape_area;
        currentCycleVariances[j]   = variance;
        currentCycleAnalyticAreas[j] = (prep1.analyticArea || 0)
                                     + (prep2 ? ratio21 * (prep2.analyticArea || 0) : 0);
    });
};
/* eslint-enable */
// --- END replacement -------------------------------------------------------


// ===========================================================================
//  3. refineParametersLM -- separable in Pawley mode
// ===========================================================================
/**
 * Levenberg-Marquardt, with the Pawley intensities projected out.
 *
 * WHAT CHANGED AND WHY
 * --------------------
 * The intensities are no longer LM parameters. For any set of non-linear
 * parameters they are the exact solution of a banded linear least-squares
 * problem (see pawley_linear.js), so they are SOLVED, once per iteration, and
 * LM searches only the ~12-30 dimensions that are genuinely non-linear.
 *
 * Consequences, all of them the point of the change:
 *
 *   * The intensities are exactly optimal at every iteration. Previously the
 *     polish was capped at 40 LM steps for ~2500 coupled linear parameters,
 *     so the reported Pawley fit was a partially converged one -- and the
 *     residual it reported was, correspondingly, not the residual of the
 *     model it claimed to have fitted.
 *   * n_params in every dense block drops from ~2525 to ~25. JtJ, A_lm, the
 *     `predicted` product and solveCholesky go from O(2500^2)-O(2500^3) per
 *     iteration to nothing. The ~125 MB of per-iteration garbage that drove
 *     the worker into continuous major GC is gone with them.
 *   * lambda and the trust region are no longer being steered by 2500
 *     parameters that did not need iterating.
 *
 * Le Bail mode is structurally unchanged -- it never had intensity parameters
 * to begin with -- but picks up the analytic background columns and the
 * hoisted allocations below.
 *
 * @param {ProfileParams} initialParams
 * @param {object} fitFlags
 * @param {number} maxIter
 * @param {Reflection[]} hklList
 * @param {string} system
 * @param {'pawley'|'le-bail'} refinementMode
 * @param {number} [leBailCycle=0]
 * @param {number} [totalLeBailCycles=1]
 * @param {{base:number, span:number}} [progressWindow] Sub-range of the
 *        progress bar this run owns. Without it a PT polish reset the bar
 *        from 98% back to 2% and crawled up again, which is what "the fit
 *        seems stuck" actually looked like.
 * @returns {object}
 */
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
    let lastBuckets = null, lastWindows = null, lastCols = null, lastDegenerate = [];
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

    const JtJ = new Array(n);
    for (let i = 0; i < n; i++) JtJ[i] = new Float64Array(n);
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


// ===========================================================================
//  4. refineParametersPT -- search only what is worth searching
// ===========================================================================
/**
 * Parallel tempering over the NON-LINEAR parameters, with the Pawley
 * intensities solved exactly at every cost evaluation.
 *
 * WHY THIS IS THE WHOLE FIX FOR PT
 * --------------------------------
 * PT proposes ONE random scalar per replica per iteration. The expected
 * number of proposals a given parameter receives is
 *
 *     maxIter * PT_NUM_REPLICAS / n_params
 *
 * With 2500 Pawley intensities, 8 replicas and 200 iterations that is 0.63.
 * Most intensities were never proposed even once, so the "global search" was
 * a very expensive way of leaving the starting point almost unchanged -- and
 * every one of those 1600 proposals paid for a FULL pattern evaluation
 * (24 ms at 2500 reflections) to test a change to a single peak's height.
 *
 * Projecting the intensities out reduces the search to the ~12-30 parameters
 * that are actually non-linear and actually multi-modal -- the cell, the
 * zero-point, the widths -- which is what PT is for. At the same budget each
 * of them now gets ~130 proposals instead of 0.63.
 *
 * @returns {object}
 */
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

    // Tell the user when the budget cannot support the search, rather than
    // burning it and reporting a result that looks like a global optimum.
    const proposalsPerParam = (maxIter * PT_NUM_REPLICAS) / paramMapping.length;
    if (proposalsPerParam < PT_MIN_PROPOSALS_PER_PARAM) {
        postMessage({ type: 'warning', message:
            `Parallel tempering has ${proposalsPerParam.toFixed(1)} proposals per parameter `
          + `(${paramMapping.length} parameters, ${maxIter} iterations). Below about `
          + `${PT_MIN_PROPOSALS_PER_PARAM} the search cannot explore; raise the iteration count `
          + `or refine fewer parameters.` });
    }

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
        ss_res: isFinite(finalCost) ? finalCost : bestOverallCost
    };
}


// ===========================================================================
//  5. polishWithLM -- and the progress bar that made this look like a hang
// ===========================================================================
/**
 * Converge from the PT optimum with LM, and produce J^T J.
 *
 * FIX (the reported symptom). polishWithLM posted progress 0.98 and then
 * called refineParametersLM with leBailCycle = 0, totalLeBailCycles = 1, so
 * the polish posted its own progress from (iter+1)/maxIter -- the bar jumped
 * from 98% back to 2% and crawled up again over tens of seconds. Nothing was
 * wrong; the run just looked wedged. The polish now owns the last 10% of the
 * bar and moves monotonically through it.
 *
 * FIX (iteration budget). The old cap of 40 was sized for a problem with
 * ~2500 parameters, where it was hopeless anyway. With the intensities
 * projected out the polish has ~12-30 parameters and converges properly, so
 * the cap can be both generous and rarely reached.
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
            return Object.assign({}, lm, { algorithm: 'pt+lm' });
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


// ===========================================================================
//  6. onmessage -- give PT and LM their own slices of the bar
// ===========================================================================
/*
 *  In the two `else` branches of self.onmessage (lines 2179 and 2208),
 *  replace
 *
 *      finalResults = refineParametersPT(currentParams, fitFlags, maxIterations,
 *                                        currentHklList, system, refinementMode, 0, 1);
 *      finalResults = polishWithLM(finalResults, fitFlags, currentHklList,
 *                                  system, refinementMode, maxIterations);
 *  with
 *
 *      finalResults = refineParametersPT(currentParams, fitFlags, maxIterations,
 *                                        currentHklList, system, refinementMode, 0, 1,
 *                                        { base: 0.05, span: 0.85 });
 *      finalResults = polishWithLM(finalResults, fitFlags, currentHklList,
 *                                  system, refinementMode, maxIterations);
 *
 *  and pass { base: 0.05, span: 0.95 } to the plain-LM branches. The bar then
 *  advances monotonically from 5% to 100% for every algorithm.
 *
 *  polishWithLM's third argument `hklList` has never been read -- it uses
 *  ptResults.hklList. Left in place so the call sites do not have to change.
 */
