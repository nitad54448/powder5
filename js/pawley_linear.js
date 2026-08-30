// pawley_linear.js
// version 1, replaces the "2500 intensities are ordinary LM parameters" model.
// ---------------------------------------------------------------------------
//  SEPARABLE (VARIABLE-PROJECTION) PAWLEY
//
//  A Pawley model is
//
//      y(x) = SUM_j I_j * phi_j(x; theta) + b(x; beta)
//
//  and it is LINEAR in every I_j. Refining the intensities as ordinary
//  Levenberg-Marquardt parameters therefore throws away the one structural
//  fact that makes the problem tractable:
//
//    * LM crawls towards an answer that a single linear solve gives exactly.
//      With ~2500 reflections and a polish capped at 40 iterations, it never
//      arrives -- the reported fit is a partially converged one.
//    * Every derived quantity (lambda, the trust region, the convergence
//      test) is dominated by 2500 parameters that did not need iterating.
//    * Parallel tempering becomes meaningless: with one random scalar per
//      replica per iteration, the expected number of proposals per parameter
//      is maxIter * nReplicas / nParams. At 200 x 8 / 2525 that is 0.63 --
//      most intensities are never proposed even once, so PT spends the whole
//      time budget doing nothing while looking busy.
//
//  Instead: for any given theta (cell, zero-point, profile) and background,
//  solve for I EXACTLY, then let LM or PT search only theta. That is
//  variable projection (Golub & Pereyra 1973). The reduced objective
//
//      C(theta) = || sqrt(W) (y - b - Phi(theta) * Ihat(theta)) ||^2
//
//  has the same minimiser in theta as the full problem and a search space of
//  ~12-30 dimensions instead of ~2500.
//
//  THE NORMAL MATRIX IS BANDED. Phi^T W Phi has a non-zero entry at (j,k)
//  exactly when reflections j and k have overlapping profile windows. With
//  the reflection list sorted by 2-theta -- which generateMasterHklList
//  guarantees -- that is a narrow skyline, so the exact solve costs
//  O(n * band^2), not O(n^3), and the whole matrix is stored in O(n * band).
//
//  Loading:  browser  script src="js/pawley_linear.js"   (after profile.js)
//            worker   importScripts('pawley_linear.js')  (after profile.js)
//            node     require('./pawley_linear.js')
// ---------------------------------------------------------------------------

'use strict';

// ===========================================================================
//  1. Contribution bookkeeping
// ===========================================================================

/**
 * Group a contribution list by the reflection it belongs to.
 *
 * PERF FIX. intensityDerivativeColumn() scanned the WHOLE contribution array
 * twice per reflection, filtering on peakIndex, and buildNormalEquations()
 * called it once per reflection -- an O(n_peaks * n_contributions) sweep to
 * produce n_peaks local columns. Measured on 4000 reflections / 8000 points:
 * 114 ms per LM iteration, against 29 ms once bucketed, and the gap widens
 * with n because one side is quadratic and the other is not.
 *
 * @param {Array<{peakIndex:number}>} contributions From buildPeakContributions.
 * @param {number} nPeaks Length of the reflection list.
 * @returns {Array<Array<object>|undefined>} buckets[j] = contributions of j.
 */
function bucketContributionsByPeak(contributions, nPeaks) {
    const buckets = new Array(nPeaks);
    for (let c = 0; c < contributions.length; c++) {
        const con = contributions[c];
        const j = con.peakIndex;
        if (j < 0 || j >= nPeaks) continue;
        const b = buckets[j];
        if (b) b.push(con); else buckets[j] = [con];
    }
    return buckets;
}

/**
 * Per-reflection data window: the union of its Ka1 and Ka2 windows, as
 * half-open index bounds into the 2-theta axis.
 *
 * A reflection with no bucket, or one whose window fell entirely outside the
 * fitted range, gets start === stop === 0, i.e. an empty window. That is the
 * honest representation: nothing in the data constrains it, and the solve
 * below reports it as undetermined rather than inventing a value.
 *
 * @param {Array<Array<object>|undefined>} buckets From bucketContributionsByPeak.
 * @param {number} nPeaks
 * @returns {{start: Int32Array, stop: Int32Array}}
 */
function peakWindows(buckets, nPeaks) {
    const start = new Int32Array(nPeaks);
    const stop = new Int32Array(nPeaks);
    for (let j = 0; j < nPeaks; j++) {
        const list = buckets[j];
        if (!list || list.length === 0) { start[j] = 0; stop[j] = 0; continue; }
        let lo = Infinity, hi = -Infinity;
        for (let c = 0; c < list.length; c++) {
            if (list[c].start < lo) lo = list[c].start;
            if (list[c].stop > hi) hi = list[c].stop;
        }
        if (!(hi > lo)) { start[j] = 0; stop[j] = 0; continue; }
        start[j] = lo; stop[j] = hi;
    }
    return { start, stop };
}

/**
 * d(pattern)/d(I_j) over reflection j's window: the sum of its contributions'
 * UNIT-HEIGHT profiles, optionally multiplied point-by-point by `weight`
 * (pass sqrt(w) to get the weighted Jacobian column directly).
 *
 * Same quantity as profile.js's intensityDerivativeColumn, but reading a
 * bucket instead of filtering the whole list.
 *
 * @param {ArrayLike<number>} tthAxis
 * @param {Array<object>|undefined} bucket Contributions of this reflection.
 * @param {number} start Window start index (from peakWindows).
 * @param {number} stop  Window stop index, exclusive.
 * @param {ArrayLike<number>|null} weight Optional per-point multiplier.
 * @param {Float64Array|null} out Optional scratch of length >= stop-start.
 * @returns {Float64Array|null} Values over [start, stop), or null if empty.
 */
function peakProfileColumn(tthAxis, bucket, start, stop, weight, out) {
    const len = stop - start;
    if (!bucket || len <= 0) return null;
    const values = (out && out.length >= len) ? out.subarray(0, len) : new Float64Array(len);
    values.fill(0);
    for (let c = 0; c < bucket.length; c++) {
        const con = bucket[c];
        const prep = con.prep, w = con.weight;
        const a = con.start, z = con.stop;
        for (let i = a; i < z; i++) values[i - start] += w * evalVoigt(tthAxis[i], prep);
    }
    if (weight) for (let k = 0; k < len; k++) values[k] *= weight[start + k];
    return values;
}

// ===========================================================================
//  2. Skyline symmetric positive-definite solver
// ===========================================================================

/**
 * Packed skyline (variable-band) storage for a symmetric matrix.
 *
 * Row i holds columns [first[i], i] contiguously at ptr[i]. This is the exact
 * storage the Cholesky factor needs: the profile/skyline theorem says L has
 * the same envelope as A, so no fill-in escapes it and no zero outside it is
 * ever touched.
 *
 * The n x n dense form of the same matrix would be 50 MB at n = 2500 -- and
 * refineParametersLM allocated one PER ITERATION. This is O(n * band):
 * ~200 kB for the same problem, allocated once.
 */
class SkylineMatrix {
    /** @param {Int32Array|number[]} first first[i] = lowest column in row i. */
    constructor(first) {
        const n = first.length;
        this.n = n;
        this.first = first;
        const ptr = new Int32Array(n + 1);
        let total = 0;
        for (let i = 0; i < n; i++) { ptr[i] = total; total += i - first[i] + 1; }
        ptr[n] = total;
        this.ptr = ptr;
        this.a = new Float64Array(total);
    }
    /** @returns {number} Flat index of (i,j), j <= i, or -1 if outside the envelope. */
    idx(i, j) {
        const f = this.first[i];
        return (j < f || j > i) ? -1 : this.ptr[i] + (j - f);
    }
    get(i, j) { const k = this.idx(i, j); return k < 0 ? 0 : this.a[k]; }
    add(i, j, v) { const k = this.idx(i, j); if (k >= 0) this.a[k] += v; }
    set(i, j, v) { const k = this.idx(i, j); if (k >= 0) this.a[k] = v; }
    zero() { this.a.fill(0); }
    /** @returns {number} Bytes of the packed store. */
    get bytes() { return this.a.byteLength; }
}

/**
 * In-place MODIFIED skyline Cholesky, A + E = L L^T, writing L over A.
 *
 * @param {Float64Array} a Packed store (mutated).
 * @param {Int32Array} first
 * @param {Int32Array} ptr
 * @param {number} n
 * @param {number[]} [repairedOut] Collects the indices of repaired pivots.
 * @returns {boolean} False only if a non-finite entry was encountered, which
 *          means the matrix is corrupt rather than merely rank-deficient.
 */
function skylineCholeskyInPlace(a, first, ptr, n, repairedOut) {
    for (let i = 0; i < n; i++) {
        const fi = first[i], pi = ptr[i];
        for (let j = fi; j <= i; j++) {
            const fj = first[j], pj = ptr[j];
            let s = a[pi + (j - fi)];
            if (!isFinite(s)) return false;
            const k0 = fi > fj ? fi : fj;
            for (let k = k0; k < j; k++) s -= a[pi + (k - fi)] * a[pj + (k - fj)];
            if (i === j) {
                const orig = Math.abs(a[pi + (j - fi)]);
                const floor = 1e-14 * Math.max(1, orig);
                if (!(s > floor)) {
                    // MODIFIED CHOLESKY (Gill-Murray). The pivot is repaired
                    // IN PLACE and the factorisation continues, rather than
                    // restarting the whole thing with a larger global ridge.
                    //
                    // WHY. Exact coincidences are normal, not exceptional:
                    // 333/511, 700/522 and their relatives in a cubic F cell,
                    // and every accidental overlap elsewhere. A global ridge
                    // big enough to factor past all of them also damps the
                    // several hundred intensities that were perfectly well
                    // determined -- measured on a cubic F pattern the global
                    // ladder had to reach 1e-6 relative, a visible bias on a
                    // weak reflection. Repairing only the pivots that fail
                    // leaves every other row exactly as the data determined
                    // it and confines the arbitrariness to the degenerate
                    // subspace, where it belongs.
                    //
                    // The resulting split within a degenerate group is
                    // arbitrary; their SUM is not. repairedOut names them so
                    // the report can say which reflections are affected
                    // instead of printing a number the data cannot support.
                    s = Math.max(floor, (orig > 0 ? orig : 1) * PAWLEY_PIVOT_REPAIR);
                    if (repairedOut) repairedOut.push(i);
                }
                a[pi + (j - fi)] = Math.sqrt(s);
            } else {
                a[pi + (j - fi)] = s / a[pj + (j - fj)];
            }
        }
    }
    return true;
}

/**
 * Forward/back substitution against a factor produced by the routine above.
 * @returns {boolean} False if the result is not finite.
 */
function skylineSolveInPlace(a, first, ptr, n, x) {
    for (let i = 0; i < n; i++) {
        const fi = first[i], pi = ptr[i];
        let s = x[i];
        for (let k = fi; k < i; k++) s -= a[pi + (k - fi)] * x[k];
        x[i] = s / a[pi + (i - fi)];
    }
    for (let i = n - 1; i >= 0; i--) {
        let s = x[i];
        for (let k = i + 1; k < n; k++) {
            const fk = first[k];
            if (fk <= i) s -= a[ptr[k] + (i - fk)] * x[k];
        }
        x[i] = s / (a[ptr[i] + (i - first[i])]);
    }
    for (let i = 0; i < n; i++) if (!isFinite(x[i])) return false;
    return true;
}

/**
 * Solve A x = b for a symmetric skyline A, repairing rank deficiency locally.
 *
 * WHY THE RIDGE. Exactly coincident reflections -- 333 and 511 in a cubic
 * cell, and every other accidental coincidence -- give identical profile
 * columns, so Phi^T W Phi is exactly singular and only their SUM is
 * determined by the data. A plain Cholesky returns null there and the whole
 * intensity solve fails, on a situation that is not an error but a property
 * of the sample.
 *
 * The ridge escalates from zero, so a well-posed system is solved EXACTLY;
 * only a genuinely degenerate one is damped, and then by the smallest amount
 * that makes it factorable. What comes back for a coincident pair is the
 * minimum-norm split of their combined intensity, which is the honest answer:
 * the sum is right, the split is arbitrary, and intensityOverlapClusters in
 * the report is what tells the user so.
 *
 * @param {SkylineMatrix} A Untouched (the factorisation works on a copy).
 * @param {Float64Array} b Right-hand side (untouched).
 * @param {Float64Array} [scratchA] Reusable factor store, length >= A.a.length.
 * @param {Float64Array} [x] Reusable solution vector, length >= n.
 * @returns {{x: Float64Array, repaired: number[]}|null} null only if the
 *          matrix holds a non-finite entry. `repaired` lists the reflection
 *          indices whose pivot had to be propped up: those reflections are
 *          not individually determined, only their group sum is.
 */
function solveSkylineSPD(A, b, scratchA, x) {
    const n = A.n, a0 = A.a, first = A.first, ptr = A.ptr;
    const fac = (scratchA && scratchA.length >= a0.length)
        ? scratchA.subarray(0, a0.length) : new Float64Array(a0.length);
    const sol = (x && x.length >= n) ? x.subarray(0, n) : new Float64Array(n);

    const repaired = [];
    fac.set(a0);
    if (!skylineCholeskyInPlace(fac, first, ptr, n, repaired)) return null;
    sol.set(b.subarray(0, n));
    if (!skylineSolveInPlace(fac, first, ptr, n, sol)) return null;
    return { x: sol, repaired };
}

// ===========================================================================
//  3. The intensity solve
// ===========================================================================

/**
 * Envelope of the intensity normal matrix: first[j] is the lowest-indexed
 * reflection whose data window overlaps j's.
 *
 * Computed exactly, not from a fixed bandwidth guess. The backward walk stops
 * only when NO earlier reflection can still reach -- `maxStopPrefix` carries
 * the running maximum of stop[], so one unusually wide low-angle reflection
 * cannot be skipped over.
 *
 * @param {Int32Array} start
 * @param {Int32Array} stop
 * @returns {Int32Array} first, with first[j] <= j.
 */
function intensityEnvelope(start, stop) {
    const n = start.length;
    const first = new Int32Array(n);
    const maxStopPrefix = new Int32Array(n);
    let run = -1;
    for (let j = 0; j < n; j++) { if (stop[j] > run) run = stop[j]; maxStopPrefix[j] = run; }
    for (let j = 0; j < n; j++) {
        let lo = j;
        if (stop[j] > start[j]) {
            for (let k = j - 1; k >= 0; k--) {
                if (stop[k] > start[j] && stop[k] > start[k]) lo = k;
                else if (maxStopPrefix[k] <= start[j]) break;
            }
        }
        first[j] = lo;
    }
    return first;
}

/**
 * Build Phi^T W Phi (skyline) and Phi^T W (y - b) for the current profile.
 *
 * @param {ArrayLike<number>} tthAxis
 * @param {Array<Array<object>|undefined>} buckets
 * @param {{start:Int32Array, stop:Int32Array}} win
 * @param {ArrayLike<number>} sqrtW Per-point sqrt(weight).
 * @param {ArrayLike<number>} target y_obs - background, UNWEIGHTED.
 * @param {object} [cache] Reusable {A, rhs, cols, colBuf} from a previous call
 *        with the same envelope. Pass it back to avoid reallocating.
 * @returns {{A: SkylineMatrix, rhs: Float64Array, cols: Array<Float64Array|null>}}
 */
function buildIntensityNormalEquations(tthAxis, buckets, win, sqrtW, target, cache) {
    const n = win.start.length;
    const first = (cache && cache.A && cache.A.n === n && cache.first) ? cache.first
                                                                      : intensityEnvelope(win.start, win.stop);
    let A = (cache && cache.A && cache.A.n === n && cache.first === first) ? cache.A : new SkylineMatrix(first);
    A.zero();
    const rhs = (cache && cache.rhs && cache.rhs.length === n) ? cache.rhs : new Float64Array(n);
    rhs.fill(0);

    // WEIGHTED columns, cached for this profile. They are needed three times
    // over -- for the normal matrix, for the right-hand side, and afterwards
    // for the pattern -- and each one costs a run of evalVoigt.
    const cols = (cache && cache.cols && cache.cols.length === n) ? cache.cols : new Array(n).fill(null);
    for (let j = 0; j < n; j++) {
        cols[j] = peakProfileColumn(tthAxis, buckets[j], win.start[j], win.stop[j], sqrtW, null);
    }

    for (let j = 0; j < n; j++) {
        const cj = cols[j];
        if (!cj) continue;
        const sj = win.start[j], ej = win.stop[j];

        // Phi^T W (y - b): the column already carries sqrt(w), so the target
        // needs the other sqrt(w) to make a full w.
        let r = 0;
        for (let i = sj; i < ej; i++) r += cj[i - sj] * sqrtW[i] * target[i];
        rhs[j] = r;

        for (let k = first[j]; k <= j; k++) {
            const ck = cols[k];
            if (!ck) continue;
            const sk = win.start[k], ek = win.stop[k];
            const lo = sj > sk ? sj : sk;
            const hi = ej < ek ? ej : ek;
            if (lo >= hi) continue;
            let s = 0;
            for (let i = lo; i < hi; i++) s += cj[i - sj] * ck[i - sk];
            A.set(j, k, s);
        }
    }
    return { A, rhs, cols, first };
}

/**
 * The exact Pawley intensities for the CURRENT profile and background.
 *
 * Reflections with an empty data window are left at zero and reported in
 * `undetermined`: no data point constrains them, so any value would be an
 * invention. That is the same convention refineParametersLM already used for
 * a zero-length Jacobian slice.
 *
 * @param {ArrayLike<number>} tthAxis
 * @param {Array<Array<object>|undefined>} buckets
 * @param {{start:Int32Array, stop:Int32Array}} win
 * @param {ArrayLike<number>} sqrtW
 * @param {ArrayLike<number>} target y_obs - background.
 * @param {object} [cache] Reusable working set; see buildIntensityNormalEquations.
 * @returns {{I: Float64Array, cols: Array, degenerate: number[],
 *            degenerateGroups: number[][], undetermined: number[],
 *            cache: object}|null}
 */
function solvePawleyIntensities(tthAxis, buckets, win, sqrtW, target, cache) {
    const n = win.start.length;
    if (n === 0) return { I: new Float64Array(0), cols: [], degenerate: [], degenerateGroups: [],
                 undetermined: [], cache: cache || {} };

    const ne = buildIntensityNormalEquations(tthAxis, buckets, win, sqrtW, target, cache);
    const c = cache || {};
    c.A = ne.A; c.rhs = ne.rhs; c.cols = ne.cols; c.first = ne.first;
    if (!c.fac || c.fac.length < ne.A.a.length) c.fac = new Float64Array(ne.A.a.length);
    if (!c.sol || c.sol.length < n) c.sol = new Float64Array(n);

    // A reflection outside the fitted range contributes an all-zero row and
    // column. Give it a unit diagonal and a zero right-hand side so the factor
    // stays positive definite and the solution comes back as an exact 0,
    // rather than poisoning the whole solve with a singular pivot.
    const undetermined = [];
    for (let j = 0; j < n; j++) {
        if (!ne.cols[j] || !(ne.A.get(j, j) > 0)) {
            // Row j and column j are ALREADY all zero: every entry that
            // involves reflection j is an inner product with its (absent or
            // identically zero) column, and buildIntensityNormalEquations
            // skips those. The earlier version zeroed them explicitly, and
            // the column loop it used -- `i < n && ne.first[i] <= j` -- would
            // stop at the first row whose envelope missed j even though a
            // later, wider row could still reach back to it. Setting the
            // diagonal is the only thing actually needed, and it cannot be
            // wrong: first[j] <= j always, so (j,j) is inside the envelope.
            ne.A.set(j, j, 1);
            ne.rhs[j] = 0;
            undetermined.push(j);
        }
    }

    const res = solveSkylineSPD(ne.A, ne.rhs, c.fac, c.sol);
    if (!res) return null;

    // -------------------------------------------------------------------
    //  DEGENERATE GROUPS, not just repaired pivots.
    //
    //  The modified Cholesky repairs the pivot of the REDUNDANT member of a
    //  coincident set -- the later index. Its partner, the one that took the
    //  whole group's intensity, factors perfectly well and is not flagged.
    //  Reporting only the repaired indices is therefore actively misleading:
    //  on a cubic F pattern, 775 came back at 2340 against a true 406 because
    //  it had absorbed 11-1-1 (49+49+25 = 121+1+1 = 123, the same |Q|), and
    //  nothing in the output said so. The number looked ordinary.
    //
    //  Each repaired pivot is therefore expanded to the full set of
    //  reflections its column is parallel to, using the normalised
    //  correlation from the normal matrix itself. Groups that touch are
    //  merged, so a triple coincidence is reported once.
    // -------------------------------------------------------------------
    const groups = [];
    for (let gi = 0; gi < res.repaired.length; gi++) {
        const i = res.repaired[gi];
        const dii = ne.A.get(i, i);
        if (!(dii > 0)) continue;
        const members = [i];
        for (let j = ne.first[i]; j < i; j++) {
            const djj = ne.A.get(j, j);
            if (!(djj > 0)) continue;
            const r = ne.A.get(i, j) / Math.sqrt(dii * djj);
            if (Math.abs(r) > PAWLEY_DEGENERACY_CORRELATION) members.push(j);
        }
        // Columns above i can be parallel to it too.
        for (let k = i + 1; k < n; k++) {
            if (ne.first[k] > i) continue;
            const dkk = ne.A.get(k, k);
            if (!(dkk > 0)) continue;
            const r = ne.A.get(k, i) / Math.sqrt(dii * dkk);
            if (Math.abs(r) > PAWLEY_DEGENERACY_CORRELATION) members.push(k);
        }
        members.sort((x, y) => x - y);
        const existing = groups.find(g => members.some(m => g.includes(m)));
        if (existing) {
            for (const m of members) if (!existing.includes(m)) existing.push(m);
            existing.sort((x, y) => x - y);
        } else {
            groups.push(members);
        }
    }

    return { I: res.x, cols: ne.cols,
             degenerate: res.repaired,        // the propped-up pivots
             degenerateGroups: groups,        // every reflection each one affects
             undetermined, cache: c };
}

/**
 * Sum the weighted residual sum of squares for a given intensity vector,
 * reusing the cached weighted columns instead of re-rendering the pattern.
 *
 * @param {Float64Array} I
 * @param {Array<Float64Array|null>} cols Weighted columns from the solve.
 * @param {{start:Int32Array}} win
 * @param {ArrayLike<number>} sqrtW
 * @param {ArrayLike<number>} target y_obs - background.
 * @param {Float64Array} scratch Length n_points; overwritten with sqrt(w)*y_calc_net.
 * @returns {number} sum_i w_i (y_i - b_i - y_calc_net_i)^2
 */
function weightedCostFromColumns(I, cols, win, sqrtW, target, scratch) {
    scratch.fill(0);
    const n = I.length, start = win.start;
    for (let j = 0; j < n; j++) {
        const cj = cols[j];
        if (!cj) continue;
        const amp = I[j];
        if (amp === 0) continue;
        const s = start[j];
        for (let k = 0; k < cj.length; k++) scratch[s + k] += amp * cj[k];
    }
    let cost = 0;
    for (let i = 0; i < scratch.length; i++) {
        const r = sqrtW[i] * target[i] - scratch[i];
        if (!isFinite(r)) return Infinity;
        cost += r * r;
    }
    return cost;
}

// ---------------------------------------------------------------------------
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        bucketContributionsByPeak, peakWindows, peakProfileColumn,
        SkylineMatrix, skylineCholeskyInPlace, skylineSolveInPlace, solveSkylineSPD,
        intensityEnvelope, buildIntensityNormalEquations,
        solvePawleyIntensities, weightedCostFromColumns
    };
}
