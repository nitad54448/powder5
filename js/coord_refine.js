// Bumped whenever this file changes, same convention as swarm_wyckoff.js.
const SHARKO_COORD_REFINE_VERSION = '2026-08-13a';

/* ------------------------------------------------------------------
   coord_refine.js - the last step: Wyckoff-constrained coordinate
   refinement of a proposed structure against the Pawley intensities.

   WHY THIS EXISTS, AND WHY THE MAP FIT CANNOT DO IT.

   The charge-flipping route scores a structure by the electron density
   underneath its atoms, sampled from a grid by trilinear interpolation.
   Trilinear interpolation is multilinear: along any one axis it is a straight
   line between two grid planes, so over a grid cell its maximum is always AT A
   VERTEX. An atom sitting on real density is therefore pulled onto a grid
   node and stays there, and no amount of further searching moves it off -
   the quench descends against the same interpolant and inherits the same
   bias. The coordinates that come out are quantised to 1/N.

       N = 32  ->  0.031 fractional  ->  about 0.26 A
       N = 64  ->  0.016 fractional  ->  about 0.13 A

   Doubling the grid halves the error; it does not remove the vertex bias,
   and the cost is 8x the memory and transform time. Neither figure is good
   enough to report a bond length to, and the signature is unmistakable once
   you look for it: sites on strong density land on exact multiples of 1/N
   while sites on flat density float anywhere.

   The fix is not a finer grid. It is to stop refining against the map at all.
   The map is an intermediate - a device for finding the right basin out of
   the phase problem - and the DATA is the Pawley intensity list. So the last
   step re-refines the proposed coordinates against |Fo|^2 directly, by
   least squares, with no grid anywhere in the objective. The map fit gets the
   structure into the right basin; this gets it to the bottom of one.

   WHAT IS AND IS NOT REFINED.

   Only the free coordinates, and only those. Every site is re-projected onto
   its Wyckoff position after every step, so a site on 4c stays on (x, 1/4, z)
   exactly rather than drifting off it by a least-squares residual, and a site
   on a fixed point such as 4a does not move at all. Cell, profile and
   background belong to the Pawley/Rietveld refinement in refinement_worker.js
   and are not touched here.

   Requires, already loaded: wyckoff_assign.js (projectOnto, wyckoffFreedom).
   ------------------------------------------------------------------ */

const CR_DEFAULTS = Object.freeze({
    maxIter: 60,
    // Levenberg-Marquardt damping. Starts genuinely damped rather than at
    // pure Gauss-Newton: the starting point is a grid-snapped structure that
    // can be a quarter of an Angstrom out, which is far enough that the
    // quadratic model is not yet trustworthy.
    lambda0: 1e-3,
    lambdaUp: 10,
    lambdaDown: 0.4,
    lambdaMax: 1e8,
    // Convergence: stop when a step moves every coordinate by less than this,
    // in fractional units. 1e-6 of a 8 A axis is 1e-5 A, comfortably below
    // anything reportable.
    tolStep: 1e-6,
    // ... or when the weighted residual stops improving by a relative amount.
    tolChi: 1e-9,
    // Reflections closer in d than this are one observation, exactly as
    // swPackReflections does it. A powder cannot resolve them and pretending
    // otherwise refines against Pawley's arbitrary split of a shared peak.
    overlapTol: 0.002,
    // Trust region on one step, fractional. A least-squares step from a bad
    // start can be enormous; capping it keeps the projection meaningful and
    // stops a site jumping into a neighbouring atom's basin.
    maxStep: 0.05
});

/**
 * Groups reflections a powder cannot resolve, and sums their contributions.
 *
 * Same rule and the same tolerance as swPackReflections in swarm_wyckoff.js,
 * so the refinement is fitting the same observation the search was scoring.
 * The observation is the multiplicity-weighted sum over the group, which is
 * what the diffractometer actually measured.
 */
function crGroupReflections(obsRows, overlapTol) {
    const sorted = [...obsRows].sort((a, b) => b.d - a.d);
    const groups = [];
    let cur = null;
    for (const r of sorted) {
        if (cur && Math.abs(r.d - cur.d) / cur.d < overlapTol) {
            cur.members.push(r);
        } else {
            cur = { members: [r], d: r.d, Iobs: 0 };
            groups.push(cur);
        }
    }
    for (const g of groups) {
        g.Iobs = g.members.reduce((s, r) => s + r.mult * r.Fo2, 0);

        // The group's VARIANCE, summed the same way its intensity is.
        //
        // The group is a sum, so its members' variances add. Treating them as
        // independent understates the true covariance -- reflections this
        // close in d are correlated in the Pawley decomposition -- but the
        // off-diagonal terms are not carried this far, and summing variances
        // is conservative rather than optimistic.
        //
        // null, not zero, when any member lacks a sigma: a group with no
        // uncertainty is not a group with zero uncertainty, and the caller
        // has to be able to tell the difference to pick a fallback weight.
        let varSum = 0, haveAll = g.members.length > 0;
        for (const r of g.members) {
            if (Number.isFinite(r.sigma) && r.sigma > 0) {
                const sm = r.mult * r.sigma;
                varSum += sm * sm;
            } else { haveAll = false; }
        }
        g.varObs = (haveAll && varSum > 0) ? varSum : null;
        // TWO resolution variables, because two conventions are in play and
        // conflating them cost a factor of two.
        //
        //   sinTh = sin(theta)/lambda = 1/(2d)   -- "stol"
        //   sInvD = 1/d = 2 * stol
        //
        // makeFormFactor takes 1/d and halves it internally to reach stol, and
        // swScatteringTable feeds it 1/d correctly. This function was passing
        // sinTh, so the form factor was evaluated at HALF the true resolution:
        // f came out 40% too high for Pb at d = 0.83 A, and the error grows
        // with angle, which is exactly the shape a wrong structure has. The
        // search and the refinement were therefore still using different
        // scattering even after both were given the same tables.
        g.sinTh = 1 / (2 * g.d);
        g.sInvD = 1 / g.d;
    }
    return groups;
}

/**
 * Expands the asymmetric unit, keeping the chain rule intact.
 *
 * Each generated atom records WHICH site it came from and the transpose of
 * the rotation that produced it. That transpose is the whole reason this is
 * a separate step: for an atom at r = R u + t, the phase is
 *
 *     2 pi h . r = 2 pi (R^T h) . u + 2 pi h . t
 *
 * so the derivative of the phase with respect to the free site coordinate u
 * is R^T h, not h. Getting that wrong gives a gradient that is right only for
 * the identity operator and quietly wrong for every other one, which shows up
 * as a refinement that will not converge in low-symmetry settings and
 * converges to the wrong place in high-symmetry ones.
 */
function crExpand(sites, symOps) {
    const atoms = [];
    for (let si = 0; si < sites.length; si++) {
        const s = sites[si];
        const seen = [];
        for (const op of symOps) {
            const r = op.r, t = op.t || [0, 0, 0];
            const q = [
                r[0] * s.x + r[1] * s.y + r[2] * s.z + t[0],
                r[3] * s.x + r[4] * s.y + r[5] * s.z + t[1],
                r[6] * s.x + r[7] * s.y + r[8] * s.z + t[2]
            ].map(v => ((v % 1) + 1) % 1);
            // Distinct positions only: on a special position several
            // operators map to the same point and counting it twice would
            // double that site's scattering.
            //
            // THE COMPARISON HAS TO WRAP. Coordinates are reduced to [0, 1),
            // so an atom at x = 1e-9 and its image at x = 0.999999999 are the
            // SAME point one lattice translation apart -- but they differ by
            // almost exactly 1, so a plain |a - b| test called them distinct.
            // Any site sitting on or near a coordinate of zero therefore had
            // its multiplicity doubled and scattered twice as strongly as it
            // should, which is a large systematic error on exactly the special
            // positions a Wyckoff search puts atoms on.
            const same = (a, b) => { const d = a - b; return Math.abs(d - Math.round(d)) < 1e-6; };
            if (seen.some(p => same(p[0], q[0]) && same(p[1], q[1]) && same(p[2], q[2]))) continue;
            seen.push(q);
            atoms.push({
                siteIdx: si, x: q[0], y: q[1], z: q[2],
                // R^T, row-major, for the chain rule above.
                rt: [r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]]
            });
        }
    }
    return atoms;
}

/**
 * |Fcalc|^2 for one reflection, and its derivatives w.r.t. every free coord.
 *
 * F = sum_j f_j exp(2 pi i h.r_j), so with A and B the real and imaginary
 * parts, |F|^2 = A^2 + B^2 and
 *
 *     d|F|^2/du = 2 (A dA/du + B dB/du)
 *     dA/du = -2 pi f (R^T h) sin(2 pi h.r)
 *     dB/def =  2 pi f (R^T h) cos(2 pi h.r)
 *
 * grad is accumulated into, 3 slots per site, and is left untouched when the
 * caller passes null (the value-only path used by the line search).
 */
function crStructureFactor(h, k, l, atoms, fOf, nSites, grad) {
    const TWO_PI = 2 * Math.PI;
    let A = 0, B = 0;
    const n = atoms.length;
    const cs = new Float64Array(n), sn = new Float64Array(n), ff = new Float64Array(n);
    for (let j = 0; j < n; j++) {
        const a = atoms[j];
        const f = fOf(a.siteIdx);
        const ph = TWO_PI * (h * a.x + k * a.y + l * a.z);
        const c = Math.cos(ph), s = Math.sin(ph);
        cs[j] = c; sn[j] = s; ff[j] = f;
        A += f * c; B += f * s;
    }
    const F2 = A * A + B * B;
    if (!grad) return F2;

    for (let j = 0; j < n; j++) {
        const a = atoms[j], rt = a.rt, f = ff[j];
        // R^T h - the reflection index seen in the site's own frame.
        const hx = rt[0] * h + rt[1] * k + rt[2] * l;
        const hy = rt[3] * h + rt[4] * k + rt[5] * l;
        const hz = rt[6] * h + rt[7] * k + rt[8] * l;
        // 2 (A dA + B dB) with dA = -2pi f hx sin, dB = 2pi f hx cos
        const w = 2 * TWO_PI * f * (B * cs[j] - A * sn[j]);
        const b = a.siteIdx * 3;
        grad[b]     += w * hx;
        grad[b + 1] += w * hy;
        grad[b + 2] += w * hz;
    }
    return F2;
}

/**
 * Wyckoff-constrained least-squares refinement of coordinates against |Fo|^2.
 *
 * @param {Object} o
 *   sites      [{ element, zn, x, y, z, w }]  zn is the ATOMIC NUMBER and
 *              x, y, z are fractional coordinates - see the note at fTab
 *              below; w is the Wyckoff position object from the space-group
 *              JSON, as carried by the search results
 *   symOps     [{ r:[9], t:[3] }]
 *   obsRows    [{ h, k, l, d, mult, Fo2 }]   Pawley observations
 *   formFactor (element, s) -> f, optional; falls back to the atomic number,
 *              i.e. point atoms, which is honest rather than silently wrong
 *   overallB   isotropic B, applied on the MODEL side as exp(-B s^2 / 4)
 *   maxIter, overlapTol, maxStep
 *
 * @returns { sites, R, Rstart, chi2, chi2Start, iterations, converged, scale,
 *            shifts }  shifts is the fractional distance each site moved.
 */
function refineCoordinatesAgainstPawley(o) {
    const opt = { ...CR_DEFAULTS, ...o };
    const symOps = o.symOps || [];
    const nSites = o.sites.length;
    if (!nSites || !symOps.length) {
        return { error: 'Nothing to refine: no sites or no symmetry operators.' };
    }
    const rows = (o.obsRows || []).filter(r =>
        Number.isFinite(r.Fo2) && r.Fo2 >= 0 && Number.isFinite(r.d) && r.d > 0 &&
        !(r.h === 0 && r.k === 0 && r.l === 0));
    if (rows.length < 3 * nSites) {
        return { error: `Only ${rows.length} usable observations for ${nSites} sites; ` +
                        `the refinement would be underdetermined.` };
    }

    const groups = crGroupReflections(rows, opt.overlapTol);
    const B = Number.isFinite(o.overallB) ? o.overallB : 0;
    const ffn = o.formFactor || (() => null);

    // Scattering factor per site per group, folded with the Debye-Waller
    // term, precomputed because it depends only on s and never on the
    // coordinates being refined.
    //
    // ATOMIC NUMBER IS `zn`, NEVER `z`. `z` is the fractional coordinate,
    // everywhere in this file and everywhere in the rest of the search. The
    // two are both small positive numbers, a swap produces a structure factor
    // that is wrong but perfectly finite, and nothing downstream would flag
    // it - the refinement would simply converge somewhere slightly wrong.
    const fTab = groups.map(g => {
        // exp(-B (sin(theta)/lambda)^2), written on 1/d as exp(-B s^2 / 4) so
        // it matches swScatteringTable term for term. It was exp(-B stol^2 / 4)
        // here, which is the same expression fed the wrong variable and comes
        // out as a QUARTER of the intended exponent -- at d = 0.83 A and
        // B = 1.5 the damping was 0.87 where it should have been 0.58.
        const dw = B > 0 ? Math.exp(-B * g.sInvD * g.sInvD / 4) : 1;
        return o.sites.map(s => {
            let f = ffn(s.element, g.sInvD);
            if (!Number.isFinite(f)) f = s.zn;
            if (!Number.isFinite(f)) {
                throw new Error(`Site ${s.element} has no scattering factor and no zn ` +
                                `(atomic number). Supply zn on every site, or a ` +
                                `formFactor callback that covers every element.`);
            }
            return f * dw;
        });
    });

    // Working copy: the caller's sites are not mutated.
    const cur = o.sites.map(s => ({ ...s }));
    const start = cur.map(s => ({ x: s.x, y: s.y, z: s.z }));
    const nP = nSites * 3;

    const project = (arr) => {
        const tmp = new Float64Array(3);
        for (let i = 0; i < nSites; i++) {
            const s = arr[i];
            if (s.w) {
                projectOnto(s.w, s.x, s.y, s.z, tmp);
                s.x = ((tmp[0] % 1) + 1) % 1;
                s.y = ((tmp[1] % 1) + 1) % 1;
                s.z = ((tmp[2] % 1) + 1) % 1;
            } else {
                s.x = ((s.x % 1) + 1) % 1;
                s.y = ((s.y % 1) + 1) % 1;
                s.z = ((s.z % 1) + 1) % 1;
            }
        }
    };
    project(cur);

    // WEIGHTS: 1/sigma^2 where the decomposition gave a sigma, 1/I otherwise.
    //
    // 1/I is counting statistics, and counting statistics are the wrong model
    // for a Pawley intensity. The uncertainty on one of these numbers comes
    // from the DECOMPOSITION, not from how many photons arrived, and for a
    // heavily overlapped reflection it is far larger than sqrt(I). A measured
    // 200 with an ESD of 500 gets weight 1/200 = 5e-3 under counting
    // statistics and 1/500^2 = 4e-6 from its actual uncertainty -- three
    // orders of magnitude of over-trust in a number consistent with anything
    // from -300 to 700.
    //
    // This is the same rule swPackReflections applies in the search, so the
    // two stages now agree about how much each observation is worth. They did
    // not before: the search down-weighted a poorly determined reflection a
    // hundredfold while the refinement still listened to it.
    //
    // The 1/I fallback keeps its floor, so a near-zero observation with no
    // sigma cannot dominate by having an enormous weight on mostly noise.
    const Imax = groups.reduce((m, g) => Math.max(m, g.Iobs), 0) || 1;
    let nSigmaWeighted = 0;
    const wt = groups.map(g => {
        if (Number.isFinite(g.varObs) && g.varObs > 0) { nSigmaWeighted++; return 1 / g.varObs; }
        return 1 / Math.max(g.Iobs, 1e-4 * Imax);
    });
    // Only the RATIOS of the weights matter to the least squares, and mixing
    // 1/sigma^2 with 1/I puts two different scales in one array. Normalising
    // to a mean of one keeps chi-square a readable number instead of one that
    // depends on which fallback happened to dominate.
    {
        const meanW = wt.reduce((a, b) => a + b, 0) / Math.max(1, wt.length);
        if (meanW > 0) for (let i = 0; i < wt.length; i++) wt[i] /= meanW;
    }

    // THE JACOBIAN HAS TO LIVE IN THE CONSTRAINED SPACE.
    //
    // A site on 4c has two free parameters, not three, and the third
    // direction is not merely unhelpful - it is not there. Building the
    // normal equations in the full 3N ambient space and projecting the
    // resulting step afterwards does not work: least squares chooses the step
    // that best reduces the residual in ambient space, projection then throws
    // part of it away, and what arrives is a direction nobody optimised. The
    // predicted decrease and the actual decrease disagree, every trial step
    // is rejected, lambda runs to its ceiling and the refinement reports
    // convergence without having moved. That is exactly what this did before.
    //
    // Since the position is r = P u + T, the chain rule gives dr/du = P, so
    // the constrained Jacobian is J_r P and its transpose-product is P^T J^T.
    // Transforming the accumulated gradient by P^T here, and the solved step
    // by P below, is the whole fix. P is idempotent, so a step already in the
    // subspace is left alone by it.
    const projMat = cur.map(s => {
        if (!s.w) return [1, 0, 0, 0, 1, 0, 0, 0, 1];
        return wyckoffProjector(s.w).P;
    });
    const applyPT = (g) => {
        for (let i = 0; i < nSites; i++) {
            const P = projMat[i], b = i * 3;
            const g0 = g[b], g1 = g[b + 1], g2 = g[b + 2];
            g[b]     = P[0] * g0 + P[3] * g1 + P[6] * g2;
            g[b + 1] = P[1] * g0 + P[4] * g1 + P[7] * g2;
            g[b + 2] = P[2] * g0 + P[5] * g1 + P[8] * g2;
        }
    };
    const applyP = (d) => {
        for (let i = 0; i < nSites; i++) {
            const P = projMat[i], b = i * 3;
            const d0 = d[b], d1 = d[b + 1], d2 = d[b + 2];
            d[b]     = P[0] * d0 + P[1] * d1 + P[2] * d2;
            d[b + 1] = P[3] * d0 + P[4] * d1 + P[5] * d2;
            d[b + 2] = P[6] * d0 + P[7] * d1 + P[8] * d2;
        }
    };

    /** Model intensities and, optionally, their Jacobian rows. */
    const evaluate = (state, J) => {
        const atoms = crExpand(state, symOps);
        const Ical = new Float64Array(groups.length);
        const g = J ? new Float64Array(nP) : null;
        for (let gi = 0; gi < groups.length; gi++) {
            const grp = groups[gi], ftab = fTab[gi];
            if (g) g.fill(0);
            let acc = 0;
            for (const r of grp.members) {
                const F2 = crStructureFactor(r.h, r.k, r.l, atoms,
                                             si => ftab[si], nSites, g);
                acc += r.mult * F2;
                if (g) {
                    // Into the constrained space before it is stored, so the
                    // normal equations below never see the frozen directions.
                    applyPT(g);
                    for (let p = 0; p < nP; p++) J[gi * nP + p] += r.mult * g[p];
                    g.fill(0);
                }
            }
            Ical[gi] = acc;
        }
        return Ical;
    };

    /**
     * The scale factor is not a refined parameter.
     *
     * For a fixed set of coordinates the least-squares optimum in s is
     * available in closed form, so solving for it each time removes one
     * parameter from the normal equations and removes the strong correlation
     * it would otherwise have with every coordinate. Refining it alongside
     * the coordinates is a common way to make this kind of fit crawl.
     */
    const bestScale = (Ical) => {
        let num = 0, den = 0;
        for (let i = 0; i < groups.length; i++) {
            num += wt[i] * groups[i].Iobs * Ical[i];
            den += wt[i] * Ical[i] * Ical[i];
        }
        return den > 0 ? num / den : 0;
    };

    const chiOf = (Ical, s) => {
        let c = 0;
        for (let i = 0; i < groups.length; i++) {
            const d = groups[i].Iobs - s * Ical[i];
            c += wt[i] * d * d;
        }
        return c;
    };
    // wR(F^2) = sqrt( sum w (Io - s Ic)^2 / sum w Io^2 ).
    //
    // The WEIGHTED residual, because that is the quantity being minimised. An
    // unweighted R reported next to a weighted refinement describes a fit
    // nobody performed: it lets a badly determined reflection, which the least
    // squares correctly ignored, drag the reported number around anyway.
    const rOf = (Ical, s) => {
        let num = 0, den = 0;
        for (let i = 0; i < groups.length; i++) {
            const d = groups[i].Iobs - s * Ical[i];
            num += wt[i] * d * d;
            den += wt[i] * groups[i].Iobs * groups[i].Iobs;
        }
        return den > 0 ? Math.sqrt(num / den) : NaN;
    };

    let Ical = evaluate(cur, null);
    let scale = bestScale(Ical);
    let chi = chiOf(Ical, scale);
    const chi2Start = chi, Rstart = rOf(Ical, scale);

    let lambda = opt.lambda0, iter = 0, converged = false;
    const J = new Float64Array(groups.length * nP);
    const JTJ = new Float64Array(nP * nP);
    const JTr = new Float64Array(nP);

    for (iter = 0; iter < opt.maxIter; iter++) {
        J.fill(0);
        Ical = evaluate(cur, J);
        scale = bestScale(Ical);
        chi = chiOf(Ical, scale);

        // Normal equations in the standard model-fitting form:
        //
        //     (J^T W J) dx = J^T W (Iobs - s Ical),   J = d(s Ical)/dp
        //
        // J is the derivative of the MODEL, not of the residual. Those differ
        // by a sign, and using the residual derivative here while keeping the
        // right-hand side as J^T W r solves for the exact negative of the
        // correct step - every trial then goes uphill, every one is rejected,
        // lambda runs to its ceiling and the refinement stops on iteration
        // one having moved nothing. The symptom is distinctive: chiTr
        // approaches chi from ABOVE as lambda grows, because the only step
        // that does no harm is the zero one.
        JTJ.fill(0); JTr.fill(0);
        for (let gi = 0; gi < groups.length; gi++) {
            const base = gi * nP, w = wt[gi];
            const res = groups[gi].Iobs - scale * Ical[gi];
            for (let p = 0; p < nP; p++) {
                const jp = scale * J[base + p];
                if (jp === 0) continue;
                JTr[p] += w * jp * res;
                for (let q = p; q < nP; q++) {
                    JTJ[p * nP + q] += w * jp * (scale * J[base + q]);
                }
            }
        }
        for (let p = 0; p < nP; p++)
            for (let q = 0; q < p; q++) JTJ[p * nP + q] = JTJ[q * nP + p];

        let stepped = false;
        while (lambda < opt.lambdaMax) {
            // Marquardt damping, scaled by the diagonal so the trust region
            // follows the curvature of each parameter rather than being the
            // same absolute size for a well-determined coordinate and a flat
            // direction. A fixed site contributes a zero row and column, and
            // the +lambda*|diag| floor is what keeps that solvable.
            const Aug = Float64Array.from(JTJ);
            let dmax = 0;
            for (let p = 0; p < nP; p++) dmax = Math.max(dmax, Aug[p * nP + p]);
            if (!(dmax > 0)) { converged = true; break; }
            for (let p = 0; p < nP; p++) {
                Aug[p * nP + p] += lambda * Math.max(Aug[p * nP + p], 1e-6 * dmax);
            }

            const dx = crSolve(Aug, JTr, nP);
            if (!dx) { lambda *= opt.lambdaUp; continue; }

            // Back out of the constrained space: the step was solved for u,
            // the position moves by P du. project() below is then only a
            // guard against accumulated rounding, not the thing doing the
            // constraining.
            applyP(dx);

            // Trust region, then projection back onto the Wyckoff subspaces.
            let big = 0;
            for (let p = 0; p < nP; p++) big = Math.max(big, Math.abs(dx[p]));
            const cap = big > opt.maxStep ? opt.maxStep / big : 1;

            const trial = cur.map((s, i) => ({
                ...s,
                x: s.x + cap * dx[i * 3],
                y: s.y + cap * dx[i * 3 + 1],
                z: s.z + cap * dx[i * 3 + 2]
            }));
            project(trial);

            const Itr = evaluate(trial, null);
            const sTr = bestScale(Itr);
            const chiTr = chiOf(Itr, sTr);

            if (chiTr < chi) {
                let moved = 0;
                for (let i = 0; i < nSites; i++) {
                    moved = Math.max(moved,
                        Math.abs(trial[i].x - cur[i].x),
                        Math.abs(trial[i].y - cur[i].y),
                        Math.abs(trial[i].z - cur[i].z));
                }
                const rel = (chi - chiTr) / Math.max(chi, 1e-300);
                for (let i = 0; i < nSites; i++) {
                    cur[i].x = trial[i].x; cur[i].y = trial[i].y; cur[i].z = trial[i].z;
                }
                scale = sTr; chi = chiTr;
                lambda = Math.max(lambda * opt.lambdaDown, 1e-12);
                stepped = true;
                if (moved < opt.tolStep || rel < opt.tolChi) converged = true;
                break;
            }
            lambda *= opt.lambdaUp;
        }
        if (converged || !stepped) break;
    }

    Ical = evaluate(cur, null);
    scale = bestScale(Ical);

    const shifts = cur.map((s, i) => {
        const dx = s.x - start[i].x, dy = s.y - start[i].y, dz = s.z - start[i].z;
        const w = v => v - Math.round(v);
        return Math.sqrt(w(dx) ** 2 + w(dy) ** 2 + w(dz) ** 2);
    });

    // THE FREE PARAMETER COUNT IS NOT 3N.
    //
    // nP is the AMBIENT dimension the normal equations are assembled in; the
    // constraint enters through the projector, not through the parameter
    // count. Reporting nP as "free coordinates" said 15 for a five-site PbSO4
    // model whose real freedom is nine -- and understating the constraint is
    // exactly backwards, because putting atoms on Wyckoff positions IS the
    // point of the method.
    //
    // For an idempotent projector the rank equals the trace, so the free
    // dimension of a site is just P[0] + P[4] + P[8]: 4c -> (x, 1/4, z) gives
    // diag(1,0,1) = 2, a general position gives 3, and a fully fixed site like
    // 4b at (0,0,1/2) gives 0. Taken from the matrix rather than a database
    // field, so it stays correct even where n_free is missing.
    const nFree = projMat.reduce(
        (acc, P) => acc + Math.round(P[0] + P[4] + P[8]), 0);

    return {
        sites: cur.map(s => ({ ...s })),
        scale,
        chi2: chi, chi2Start,
        R: rOf(Ical, scale), Rstart,
        // So the caller can say wR rather than R, and say what it was weighted
        // on. A residual whose weighting is unstated is not interpretable.
        weighted: true,
        nSigmaWeighted, nGroupsTotal: groups.length,
        iterations: iter + 1, converged,
        nObs: groups.length,
        nParams: nFree,        // free coordinates, after the Wyckoff constraints
        nAmbient: nP,          // 3N, the space the normal equations live in
        shifts
    };
}

/** Cholesky with a partial-pivot LU fallback. Returns null if singular. */
function crSolve(A, b, n) {
    const L = new Float64Array(n * n);
    let ok = true;
    for (let i = 0; i < n && ok; i++) {
        for (let j = 0; j <= i; j++) {
            let s = A[i * n + j];
            for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
            if (i === j) {
                if (!(s > 0)) { ok = false; break; }
                L[i * n + i] = Math.sqrt(s);
            } else {
                L[i * n + j] = s / L[j * n + j];
            }
        }
    }
    if (!ok) return null;
    const y = new Float64Array(n), x = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        let s = b[i];
        for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k];
        y[i] = s / L[i * n + i];
    }
    for (let i = n - 1; i >= 0; i--) {
        let s = y[i];
        for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
        x[i] = s / L[i * n + i];
    }
    for (let i = 0; i < n; i++) if (!Number.isFinite(x[i])) return null;
    return x;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { SHARKO_COORD_REFINE_VERSION, CR_DEFAULTS,
                       refineCoordinatesAgainstPawley,
                       crGroupReflections, crExpand, crStructureFactor, crSolve };
}
