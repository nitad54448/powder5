/* ------------------------------------------------------------------
   Wyckoff-space search for sHarko.

   The swarm currently parametrises every independent atom as three free
   coordinates in a general position, so the cell contents it produces are
   whatever the symmetry happens to generate - Pb8S8O24 rather than PbSO4.
   Composition is checked after the fact and cannot be imposed.

   This module inverts that. A trial structure is described by

       (Wyckoff assignment) x (free parameters of that assignment)

   where an assignment picks, for each element, a multiset of Wyckoff
   positions whose multiplicities sum to the required number of atoms per
   cell. Every structure the swarm can then express is stoichiometric by
   construction, and the dimensionality drops: a site on 4c (x,1/4,z) has two
   parameters, one on 4a (0,0,0) has none.

   Three things are exported:

     enumerateAssignments()  the combinatorics, with pruning
     rankAssignments()       a cheap prior, using Harker evidence when present
     buildAssignmentTables() flattening for the multi-assignment GPU kernel

   plus projectSites(), which is what keeps the PSO inside the subspace.
   ------------------------------------------------------------------ */

/**
 * Module version.
 *
 * It was in module.exports and defined nowhere. The browser never runs that
 * block so nothing showed; under Node it is a ReferenceError at load, which
 * means the unit tests these exports exist for could not require the module.
 */
const SHARKO_WYCKOFF_ASSIGN_VERSION = '1.1.0';

/* ------------------------------------------------------------------
   KERNEL CAPACITIES.

   These are compile-time array bounds in swarm_reflection.wgsl, not
   preferences, and this file is where they are
   enforced because this file is what builds the tables that would overrun
   them. WGSL does not trap an out-of-range index - it clamps it - so every
   one of these overflows produces a plausible-looking structure scored
   against somebody else's data rather than an error anyone can see.

     WY_MAX_SITES  width of prop_sites: independent Wyckoff sites in ONE
                   assignment, summed over all elements. Must equal MAX_SITES
                   in swarm_reflection.wgsl -- the two are checked against each
                   other only by this file throwing before a dispatch, so
                   raising one without the other silently truncates instead.
     WY_MAX_ELEM   local_rMin is WY_MAX_ELEM^2 floats of workgroup memory.
     WY_MAX_OPS    genOp is packed into 12 bits. Group orders top out at 192.
     WY_MAX_RULES  MAX_BOND_RULES.
     WY_MAX_COORD_SLOTS  width of the private `slot` array, shared across every
                   counted rule. Private, not workgroup - it does not compete
                   with MAX_GEN_ATOMS for the 16 kB, and the kernel initialises
                   only the slots the rules claimed - so 64 (= WY_MAX_RULES x 8)
                   is affordable and stops this bound binding on any real
                   coordination shell. PbSO4 with S-O 4, Pb-O 8 and Pb-S 8
                   needs 20, which the old 16 refused.
   ------------------------------------------------------------------ */
const WY_MAX_SITES = 32;
const WY_MAX_ELEM = 8;
const WY_MAX_OPS = 4096;
const WY_MAX_RULES = 8;
const WY_MAX_COORD_SLOTS = 64;

/* ------------------------------------------------------------------ */
/*  Wyckoff position helpers                                           */
/* ------------------------------------------------------------------ */

/**
 * Degrees of freedom of a position.
 *
 * Uses the generator's exact rank(P). There is deliberately NO string fallback:
 * cctbx does not always print a special operator in the tidy form one expects.
 * The projector for (x,x,x) in Ia-3 comes out as
 *
 *     1/3*x+1/3*y+1/3*z, 1/3*x+1/3*y+1/3*z, 1/3*x+1/3*y+1/3*z
 *
 * which mentions all three symbols but has rank 1. Counting symbols would call
 * that three free parameters, giving the site three search dimensions where it
 * has one - and the excess two would be a flat direction the swarm wastes its
 * time in, with nothing to signal the mistake. Better to fail loudly: if
 * n_free is missing the database was built by an older generator, and the
 * answer is to re-run it rather than to guess.
 */
function wyckoffFreedom(w) {
    if (Number.isFinite(w.n_free)) return w.n_free;
    throw new Error(
        `Wyckoff position ${w.multiplicity || ''}${w.letter || '?'} has no n_free. ` +
        `Regenerate the space-group database with cctbx_Harko_v1.py - the degrees ` +
        `of freedom cannot be recovered reliably from the special_op string.`);
}

/** Projection matrix P and offset T as plain floats, cached on the position. */
function wyckoffProjector(w) {
    if (w._proj) return w._proj;
    let P, T;
    if (Array.isArray(w.P_num) && w.P_num.length === 9 && w.P_den) {
        P = w.P_num.map(v => v / w.P_den);
        T = (w.T_num || [0, 0, 0]).map(v => v / (w.T_den || 1));
    } else {
        // No matrix in the database: treat as a general position. This is the
        // honest fallback - it constrains nothing rather than constraining
        // something wrong - but the multiplicity will then be whatever the
        // group generates, so the caller should warn.
        P = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        T = [0, 0, 0];
        w._degraded = true;
    }
    w._proj = { P, T };
    return w._proj;
}

/**
 * Projects a fractional point onto a Wyckoff position, writing into `out`.
 *
 * special_op is idempotent, so this is a true projection: applying it to a
 * point already on the position leaves it alone. That is what allows the PSO
 * update to run untouched - see projectSites().
 */
function projectOnto(w, x, y, z, out) {
    const { P, T } = wyckoffProjector(w);
    out[0] = P[0] * x + P[1] * y + P[2] * z + T[0];
    out[1] = P[3] * x + P[4] * y + P[5] * z + T[1];
    out[2] = P[6] * x + P[7] * y + P[8] * z + T[2];
    return out;
}

const _tmp3 = new Float64Array(3);

/**
 * How far a point has to move, in Angstrom, to sit on a Wyckoff position.
 *
 * This is the quantity that turns Harker output into evidence about
 * assignments: a consolidated site 0.04 A from 4c but 1.8 A from 4a says the
 * heavy atom is on 4c, before any swarm has run.
 */
function wyckoffShift(w, p, metric) {
    projectOnto(w, p.x, p.y, p.z, _tmp3);
    let dx = _tmp3[0] - p.x, dy = _tmp3[1] - p.y, dz = _tmp3[2] - p.z;

    // `metric` is either a bare fractional->Cartesian matrix or a reduction
    // from sharkoReducedCell(). With the reduction, the rounding below happens
    // in the reduced basis, where it really is the nearest image; with a bare
    // matrix it is the old convention, which is only right for a near-
    // orthogonal cell. The shift being measured is small, so the difference
    // rarely bites - but this feeds the pre-search ranking, and an assignment
    // wrongly reported as far from a Harker peak is one that gets searched last
    // or not at all.
    const red = (metric && metric.xform) ? metric : null;
    const orth = red ? red.orth : (metric || [1, 0, 0, 0, 1, 0, 0, 0, 1]);
    if (red) {
        const X = red.xform;
        const u = X[0] * dx + X[1] * dy + X[2] * dz;
        const v = X[3] * dx + X[4] * dy + X[5] * dz;
        const t = X[6] * dx + X[7] * dy + X[8] * dz;
        dx = u; dy = v; dz = t;
    }
    dx -= Math.round(dx); dy -= Math.round(dy); dz -= Math.round(dz);

    const cx = orth[0] * dx + orth[1] * dy + orth[2] * dz;
    const cy = orth[3] * dx + orth[4] * dy + orth[5] * dz;
    const cz = orth[6] * dx + orth[7] * dy + orth[8] * dz;
    return Math.sqrt(cx * cx + cy * cy + cz * cz);
}

/* ------------------------------------------------------------------ */
/*  Enumeration                                                        */
/* ------------------------------------------------------------------ */

/**
 * Multisets of positions whose multiplicities sum to `target`.
 *
 * Indices are non-decreasing, which enumerates each multiset once instead of
 * once per ordering. A position with no free parameter is a single point in
 * the cell, so it can be used at most once; one with parameters can hold
 * several independent atoms at different x (PbSO4 has two distinct oxygens on
 * 4c), capped by maxRepeat to keep the space finite.
 */
function multisetsSummingTo(positions, target, opts) {
    const out = [];
    const counts = new Int32Array(positions.length);
    const chosen = [];

    const rec = (start, remaining) => {
        if (remaining === 0) { out.push(chosen.slice()); return; }
        if (chosen.length >= opts.maxSites) return;
        if (out.length >= opts.maxPerElement) return;
        for (let i = start; i < positions.length; i++) {
            const m = positions[i].multiplicity;
            if (!m || m > remaining) continue;
            const cap = wyckoffFreedom(positions[i]) > 0 ? opts.maxRepeat : 1;
            if (counts[i] >= cap) continue;
            counts[i]++; chosen.push(i);
            rec(i, remaining - m);
            chosen.pop(); counts[i]--;
            if (out.length >= opts.maxPerElement) return;
        }
    };
    rec(0, target);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Site symmetry vs. coordination geometry                            */
/* ------------------------------------------------------------------ */

/**
 * The operators of the space group that leave a Wyckoff position fixed.
 *
 * Derived from the operators themselves rather than parsed out of the
 * site_symmetry string, because the string is a Hermann-Mauguin label whose
 * axis order depends on the setting and which no two databases punctuate the
 * same way. A representative point of the position is obtained by projecting
 * a deliberately generic point - one whose coordinates are not simple
 * fractions, so it cannot land on a higher-symmetry sub-position by accident -
 * and the stabiliser is then just the operators that map it back onto itself
 * modulo a lattice translation.
 */
const _SS_GENERIC = [0.1357913, 0.2468135, 0.3579246];

function siteStabilizer(w, symOps) {
    // Cached on the position object, keyed on the IDENTITY of the operator
    // list rather than on its length. Two settings of the same group have the
    // same number of operators and different operators, so a length key would
    // hand the second setting the first one's stabiliser - and the stabiliser
    // is what decides whether a site can host a tetrahedron.
    if (w._stab && w._stabOps === symOps) return w._stab;
    const p = projectOnto(w, _SS_GENERIC[0], _SS_GENERIC[1], _SS_GENERIC[2],
                          new Float64Array(3));
    const near = (a, b) => {
        const d = a - b;
        return Math.abs(d - Math.round(d)) < 1e-4;
    };
    const stab = [];
    for (const op of (symOps || [])) {
        const r = op.r, t = op.t || [0, 0, 0];
        if (!r || r.length < 9) continue;
        const q = [
            r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + t[0],
            r[3] * p[0] + r[4] * p[1] + r[5] * p[2] + t[1],
            r[6] * p[0] + r[7] * p[1] + r[8] * p[2] + t[2]
        ];
        if (near(q[0], p[0]) && near(q[1], p[1]) && near(q[2], p[2])) stab.push(op);
    }
    w._stab = stab; w._stabOps = symOps;
    return stab;
}

/** Does the site's own symmetry group contain an inversion centre? */
function siteHasInversion(w, symOps) {
    return siteStabilizer(w, symOps).some(op => {
        const r = op.r;
        for (let i = 0; i < 9; i++) {
            const want = (i % 4 === 0) ? -1 : 0;   // -I, row-major
            if (Math.abs(r[i] - want) > 1e-6) return false;
        }
        return true;
    });
}

/**
 * Coordination polyhedra that possess no inversion centre.
 *
 * An atom whose neighbours form one of these cannot occupy a site that has
 * one, because the site's inversion would have to map the coordination shell
 * onto itself - which for a tetrahedron means pairing each vertex with an
 * antipode it does not have.
 */
const NON_CENTROSYMMETRIC_GEOMETRY = new Set([
    'tetrahedral', 'trigonal-planar', 'trigonal-pyramidal',
    'square-pyramidal', 'trigonal-prismatic', 'see-saw', 'bent'
]);

/**
 * Central atoms for which a coordination number of 4 means a TETRAHEDRON.
 *
 * Deliberately a list of the oxoanion and framework formers rather than
 * "anything with four neighbours": four-coordinate d8 cations (Cu(II), Ni(II),
 * Pd, Pt, Au) are commonly SQUARE PLANAR, which is centrosymmetric and
 * perfectly happy on an inversion centre. Guessing tetrahedral for those would
 * throw away correct assignments, which is a worse failure than admitting a
 * few wrong ones. A caller that knows better can always state `geometry`
 * explicitly on the distance constraint.
 */
const TETRAHEDRAL_AT_FOUR = new Set([
    'B', 'C', 'N', 'O', 'F', 'AL', 'SI', 'P', 'S', 'CL',
    'GA', 'GE', 'AS', 'SE', 'BR', 'V', 'CR', 'MN', 'FE', 'CO', 'ZN',
    'MO', 'W', 'TC', 'RE', 'BE', 'I', 'XE'
]);

/**
 * Can an atom of `element` needing `req` coordination sit on position `w`?
 *
 * Two independent tests, both necessary conditions rather than guesses:
 *
 *   PARITY. If the site has an inversion centre, the only point it fixes is
 *   the centre itself, so every neighbour lies in an orbit of exactly two.
 *   The coordination number must therefore be even. This is exact - it
 *   follows from the group action alone and needs no chemistry.
 *
 *   GEOMETRY. A declared (or, for the classic tetrahedral formers, inferred)
 *   non-centrosymmetric polyhedron cannot survive an inversion centre. This
 *   is the test that catches sulfur on Pnma 4a: 4a has site symmetry -1, and
 *   a sulfate tetrahedron has no inversion, so the assignment is impossible
 *   however well it happens to score.
 *
 * Returns null when compatible, or a short reason string when not.
 */
function coordinationCompatible(w, element, req, symOps) {
    if (!req || !Number.isFinite(req.count) || req.count <= 0) return null;
    if (!symOps || !symOps.length) return null;
    if (!siteHasInversion(w, symOps)) return null;

    const label = `${w.multiplicity}${w.letter}`;
    if (req.count % 2 === 1) {
        return `${element} needs ${req.count} neighbours, an odd number, but ${label} ` +
               `has an inversion centre, which pairs every neighbour with an antipode.`;
    }

    let geom = req.geometry ? String(req.geometry).toLowerCase() : null;
    if (!geom && req.count === 4 &&
        TETRAHEDRAL_AT_FOUR.has(String(element).toUpperCase())) {
        geom = 'tetrahedral';
    }
    if (geom && NON_CENTROSYMMETRIC_GEOMETRY.has(geom)) {
        return `${element} on ${label} would put a ${geom} coordination on an inversion ` +
               `centre, which no ${geom} arrangement has.`;
    }
    return null;
}

/**
 * All Wyckoff assignments consistent with a composition.
 *
 * @param {Array}  wyckoff  the setting's wyckoff array from the JSON
 * @param {Array}  demand   [{ element, z, r, count }] - count is atoms PER
 *                          CELL, i.e. Z x (atoms per formula unit)
 * @param {Object} options  maxSites, maxRepeat, maxPerElement, maxTotal
 * @returns {Array} [{ sites: [{ element, elementIdx, z, r, w, wIdx }], nSites }]
 *
 * The only cross-element constraint is that a zero-parameter position is one
 * point in the cell and cannot be shared: if Pb takes 4a, no other element
 * may. Positions with parameters are freely reusable.
 */
function enumerateAssignments(wyckoff, demand, options = {}) {
    // Caps that cannot possibly be met are not caps, they are a wall.
    //
    // Both defaults used to be a flat 4, and 4 is generous in a high-symmetry
    // group and impossible in a low one. PbSO4 as Pb4S4O16 in P-1 needs eight
    // oxygens on 2i - eight distinct sites, eight reuses of one position - and
    // the enumeration returned nothing with a message blaming Z, which was the
    // one thing that was right. The floor below is what the composition
    // arithmetically requires; the caller's value is honoured when it is
    // larger, and `ceiling` is what stops a genuinely absurd request.
    const ceiling = options.ceiling ?? 24;
    const mults = wyckoff.map(w => w.multiplicity).filter(m => m > 0);
    const mMax = mults.length ? Math.max(...mults) : 1;
    const maxCount = demand.reduce((n, d) => Math.max(n, d.count), 0);
    const needed = Math.min(ceiling, Math.ceil(maxCount / mMax));

    const opts = {
        ...options,
        // distinct sites one element may occupy
        maxSites: Math.min(ceiling, Math.max(options.maxSites ?? 4, needed)),
        // times one position may be reused
        maxRepeat: Math.min(ceiling, Math.max(options.maxRepeat ?? 4, needed)),
        maxPerElement: options.maxPerElement ?? 500,
        maxTotal: options.maxTotal ?? 4000
    };

    // Sites in ONE assignment, across every element. maxSites above is the
    // per-element cap and says nothing about the total: PbSO4 in P1 needs
    // 4 + 4 + 16 = 24 sites with a per-element cap of 16, which is inside
    // every cap here and outside what the kernel can hold.
    //
    // Filtered during enumeration rather than rejected afterwards, because an
    // over-wide assignment is one candidate among many and dropping it should
    // cost that candidate, not the run. How many were dropped is returned so
    // the caller can say so instead of the list quietly being shorter.
    const maxTotalSites = Number.isFinite(options.maxTotalSites)
        ? options.maxTotalSites : WY_MAX_SITES;

    // Largest multiplicity first: the recursion then hits complete solutions
    // early and the maxPerElement cap truncates the tail of many-small-sites
    // combinations rather than the head of plausible ones.
    const idx = wyckoff.map((w, i) => i)
                       .sort((a, b) => wyckoff[b].multiplicity - wyckoff[a].multiplicity);
    const sorted = idx.map(i => wyckoff[i]);

    // Rule out, per element, the positions whose SITE SYMMETRY cannot host
    // that element's coordination - before the combinatorics rather than
    // after, so an impossible assignment never reaches the GPU at all.
    //
    // This is what stops sulfur being offered Pnma 4a. The composition
    // arithmetic is perfectly happy with S on 4a (four sulfurs, exactly what
    // PbSO4 needs) and so is the correlation, which knows nothing about
    // tetrahedra; only the site symmetry says no, and it says so for free.
    const coordination = options.coordination || null;
    const symOps = options.symOps || null;
    const excluded = [];
    const allowedFor = demand.map(d => {
        const req = coordination ? coordination[String(d.element).toUpperCase()] : null;
        if (!req || !symOps) return null;                 // nothing to test against
        const ok = new Set();
        for (let i = 0; i < wyckoff.length; i++) {
            const why = coordinationCompatible(wyckoff[i], d.element, req, symOps);
            if (why) excluded.push(why); else ok.add(i);
        }
        return ok;
    });

    const perElement = demand.map((d, e) => {
        const ok = allowedFor[e];
        // multisetsSummingTo works on the sorted array, so filter there and
        // map back through idx exactly as before.
        const pool = ok ? sorted.filter((w, k) => ok.has(idx[k])) : sorted;
        const poolIdx = ok ? idx.filter((_, k) => ok.has(idx[k])) : idx;
        const sets = multisetsSummingTo(pool, d.count, opts);
        return sets.map(s => s.map(k => poolIdx[k]));   // back to original indices
    });

    for (let e = 0; e < demand.length; e++) {
        if (perElement[e].length === 0) {
            // Say WHICH constraint failed. "Check Z" was the old advice and it
            // was usually wrong: the multiplicities available almost always can
            // sum to the demand, and what actually stopped the recursion was a
            // cap on how many sites were allowed to do it.
            const avail = [...new Set(mults)].sort((a, b) => a - b);
            const g = avail.reduce((x, y) => { while (y) { [x, y] = [y, x % y]; } return x; }, 0);
            const reachable = demand[e].count % g === 0;
            const minSites = Math.ceil(demand[e].count / mMax);
            return { assignments: [], error:
                `No Wyckoff assignment gives ${demand[e].count} ${demand[e].element} atoms ` +
                `per cell. ` +
                (reachable
                    ? `The multiplicities available (${avail.join(', ')}) can reach that total, ` +
                      `but only with at least ${minSites} sites, and the search is capped at ` +
                      `${opts.maxSites} distinct site(s) per element and ${opts.maxRepeat} reuse(s) ` +
                      `of any one position. Raise the caps, or check Z - the formula is multiplied ` +
                      `by Z, so writing the cell contents as the formula AND setting Z above 1 ` +
                      `counts every atom twice.`
                    : `The multiplicities available (${avail.join(', ')}) share a common factor of ` +
                      `${g}, so no combination of them can sum to ${demand[e].count}. Check Z.`) };
        }
    }

    const fixedOnly = new Set();
    wyckoff.forEach((w, i) => { if (wyckoffFreedom(w) === 0) fixedOnly.add(i); });

    const assignments = [];
    let droppedTooManySites = 0;
    const cur = new Array(demand.length);
    const rec = (e, usedFixed, used) => {
        if (assignments.length >= opts.maxTotal) return;
        if (used > maxTotalSites) { droppedTooManySites++; return; }
        if (e === demand.length) {
            const sites = [];
            for (let k = 0; k < demand.length; k++) {
                for (const wIdx of cur[k]) {
                    sites.push({
                        element: demand[k].element, elementIdx: k,
                        z: demand[k].z, r: demand[k].r,
                        wIdx, w: wyckoff[wIdx]
                    });
                }
            }
            assignments.push({ sites, nSites: sites.length });
            return;
        }
        for (const combo of perElement[e]) {
            // A fixed point can host one atom, full stop.
            let clash = false;
            const added = [];
            for (const wIdx of combo) {
                if (!fixedOnly.has(wIdx)) continue;
                if (usedFixed.has(wIdx)) { clash = true; break; }
                usedFixed.add(wIdx); added.push(wIdx);
            }
            if (!clash) { cur[e] = combo; rec(e + 1, usedFixed, used + combo.length); }
            for (const wIdx of added) usedFixed.delete(wIdx);
            if (assignments.length >= opts.maxTotal) return;
        }
    };
    rec(0, new Set(), 0);

    if (!assignments.length && droppedTooManySites) {
        return { assignments: [], siteSymmetryExcluded: excluded,
                 droppedTooManySites, error:
            `Every Wyckoff assignment for this composition needs more than ` +
            `${maxTotalSites} independent sites, which is the most the search ` +
            `kernel holds. This normally means Z is too large for the group: ` +
            `the multiplicities available cannot cover ` +
            `${demand.map(d => d.element + d.count).join(' ')} in fewer sites. ` +
            `Lower Z, or choose a higher-symmetry setting.` };
    }

    return { assignments, truncated: assignments.length >= opts.maxTotal,
             // How many valid assignments were skipped for exceeding the
             // kernel's site limit. Reported rather than silent: an
             // enumeration that returns 19 of 35 candidates and says nothing
             // is indistinguishable from one that found 19.
             droppedTooManySites,
             // The site-symmetry rejections, so the caller can say what was
             // ruled out and why. Silence here would be indistinguishable
             // from the filter not running.
             siteSymmetryExcluded: excluded };
}

/* ------------------------------------------------------------------ */
/*  Ranking                                                            */
/* ------------------------------------------------------------------ */

/**
 * Orders assignments best-first on evidence available before any search.
 *
 * Two terms:
 *
 *   Harker residual. Every consolidated Harker site is a candidate position
 *   for a heavy atom. For each heavy site in the assignment, the distance to
 *   the nearest candidate that also lies ON the assigned position is a direct
 *   measure of whether that assignment is compatible with the Patterson we
 *   already computed. This is far more informative than any chemical rule of
 *   thumb and it is free - no GPU work at all.
 *
 *   Parsimony. Among assignments the data cannot distinguish, prefer fewer
 *   distinct sites. This is a tiebreaker, weighted low enough that it never
 *   overrides evidence.
 */
function rankAssignments(assignments, ctx) {
    const { orth, harkerSites = [], heavyZ = 0, parsimony = 0.15 } = ctx;

    for (const A of assignments) {
        let evidence = 0, n = 0;
        for (const s of A.sites) {
            if (s.z < heavyZ) continue;
            let best = Infinity;
            for (const h of harkerSites) {
                const d = wyckoffShift(s.w, h, orth);
                if (d < best) best = d;
            }
            if (best !== Infinity) { evidence += Math.min(best, 3.0); n++; }
        }
        A.harkerResidual = n ? evidence / n : null;
        A.priorScore = (n ? evidence / n : 1.5) + parsimony * A.nSites;
    }
    assignments.sort((a, b) => a.priorScore - b.priorScore);
    return assignments;
}

/* ------------------------------------------------------------------ */
/*  Flattening for the GPU                                             */
/* ------------------------------------------------------------------ */

/**
 * Packs every assignment into the flat buffers swarm_reflection.wgsl expects.
 *
 * The crucial invariant: the total number of atoms generated in the cell is
 * the same for every assignment, because the composition is fixed. That makes
 * the generated-atom arrays a constant size, removes the need for the shader's
 * distance-based special-position collapse entirely, and - because the pair
 * weight sum is then identical everywhere - makes raw fitness directly
 * comparable BETWEEN assignments. Comparing them is the whole point.
 *
 * Per generated atom g of assignment A, at index A*nTot + g:
 *   genSite  which site of A it belongs to (indexes the particle's coords)
 *   genOp    which symmetry operator produced it
 *   genZ     atomic number
 *   genType  element index, for the distance restraints
 *
 * genOp comes from the position's precomputed coset operators, so exactly
 * `multiplicity` distinct atoms are produced with no tolerance to tune.
 */
function buildAssignmentTables(assignments, symOpCount, demand) {
    const nTot = demand.reduce((n, d) => n + d.count, 0);
    const nA = assignments.length;
    const maxSites = assignments.reduce((m, A) => Math.max(m, A.nSites), 0);

    const genSite = new Uint32Array(nA * nTot);
    const genOp = new Uint32Array(nA * nTot);
    const genZ = new Float32Array(nA * nTot);
    const genType = new Uint32Array(nA * nTot);
    const siteProj = new Float32Array(nA * maxSites * 12);   // P(9) + T(3)
    const siteCount = new Uint32Array(nA);
    const warnings = [];

    assignments.forEach((A, ai) => {
        let g = 0;
        A.sites.forEach((s, si) => {
            const cosets = Array.isArray(s.w.coset_ops) && s.w.coset_ops.length
                ? s.w.coset_ops
                : null;
            // THE FALLBACK IS ONLY VALID FOR A GENERAL POSITION.
            //
            // `k % symOpCount` for k = 0..multiplicity-1 enumerates operators
            // 0, 1, 2, ... in database order. For a general position that IS
            // the whole group and the orbit is right. For a SPECIAL position
            // the multiplicity is smaller than the group order, so this takes
            // the first few operators as listed - which are not a transversal
            // of the stabiliser and do not generate the position's orbit. The
            // atoms come out in the wrong places, several of them coincident,
            // and the only trace was a warning in a log.
            //
            // A structure that cannot be represented is a failure, not a
            // caveat - the same rule the genPack width check below follows.
            if (!cosets && s.w.multiplicity !== symOpCount) {
                throw new Error(
                    `Wyckoff position ${s.w.multiplicity}${s.w.letter} has no coset_ops in the ` +
                    `space-group database, and it is a special position (multiplicity ` +
                    `${s.w.multiplicity} against a group order of ${symOpCount}), so its orbit ` +
                    `cannot be reconstructed from the operator order. Regenerate the database ` +
                    `with cctbx_Harko_v1.py.`);
            }
            if (!cosets) {
                warnings.push(`Wyckoff ${s.w.multiplicity}${s.w.letter} has no coset_ops in the ` +
                              `database; using the full operator list, which is correct for a ` +
                              `general position only. Re-run the generator with the Wyckoff patch.`);
            }
            const ops = cosets || Array.from({ length: s.w.multiplicity }, (_, k) => k % symOpCount);
            for (const op of ops) {
                if (g >= nTot) break;
                genSite[ai * nTot + g] = si;
                genOp[ai * nTot + g] = op;
                genZ[ai * nTot + g] = s.z;
                genType[ai * nTot + g] = s.elementIdx;
                g++;
            }
            const { P, T } = wyckoffProjector(s.w);
            const b = (ai * maxSites + si) * 12;
            for (let k = 0; k < 9; k++) siteProj[b + k] = P[k];
            for (let k = 0; k < 3; k++) siteProj[b + 9 + k] = T[k];
        });
        // Unused site slots project to the origin so a stale coordinate can
        // never leak into a shorter assignment's atom list.
        for (let si = A.nSites; si < maxSites; si++) {
            const b = (ai * maxSites + si) * 12;
            for (let k = 0; k < 12; k++) siteProj[b + k] = 0;
        }
        siteCount[ai] = A.nSites;
        if (g !== nTot) {
            warnings.push(`Assignment ${ai} generated ${g} atoms, expected ${nTot}.`);
        }
    });

    // Constant pair weight: (sum z)^2 - sum z^2, halved. Identical for every
    // assignment because the cell contents are.
    let sz = 0, sz2 = 0;
    for (const d of demand) { sz += d.z * d.count; sz2 += d.z * d.z * d.count; }
    const normW = Math.max(1e-6, 0.5 * (sz * sz - sz2));

    // CAPACITIES FIRST, PACKING SECOND.
    //
    // The bit layout (site | op << 8 | type << 20) has room for 256 sites and
    // 16 elements, and the field widths are NOT the binding constraint - the
    // kernel's workgroup arrays are. prop_sites holds WY_MAX_SITES entries and
    // local_rMin holds WY_MAX_ELEM^2, and exceeding either is not a truncation
    // the packing can detect: the values fit the bits perfectly and go on to
    // index past the end of a workgroup array, which WGSL clamps in silence.
    //
    // So the limits checked here are the kernel's, not the packing's, and they
    // are checked ONCE up front rather than per atom - the old loop tested
    // 256/4096/16 on every one of nA * nTot entries to catch a condition that
    // is a property of the tables as a whole.
    if (maxSites > WY_MAX_SITES) {
        throw new Error(
            `The widest assignment needs ${maxSites} independent Wyckoff sites; the search ` +
            `kernel holds ${WY_MAX_SITES} (the width of prop_sites). Lower Z, lower the ` +
            `per-element site cap, or raise MAX_SITES in swarm_reflection.wgsl ` +
            `together with the prop_sites array.`);
    }
    if (demand.length > WY_MAX_ELEM) {
        throw new Error(
            `${demand.length} distinct elements; the search kernel holds ${WY_MAX_ELEM} ` +
            `(local_rMin is MAX_ELEM^2). Simplify the composition, or raise MAX_ELEM in ` +
            `swarm_reflection.wgsl together with local_rMin.`);
    }
    if (symOpCount > WY_MAX_OPS) {
        throw new Error(
            `${symOpCount} symmetry operators; genOp is packed 12 bits wide, so ${WY_MAX_OPS} ` +
            `is the ceiling. No crystallographic group reaches this.`);
    }

    // Bit-packed for the kernel: site | op << 8 | type << 20. WebGPU guarantees
    // only eight storage buffers per stage and this kernel uses all eight, so
    // three per-atom u32 fields have to travel as one.
    const genPack = new Uint32Array(nA * nTot);
    for (let i = 0; i < genPack.length; i++) {
        genPack[i] = genSite[i] | (genOp[i] << 8) | (genType[i] << 20);
    }

    return { nTot, maxSites, genSite, genOp, genZ, genType, genPack,
             siteProj, siteCount, normW, warnings: [...new Set(warnings)] };
}

/**
 * The single f32 lookup buffer: [ zByType | rMin matrix | bond rules ].
 *
 * Concatenated for the same binding-budget reason as genPack. Offsets go in
 * the uniform block.
 *
 * rMin is a full nElem x nElem lower-bound matrix. Every pair starts at the
 * minimum-contact floor (`options.minContact`, a hard floor on every pair);
 * `windows` then overrides individual pairs:
 *
 *     [{ a: 'S', b: 'O', dmin: 1.35, dmax: 1.65 }]
 *
 * dmin is symmetric and goes into the matrix. dmax becomes a directed
 * nearest-neighbour rule, and is emitted BOTH ways round unless the caller
 * asks otherwise: "every S has an O within 1.65 A" and "every O has an S
 * within 1.65 A" are different statements, and for a sulfate both are true.
 * For a rule like Pb-O where the partner counts differ wildly, one direction
 * is usually what you want.
 */
/**
 * (elementA, elementB) -> the clash floor the search used, in Angstrom.
 */
function floorLookup(demand, rMin, nElem) {
    const idx = {};
    demand.forEach((d, i) => { idx[String(d.element).toUpperCase()] = i; });
    return (a, b) => {
        const i = idx[String(a).toUpperCase()], j = idx[String(b).toUpperCase()];
        return (i === undefined || j === undefined) ? 0 : rMin[i * nElem + j];
    };
}

function buildRestraintTables(demand, windows = [], options = {}) {
    const nElem = demand.length;
    const idxOf = {};
    demand.forEach((d, i) => { idxOf[d.element.toUpperCase()] = i; });

    // THE CLASH FLOOR, IN ONE RULE:
    //
    //     rMin(A,B) = the dmin the user stated for that pair,   if they stated one
    //               = the minimum-contact slider,               otherwise
    //
    // Nothing else feeds it. dmax does NOT - an upper bound on where partners
    // may sit says nothing about how close atoms may approach, and letting it
    // move the floor made the two ends of a constraint interact in a way no one
    // could follow from the input.
    //
    // WHY THE RADIUS RULE IS GONE. The floor used to be
    // min(slider, 0.65 x (rc_A + rc_B)), which meant the slider could only ever
    // LOOSEN a pair: raised above a pair's radius sum it did nothing at all.
    // On PbSO4 a slider set to 1.10 A still returned candidates with O-O at
    // 0.94, because the O-O radius floor is 0.86 and min() kept it. That is not
    // a defensible reading of a control called "minimum contact distance", and
    // no amount of documentation fixes a control whose name is the opposite of
    // its behaviour. One number, applied to every pair, is also the only
    // version a user can check: the Min d column either respects it or the code
    // is wrong.
    //
    // The old rule's motive was real and is preserved in the DEFAULT rather
    // than in the arithmetic. A floor is a guess and the two ways of being
    // wrong are not symmetric: too low merely fails to catch nonsense the
    // correlation would have rejected anyway, while too high makes the right
    // answer UNREACHABLE - and charges it more than any wrong one, since the
    // clash term scales with overlap depth. So the default is deliberately
    // permissive, and raising it is the user's decision to make.
    //
    // 1.0 A: below every bond that actually occurs (the shortest in ordinary
    // inorganics is a little over 1.1 A, B-O and C-O), and far above what a
    // wrong Wyckoff assignment produces - a Pb-O contact of 0.32 A correlated
    // at 0.9428 and ranked fourth in a real run when there was no floor at all.
    // An explicit 0 turns the floor off; only an ABSENT option falls back to
    // the default, so a caller that forgets to pass one is not silently left
    // with no clash term.
    const DEFAULT_CONTACT_FLOOR_A = 1.0;
    const floor = (options.minContact === undefined || options.minContact === null)
        ? DEFAULT_CONTACT_FLOOR_A
        : (Number.isFinite(options.minContact) && options.minContact > 0 ? options.minContact : 0);
    const zByType = demand.map(d => d.z);
    const rMin = new Float64Array(nElem * nElem);
    rMin.fill(floor);

    // Six floats per rule: [aType, bType, dmin, dmax, count, mode].
    //   mode 0 = no count requirement (nearest-neighbour rule only)
    //        1 = exactly `count` neighbours in the window
    //        2 = at least `count`
    // dmax is +inf when the constraint left the upper side open.
    const RULE_STRIDE = 6;
    const MODE = { none: 0, exact: 1, atleast: 2 };
    const problems = [];
    const rules = [];
    for (const w of windows) {
        const a = idxOf[String(w.a).toUpperCase()];
        const b = idxOf[String(w.b).toUpperCase()];
        if (a === undefined || b === undefined) {
            problems.push(`Constraint ${w.a}-${w.b} names an element that is not in the composition.`);
            continue;
        }
        const dmin = Number.isFinite(w.dmin) ? w.dmin : null;
        const dmax = Number.isFinite(w.dmax) ? w.dmax : null;
        if (dmin !== null && dmax !== null && dmax < dmin) {
            problems.push(`Constraint ${w.a}-${w.b} has dmax (${dmax}) below dmin (${dmin}).`);
            continue;
        }

        // dmin is a HARD FLOOR, always, whether or not a count was asked for.
        // It goes into the symmetric rMin matrix, which is what the kernel
        // charges as a clash - so "S O - 1.4/" is expressed entirely here and
        // needs no rule at all.
        // A stated dmin REPLACES the global floor for that pair - not max(),
        // not min(). Both numbers come from the user, and the more specific one
        // wins, in either direction: a constraint may hold one pair further
        // apart than the slider, or admit an approach the slider would forbid.
        // Anything else would make a constraint mean different things depending
        // on where the slider happened to be sitting.
        if (dmin !== null) {
            rMin[a * nElem + b] = dmin;
            rMin[b * nElem + a] = dmin;
        }

        const hasCount = Number.isFinite(w.count) && w.count >= 0;
        const mode = hasCount ? (w.countMode === 'atleast' ? MODE.atleast : MODE.exact) : MODE.none;

        // A rule is only needed when something beyond the floor is being asked:
        // a coordination number, or an upper distance.
        if (!hasCount && dmax === null) continue;

        const rec = [a, b, dmin ?? 0, dmax ?? Infinity, hasCount ? w.count : 0, mode];
        rules.push(rec);
        // A coordination number is DIRECTED - "four O around every S" says
        // nothing about how many S surround an O - so it is never mirrored.
        // A bare upper bound is a nearest-neighbour condition and is symmetric.
        if (!hasCount && w.bothWays !== false && a !== b) {
            rules.push([b, a, dmin ?? 0, dmax ?? Infinity, 0, MODE.none]);
        }
    }
    // BOTH OF THESE THROW. They used to be advisory strings on `problems`,
    // which the caller passes to say() - a line in a log - while the return
    // value quietly truncated `nRules` to 8 and left the slot overflow to the
    // kernel, where it wraps onto another rule's distances.
    //
    // A restraint the user typed and the search did not apply is the worst
    // possible failure mode here: it costs a full run, the answer looks
    // reasonable, and the structure violates a constraint the report will
    // then correctly say is violated - which reads as a bad structure rather
    // than a dropped rule.
    if (rules.length > WY_MAX_RULES) {
        throw new Error(
            `${rules.length} distance rules were generated from ${windows.length} constraint(s); ` +
            `the search kernel holds ${WY_MAX_RULES}. Note that a bare upper bound is emitted ` +
            `BOTH ways round unless bothWays is false, so it costs two rules. Drop a constraint, ` +
            `or set bothWays: false on a directed one.`);
    }
    // The kernel keeps the N nearest partners of every counted rule so the
    // coordination penalty can be charged by distance rather than by count.
    // Those slots are a fixed private array shared across all rules.
    const slotsNeeded = rules.reduce((n, r) => n + (r[5] !== MODE.none ? r[4] : 0), 0);
    if (slotsNeeded > WY_MAX_COORD_SLOTS) {
        throw new Error(
            `The coordination numbers add up to ${slotsNeeded} neighbour slots; the search ` +
            `kernel tracks ${WY_MAX_COORD_SLOTS} in total across every rule. Drop a ` +
            `constraint, or lower a coordination number.`);
    }

    const rMinOff = zByType.length;
    const ruleOff = rMinOff + rMin.length;
    const tables = new Float32Array(ruleOff + rules.length * RULE_STRIDE);
    tables.set(zByType, 0);
    tables.set(rMin, rMinOff);
    // Float32Array turns Infinity into Infinity, which compares correctly in
    // WGSL; no sentinel needed.
    rules.forEach((r, k) => tables.set(r, ruleOff + k * RULE_STRIDE));

    // The longest distance any comparison in the kernel depends on. The
    // minimum-image convention only has to be exact out to here, so this is
    // what decides whether the reduced basis alone is enough or the kernel must
    // also search the neighbouring translations. See sharkoReducedCell().
    let maxDistance = 0;
    for (let i = 0; i < rMin.length; i++) maxDistance = Math.max(maxDistance, rMin[i]);
    for (const r of rules) if (Number.isFinite(r[3])) maxDistance = Math.max(maxDistance, r[3]);

    return { tables, rMinOff, ruleOff, nRules: rules.length, nElem,
             ruleStride: RULE_STRIDE, maxDistance, problems,
             // The resolved floors, so the post-search filter can reject on
             // exactly what the search enforced. They used to disagree: the
             // search applied per-pair floors while the filter applied a flat
             // 0.9 A to everything, which let the search spend a run on
             // contacts the filter then threw out. Anything reading this gets
             // the per-pair answer including any constraint override, not the
             // slider value, which is why it stays a lookup and not a number.
             floors: floorLookup(demand, rMin, nElem) };
}

/** Symmetry operators as 12 floats each: r[9] then t[3]. */
function packSymOps(symOps) {
    const out = new Float32Array(symOps.length * 12);
    symOps.forEach((op, i) => { out.set(op.r, i * 12); out.set(op.t, i * 12 + 9); });
    return out;
}

/**
 * Projects one particle's coordinates onto its assignment's subspaces.
 *
 * Call this after initialising a chain and on EVERY proposal. It is
 * idempotent, so the cost is one 3x3 multiply per site and calling it too
 * often is harmless.
 *
 * Under the Monte Carlo sampler this is not a correction for drift - it is the
 * only thing keeping a structure on its subspace. A particle swarm stayed there
 * largely for free, because its update is a linear combination of points
 * already in the affine subspace; a Gaussian proposal is not, and lands off the
 * subspace almost surely. A move is therefore a step in the ambient cell
 * followed by a projection, which is still a symmetric proposal on the subspace
 * and leaves the Metropolis acceptance rule valid.
 */
function projectSites(positions, particleIdx, assignIdx, tables) {
    const { maxSites, siteProj, siteCount } = tables;
    const nS = siteCount[assignIdx];
    const stride = maxSites * 3;
    for (let s = 0; s < nS; s++) {
        const b = (assignIdx * maxSites + s) * 12;
        const o = particleIdx * stride + s * 3;
        const x = positions[o], y = positions[o + 1], z = positions[o + 2];
        positions[o]     = siteProj[b]     * x + siteProj[b + 1] * y + siteProj[b + 2] * z + siteProj[b + 9];
        positions[o + 1] = siteProj[b + 3] * x + siteProj[b + 4] * y + siteProj[b + 5] * z + siteProj[b + 10];
        positions[o + 2] = siteProj[b + 6] * x + siteProj[b + 7] * y + siteProj[b + 8] * z + siteProj[b + 11];
        for (let k = 0; k < 3; k++) {
            const i = o + k;
            positions[i] -= Math.floor(positions[i]);
        }
    }
    // Slots this assignment does not use are parked at the origin so they can
    // never be read as coordinates.
    for (let s = nS; s < maxSites; s++) {
        const o = particleIdx * stride + s * 3;
        positions[o] = positions[o + 1] = positions[o + 2] = 0;
    }
}

/**
 * Shares a particle budget over a set of assignments, weighted by dimension.
 *
 * Split out because it is needed twice: once when a wave starts, and again
 * every time successive halving retires the losers and the survivors inherit
 * their particles. Redistributing uniformly at that second point would undo
 * the weighting exactly when the budget per assignment finally becomes large
 * enough for it to matter.
 */
function weightedAssign(ids, numParticles, minPer, weights) {
    const out = new Uint32Array(numParticles);
    if (!ids.length) return out;

    const floorEach = Math.min(minPer, Math.floor(numParticles / ids.length));
    const w = ids.map(a => Math.max(1e-6, (weights && weights[a]) || 1));
    const total = w.reduce((x, y) => x + y, 0);
    const surplus = Math.max(0, numParticles - floorEach * ids.length);
    const counts = w.map(x => floorEach + Math.floor(surplus * x / total));

    // Rounding leaves a few over; they go to the largest problems.
    const order = ids.map((_, k) => k).sort((a, b) => w[b] - w[a]);
    let assigned = counts.reduce((x, y) => x + y, 0), oi = 0;
    while (assigned < numParticles) { counts[order[oi++ % order.length]]++; assigned++; }
    // The bounded form. The arithmetic above cannot overshoot today
    // (floorEach * n + surplus == numParticles exactly), but the loop's exit
    // depends on some count being greater than 1 - so if every count were 1 and
    // the sum were still too large, this spun forever and the tab froze with no
    // error. `guard` makes that a wrong allocation instead of a hang.
    let guard = order.length * 4 + 8;
    while (assigned > numParticles && guard-- > 0) {
        const k = order[order.length - 1 - (oi++ % order.length)];
        if (counts[k] > 1) { counts[k]--; assigned--; }
    }

    let p = 0;
    ids.forEach((a, k) => { for (let n = 0; n < counts[k] && p < numParticles; n++) out[p++] = a; });
    while (p < numParticles) out[p++] = ids[p % ids.length];
    return out;
}

/**
 * Splits the particle budget across assignments, best-first.
 *
 * Every assignment is a separate optimisation problem with its own global
 * best - particles on different assignments must not attract each other, or
 * the swarm tears itself between incompatible structures. So this is N
 * independent sub-swarms sharing one dispatch.
 *
 * `weights` should be the free-parameter count. An equal split quietly favours
 * small models: an assignment with no free parameters is solved by its first
 * particle, while an eleven-parameter one is still searching when the run ends
 * and reports a correlation below what it could reach. The ranking then
 * measures how hard each assignment was to search rather than how well it fits.
 *
 * A wave holds only as many assignments as leaves room to weight. Packing in
 * `numParticles / minPer` of them spends the entire budget on the floor - with
 * 512 particles, a floor of 24 and 21 assignments, eight particles remain to
 * distribute and the weighting does nothing. Halving that keeps a real surplus,
 * at the cost of more waves.
 */
function allocateParticles(nAssign, numParticles, minPer = 24, weights = null) {
    const headroom = weights ? 2 : 1;
    const perWave = Math.max(1, Math.floor(numParticles / (minPer * headroom)));
    const waves = [];
    for (let i = 0; i < nAssign; i += perWave) {
        const ids = [];
        for (let k = i; k < Math.min(i + perWave, nAssign); k++) ids.push(k);
        waves.push({ ids, assignOf: weights
            ? weightedAssign(ids, numParticles, minPer, weights)
            : Uint32Array.from({ length: numParticles }, (_, p) => ids[p % ids.length]) });
    }
    return waves;
}

/**
 * Drops the worst assignments and re-runs the survivors with more particles.
 *
 * Successive halving: most assignments are visibly hopeless within a few tens
 * of generations, and the search is far better spent on the few that are not.
 * `best` is the best fitness seen per assignment so far.
 */
function pruneAssignments(ids, best, keepFraction = 0.5, minKeep = 4) {
    const ranked = [...ids].sort((a, b) => (best[b] ?? -Infinity) - (best[a] ?? -Infinity));
    const keep = Math.max(minKeep, Math.ceil(ranked.length * keepFraction));
    return ranked.slice(0, Math.min(keep, ranked.length));
}

/* ------------------------------------------------------------------ */

/**
 * "PbSO4" + Z -> [{ element, count }]. Handles nested groups and counts.
 */
function parseFormula(formula, Z = 1) {
    const counts = {};
    // UNPARSEABLE CHARACTERS ARE AN ERROR, NOT A SKIP.
    //
    // The old loop did `if (!m) { i++; continue; }`, so anything the element
    // pattern did not match was silently discarded. "pbso4" - a perfectly
    // ordinary thing to type - matched nothing until the O4 and returned a
    // composition of four oxygens, which then enumerated, searched, and
    // reported a structure. Nothing anywhere said the formula had not been
    // understood.
    const junk = [];
    const parse = (s, mult) => {
        let i = 0;
        while (i < s.length) {
            if (s[i] === '(') {
                let depth = 1, j = i + 1;
                while (j < s.length && depth) { if (s[j] === '(') depth++; if (s[j] === ')') depth--; j++; }
                if (depth) { junk.push('unclosed "("'); return; }
                const inner = s.slice(i + 1, j - 1);
                const m = /^(\d+(?:\.\d+)?)/.exec(s.slice(j));
                const n = m ? parseFloat(m[1]) : 1;
                parse(inner, mult * n);
                i = j + (m ? m[1].length : 0);
                continue;
            }
            const m = /^([A-Z][a-z]?)(\d+(?:\.\d+)?)?/.exec(s.slice(i));
            if (!m) { junk.push(s[i]); i++; continue; }
            const el = m[1].toUpperCase();
            counts[el] = (counts[el] || 0) + mult * (m[2] ? parseFloat(m[2]) : 1);
            i += m[0].length;
        }
    };
    parse(String(formula || '').replace(/\s+/g, ''), Z);
    if (junk.length) {
        throw new Error(
            `Could not read the formula "${formula}": ${[...new Set(junk)].map(c => `"${c}"`).join(', ')} ` +
            `is not an element symbol. Element symbols start with a capital letter - ` +
            `"PbSO4", not "pbso4" or "PBSO4".`);
    }
    const out = Object.entries(counts).map(([element, count]) => ({ element, count }));
    if (!out.length) throw new Error(`The formula "${formula}" names no elements.`);
    return out;
}

/**
 * Plausible Z values from the cell volume.
 *
 * A RANGE, not a number. The familiar 18 A^3 per non-hydrogen atom is an
 * organic-chemistry rule; dense inorganics run far tighter - PbSO4 sits near
 * 13 A^3 per atom, so the 18 rule alone returns Z=3 for a structure whose Z is
 * 4. Spanning 11 to 22 A^3 covers oxides, sulfates and intermetallics without
 * pretending to a precision the estimate does not have.
 *
 * Returns { min, max, likely } with `likely` the midpoint, and the caller is
 * expected to offer all of them: enumerating assignments for two or three Z is
 * cheap, and a wrong Z is otherwise a silent, total failure.
 */
function suggestZ(cellVolume, atomsPerFormula, symOpCount) {
    if (!cellVolume || !atomsPerFormula) return null;
    const lo = Math.max(1, Math.floor(cellVolume / (22 * atomsPerFormula)));
    const hi = Math.max(1, Math.ceil(cellVolume / (11 * atomsPerFormula)));
    const all = [];
    for (let z = lo; z <= hi; z++) {
        // Z is very often the group order or a divisor of it, because the
        // asymmetric unit usually holds a whole number of formula units.
        const divides = !symOpCount || symOpCount % z === 0 || z % symOpCount === 0;
        all.push({ Z: z, favoured: divides });
    }
    const fav = all.filter(a => a.favoured).map(a => a.Z);
    return { min: lo, max: hi, candidates: all.map(a => a.Z), favoured: fav,
             likely: fav.length ? fav[Math.floor(fav.length / 2)]
                                : Math.round((lo + hi) / 2) };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { SHARKO_WYCKOFF_ASSIGN_VERSION,
        WY_MAX_SITES, WY_MAX_ELEM, WY_MAX_OPS, WY_MAX_RULES, WY_MAX_COORD_SLOTS,
        wyckoffFreedom, wyckoffProjector, projectOnto, wyckoffShift,
        siteStabilizer, siteHasInversion, coordinationCompatible,
        enumerateAssignments, rankAssignments, buildAssignmentTables,
        projectSites, allocateParticles, weightedAssign, pruneAssignments,
        buildRestraintTables, packSymOps,
        parseFormula, suggestZ
    };
}
