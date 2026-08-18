/* ===========================================================================
 *  wyckoff_report.js  --  structure plot, contact list and bond-valence
 *  analysis for the Wyckoff results panel.
 *
 *  Self-contained on purpose. It needs nothing from the rest of Powder 5
 *  except the structure object the Wyckoff search already produces:
 *
 *      st.cell          {a,b,c,alpha,beta,gamma}
 *      st.sites         [{rank, element, x, y, z, wyckoff, multiplicity}]
 *      st.spaceGroup    {number, symbol}
 *      st.opMatrices    [{r:[9], t:[3], xyz}]        <-- optional
 *
 *  If opMatrices is absent it falls back to getSymopsForSetting(st.spaceGroup),
 *  which is what the charge-flipping path already uses. With no operators at
 *  all the panel says so rather than drawing a P1 cell and calling it a
 *  structure -- a coordination number computed without symmetry equivalents is
 *  not a small error, it is a different number.
 *
 *  Entry point:  WyckoffReport.render(hostElement, structure)
 * =========================================================================== */
(function (global) {
'use strict';

/* ---------------------------------------------------------------------------
 *  Element data. Masses (IUPAC 2021), covalent radii (Cordero 2008, A) and
 *  CPK-ish colours. Radii are used ONLY to decide what to draw as a bond and
 *  how big to draw a sphere; nothing in the numbers below feeds the analysis.
 * ------------------------------------------------------------------------- */
const EL = {
 H :[1.008 ,0.31,'#e8e8e8'], He:[4.003 ,0.28,'#d9ffff'], Li:[6.94  ,1.28,'#cc80ff'],
 Be:[9.012 ,0.96,'#c2ff00'], B :[10.81 ,0.84,'#ffb5b5'], C :[12.011,0.76,'#4a4a4a'],
 N :[14.007,0.71,'#3050f8'], O :[15.999,0.66,'#d03020'], F :[18.998,0.57,'#90e050'],
 Ne:[20.180,0.58,'#b3e3f5'], Na:[22.990,1.66,'#ab5cf2'], Mg:[24.305,1.41,'#8aff00'],
 Al:[26.982,1.21,'#bfa6a6'], Si:[28.085,1.11,'#f0c8a0'], P :[30.974,1.07,'#ff8000'],
 S :[32.06 ,1.05,'#e6c619'], Cl:[35.45 ,1.02,'#1ff01f'], Ar:[39.95 ,1.06,'#80d1e3'],
 K :[39.098,2.03,'#8f40d4'], Ca:[40.078,1.76,'#3dff00'], Sc:[44.956,1.70,'#e6e6e6'],
 Ti:[47.867,1.60,'#bfc2c7'], V :[50.942,1.53,'#a6a6ab'], Cr:[51.996,1.39,'#8a99c7'],
 Mn:[54.938,1.39,'#9c7ac7'], Fe:[55.845,1.32,'#e06633'], Co:[58.933,1.26,'#f090a0'],
 Ni:[58.693,1.24,'#50d050'], Cu:[63.546,1.32,'#c88033'], Zn:[65.38 ,1.22,'#7d80b0'],
 Ga:[69.723,1.22,'#c28f8f'], Ge:[72.630,1.20,'#668f8f'], As:[74.922,1.19,'#bd80e3'],
 Se:[78.971,1.20,'#ffa100'], Br:[79.904,1.20,'#a62929'], Kr:[83.798,1.16,'#5cb8d1'],
 Rb:[85.468,2.20,'#702eb0'], Sr:[87.62 ,1.95,'#00ff00'], Y :[88.906,1.90,'#94ffff'],
 Zr:[91.224,1.75,'#94e0e0'], Nb:[92.906,1.64,'#73c2c9'], Mo:[95.95 ,1.54,'#54b5b5'],
 Tc:[98.0  ,1.47,'#3b9e9e'], Ru:[101.07,1.46,'#248f8f'], Rh:[102.91,1.42,'#0a7d8c'],
 Pd:[106.42,1.39,'#006985'], Ag:[107.87,1.45,'#b8b8b8'], Cd:[112.41,1.44,'#ffd98f'],
 In:[114.82,1.42,'#a67573'], Sn:[118.71,1.39,'#668080'], Sb:[121.76,1.39,'#9e63b5'],
 Te:[127.60,1.38,'#d47a00'], I :[126.90,1.39,'#940094'], Xe:[131.29,1.40,'#429eb0'],
 Cs:[132.91,2.44,'#57178f'], Ba:[137.33,2.15,'#00c900'], La:[138.91,2.07,'#70d4ff'],
 Ce:[140.12,2.04,'#ffffc7'], Pr:[140.91,2.03,'#d9ffc7'], Nd:[144.24,2.01,'#c7ffc7'],
 Pm:[145.0 ,1.99,'#a3ffc7'], Sm:[150.36,1.98,'#8fffc7'], Eu:[151.96,1.98,'#61ffc7'],
 Gd:[157.25,1.96,'#45ffc7'], Tb:[158.93,1.94,'#30ffc7'], Dy:[162.50,1.92,'#1fffc7'],
 Ho:[164.93,1.92,'#00ff9c'], Er:[167.26,1.89,'#00e675'], Tm:[168.93,1.90,'#00d452'],
 Yb:[173.05,1.87,'#00bf38'], Lu:[174.97,1.87,'#00ab24'], Hf:[178.49,1.75,'#4dc2ff'],
 Ta:[180.95,1.70,'#4da6ff'], W :[183.84,1.62,'#2194d6'], Re:[186.21,1.51,'#267dab'],
 Os:[190.23,1.44,'#266696'], Ir:[192.22,1.41,'#175487'], Pt:[195.08,1.36,'#d0d0e0'],
 Au:[196.97,1.36,'#ffd123'], Hg:[200.59,1.32,'#b8b8d0'], Tl:[204.38,1.45,'#a6544d'],
 Pb:[207.2 ,1.46,'#575961'], Bi:[208.98,1.48,'#9e4fb5'], Po:[209.0 ,1.40,'#ab5c00'],
 Th:[232.04,2.06,'#00baff'], Pa:[231.04,2.00,'#00a1ff'], U :[238.03,1.96,'#008fff'],
 Np:[237.0 ,1.90,'#0080ff'], Pu:[244.0 ,1.87,'#006bff']
};
const AVOGADRO = 6.02214076e23;

/* Elements that CAN be the anion. Which of them actually is one depends on
 * what else is in the cell: sulfur is the anion in PbS and the cation in
 * PbSO4, and a fixed list gets one of those two wrong every time. The rule
 * below is "most electronegative species present, and anything close to it". */
const ANION_CANDIDATES = new Set(['O', 'F', 'Cl', 'Br', 'I', 'S', 'Se', 'Te', 'N', 'P', 'As', 'H', 'C']);
const EN = {
 H:2.20, B:2.04, C:2.55, N:3.04, O:3.44, F:3.98, Si:1.90, P:2.19, S:2.58, Cl:3.16,
 Ge:2.01, As:2.18, Se:2.55, Br:2.96, Sb:2.05, Te:2.10, I:2.66, Li:0.98, Na:0.93,
 K:0.82, Rb:0.82, Cs:0.79, Be:1.57, Mg:1.31, Ca:1.00, Sr:0.95, Ba:0.89, Al:1.61,
 Ga:1.81, In:1.78, Tl:1.62, Sn:1.96, Pb:2.33, Bi:2.02
};
const enOf = el => EN[el] ?? 1.7;

/* Formal charge assumed for each anion. */
const ANION_CHARGE = { O:-2, S:-2, Se:-2, Te:-2, F:-1, Cl:-1, Br:-1, I:-1,
                       N:-3, P:-3, As:-3, C:-4, H:-1 };

/** Which elements of this structure act as anions. */
function anionSet(elements) {
    const cand = elements.filter(e => ANION_CANDIDATES.has(e));
    if (!cand.length) return new Set();
    const maxEN = Math.max(...elements.map(enOf));
    const out = new Set(cand.filter(e => enOf(e) >= maxEN - 0.5));
    return out.size ? out : new Set([cand.reduce((a, b) => enOf(a) >= enOf(b) ? a : b)]);
}

/* Default cation oxidation state. The site table lets the reader change it;
 * this is only the starting guess. */
const DEFAULT_OX = {
 H:1, Li:1, Na:1, K:1, Rb:1, Cs:1, Ag:1, Cu:2, Au:1, Tl:1,
 Be:2, Mg:2, Ca:2, Sr:2, Ba:2, Zn:2, Cd:2, Hg:2, Pb:2, Sn:2, Ni:2, Co:2, Mn:2, Fe:3,
 B:3, Al:3, Ga:3, In:3, Sc:3, Y:3, La:3, Ce:3, Pr:3, Nd:3, Sm:3, Eu:3, Gd:3, Tb:3,
 Dy:3, Ho:3, Er:3, Tm:3, Yb:3, Lu:3, Bi:3, Sb:3, As:5, Cr:3, V:5, Ti:4, Zr:4, Hf:4,
 Si:4, Ge:4, C:4, Th:4, U:6, Nb:5, Ta:5, Mo:6, W:6, Re:7, Ru:4, Rh:3, Pd:2, Pt:2,
 Ir:4, Os:6, P:5, S:6, Se:4, Te:4, N:5, Cl:7, Br:7, I:5
};

/* Bond-valence parameters R0 (A), Brese & O'Keeffe (1991) and Brown's
 * accumulated table. b = 0.37 throughout, which is the value those R0 go with.
 * Key: "<element><oxidation>-<anion>". A pair that is not here is reported as
 * "no parameters" -- an invented R0 would produce a valence that looks like a
 * measurement and is not one. */
const BV_B = 0.37;
const BV = {
/* --- oxides ------------------------------------------------------------ */
'Li1-O':1.466,'Na1-O':1.803,'K1-O':2.132,'Rb1-O':2.263,'Cs1-O':2.417,
'Be2-O':1.381,'Mg2-O':1.693,'Ca2-O':1.967,'Sr2-O':2.118,'Ba2-O':2.285,
'B3-O':1.371,'Al3-O':1.651,'Ga3-O':1.730,'In3-O':1.902,'Tl1-O':2.172,'Tl3-O':2.003,
'C4-O':1.390,'Si4-O':1.624,'Ge4-O':1.748,'Sn2-O':1.984,'Sn4-O':1.905,
'Pb2-O':2.112,'Pb4-O':2.042,
'N5-O':1.432,'P5-O':1.617,'As3-O':1.789,'As5-O':1.767,'Sb3-O':1.973,'Sb5-O':1.942,
'Bi3-O':2.094,'Bi5-O':2.060,
'S4-O':1.644,'S6-O':1.624,'Se4-O':1.811,'Se6-O':1.788,'Te4-O':1.977,'Te6-O':1.917,
'Cl3-O':1.720,'Cl7-O':1.632,'Br7-O':1.810,'I5-O':2.000,'I7-O':1.930,
'Sc3-O':1.849,'Y3-O':2.019,'Ti3-O':1.791,'Ti4-O':1.815,'Zr4-O':1.937,'Hf4-O':1.923,
'V3-O':1.743,'V4-O':1.784,'V5-O':1.803,'Nb5-O':1.911,'Ta5-O':1.920,
'Cr3-O':1.724,'Cr6-O':1.794,'Mo6-O':1.907,'W6-O':1.921,'Re7-O':1.970,
'Mn2-O':1.790,'Mn3-O':1.760,'Mn4-O':1.753,'Mn7-O':1.827,
'Fe2-O':1.734,'Fe3-O':1.759,'Co2-O':1.692,'Co3-O':1.700,'Ni2-O':1.654,
'Cu1-O':1.593,'Cu2-O':1.679,'Zn2-O':1.704,'Cd2-O':1.904,'Hg2-O':1.972,
'Ru4-O':1.834,'Rh3-O':1.791,'Pd2-O':1.792,'Ag1-O':1.805,'Pt2-O':1.768,'Pt4-O':1.879,
'Au1-O':1.830,'Au3-O':1.890,'Os6-O':1.925,'Ir4-O':1.870,
'La3-O':2.172,'Ce3-O':2.151,'Ce4-O':2.028,'Pr3-O':2.138,'Nd3-O':2.105,'Sm3-O':2.088,
'Eu3-O':2.074,'Gd3-O':2.065,'Tb3-O':2.049,'Dy3-O':2.036,'Ho3-O':2.023,'Er3-O':2.010,
'Tm3-O':1.995,'Yb3-O':1.985,'Lu3-O':1.971,
'Th4-O':2.167,'U4-O':2.112,'U6-O':2.075,
/* --- fluorides --------------------------------------------------------- */
'Li1-F':1.360,'Na1-F':1.677,'K1-F':1.992,'Rb1-F':2.160,'Cs1-F':2.330,
'Mg2-F':1.581,'Ca2-F':1.842,'Sr2-F':2.019,'Ba2-F':2.188,'Al3-F':1.545,'Si4-F':1.580,
'Sc3-F':1.760,'Y3-F':1.904,'La3-F':2.019,'Ti4-F':1.723,'Zr4-F':1.854,'Nb5-F':1.870,
'Ta5-F':1.880,'Cr3-F':1.657,'Mn2-F':1.698,'Fe2-F':1.650,'Fe3-F':1.679,'Co2-F':1.640,
'Ni2-F':1.596,'Cu2-F':1.594,'Zn2-F':1.620,'Cd2-F':1.811,'Pb2-F':2.030,'Ag1-F':1.800,
/* --- chlorides --------------------------------------------------------- */
'Li1-Cl':1.910,'Na1-Cl':2.150,'K1-Cl':2.519,'Rb1-Cl':2.652,'Cs1-Cl':2.791,
'Mg2-Cl':2.080,'Ca2-Cl':2.370,'Sr2-Cl':2.510,'Ba2-Cl':2.690,'Al3-Cl':2.030,
'Ti4-Cl':2.190,'Zr4-Cl':2.330,'Mn2-Cl':2.133,'Fe2-Cl':2.060,'Fe3-Cl':2.090,
'Co2-Cl':2.033,'Ni2-Cl':2.020,'Cu1-Cl':1.850,'Cu2-Cl':2.000,'Zn2-Cl':2.010,
'Ag1-Cl':2.090,'Cd2-Cl':2.212,'Hg2-Cl':2.280,'Sn2-Cl':2.360,'Pb2-Cl':2.530,
'Bi3-Cl':2.480,'La3-Cl':2.545,'Y3-Cl':2.400,
/* --- sulfides ---------------------------------------------------------- */
'Li1-S':1.940,'Na1-S':2.300,'K1-S':2.590,'Mg2-S':2.160,'Ca2-S':2.450,'Sr2-S':2.640,
'Ba2-S':2.880,'Al3-S':2.130,'Ga3-S':2.163,'In3-S':2.370,'Ge4-S':2.217,'Sn2-S':2.440,
'Sn4-S':2.400,'Pb2-S':2.541,'As3-S':2.272,'Sb3-S':2.474,'Bi3-S':2.550,
'Ti4-S':2.290,'V3-S':2.190,'Cr3-S':2.160,'Mn2-S':2.220,'Fe2-S':2.125,'Fe3-S':2.149,
'Co2-S':2.060,'Ni2-S':2.040,'Cu1-S':1.898,'Cu2-S':2.054,'Zn2-S':2.090,'Ag1-S':2.119,
'Cd2-S':2.304,'Hg2-S':2.320,'Mo6-S':2.350,'W6-S':2.390,'La3-S':2.643,
/* --- nitrides ---------------------------------------------------------- */
'Li1-N':1.610,'Na1-N':1.930,'Mg2-N':1.850,'Ca2-N':2.140,'Al3-N':1.790,'Si4-N':1.724,
'Ga3-N':1.840,'Ge4-N':1.790,'B3-N':1.470,'C4-N':1.442,'Ti4-N':1.930,'P5-N':1.704,
/* --- hydrides / hydrogen bonds ----------------------------------------- */
'H1-O':0.882,'H1-F':0.870,'H1-Cl':1.280,'H1-N':0.880
};

/* ------------------------------------------------------------------ maths */
function cellMatrix(c) {
    const a = +c.a, b = +c.b, cc = +c.c;
    const al = (c.alpha ?? 90) * Math.PI / 180,
          be = (c.beta  ?? 90) * Math.PI / 180,
          ga = (c.gamma ?? 90) * Math.PI / 180;
    const ca = Math.cos(al), cb = Math.cos(be), cg = Math.cos(ga), sg = Math.sin(ga);
    const vFac = Math.sqrt(Math.max(0,
        1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg));
    // Columns are a, b, c in Cartesian. frac -> cart is M * [x,y,z].
    const M = [
        a, b * cg,  cc * cb,
        0, b * sg,  cc * (ca - cb * cg) / sg,
        0, 0,       cc * vFac / sg
    ];
    return { M, volume: a * b * cc * vFac };
}
const toCart = (M, f) => [
    M[0] * f[0] + M[1] * f[1] + M[2] * f[2],
    M[3] * f[0] + M[4] * f[1] + M[5] * f[2],
    M[6] * f[0] + M[7] * f[1] + M[8] * f[2]
];
const wrap01 = v => { const f = v - Math.floor(v); return (f >= 1 - 1e-8) ? 0 : f; };
const sym = (el) => String(el || '').replace(/[^A-Za-z]/g, '')
                        .replace(/^(.)(.*)$/, (m, a, b) => a.toUpperCase() + b.toLowerCase());

/* ------------------------------------------------------- symmetry expansion */
function getOps(st) {
    if (Array.isArray(st.opMatrices) && st.opMatrices.length) return st.opMatrices;
    if (typeof global.getSymopsForSetting === 'function' && st.spaceGroup) {
        const got = global.getSymopsForSetting(st.spaceGroup);
        if (got && Array.isArray(got.symops) && got.symops.length) return got.symops;
    }
    return null;
}

/** One entry per atom in the unit cell, tagged with the asymmetric site it
 *  came from. Multiplicity is COUNTED here rather than trusted from the
 *  Wyckoff label: a site that has drifted onto a special position during
 *  refinement generates fewer distinct images, and the count is what the
 *  contacts and the cell content are actually built from. */
function expandCell(sites, ops) {
    const atoms = [];
    const mult = [];
    sites.forEach((s, si) => {
        const seen = new Set();
        for (const op of ops) {
            const r = op.r, t = op.t || [0, 0, 0];
            const f = [
                wrap01(r[0] * s.x + r[1] * s.y + r[2] * s.z + t[0]),
                wrap01(r[3] * s.x + r[4] * s.y + r[5] * s.z + t[1]),
                wrap01(r[6] * s.x + r[7] * s.y + r[8] * s.z + t[2])
            ];
            const key = f.map(v => Math.round(v * 2000)).join(',');
            if (seen.has(key)) continue;
            seen.add(key);
            atoms.push({ site: si, element: sym(s.element), f });
        }
        mult.push(seen.size);
    });
    return { atoms, mult };
}

/* Neighbours of one fractional position, over the 27 nearest cells. */
function neighboursOf(f0, atoms, M, cutoff, selfTol) {
    const p0 = toCart(M, f0);
    const out = [];
    for (const at of atoms) {
        for (let i = -1; i <= 1; i++)
        for (let j = -1; j <= 1; j++)
        for (let k = -1; k <= 1; k++) {
            const f = [at.f[0] + i, at.f[1] + j, at.f[2] + k];
            const p = toCart(M, f);
            const dx = p[0] - p0[0], dy = p[1] - p0[1], dz = p[2] - p0[2];
            const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (d < (selfTol ?? 0.05) || d > cutoff) continue;
            out.push({ d, element: at.element, site: at.site, f, p });
        }
    }
    out.sort((u, v) => u.d - v.d);
    return out;
}

/** Contacts up to the first break larger than `gap` in the sorted distances,
 *  taken SEPARATELY for each counter-ion element.
 *
 *  An oxygen in a sulfate sits 1.5 A from its sulfur and 2.4 A from the
 *  nearest lead. Cutting the merged list at the first break would stop at the
 *  S-O bond and report CN 1, which describes the sulfate group rather than the
 *  oxygen's coordination. Each element gets its own break, and the results are
 *  merged. */
function firstShell(sorted, gap) {
    if (!sorted.length) return [];
    const byEl = new Map();
    for (const n of sorted) {
        if (!byEl.has(n.element)) byEl.set(n.element, []);
        byEl.get(n.element).push(n);
    }
    const out = [];
    for (const list of byEl.values()) {
        out.push(list[0]);
        for (let i = 1; i < list.length; i++) {
            if (list[i].d - list[i - 1].d > gap) break;
            out.push(list[i]);
        }
    }
    out.sort((u, v) => u.d - v.d);
    return out;
}

/* ------------------------------------------------------------- the analysis */
/**
 * @param {object} st       structure from the Wyckoff search
 * @param {object} opts     {cutoff, ox:{elementOrRank: state}}
 */
function analyse(st, opts) {
    opts = opts || {};
    const cutoff = Number.isFinite(opts.cutoff) ? opts.cutoff : 3.5;
    const cell = st.cell || {};
    const { M, volume } = cellMatrix(cell);
    const ops = getOps(st);
    if (!ops) return { error: 'No symmetry operators available for this space group, ' +
                              'so the cell cannot be expanded. Contacts and coordination ' +
                              'numbers would be those of P1 and are not shown.' };

    const sites = st.sites.map((s, i) => Object.assign({}, s, { element: sym(s.element), __i: i }));
    const anions = anionSet(sites.map(s => s.element));
    const { atoms, mult } = expandCell(sites, ops);

    /* Cell content and density from the counted multiplicities. */
    let mass = 0;
    const counts = {};
    sites.forEach((s, i) => {
        const occ = Number.isFinite(s.occupancy) ? s.occupancy : 1;
        const n = mult[i] * occ;
        counts[s.element] = (counts[s.element] || 0) + n;
        mass += n * ((EL[s.element] || [0])[0] || 0);
    });
    const density = volume > 0 ? mass / (AVOGADRO * volume * 1e-24) : NaN;

    /* Per-site contacts, coordination and bond valence. */
    const report = sites.map((s, i) => {
        const nb = neighboursOf([s.x, s.y, s.z], atoms, M, cutoff, 0.2);
        const isAnion = anions.has(s.element);
        const ox = oxidationOf(s, opts.ox, anions);
        // Coordination shell: the counter-ions. For an anion that means the
        // cations around it, which is what makes the anion's own valence sum
        // meaningful.
        const counter = nb.filter(n => anions.has(n.element) !== isAnion);
        // THE COORDINATION SHELL IS CUT AT THE GAP, NOT AT THE CUTOFF.
        //
        // A single distance cutoff cannot serve both a sulfate sulfur and a
        // twelve-coordinate lead: at 3.5 A the S would be reported as CN 5 or
        // 6, with a mean bond length halfway between the SO4 tetrahedron and
        // the next O in the lattice, and a spread that describes nothing.
        // Sorting the contacts and cutting at the first clear break gives the
        // first coordination sphere in both cases. The bond-valence sum below
        // still runs over everything inside the cutoff, where the exponential
        // makes the far contacts irrelevant on their own terms.
        const shell = firstShell(counter, 0.40);
        const ds = shell.map(n => n.d);
        const mean = ds.length ? ds.reduce((a, b) => a + b, 0) / ds.length : NaN;
        const spread = ds.length
            ? Math.sqrt(ds.reduce((a, b) => a + (b - mean) * (b - mean), 0) / ds.length)
            : NaN;

        let bvs = 0, missing = new Set(), used = 0;
        for (const n of counter) {
            const cat = isAnion ? n.element : s.element;
            const an  = isAnion ? s.element : n.element;
            const cOx = isAnion
                ? oxidationOf(sites[n.site], opts.ox, anions)
                : ox;
            const key = `${cat}${Math.abs(cOx)}-${an}`;
            const r0 = BV[key];
            if (!Number.isFinite(r0)) { missing.add(key); continue; }
            bvs += Math.exp((r0 - n.d) / BV_B);
            used++;
        }
        const expected = Math.abs(ox);
        const haveBV = used > 0 && used === counter.length && expected > 0;
        const dev = haveBV ? Math.abs(bvs - expected) / expected : NaN;

        return {
            index: i, rank: s.rank, element: s.element, wyckoff: s.wyckoff || String(mult[i]),
            multiplicity: mult[i], x: s.x, y: s.y, z: s.z, isAnion, ox, expected,
            contacts: nb, shell, cn: shell.length,
            shellComp: composition(shell), mean, spread, shellMax: ds.length ? ds[ds.length - 1] : NaN,
            bvs: haveBV ? bvs : NaN, bvsPartial: used > 0 && used < counter.length,
            missingParams: Array.from(missing),
            nBvs: used,
            verdict: !haveBV ? 'no parameters'
                   : dev < 0.15 ? 'plausible'
                   : dev < 0.30 ? 'borderline' : 'check',
            deviation: dev
        };
    });

    // Charge balance over the whole cell. A structure whose formal charges do
    // not cancel is either mis-assigned above or is not the composition the
    // search was given, and every bond-valence verdict below inherits that.
    const netCharge = report.reduce((a, r) => a + r.ox * mult[r.index], 0);

    return { M, volume, mass, density, counts, atoms, mult, sites, report,
             cutoff, ops, anions, netCharge };
}

/** "3 O + 1 S", in order of decreasing count. */
function composition(shell) {
    const c = {};
    shell.forEach(n => { c[n.element] = (c[n.element] || 0) + 1; });
    return Object.keys(c).sort((a, b) => c[b] - c[a])
                 .map(e => `${c[e]}\u00d7${e}`).join(' + ');
}

function oxidationOf(site, overrides, anions) {
    if (!site) return 0;
    const el = sym(site.element);
    if (overrides && Number.isFinite(overrides[site.__i])) return overrides[site.__i];
    if (anions && anions.has(el)) return ANION_CHARGE[el] ?? -1;
    return DEFAULT_OX[el] ?? 2;
}

/* ==========================================================================
 *  Structure plot. Canvas 2D, orthographic, depth-sorted spheres and split
 *  bonds. Three.js is loaded for the density viewer but that viewer owns a
 *  WebGL context and a whole map; a report figure that has to survive being
 *  re-rendered on every run selection is better off without one.
 * ======================================================================== */
function Plot(canvas, data, state) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.data = data;          // {M, atoms, report, cutoff}
    this.state = state;        // {yaw, pitch, zoom, hidden:Set, selected, showCell, showLabels, pack}
    this.build();
}

Plot.prototype.build = function () {
    const { M, atoms } = this.data;
    const pack = this.state.pack ? 0.08 : 0.001;   // fractional margin for images
    const list = [];
    for (const at of atoms) {
        for (let i = -1; i <= 1; i++)
        for (let j = -1; j <= 1; j++)
        for (let k = -1; k <= 1; k++) {
            const f = [at.f[0] + i, at.f[1] + j, at.f[2] + k];
            if (f.some(v => v < -pack || v > 1 + pack)) continue;
            list.push({ site: at.site, element: at.element, f, p: toCart(M, f) });
        }
    }
    // Bonds: cation-anion pairs inside a radius-scaled cutoff. Drawn between
    // the atoms actually on screen, so a polyhedron cut by the box edge shows
    // the bonds it has and no more.
    const bonds = [];
    for (let a = 0; a < list.length; a++) {
        for (let b = a + 1; b < list.length; b++) {
            const A = list[a], B = list[b];
            if (this.data.anions.has(A.element) === this.data.anions.has(B.element)) continue;
            const ra = (EL[A.element] || [0, 1.2])[1], rb = (EL[B.element] || [0, 1.2])[1];
            const lim = Math.min(this.data.cutoff, 1.40 * (ra + rb));
            const dx = A.p[0] - B.p[0], dy = A.p[1] - B.p[1], dz = A.p[2] - B.p[2];
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 <= lim * lim && d2 > 0.04) bonds.push([a, b]);
        }
    }
    this.list = list;
    this.bonds = bonds;

    // Cell frame.
    const corners = [];
    for (let i = 0; i < 8; i++)
        corners.push(toCart(M, [(i >> 2) & 1, (i >> 1) & 1, i & 1]));
    this.corners = corners;
    this.edges = [[0,1],[0,2],[0,4],[1,3],[1,5],[2,3],[2,6],[3,7],[4,5],[4,6],[5,7],[6,7]];
    this.centre = toCart(M, [0.5, 0.5, 0.5]);
    this.extent = Math.max(...corners.map(c => Math.hypot(
        c[0] - this.centre[0], c[1] - this.centre[1], c[2] - this.centre[2]))) || 1;
};

Plot.prototype.rot = function (p) {
    const s = this.state;
    const x = p[0] - this.centre[0], y = p[1] - this.centre[1], z = p[2] - this.centre[2];
    const cy = Math.cos(s.yaw), sy = Math.sin(s.yaw);
    const x1 = cy * x + sy * z, z1 = -sy * x + cy * z;
    const cp = Math.cos(s.pitch), sp = Math.sin(s.pitch);
    const y1 = cp * y - sp * z1, z2 = sp * y + cp * z1;
    return [x1, y1, z2];
};

Plot.prototype.draw = function () {
    const cv = this.canvas, ctx = this.ctx, s = this.state;
    const dpr = global.devicePixelRatio || 1;
    const w = cv.clientWidth || 480, h = cv.clientHeight || 320;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
        cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);

    const scale = (Math.min(w, h) * 0.42 / this.extent) * s.zoom;
    const px = p => { const r = this.rot(p); return [w / 2 + r[0] * scale, h / 2 - r[1] * scale, r[2]]; };

    // Cell frame, behind everything.
    if (s.showCell) {
        ctx.strokeStyle = cssVar('--rule-strong', '#8a8f96');
        ctx.lineWidth = 1;
        ctx.globalAlpha = 0.75;
        for (const [i, j] of this.edges) {
            const A = px(this.corners[i]), B = px(this.corners[j]);
            ctx.beginPath(); ctx.moveTo(A[0], A[1]); ctx.lineTo(B[0], B[1]); ctx.stroke();
        }
        ctx.globalAlpha = 1;
        // a, b, c labels at the three edges from the origin.
        const O = px(this.corners[0]);
        [['c', 1], ['b', 2], ['a', 4]].forEach(([name, idx]) => {
            const P = px(this.corners[idx]);
            ctx.fillStyle = cssVar('--fg-3', '#7a7f86');
            ctx.font = '10px ' + cssVar('--font-mono', 'monospace');
            ctx.fillText(name, O[0] + (P[0] - O[0]) * 1.06, O[1] + (P[1] - O[1]) * 1.06);
        });
    }

    const visible = i => !s.hidden.has(this.list[i].site);

    // Painter's algorithm over bonds and atoms together.
    const items = [];
    this.bonds.forEach(([a, b]) => {
        if (!visible(a) || !visible(b)) return;
        const A = px(this.list[a].p), B = px(this.list[b].p);
        items.push({ z: (A[2] + B[2]) / 2, kind: 'bond', A, B, a, b });
    });
    this.list.forEach((at, i) => {
        if (!visible(i)) return;
        const P = px(at.p);
        items.push({ z: P[2], kind: 'atom', P, i });
    });
    items.sort((u, v) => u.z - v.z);

    const rOf = el => 0.34 * ((EL[el] || [0, 1.2])[1]) * scale + 2;

    for (const it of items) {
        if (it.kind === 'bond') {
            const cA = (EL[this.list[it.a].element] || [0, 0, '#888'])[2];
            const cB = (EL[this.list[it.b].element] || [0, 0, '#888'])[2];
            const mx = (it.A[0] + it.B[0]) / 2, my = (it.A[1] + it.B[1]) / 2;
            ctx.lineWidth = Math.max(1.5, scale * 0.045);
            ctx.lineCap = 'round';
            ctx.strokeStyle = cA;
            ctx.beginPath(); ctx.moveTo(it.A[0], it.A[1]); ctx.lineTo(mx, my); ctx.stroke();
            ctx.strokeStyle = cB;
            ctx.beginPath(); ctx.moveTo(mx, my); ctx.lineTo(it.B[0], it.B[1]); ctx.stroke();
        } else {
            const at = this.list[it.i];
            const r = rOf(at.element);
            const col = (EL[at.element] || [0, 0, '#888'])[2];
            const g = ctx.createRadialGradient(
                it.P[0] - r * 0.35, it.P[1] - r * 0.35, r * 0.1, it.P[0], it.P[1], r);
            g.addColorStop(0, mix(col, '#ffffff', 0.55));
            g.addColorStop(1, mix(col, '#000000', 0.25));
            ctx.fillStyle = g;
            ctx.beginPath(); ctx.arc(it.P[0], it.P[1], r, 0, 6.2832); ctx.fill();
            ctx.lineWidth = 0.8;
            ctx.strokeStyle = 'rgba(0,0,0,0.35)';
            ctx.stroke();
            if (at.site === s.selected) {
                ctx.strokeStyle = cssVar('--accent', '#4a7fb5');
                ctx.lineWidth = 2;
                ctx.beginPath(); ctx.arc(it.P[0], it.P[1], r + 3, 0, 6.2832); ctx.stroke();
            }
            if (s.showLabels) {
                ctx.fillStyle = cssVar('--fg', '#1a1c1f');
                ctx.font = '10px ' + cssVar('--font-mono', 'monospace');
                ctx.fillText(at.element, it.P[0] + r + 2, it.P[1] - r * 0.2);
            }
        }
    }
    this._px = px;
};

/** Nearest drawn atom to a canvas point, for click-to-select. */
Plot.prototype.pick = function (cx, cy) {
    if (!this._px) return -1;
    let best = -1, bestD = 18 * 18, bestZ = -Infinity;
    this.list.forEach((at, i) => {
        if (this.state.hidden.has(at.site)) return;
        const P = this._px(at.p);
        const d = (P[0] - cx) ** 2 + (P[1] - cy) ** 2;
        if (d <= bestD && P[2] > bestZ) { best = i; bestD = d; bestZ = P[2]; }
    });
    return best < 0 ? -1 : this.list[best].site;
};

function cssVar(name, fallback) {
    try {
        const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
        return v || fallback;
    } catch (e) { return fallback; }
}
function mix(hex, other, t) {
    const p = h => {
        h = String(h).replace('#', '');
        if (h.length === 3) h = h.split('').map(c => c + c).join('');
        return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
    };
    const A = p(hex), B = p(other);
    const c = A.map((v, i) => Math.round(v + (B[i] - v) * t));
    return `rgb(${c[0]},${c[1]},${c[2]})`;
}

/* ============================================================== the markup */
const CSS = `
.wy-an { margin-top: 18px; }
.wy-an h3 { font-size: 13px; margin: 16px 0 6px 0; }
.wy-an .wy-head { display:flex; justify-content:space-between; align-items:baseline;
  gap:12px; flex-wrap:wrap; font-size:11px; color:var(--fg-3,#7a7f86);
  font-family:var(--font-mono,monospace); }
.wy-an .wy-controls { display:flex; gap:14px; align-items:center; flex-wrap:wrap;
  font-size:11px; color:var(--fg-2,#4a4e54); margin:8px 0; }
.wy-an .wy-controls label { display:flex; gap:5px; align-items:center; }
.wy-an .wy-controls input[type=number] { width:64px; padding:2px 4px; font-size:11px;
  font-family:var(--font-mono,monospace); }
.wy-plot-wrap { position:relative; border:1px solid var(--rule,#c4c7cb);
  border-radius:var(--r-3,4px); background:var(--bg-stage,#fff); overflow:hidden; }
.wy-plot { display:block; width:100%; height:680px; cursor:grab; touch-action:none; }
.wy-plot.dragging { cursor:grabbing; }
.wy-plot-hint { position:absolute; left:8px; bottom:6px; font-size:10px;
  color:var(--fg-4,#a0a4ab); font-family:var(--font-mono,monospace); pointer-events:none; }
.wy-plot-btns { position:absolute; right:6px; top:6px; display:flex; gap:4px; }
.wy-plot-btns button { font-size:10px; padding:2px 6px; }
.wy-an table.wy-tab { width:100%; border-collapse:collapse; font-size:11.5px;
  font-family:var(--font-mono,monospace); }
.wy-an table.wy-tab th { text-align:right; font-weight:500; color:var(--fg-3,#7a7f86);
  border-bottom:1px solid var(--rule,#c4c7cb); padding:4px 6px; font-size:10px;
  text-transform:uppercase; letter-spacing:0.04em; }
.wy-an table.wy-tab td { text-align:right; padding:3px 6px;
  border-bottom:1px solid var(--rule-soft,#d8dadd); }
.wy-an table.wy-tab th:first-child, .wy-an table.wy-tab td:first-child { text-align:left; }
.wy-an tr.wy-site { cursor:pointer; }
.wy-an tr.wy-site:hover { background:var(--bg-hover,rgba(0,0,0,0.045)); }
.wy-an tr.wy-site.sel { background:var(--accent-tint,rgba(74,127,181,0.14)); }
.wy-an tr.wy-site.off td { opacity:0.4; }
.wy-eye { border:none; background:none; cursor:pointer; font-size:12px; padding:0 2px;
  color:var(--fg-3,#7a7f86); }
.wy-ox { width:46px; font-size:11px; padding:1px 3px; font-family:var(--font-mono,monospace); }
.wy-contacts { max-height:190px; overflow:auto; border:1px solid var(--rule,#c4c7cb);
  border-radius:var(--r-3,4px); padding:8px 10px; font-family:var(--font-mono,monospace);
  font-size:11.5px; }
.wy-contacts .wy-ctitle { color:var(--fg-3,#7a7f86); margin-bottom:6px; }
.wy-contacts .wy-crow { display:grid; grid-template-columns:2.5em 6em 1fr; gap:8px;
  padding:1px 0; }
.wy-contacts .wy-gap { border-top:1px dashed var(--rule,#c4c7cb); margin:4px 0; }
.wy-v-plausible { color:var(--success,#3a7a4a); }
.wy-v-borderline { color:var(--warning,#b88318); }
.wy-v-check { color:var(--danger,#b03a3a); }
.wy-v-none { color:var(--fg-3,#7a7f86); }
.wy-note { font-size:11px; color:var(--fg-3,#7a7f86); margin:6px 0 0 0; line-height:1.5; }
`;

function ensureCSS() {
    if (document.getElementById('wy-report-css')) return;
    const st = document.createElement('style');
    st.id = 'wy-report-css';
    st.textContent = CSS;
    document.head.appendChild(st);
}

const esc = s => String(s).replace(/[&<>"]/g, c =>
    ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const f2 = (v, n) => Number.isFinite(v) ? v.toFixed(n) : '\u2013';

/* Per-structure UI state, keyed so re-rendering the same run keeps the view. */
const uiState = new WeakMap();

/**
 * Draw the whole analysis block into `host` for structure `st`.
 * Safe to call repeatedly; it replaces its own contents.
 */
function render(host, st) {
    if (!host) return;
    ensureCSS();
    if (!st || !Array.isArray(st.sites) || !st.sites.length || !st.cell) {
        host.innerHTML = '';
        return;
    }

    let S = uiState.get(st);
    if (!S) {
        S = { cutoff: 3.5, ox: {}, hidden: new Set(), selected: 0,
              yaw: -0.6, pitch: 0.35, zoom: 1, showCell: true, showLabels: false, pack: true };
        uiState.set(st, S);
    }

    const A = analyse(st, { cutoff: S.cutoff, ox: S.ox });
    if (A.error) {
        host.innerHTML = `<div class="wy-an"><p class="wy-note">${esc(A.error)}</p></div>`;
        return;
    }

    const formula = Object.keys(A.counts).sort().map(e =>
        `${e}${A.counts[e] === 1 ? '' : fmtCount(A.counts[e])}`).join(' ');

    host.innerHTML = `
<div class="wy-an">
  <div class="wy-head">
    <span>ASYMMETRIC UNIT &nbsp;\u00b7&nbsp; ${esc(formula)} per cell</span>
    <span>Cell content: ${f2(A.mass, 2)} g/mol &nbsp;|&nbsp; V = ${f2(A.volume, 2)} \u00c5\u00b3
      &nbsp;|&nbsp; Calculated density: ${f2(A.density, 3)} g/cm\u00b3${
        Math.abs(A.netCharge) > 0.01
          ? ` &nbsp;|&nbsp; <span class="wy-v-borderline">net formal charge ${A.netCharge > 0 ? '+' : ''}${f2(A.netCharge, 1)}</span>`
          : ''}</span>
  </div>

  <div class="wy-plot-wrap">
    <canvas class="wy-plot" id="wy-canvas"></canvas>
    <div class="wy-plot-btns">
      <button class="btn" id="wy-reset">Reset view</button>
      <button class="btn" id="wy-png">PNG</button>
    </div>
    <div class="wy-plot-hint">drag to rotate \u00b7 wheel to zoom \u00b7 click an atom to list its contacts</div>
  </div>

  <div class="wy-controls">
    <label>Contact cutoff
      <input type="number" id="wy-cutoff" min="1" max="8" step="0.1" value="${S.cutoff}"> \u00c5</label>
    <label><input type="checkbox" id="wy-cell" ${S.showCell ? 'checked' : ''}> unit cell</label>
    <label><input type="checkbox" id="wy-pack" ${S.pack ? 'checked' : ''}> boundary images</label>
    <label><input type="checkbox" id="wy-lab" ${S.showLabels ? 'checked' : ''}> labels</label>
  </div>

  <table class="wy-tab" id="wy-sites">
    <thead><tr>
      <th>El</th><th>Wyck</th><th>x</th><th>y</th><th>z</th>
      <th title="Number of symmetry-equivalent atoms actually generated in the cell.">Mult</th>
      <th title="Formal oxidation state used for the bond-valence sum. Editable.">Ox</th>
      <th title="Show or hide this site in the plot.">Plot</th>
    </tr></thead>
    <tbody>
      ${A.report.map(r => `
      <tr class="wy-site ${r.index === S.selected ? 'sel' : ''} ${S.hidden.has(r.index) ? 'off' : ''}"
          data-row="${r.index}">
        <td>${esc(r.element)}</td><td>${esc(r.wyckoff)}</td>
        <td>${f2(r.x, 4)}</td><td>${f2(r.y, 4)}</td><td>${f2(r.z, 4)}</td>
        <td>${r.multiplicity}</td>
        <td><input class="wy-ox" type="number" step="1" data-i="${r.index}" value="${r.ox}"></td>
        <td><button class="wy-eye" data-eye="${r.index}"
             title="${S.hidden.has(r.index) ? 'Show' : 'Hide'} in plot">${S.hidden.has(r.index) ? '\u25cb' : '\u25c9'}</button></td>
      </tr>`).join('')}
    </tbody>
  </table>

  <h3>Contacts</h3>
  <div class="wy-contacts" id="wy-contacts"></div>

  <h3>Structure quality</h3>
  <table class="wy-tab" id="wy-quality">
    <thead><tr>
      <th>Site</th>
      <th title="Counter-ions within the contact cutoff.">CN</th>
      <th>Mean (\u00c5)</th>
      <th title="Standard deviation of the distances in the coordination shell.">Spread</th>
      <th title="Sum of exp((R0-d)/0.37) over the coordination shell, against the formal valence.">Bond valence</th>
      <th>Verdict</th>
    </tr></thead>
    <tbody>${A.report.map(r => qualityRow(r)).join('')}</tbody>
  </table>
  <p class="wy-note">
    Bond-valence parameters are Brese &amp; O'Keeffe (1991), b = 0.37 \u00c5. A sum within
    15% of the formal valence is reported as plausible, within 30% as borderline.
    Both the sum and the coordination number depend on the ${f2(A.cutoff, 1)} \u00c5 cutoff
    above; distant contacts contribute little to the valence but do raise CN.
    ${A.report.some(r => r.missingParams.length)
      ? 'Some pairs in this structure have no tabulated R<sub>0</sub> and are left out of the sum rather than estimated.'
      : ''}
  </p>
</div>`;

    /* ---- plot ---- */
    const canvas = host.querySelector('#wy-canvas');
    const plot = new Plot(canvas, { M: A.M, atoms: A.atoms, cutoff: A.cutoff, anions: A.anions }, S);
    const redraw = () => plot.draw();
    requestAnimationFrame(redraw);

    let drag = null;
    canvas.addEventListener('pointerdown', e => {
        drag = { x: e.clientX, y: e.clientY, moved: 0 };
        canvas.classList.add('dragging');
        canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener('pointermove', e => {
        if (!drag) return;
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
        drag.moved += Math.abs(dx) + Math.abs(dy);
        S.yaw += dx * 0.01;
        S.pitch = Math.max(-1.5, Math.min(1.5, S.pitch + dy * 0.01));
        drag.x = e.clientX; drag.y = e.clientY;
        redraw();
    });
    const endDrag = (e) => {
        if (!drag) return;
        canvas.classList.remove('dragging');
        if (drag.moved < 4) {
            const r = canvas.getBoundingClientRect();
            const hit = plot.pick(e.clientX - r.left, e.clientY - r.top);
            if (hit >= 0) { S.selected = hit; render(host, st); return; }
        }
        drag = null;
    };
    canvas.addEventListener('pointerup', endDrag);
    canvas.addEventListener('pointercancel', () => { drag = null; canvas.classList.remove('dragging'); });
    canvas.addEventListener('wheel', e => {
        e.preventDefault();
        S.zoom = Math.max(0.3, Math.min(6, S.zoom * (e.deltaY > 0 ? 0.9 : 1.1)));
        redraw();
    }, { passive: false });
    // The panel is often hidden when this runs -- a hidden canvas measures
    // zero, and nothing would ever ask it to draw again once the tab is shown.
    // The observer covers that and window resizing in one go. Disconnected on
    // the next render so repeated run selections do not stack observers.
    if (host.__wyRO) { try { host.__wyRO.disconnect(); } catch (e) {} }
    if (typeof ResizeObserver === 'function') {
        host.__wyRO = new ResizeObserver(() => redraw());
        host.__wyRO.observe(canvas);
    } else {
        host.__wyRO = null;
        global.addEventListener('resize', redraw);
    }

    host.querySelector('#wy-reset').addEventListener('click', () => {
        S.yaw = -0.6; S.pitch = 0.35; S.zoom = 1; redraw();
    });
    host.querySelector('#wy-png').addEventListener('click', () => {
        const a = document.createElement('a');
        a.download = `structure_${(st.spaceGroup && st.spaceGroup.symbol || 'cell').replace(/[^\w]/g, '')}.png`;
        a.href = canvas.toDataURL('image/png');
        a.click();
    });

    /* ---- controls ---- */
    const rerender = () => render(host, st);
    host.querySelector('#wy-cutoff').addEventListener('change', e => {
        const v = parseFloat(e.target.value);
        if (Number.isFinite(v) && v > 0.5) { S.cutoff = v; rerender(); }
    });
    host.querySelector('#wy-cell').addEventListener('change', e => {
        S.showCell = e.target.checked; redraw();
    });
    host.querySelector('#wy-lab').addEventListener('change', e => {
        S.showLabels = e.target.checked; redraw();
    });
    host.querySelector('#wy-pack').addEventListener('change', e => {
        S.pack = e.target.checked; plot.build(); redraw();
    });

    host.querySelectorAll('.wy-ox').forEach(inp => {
        inp.addEventListener('change', e => {
            e.stopPropagation();
            const v = parseInt(e.target.value, 10);
            if (Number.isFinite(v)) { S.ox[+e.target.dataset.i] = v; rerender(); }
        });
        inp.addEventListener('click', e => e.stopPropagation());
    });
    host.querySelectorAll('[data-eye]').forEach(b => {
        b.addEventListener('click', e => {
            e.stopPropagation();
            const i = +b.dataset.eye;
            if (S.hidden.has(i)) S.hidden.delete(i); else S.hidden.add(i);
            rerender();
        });
    });
    host.querySelectorAll('tr.wy-site').forEach(tr => {
        tr.addEventListener('click', () => { S.selected = +tr.dataset.row; rerender(); });
    });

    /* ---- contacts panel ---- */
    renderContacts(host.querySelector('#wy-contacts'), A, S.selected);
}

function fmtCount(n) {
    return Number.isInteger(n) ? String(n) : n.toFixed(2);
}

function qualityRow(r) {
    const cls = r.verdict === 'plausible' ? 'wy-v-plausible'
              : r.verdict === 'borderline' ? 'wy-v-borderline'
              : r.verdict === 'check' ? 'wy-v-check' : 'wy-v-none';
    const bv = Number.isFinite(r.bvs)
        ? `${r.bvs.toFixed(2)} of ${r.expected} expected`
        : (r.missingParams.length
            ? `no R\u2080 for ${esc(r.missingParams[0])}`
            : (r.cn ? '\u2013' : 'no coordination shell'));
    const spreadCls = Number.isFinite(r.spread) && r.spread > 0.15 ? ' class="wy-v-borderline"' : '';
    return `<tr>
      <td>${esc(r.element)} ${esc(r.wyckoff)}</td>
      <td title="${esc(r.shellComp)}">${r.cn}${r.shellComp.indexOf('+') >= 0 ? ` <span class="wy-v-none">(${esc(r.shellComp)})</span>` : ''}</td>
      <td>${f2(r.mean, 3)}</td>
      <td${spreadCls}>${f2(r.spread, 3)}</td>
      <td>${bv}</td>
      <td class="${cls}">${r.verdict}</td>
    </tr>`;
}

function renderContacts(box, A, idx) {
    if (!box) return;
    const r = A.report[idx];
    if (!r) { box.innerHTML = ''; return; }
    const rows = r.contacts.map((c, i) => {
        // A visible break where the distance jumps: the gap after the first
        // coordination sphere is the thing a reader is looking for, and a flat
        // list hides it.
        const prev = r.contacts[i - 1];
        const gap = (prev && c.d - prev.d > 0.45) ? '<div class="wy-gap"></div>' : '';
        const tag = A.report[c.site]
            ? `${A.report[c.site].element} ${A.report[c.site].wyckoff}`
            : `site ${c.site + 1}`;
        return gap + `<div class="wy-crow"><span>${esc(c.element)}</span>` +
               `<span>${c.d.toFixed(3)} \u00c5</span><span>${esc(tag)}</span></div>`;
    }).join('');
    box.innerHTML =
        `<div class="wy-ctitle">Contacts around ${esc(r.element)} (${esc(r.wyckoff)}) at ` +
        `${f2(r.x, 4)}, ${f2(r.y, 4)}, ${f2(r.z, 4)} \u2014 within ${f2(A.cutoff, 1)} \u00c5</div>` +
        (rows || '<div class="wy-crow"><span>\u2013</span><span>no contact within the cutoff</span><span></span></div>');
}

global.WyckoffReport = { render, analyse, Plot, EL, BV };

})(typeof window !== 'undefined' ? window : globalThis);
