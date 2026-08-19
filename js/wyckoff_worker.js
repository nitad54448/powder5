// wyckoff_worker.js
// version 1, 19 august 2026
//
// ---------------------------------------------------------------------------
//  THE WYCKOFF SEARCH, ON ITS OWN.
//
//  This route was living inside charge_flipping_worker.js. Nothing about it
//  needed to be: it takes a space group, a cell, a composition and |Fobs|^2
//  from the Pawley fit, and it never touches an electron-density map. The only
//  thing the two shared was the file, which is why a Wyckoff failure reported
//  itself from a charge-flipping stack ("charge_flipping_worker.js:2782") and
//  read as a coupling that no longer exists.
//
//  Everything the search needs that used to be defined next to the flipping
//  kernel -- the polarisation helpers, applyOp, the atom table and the two
//  observation-row builders -- is duplicated here rather than imported from
//  there, because importing it would put the coupling straight back.
//
//  Messages in:   'wyckoff-search', 'cancel-wyckoff'
//  Messages out:  'wy-start', 'wy-progress', 'wy-log', 'wy-info', 'wy-error',
//                 'wy-structure'  -- powder5.html listens for these names, not
//                 the old 'cf-*' ones.
// ---------------------------------------------------------------------------

// NOTE ON LOAD ORDER. crystal.js must come before symmetry_utils.js: the
// duplicate sharkoLorentzPolarization() that symmetry_utils.js used to carry
// has been removed, and it now uses crystal.js's. Loading it first only costs
// a console warning, but the order below is the supported one.
if (typeof importScripts === 'function') {
    importScripts('constants.js', 'crystal.js',
                  'symmetry_utils.js', 'observations.js', 'wyckoff_assign.js',
                  'scatterers.js', 'contacts.js', 'swarm_wyckoff.js', 'coord_refine.js');
}

// ---------------------------------------------------------------------------
//  Polarisation. Thin names over crystal.js, kept because the observation
//  builders below and the host both call them.
// ---------------------------------------------------------------------------
function polarizationK(pol) {
    // Delegates to crystal.js, which every entry point loads. Kept as a name
    // here because callers in this file and in the host use it; the body is
    // not duplicated, because three copies of this arithmetic is how the
    // K-versus-descriptor confusion started.
    return sharkoPolarizationK(pol);
}

/**
 * Lp from 2-theta and a polarisation DESCRIPTOR (object, or a bare number read
 * as a legacy monochromator 2-theta).
 *
 * If you already hold the ratio K, call sharkoLorentzPolarization(tth, K)
 * instead. Passing K here would be wrong in a way that is easy to miss: K = 0,
 * which is what a neutron pattern needs, would be read as a monochromator
 * angle of zero degrees and fall back to K = 1.
 */
function lorentzPolarization(tthDeg, pol) {
    return sharkoLp(tthDeg, pol);
}

// ---------------------------------------------------------------------------
//  Apply a symmetry operator to a FRACTIONAL COORDINATE (rotation + translation).
//  Not the reciprocal transform used for hkl.
// ---------------------------------------------------------------------------
function applyOp(op, p) {
    const r = op.r, t = op.t;
    return [
        r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + t[0],
        r[3] * p[0] + r[4] * p[1] + r[5] * p[2] + t[1],
        r[6] * p[0] + r[7] * p[1] + r[8] * p[2] + t[2]
    ];
}

// ===========================================================================
//  ATOM DATA. z = atomic number, m = atomic mass (needed for the density-based
//  Z scan), r = van der Waals radius, rc = covalent radius.
// ===========================================================================
const WY_ATOM_DATA = {
    "H":  { z: 1,  r: 1.20, m: 1.008, rc: 0.31   }, "HE": { z: 2,  r: 1.40, m: 4.003, rc: 0.28   }, "LI": { z: 3,  r: 1.82, m: 6.941, rc: 1.28   }, "BE": { z: 4,  r: 1.53, m: 9.012, rc: 0.96   },
    "B":  { z: 5,  r: 1.92, m: 10.811, rc: 0.84  }, "C":  { z: 6,  r: 1.70, m: 12.011, rc: 0.76  }, "N":  { z: 7,  r: 1.55, m: 14.007, rc: 0.71  }, "O":  { z: 8,  r: 1.52, m: 15.999, rc: 0.66  },
    "F":  { z: 9,  r: 1.47, m: 18.998, rc: 0.57  }, "NE": { z: 10, r: 1.54, m: 20.180, rc: 0.58  }, "NA": { z: 11, r: 2.27, m: 22.990, rc: 1.66  }, "MG": { z: 12, r: 1.73, m: 24.305, rc: 1.41  },
    "AL": { z: 13, r: 1.18, m: 26.982, rc: 1.21  }, "SI": { z: 14, r: 2.10, m: 28.085, rc: 1.11  }, "P":  { z: 15, r: 1.80, m: 30.974, rc: 1.07  }, "S":  { z: 16, r: 1.80, m: 32.06, rc: 1.05   },
    "CL": { z: 17, r: 1.75, m: 35.45, rc: 1.02   }, "AR": { z: 18, r: 1.88, m: 39.948, rc: 1.06  }, "K":  { z: 19, r: 2.75, m: 39.098, rc: 2.03  }, "CA": { z: 20, r: 2.31, m: 40.078, rc: 1.76  },
    "SC": { z: 21, r: 2.11, m: 44.956, rc: 1.70  }, "TI": { z: 22, r: 2.00, m: 47.867, rc: 1.60  }, "V":  { z: 23, r: 1.90, m: 50.942, rc: 1.53  }, "CR": { z: 24, r: 1.90, m: 51.996, rc: 1.39  },
    "MN": { z: 25, r: 1.90, m: 54.938, rc: 1.50  }, "FE": { z: 26, r: 1.26, m: 55.845, rc: 1.42  }, "CO": { z: 27, r: 1.90, m: 58.933, rc: 1.38  }, "NI": { z: 28, r: 1.63, m: 58.693, rc: 1.24  },
    "CU": { z: 29, r: 1.40, m: 63.546, rc: 1.32  }, "ZN": { z: 30, r: 1.39, m: 65.38, rc: 1.22   }, "GA": { z: 31, r: 1.87, m: 69.723, rc: 1.22  }, "GE": { z: 32, r: 2.11, m: 72.630, rc: 1.20  },
    "AS": { z: 33, r: 1.85, m: 74.922, rc: 1.19  }, "SE": { z: 34, r: 1.90, m: 78.971, rc: 1.20  }, "BR": { z: 35, r: 1.85, m: 79.904, rc: 1.20  }, "KR": { z: 36, r: 2.02, m: 83.798, rc: 1.16  },
    "RB": { z: 37, r: 3.03, m: 85.468, rc: 2.20  }, "SR": { z: 38, r: 2.49, m: 87.62, rc: 1.95   }, "Y":  { z: 39, r: 2.19, m: 88.906, rc: 1.90  }, "ZR": { z: 40, r: 2.06, m: 91.224, rc: 1.75  },
    "NB": { z: 41, r: 1.98, m: 92.906, rc: 1.64  }, "MO": { z: 42, r: 1.90, m: 95.95, rc: 1.54   }, "TC": { z: 43, r: 1.83, m: 98.0, rc: 1.47    }, "RU": { z: 44, r: 1.78, m: 101.07, rc: 1.46  },
    "RH": { z: 45, r: 1.73, m: 102.906, rc: 1.42 }, "PD": { z: 46, r: 1.63, m: 106.42, rc: 1.39  }, "AG": { z: 47, r: 1.72, m: 107.868, rc: 1.45 }, "CD": { z: 48, r: 1.58, m: 112.411, rc: 1.44 },
    "IN": { z: 49, r: 1.93, m: 114.818, rc: 1.42 }, "SN": { z: 50, r: 2.17, m: 118.71, rc: 1.39  }, "SB": { z: 51, r: 2.06, m: 121.760, rc: 1.39 }, "TE": { z: 52, r: 2.06, m: 127.6, rc: 1.38   },
    "I":  { z: 53, r: 1.98, m: 126.904, rc: 1.39 }, "XE": { z: 54, r: 2.16, m: 131.293, rc: 1.40 }, "CS": { z: 55, r: 3.43, m: 132.905, rc: 2.44 }, "BA": { z: 56, r: 2.68, m: 137.327, rc: 2.15 },
    "LA": { z: 57, r: 2.43, m: 138.905, rc: 2.07 }, "CE": { z: 58, r: 2.42, m: 140.116, rc: 2.04 }, "PR": { z: 59, r: 2.40, m: 140.908, rc: 2.03 }, "ND": { z: 60, r: 2.39, m: 144.242, rc: 2.01 },
    "PM": { z: 61, r: 2.38, m: 145.0, rc: 1.99   }, "SM": { z: 62, r: 2.36, m: 150.36, rc: 1.98  }, "EU": { z: 63, r: 2.35, m: 151.964, rc: 1.98 }, "GD": { z: 64, r: 2.34, m: 157.25, rc: 1.96  }, 
    "TB": { z: 65, r: 2.33, m: 158.925, rc: 1.94 }, "DY": { z: 66, r: 2.31, m: 162.500, rc: 1.92 }, "HO": { z: 67, r: 2.30, m: 164.930, rc: 1.92 }, "ER": { z: 68, r: 2.29, m: 167.259, rc: 1.89 }, 
    "TM": { z: 69, r: 2.27, m: 168.934, rc: 1.90 }, "YB": { z: 70, r: 2.26, m: 173.054, rc: 1.87 }, "LU": { z: 71, r: 2.24, m: 174.967, rc: 1.87 }, "HF": { z: 72, r: 2.12, m: 178.49, rc: 1.75  }, 
    "TA": { z: 73, r: 2.06, m: 180.948, rc: 1.70 }, "W":  { z: 74, r: 2.02, m: 183.84, rc: 1.62  }, "RE": { z: 75, r: 1.97, m: 186.207, rc: 1.51 }, "OS": { z: 76, r: 1.92, m: 190.23, rc: 1.44  }, 
    "IR": { z: 77, r: 1.87, m: 192.217, rc: 1.41 }, "PT": { z: 78, r: 1.75, m: 195.084, rc: 1.36 }, "AU": { z: 79, r: 1.66, m: 196.967, rc: 1.36 }, "HG": { z: 80, r: 1.55, m: 200.592, rc: 1.32 }, 
    "TL": { z: 81, r: 1.96, m: 204.38, rc: 1.45  }, "PB": { z: 82, r: 2.02, m: 207.2, rc: 1.46   }, "BI": { z: 83, r: 2.07, m: 208.980, rc: 1.48 }, "PO": { z: 84, r: 1.97, m: 209.0, rc: 1.40   }, 
    "AT": { z: 85, r: 2.02, m: 210.0, rc: 1.50   }, "RN": { z: 86, r: 2.20, m: 222.018, rc: 1.50 }, "FR": { z: 87, r: 3.48, m: 223.0, rc: 2.60   }, "RA": { z: 88, r: 2.83, m: 226.0, rc: 2.21   }, 
    "AC": { z: 89, r: 2.47, m: 227.0, rc: 2.15   }, "TH": { z: 90, r: 2.37, m: 232.038, rc: 2.06 }, "PA": { z: 91, r: 2.43, m: 231.036, rc: 2.00 }, "U":  { z: 92, r: 1.86, m: 238.029, rc: 1.96 }
};

/**
 * The raw observation rows, before normalisation.
 *
 * wyObservationsFromHkl builds these and then hands them to
 * normaliseObservations. runWyckoffSearch wants them UNNORMALISED, because it
 * runs the normalisation itself and derives the Wilson B from the result --
 * so this shares the row construction rather than duplicating it, and the map
 * and the search cannot end up disagreeing about Lp.
 *
 * @param {object} job
 * @returns {object[]} rows carrying h, k, l, Ihkl, lp, multiplicity, d, tth
 */
function wyRawObservationRows(job) {
    const applyLP = job.applyLP !== false;
    const polModel = (job.polarization !== undefined && job.polarization !== null)
        ? job.polarization : job.monochromatorTth;

    const rows = [];
    for (const peak of job.hklList || []) {
        if (!peak) continue;
        const I = Number(peak.intensity);
        const h = peak.h_orig, k = peak.k_orig, l = peak.l_orig;
        if (!Number.isFinite(I) || I < 0) continue;
        if (!Number.isFinite(h) || !Number.isFinite(k) || !Number.isFinite(l)) continue;
        if (h === 0 && k === 0 && l === 0) continue;
        // Ten bits per index in the kernel's packing, biased by 512.
        if (Math.abs(h) > 511 || Math.abs(k) > 511 || Math.abs(l) > 511) continue;

        rows.push({
            h, k, l, Ihkl: I,
            multiplicity: peak.multiplicity,
            d: peak.d, twoTheta: peak.tth,
            sigma: peak.sigma,
            lp: applyLP
                ? ((Number.isFinite(peak.lp) && peak.lp > 0) ? peak.lp
                                                            : lorentzPolarization(peak.tth, polModel))
                : 1
        });
    }
    return rows;
}

function wyObservationsFromHkl(job, symops) {
    const applyLP = job.applyLP !== false;
    const polModel = (job.polarization !== undefined && job.polarization !== null)
        ? job.polarization : job.monochromatorTth;

    const src = [];
    for (const peak of job.hklList || []) {
        if (!peak) continue;
        const I = Number(peak.intensity);
        const h = peak.h_orig, k = peak.k_orig, l = peak.l_orig;
        if (!Number.isFinite(I) || I < 0) continue;
        if (!Number.isFinite(h) || !Number.isFinite(k) || !Number.isFinite(l)) continue;
        if (h === 0 && k === 0 && l === 0) continue;

        const row = { h, k, l, Ihkl: I, multiplicity: peak.multiplicity,
                      d: peak.d, twoTheta: peak.tth };
        if (applyLP) {
            // Prefer the host's Lp - it was evaluated at the refined 2-theta
            // including the zero shift, and is the same number the map was
            // built from. Fall back to computing it HERE, with this file's own
            // lorentzPolarization(), rather than leaving the field absent and
            // letting normaliseObservations reconstruct it. Two reasons, and
            // both matter:
            //
            //   It is the identical call buildReflectionModel() makes in
            //   charge_flipping_worker.js, so a map built there and a Wyckoff
            //   search run here cannot end up with different Lp for the same
            //   reflection. The two workers no longer share a file; they still
            //   have to share this number.
            //
            //   It also keeps the rows on normaliseObservations' m_lp route
            //   rather than its reconstruct-from-2-theta fallback, so there is
            //   exactly one place in this program where Lp is decided.
            row.lp = (Number.isFinite(peak.lp) && peak.lp > 0)
                ? peak.lp
                : lorentzPolarization(peak.tth, polModel);
        } else {
            // Already Lp-free: divide out the multiplicity and nothing else.
            row.lp = 1;
        }
        src.push(row);
    }
    if (!src.length) return { rows: [], route: null, warnings: [], errors: ['No usable Pawley intensities.'] };

    return normaliseObservations(src, {
        rotations: (symops || []).map(op => op.r),
        wavelength: job.wavelength,
        radiation: job.radiation || 'xray',
        polarisationK: polarizationK(polModel)
    });
}

// ===========================================================================
//  AUTOMATIC Z: A SCAN, NOT A GUESS.
//
//  The old auto-Z called suggestZ(), took its single `likely` value, and ran
//  with it. suggestZ works from a volume-per-atom rule of thumb (11-22 A^3)
//  and then prefers a Z that divides the group order, which for a pyrochlore
//  in Fd-3m returned Z = 6 -- and Z = 6 puts 12 Gd in a cell whose Wyckoff
//  multiplicities are all multiples of 8. The search then died on an
//  arithmetic impossibility that was never the user's to fix, while the
//  correct answer, Z = 8, sat two entries away in the same list.
//
//  What replaces it:
//
//   1. A WIDE SCAN. Every integer Z whose crystallographic density falls
//      between WY_RHO_MIN and WY_RHO_MAX g/cm^3 is a candidate. That spans
//      framework silicates and MOFs at the light end through the densest
//      intermetallics and actinides at the heavy end, which is far wider than
//      any volume-per-atom rule and cannot silently exclude the true Z.
//
//   2. AN ARITHMETIC FILTER, APPLIED BEFORE THE SEARCH. Each candidate is
//      tested by the same enumerateAssignments() the search itself uses, so a
//      Z that cannot produce a single legal assignment is discarded here --
//      cheaply, in a few milliseconds -- instead of ending the run half a
//      second later with a stack trace.
//
//   3. A PLAUSIBILITY RANKING among the survivors: volume per atom nearest
//      WY_VPA_TYPICAL, with a mild preference for a Z commensurate with the
//      group order. The whole ranked list is kept, so if the chosen Z fails
//      later for a reason the filter cannot see (a GPU atom budget, a
//      composition the coordination windows reject) the next one is tried
//      rather than the run being abandoned.
//
//  A Z the user typed is still honoured exactly as typed. This is what "Auto"
//  means, not what Z means.
// ===========================================================================
const WY_RHO_MIN = 1.0;    // g/cm^3
const WY_RHO_MAX = 20.0;   // g/cm^3
const WY_VPA_TYPICAL = 16.0;  // A^3 per atom, the centre of the plausible band
const WY_Z_HARD_MAX = 192;    // no Z above the largest general multiplicity
const WY_MAX_ATOMS_PER_CELL = 512;  // beyond this no GPU in this program will run

/** Formula mass in g/mol for one formula unit, or null if an element is unknown. */
function wyFormulaMass(formula) {
    let M = 0;
    for (const c of parseFormula(formula, 1)) {
        const ad = WY_ATOM_DATA[String(c.element).toUpperCase()];
        if (!ad || !Number.isFinite(ad.m)) return null;
        M += ad.m * c.count;
    }
    return M > 0 ? M : null;
}

/** Crystallographic density in g/cm^3. 1 amu / A^3 = 1.66054 g/cm^3. */
function wyDensity(Z, formulaMass, volume) {
    return (volume > 0) ? 1.66053906660 * Z * formulaMass / volume : NaN;
}

/**
 * Coordination requirements read off the distance windows, in the same way
 * runWyckoffSearch reads them -- so the pre-filter and the search agree about
 * which assignments exist.
 */
function wyCoordinationFromWindows(windows) {
    const coordination = {};
    for (const w of (windows || [])) {
        if (!w || !Number.isFinite(w.count) || w.count <= 0) continue;
        const key = String(w.a).toUpperCase();
        if (!coordination[key] || w.count > coordination[key].count) {
            coordination[key] = { count: w.count, geometry: w.geometry || null };
        }
    }
    return coordination;
}

/**
 * Does this Z admit at least one Wyckoff assignment?
 *
 * The same call the search makes, with the same caps, so a candidate that
 * passes here cannot fail there for a combinatorial reason.
 *
 * @returns {{ok: boolean, why: string, nAtoms: number}}
 */
function wyZFeasible(Z, job, symops) {
    let comp;
    try { comp = parseFormula(job.targetComposition, Z); }
    catch (e) { return { ok: false, why: (e && e.message) || String(e), nAtoms: 0 }; }

    const demand = [];
    let nAtoms = 0;
    for (const c of comp) {
        if (!Number.isInteger(c.count)) {
            return { ok: false, why: `${c.element} would need ${c.count} atoms, not a whole number`,
                     nAtoms: 0 };
        }
        const ad = WY_ATOM_DATA[String(c.element).toUpperCase()];
        if (!ad) return { ok: false, why: `no atomic data for ${c.element}`, nAtoms: 0 };
        demand.push({ element: c.element, count: c.count, z: ad.z, r: ad.rc ?? ad.r });
        nAtoms += c.count;
    }
    if (nAtoms > WY_MAX_ATOMS_PER_CELL) {
        return { ok: false, why: `${nAtoms} atoms per cell is beyond what the kernel can hold`,
                 nAtoms };
    }

    let en;
    try {
        en = enumerateAssignments(job.wyckoffTable || [], demand, {
            // The caps the SEARCH will use, so a Z that passes here cannot
            // fail there for a combinatorial reason. Passing undefined is
            // deliberate: enumerateAssignments then applies its own defaults,
            // which is exactly what runWyckoffSearch does with them.
            maxSites: job.wyckoffMaxSites || undefined,
            maxRepeat: job.wyckoffMaxRepeat || undefined,
            ceiling: 24,
            coordination: wyCoordinationFromWindows(job.distanceConstraints),
            symOps: symops,
            // ONE ASSIGNMENT IS THE WHOLE QUESTION. This is a feasibility
            // test, not the search's enumeration, and the recursion stops as
            // soon as the cap is reached -- so a scan over twenty-odd
            // candidate Z costs about what one enumeration costs, instead of
            // building up to 4000 assignments per candidate and discarding
            // all but the count.
            maxTotal: 1
        });
    } catch (e) {
        return { ok: false, why: (e && e.message) || String(e), nAtoms };
    }
    if (en.error || !en.assignments || !en.assignments.length) {
        // First clause only: the enumeration's own message is a paragraph
        // aimed at a user who has already chosen a Z, and the scan prints one
        // line per candidate. The full message is still what a failed RUN
        // reports.
        const why = (en.error || 'no assignment matches the composition').split('. ')[0];
        return { ok: false, why, nAtoms };
    }
    return { ok: true, why: '', nAtoms };
}

/**
 * The ranked list of Z values worth trying, best first.
 *
 * @returns {{list: Array, tried: Array, mass: number|null}}
 *   list  [{Z, rho, vpa, nAtoms, favoured, score}]  feasible only
 *   tried every candidate with its verdict, for the log
 */
function wyAutoZ(cell, symops, job) {
    const vol = cellVolume(cell);
    const mass = wyFormulaMass(job.targetComposition);
    const nOps = (symops || []).length;

    // The Z window from the density bounds. Without a mass -- an element the
    // table does not carry -- fall back to a volume-per-atom window, which is
    // the old rule but used as a SCAN RANGE rather than as an answer.
    let zLo, zHi;
    if (mass && vol > 0) {
        zLo = Math.max(1, Math.floor(vol * WY_RHO_MIN / (1.66053906660 * mass)));
        zHi = Math.min(WY_Z_HARD_MAX, Math.ceil(vol * WY_RHO_MAX / (1.66053906660 * mass)));
    } else {
        const perFormula = parseFormula(job.targetComposition, 1)
            .reduce((n, c) => n + c.count, 0) || 1;
        zLo = 1;
        zHi = Math.min(WY_Z_HARD_MAX, Math.ceil(vol / (6 * perFormula)));
    }
    if (!(zHi >= zLo)) { zLo = 1; zHi = Math.max(1, zHi); }

    const perFormulaAtoms = parseFormula(job.targetComposition, 1)
        .reduce((n, c) => n + c.count, 0) || 1;

    const tried = [], list = [];
    for (let Z = Math.max(1, zLo); Z <= zHi; Z++) {
        const rho = mass ? wyDensity(Z, mass, vol) : NaN;
        const vpa = vol / (Z * perFormulaAtoms);
        const f = wyZFeasible(Z, job, symops);
        tried.push({ Z, rho, vpa, ok: f.ok, why: f.why });
        if (!f.ok) continue;
        // Plausibility. Volume per atom is the transferable quantity -- a
        // density of 6.6 is unremarkable for a pyrochlore and impossible for a
        // zeolite, while 12 A^3 per atom is dense-but-normal for both. The
        // symmetry term is a nudge, not a rule: it is worth about a 10%
        // difference in volume per atom, so it breaks ties and nothing else.
        const favoured = !nOps || nOps % Z === 0 || Z % nOps === 0;
        const score = Math.abs(Math.log(vpa / WY_VPA_TYPICAL)) - (favoured ? 0.10 : 0);
        list.push({ Z, rho, vpa, nAtoms: f.nAtoms, favoured, score });
    }
    list.sort((a, b) => a.score - b.score);
    return { list, tried, mass, zLo, zHi, volume: vol };
}

/**
 * The Wyckoff search. Space group + cell + composition + |Fobs|^2, and nothing
 * else - in particular no electron-density map, by design.
 *
 * There is deliberately no map parameter. Its output is a site list and a
 * CIF; it is not a stage of charge flipping and does not consume one.
 */
async function runWyckoffFromIntensities(cell, symops, job) {

    const formula = job.targetComposition;
    if (!formula) return { sites: [] };

    let device;
    try {
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) throw new Error('WebGPU unavailable in worker.');
        device = await adapter.requestDevice();
    } catch (e) {
        // Returned, not just logged. This is the same swallowing bug as the one
        // around the search below: an empty site list reached the host, which
        // reported "produced no sites", and the actual cause -- no WebGPU --
        // was visible only in the console.
        console.error('[Wyckoff] WebGPU error:', e);
        return { sites: [], error:
            `The Wyckoff search needs WebGPU and could not get a device: ` +
            `${(e && e.message) || e}. In Firefox and older Safari it may need ` +
            `enabling; on Linux it may need a supported driver.` };
    }

    // ONE OBJECTIVE: |Fcalc|^2 against the Pawley intensities. There is no
    // density alternative: scoring a trial by point-sampling a map is both the
    // coupling this method must not have and a biased functional in its own
    // right, since a Z-weighted point sample ranks a heavy atom BELOW a light
    // one -- the heavy atom's own series-termination ripples dig a trough under
    // it (the effect integrateSphere(), over in the charge-flipping worker,
    // was written to work around).
    const shaderFile = '../swarm_reflection.wgsl';
    const shaderResp = await fetch(shaderFile);
    if (!shaderResp.ok) throw new Error(shaderFile + ' could not be loaded.');
    const swarmShaderSrc = await shaderResp.text();

    // Raw observation rows for the reflection objective. runWyckoffSearch
    // normalises them itself and estimates the Wilson B from them, so they are
    // handed over unprocessed rather than pre-reduced.
    let rawRows = [];
    {
        rawRows = wyRawObservationRows(job);
        if (!rawRows.length) {
            return { sites: [], error: 'No usable Pawley intensities for the Wyckoff search.' };
        }
    }

    // Loaded for BOTH objectives, not just the mapless one, because the FINAL
    // REFINEMENT needs them too and could not get them: job.formFactor is a
    // function, and a function cannot cross a postMessage boundary, so
    // coord_refine's `o.formFactor || (() => null)` fell through to f = Z on
    // every single run.
    //
    // That is the discrepancy behind a search reporting CC = 0.99 and the
    // refinement then reporting R = 50%: THEY WERE NOT SCORING THE SAME MODEL.
    // The search used tabulated f(s) with its angular fall-off; the refinement
    // used a flat atomic number, which overweights every heavy atom at high
    // angle. No arrangement of atoms satisfies both descriptions at once, so
    // the refinement moved the structure away from the search's answer and
    // still could not fit.
    //
    // '../scatters', not the default 'scatters': a relative fetch in a worker
    // resolves against the WORKER script's URL, and this worker lives in js/.
    let scatterTables = null;
    try {
        scatterTables = await loadScatteringTables('../scatters');
    } catch (e) {
        console.warn('[Wyckoff] No scattering tables; falling back to f = Z.', e);
    }
    const refineFormFactor = (scatterTables && typeof makeFormFactor === 'function')
        ? makeFormFactor(scatterTables, { radiation: job.radiation || 'xray', ions: job.ions })
        : null;

    const stg = { sym_ops: symops, wyckoff: job.wyckoffTable || [] };

    // ----------------------------------------------------------------------
    //  Z. Either the user's number, or the scan -- see wyAutoZ above for why
    //  the old single-guess route could not work for a pyrochlore.
    // ----------------------------------------------------------------------
    const userZ = Number(job.wyckoffZ);
    const autoZ = !(Number.isFinite(userZ) && userZ >= 1);
    let zQueue;
    if (!autoZ) {
        zQueue = [Math.round(userZ)];
    } else {
        const scan = wyAutoZ(cell, symops, job);
        const vLine = `Cell volume ${scan.volume.toFixed(1)} A^3` +
                      (scan.mass ? `, formula mass ${scan.mass.toFixed(2)} g/mol` : '');
        wySay(`Z = Auto. ${vLine}. Scanning Z = ${scan.zLo}-${scan.zHi} ` +
              `(the range that keeps the density between ${WY_RHO_MIN} and ` +
              `${WY_RHO_MAX} g/cm^3).`);
        // Every candidate, with its verdict. A rejected Z is as informative as
        // an accepted one -- it is usually the arithmetic that decides, and
        // seeing WHICH arithmetic is what tells the user their formula, not
        // their Z, is the thing to change.
        for (const t of scan.tried) {
            wySay(`  Z = ${String(t.Z).padStart(3)}  ` +
                  `rho = ${Number.isFinite(t.rho) ? t.rho.toFixed(2).padStart(6) : '   ? '} g/cm^3  ` +
                  `V/atom = ${t.vpa.toFixed(1).padStart(5)} A^3  ` +
                  (t.ok ? 'assignments exist' : `rejected: ${t.why}`));
        }
        if (!scan.list.length) {
            // WHEN EVERY Z FAILS FOR THE SAME NON-ARITHMETIC REASON, THAT
            // REASON IS THE ANSWER.
            //
            // A broken Wyckoff table -- a position with no n_free, say -- gets
            // rejected at every single candidate, and blaming the composition
            // for it would send the user off to re-check a formula that was
            // never the problem. Only a reason that varies with Z is about Z.
            const reasons = [...new Set(scan.tried.filter(t => !t.ok).map(t => t.why))];
            const arithmetic = /atoms per cell|whole number|assignment|multiplicit/i;
            const arith = reasons.filter(r => arithmetic.test(r));
            const other = reasons.filter(r => !arithmetic.test(r));
            // NOT "are the other reasons all identical" -- a broken table names
            // a different position at each Z ("8a has no n_free", "16c has no
            // n_free"), so counting distinct strings would miss exactly the
            // case this is for. What matters is that NO candidate failed on the
            // arithmetic: if none did, Z was never the variable under test.
            if (!arith.length && other.length) {
                return { sites: [], error:
                    `The Wyckoff search could not start, and Z is not the reason -- ` +
                    `every candidate from ${scan.zLo} to ${scan.zHi} failed for the same ` +
                    `kind of reason: ${other.slice(0, 2).join(' / ')}` +
                    (other.length > 2 ? ` (and ${other.length - 2} more like it).` : '') };
            }
            return { sites: [], error:
                `No value of Z between ${scan.zLo} and ${scan.zHi} -- the whole range that ` +
                `gives a density of ${WY_RHO_MIN} to ${WY_RHO_MAX} g/cm^3 for this cell -- ` +
                `produces a legal Wyckoff assignment for ${formula}. That is a statement ` +
                `about the FORMULA and the SPACE GROUP, not about Z: check that the formula ` +
                `is the formula unit (not the cell contents), that the space group is right, ` +
                `and that any coordination counts in the distance windows are achievable.` };
        }
        // Best first, then the rest as fallbacks. Three is enough: past that
        // the density is implausible and a further attempt costs a full search
        // to tell the user the same thing.
        zQueue = scan.list.slice(0, 3).map(c => c.Z);
        const b = scan.list[0];
        wySay(`Chosen Z = ${b.Z} (density ${Number.isFinite(b.rho) ? b.rho.toFixed(2) : '?'} g/cm^3, ` +
              `${b.vpa.toFixed(1)} A^3 per atom, ${b.nAtoms} atoms per cell)` +
              (zQueue.length > 1 ? `; fallbacks ${zQueue.slice(1).join(', ')}.` : '.'));
        postMessage({ type: 'wy-info', message:
            `Z = ${b.Z} chosen automatically (density ${Number.isFinite(b.rho) ? b.rho.toFixed(2) : '?'} g/cm^3).` });
    }
    let Z = zQueue[0];


    // "Compiling WebGPU Swarm..." described an internal implementation step to
    // a user who has no way to act on it, and it stayed on screen as the last
    // thing said until the first progress message arrived -- so an idle panel
    // read as if a shader compile were stuck. The phase is now named for what
    // the program is doing on the user's behalf.
    postBoth({ type: 'wy-start', message: 'Preparing\u2026' });

    let out = null, lastErr = null;
    try {
      for (let attempt = 0; attempt < zQueue.length; attempt++) {
        Z = zQueue[attempt];
        if (attempt > 0) {
            wySay(`Z = ${zQueue[attempt - 1]} failed in the search ` +
                  `(${(lastErr && lastErr.message) || 'unknown reason'}). Retrying with Z = ${Z}.`);
            postMessage({ type: 'wy-info', message: `Retrying with Z = ${Z}.` });
        }
        try {
        out = await runWyckoffSearch({
            device, 
            ccShaderSource: swarmShaderSrc,
            // No groupStride here: swPackReflections decides it, and
            // hard-coding 3 would have overridden the weighted layout
            // it now emits, leaving the kernel reading Iobs where the
            // weight should be.

            setting: stg,
            cell,
            reflections: rawRows,
            wavelength: job.wavelength,
            radiation: job.radiation || 'xray',
            polarisationK: polarizationK((job.polarization !== undefined && job.polarization !== null)
                                          ? job.polarization : job.monochromatorTth),
            // No gridSize. It described the edge of the density grid, which
            // runWyckoffSearch read only in its density mode; with that mode
            // gone the option has no meaning, and the local N it referred to
            // went with the map parameter -- so this line was a ReferenceError
            // waiting for the first run.
            formula, 
            Z: Z,
            maxSitesPerElement: job.wyckoffMaxSites || undefined,
            maxRepeat: job.wyckoffMaxRepeat || undefined,
            wyckoffCapCeiling: 24,
            atomData: WY_ATOM_DATA, 
            scatterTables: scatterTables,
            // Accept the pre-parsed array directly from the UI
            windows: job.distanceConstraints || [],
            harkerSites: [],
            numParticles: job.swarmParticles || 512,
            generations: job.swarmIterations || 1000,
            restarts: job.swarmRestarts || 4,
            minContact: job.minContact || 1.0,
onLog: m => wySay(m),
            onProgress: p => {
                postBoth({ 
                    type: 'wy-progress', 
                    wave: p.wave, waves: p.waves, 
                    generation: p.generation, generations: p.generations, 
                    restart: p.restart, restarts: p.restarts, 
                    best: p.best 
                });
            },
            shouldStop: () => wyckoffStopRequested
        });
        break;                       // this Z worked
        } catch (e) {
            lastErr = e;
            // A COMBINATORIAL FAILURE IS A STATEMENT ABOUT Z, AND NOTHING
            // ELSE. On auto it is a reason to try the next candidate; on a Z
            // the user typed it is the answer to their question, and swapping
            // in a different Z behind their back would be worse than failing.
            if (autoZ && attempt < zQueue.length - 1 &&
                !wyckoffStopRequested && wyRetryableZError(e)) continue;
            throw e;
        }
      }
    } catch (e) {
        console.error('[Wyckoff] the search failed:', e);
        // The message was going to the console and NOWHERE ELSE: the caller
        // saw an empty site list and reported "produced no sites", which is
        // how a plain arithmetic complaint about Z ended up looking like a
        // crash. It is returned now, and the host shows it.
        return { sites: [], error: (e && e.message) || String(e) };
    } finally {
        if (device) device.destroy();
    }
    if (!out) return { sites: [], error: (lastErr && lastErr.message) || 'The Wyckoff search returned nothing.' };

if (out && out.candidates && out.candidates.length > 0) {
        const G = metricTensor(cell);

        // ------------------------------------------------------------------
        //  A CONTACT FLOOR IS A CONSTRAINT, NOT A PREFERENCE.
        //
        //  The kernel charges a clash as a soft penalty, scaled by a ramp that
        //  starts gentle so the swarm can move through crowded configurations
        //  early. That is the right behaviour DURING the search. It is the
        //  wrong behaviour for the answer: a solution containing a 0.47 A Pb-O
        //  contact is not a slightly worse structure, it is not a structure,
        //  and no correlation coefficient should be able to buy it.
        //
        //  runWyckoffSearch returns `floors` precisely so the caller can reject
        //  on the same numbers the search enforced -- and nothing called it.
        //  Candidates are now walked in score order and the first one whose
        //  contacts are all legal is taken. Rejections are reported rather than
        //  silently dropped: if every candidate fails, that is a statement
        //  about the composition or the floor, and the user needs to see it.
        // ------------------------------------------------------------------
        let cfWorstContact = null;
        const floorOf = (typeof out.floors === 'function')
            ? out.floors
            : () => (Number.isFinite(job.minContact) ? job.minContact : 1.0);

        const worstContact = (sites) => {
            // Expand to the full cell, then take the shortest unlike-pair
            // distance as a fraction of that pair's floor. < 1 means illegal.
            const exp = [];
            for (const st of sites || []) {
                for (const op of symops) {
                    // applyOp, NOT op.apply. The latter is the RECIPROCAL
                    // transform used for hkl, which is the transpose of the
                    // rotation and carries no translation -- correct for
                    // indices, silently wrong for coordinates.
                    const p = applyOp(op, [st.x, st.y, st.z])
                        .map(v => ((v % 1) + 1) % 1);
                    if (!exp.some(q => q.el === st.element &&
                                       fracDistance(G, q.p[0] - p[0], q.p[1] - p[1], q.p[2] - p[2]) < 1e-3)) {
                        exp.push({ el: st.element, p });
                    }
                }
            }
            // fracDistance handles the minimum image itself, including the
            // 27-point search an oblique cell needs. An earlier version did
            // that search HERE, passing dx+1 into fracDistance -- which
            // immediately rounded it back, making the whole loop a no-op and
            // the test that 'verified' it meaningless, since both sides were
            // wrapping.
            let worst = Infinity, pair = null;
            for (let i = 0; i < exp.length; i++) {
                for (let j = i + 1; j < exp.length; j++) {
                    let dx = exp[i].p[0] - exp[j].p[0];
                    let dy = exp[i].p[1] - exp[j].p[1];
                    let dz = exp[i].p[2] - exp[j].p[2];
                    dx -= Math.round(dx); dy -= Math.round(dy); dz -= Math.round(dz);
                    const d = fracDistance(G, dx, dy, dz);
                    const fl = floorOf(exp[i].el, exp[j].el);
                    if (!(fl > 0)) continue;
                    const ratio = d / fl;
                    if (ratio < worst) {
                        worst = ratio;
                        pair = { a: exp[i].el, b: exp[j].el, d, floor: fl };
                    }
                }
            }
            return { worst, pair };
        };

        // Kept in scope: the refinement below has to be judged by the SAME
        // test the candidates were, or the filter only guards the seed.
        cfWorstContact = worstContact;

        let best = null;
        const contactRejects = [];
        for (const cand of out.candidates) {
            const chk = worstContact(cand.sites);
            if (chk.worst >= 1.0 || !chk.pair) { best = cand; break; }
            contactRejects.push(
                `${cand.assignment}: ${chk.pair.a}-${chk.pair.b} at ${chk.pair.d.toFixed(2)} A ` +
                `against a floor of ${chk.pair.floor.toFixed(2)} A`);
        }
        if (!best) {
            return {
                sites: [],
                error: `Every candidate breaks a minimum-contact distance. ` +
                       `Closest offenders: ${contactRejects.slice(0, 3).join('; ')}. ` +
                       `Either the composition or Z is wrong for this cell, or the contact ` +
                       `floor is set higher than the structure allows.`
            };
        }
        if (contactRejects.length) {
            console.warn(`[Wyckoff] ${contactRejects.length} higher-scoring candidate(s) ` +
                         `rejected on contact distance: ${contactRejects.slice(0, 3).join('; ')}`);
        }

        // ------------------------------------------------------------------
        // FINAL COORDINATE REFINEMENT.
        //
        // The search above scores |Fcalc|^2 against the Pawley intensities as
        // a correlation, which is scale-free and insensitive to a uniform
        // error but says nothing about how well any INDIVIDUAL reflection is
        // reproduced. This step refines the free coordinates against |Fo|^2 by
        // least squares, with every site re-projected onto its Wyckoff
        // position at each step, and reports wR(F^2) -- a figure the
        // correlation cannot give.
        //
        // (This preamble used to describe grid quantisation from trilinear
        // sampling of a density map. That objective is gone: this route never
        // sees a map, so its coordinates were never quantised to 1/N.)
        //
        // Failure here is not fatal: the unrefined coordinates are still the
        // answer the search gave, so the reason is recorded and the pipeline
        // carries on.
        // ------------------------------------------------------------------
        // THE SEARCH RESULT, before the refinement is allowed to touch it.
        //
        // The refinement mutates best.sites in place, so without this copy the
        // Monte Carlo answer is unrecoverable by the time anything is
        // reported -- and the two are worth comparing: they optimise different
        // quantities, and when they disagree that disagreement is the
        // interesting part.
        const searchSites = best.sites.map(x => ({ ...x }));

        let refinement = null;
        try {
            const obs = wyObservationsFromHkl(job, symops);
            const obsRows = obs.rows || [];
            for (const w of (obs.warnings || [])) console.warn('[Wyckoff refine] ' + w);
            for (const e of (obs.errors || [])) console.warn('[Wyckoff refine] ' + e);
            if (obsRows.length) {
                const seed = best.sites.map(s => ({
                    element: s.element,
                    zn: Number.isFinite(s.zn)
                        ? s.zn
                        : (WY_ATOM_DATA[String(s.element).toUpperCase()] || {}).z,
                    x: s.x, y: s.y, z: s.z, w: s.w
                }));
                if (seed.every(s => Number.isFinite(s.zn) && s.w)) {
                    // out.overallB is the WILSON B the search itself used to
                    // build its scattering table. job.overallB is only set if
                    // the user typed one, so this was refining with B = 0 on
                    // every ordinary run -- fitting sharp point atoms to data
                    // that falls off, which loads the error onto the light
                    // atoms first. The search and the refinement must use the
                    // same B or they are not fitting the same model.
                    const refB = Number.isFinite(job.overallB) ? job.overallB
                               : (Number.isFinite(out.overallB) ? out.overallB : 0);
                    const r = refineCoordinatesAgainstPawley({
                        sites: seed, symOps: symops, obsRows,
                        overallB: refB, formFactor: job.formFactor || refineFormFactor
                    });

                    // ------------------------------------------------------
                    //  THE REFINEMENT IS SUBJECT TO THE CONTACT FLOOR TOO.
                    //
                    //  The candidate filter runs on the SEARCH output, and the
                    //  refinement then moves every free coordinate with no
                    //  contact term at all -- so it could and did walk a
                    //  legal seed into an illegal structure: a candidate that
                    //  passed at 1.4 A came out with S-O at 0.85 A against a
                    //  floor of 1.35, because shortening that contact happened
                    //  to lower R. Guarding only the input to a step that
                    //  moves the atoms guards nothing.
                    //
                    //  A refinement that breaks a floor is discarded rather
                    //  than clamped. Clamping would leave the coordinates at
                    //  an arbitrary point on the constraint surface and report
                    //  an R that belongs to neither position; the seed is at
                    //  least a structure the search actually scored.
                    // ------------------------------------------------------
                    let contactBreak = null;
                    if (r && r.sites && !r.error && cfWorstContact) {
                        const chk = cfWorstContact(r.sites);
                        if (chk.pair && chk.worst < 1.0) contactBreak = chk.pair;
                    }

                    if (contactBreak) {
                        // Not a `return`: this runs inside the try block of
                        // runWyckoffFromIntensities, and returning here would skip
                        // the site formatting and hand back an empty result.
                        console.warn(`[Wyckoff refine] Refinement rejected: it produced ` +
                            `${contactBreak.a}-${contactBreak.b} at ${contactBreak.d.toFixed(2)} A ` +
                            `against a floor of ${contactBreak.floor.toFixed(2)} A. ` +
                            `Keeping the unrefined positions.`);
                        refinement = {
                            skipped: `the refined coordinates broke the ${contactBreak.a}-${contactBreak.b} ` +
                                     `contact floor (${contactBreak.d.toFixed(2)} A against ` +
                                     `${contactBreak.floor.toFixed(2)} A), so the search positions were kept`,
                            rejectedR: r.R
                        };
                    } else if (!r.error) {
                        r.sites.forEach((rs, i) => {
                            best.sites[i].x = rs.x;
                            best.sites[i].y = rs.y;
                            best.sites[i].z = rs.z;
                        });
                        refinement = {
                            // Which normalisation route the observations took,
                            // so a surprising R can be traced to how |Fo|^2 was
                            // obtained rather than guessed at.
                            route: obs.route, lpApplied: job.applyLP !== false,
                            R: r.R, Rstart: r.Rstart,
                            iterations: r.iterations, converged: r.converged,
                            nObs: r.nObs, nParams: r.nParams,
                            // How far each site actually moved, in fractional
                            // units. A shift of about 1/N is the grid bias
                            // being taken out; a much larger one means the map
                            // and the intensities disagree about that site and
                            // is worth looking at rather than accepting.
                            shifts: r.shifts
                        };
                    } else {
                        refinement = { skipped: r.error };
                    }
                } else {
                    refinement = { skipped: 'Sites lack an atomic number or a Wyckoff position.' };
                }
            } else {
                refinement = { skipped: (obs.errors && obs.errors[0]) ||
                                        'No usable Pawley intensities were supplied.' };
            }
        } catch (err) {
            refinement = { skipped: (err && err.message) || String(err) };
        }

        // Minimum image over the 27 neighbouring translations. The delta is
        // wrapped to the nearest cell first, which is enough on its own for a
        // near-orthogonal cell; the shell search covers the oblique ones,
        // where the nearest image is not always the wrapped one.
        const minImage = (dx, dy, dz) => {
            dx -= Math.round(dx); dy -= Math.round(dy); dz -= Math.round(dz);
            let best = Infinity;
            for (let i = -1; i <= 1; i++)
                for (let j = -1; j <= 1; j++)
                    for (let k = -1; k <= 1; k++) {
                        const d = fracDistance(G, dx + i, dy + j, dz + k);
                        if (d < best) best = d;
                    }
            return best;
        };

        // Expand the asymmetric unit to the whole cell once, then measure
        // every unique site against it.
        const expanded = [];
        best.sites.forEach((s, si) => {
            const imgs = [];
            for (const op of symops) {
                const q = applyOp(op, [s.x, s.y, s.z]).map(v => ((v % 1) + 1) % 1);
                if (!imgs.some(im => minImage(im[0] - q[0], im[1] - q[1], im[2] - q[2]) < 1e-3)) {
                    imgs.push(q);
                }
            }
            for (const q of imgs) {
                expanded.push({ element: s.element, siteIdx: si, x: q[0], y: q[1], z: q[2] });
            }
        });

        const formattedSites = best.sites.map((s, i) => {
            // THE DENSITY ACTUALLY UNDER THE SITE.
            //
            // This was hardcoded to 1.0, so every site claimed a full peak
            // whether it sat on one or on vacuum, and the one column that
            // would have shown a site fitted to empty map read the same as a
            // site fitted to the heaviest atom in the cell. The map arrives
            // scaled so its maximum is 1, so this is directly comparable
            // between sites and between runs.
            // No map, no map density. Reporting 0 would read as "sits on
            // nothing" rather than "was never measured".
            // NO MAP ON THIS ROUTE, SO NO MAP DENSITY. Explicitly null rather
            // than absent: the site table can then show "not applicable"
            // instead of a number the method never had access to. Filling this
            // in would mean handing the Wyckoff search the charge-flipping map
            // to score itself against, which is the coupling being removed.
            const mapDensity = null;

            // Nearest neighbours, as measured - not as assumed. Nothing here
            // asks what coordination the site OUGHT to have; it reports what
            // the fitted structure has, so a missing bond or an impossible
            // contact is visible in the results table instead of having to be
            // computed by hand afterwards.
            const self = expanded.find(a => a.siteIdx === i);
            const near = [];
            for (const b of expanded) {
                if (b === self) continue;
                const d = minImage(self.x - b.x, self.y - b.y, self.z - b.z);
                if (d > 1e-3) near.push({ element: b.element, distance: d });
            }
            near.sort((p, q) => p.distance - q.distance);

            return {
                rank: i + 1,
                element: s.element,
                x: s.x,
                y: s.y,
                z: s.z,
                multiplicity: s.multiplicity,
                wyckoff: s.wyckoff,
                special: s.multiplicity < symops.length,
                // TWO DIFFERENT NUMBERS, BOTH USEFUL, NEITHER A PEAK SEARCH.
                //
                //   mapDensity  the density under this site as a fraction of
                //               the MAP's maximum. An absolute statement: 0.93
                //               means this site sits on something nearly as
                //               strong as the strongest feature in the cell,
                //               and a value at or below zero means it sits on
                //               nothing at all - the map had no opinion about
                //               it and it was placed by composition and
                //               penalties alone.
                //
                //   relative    the same sample against the strongest FITTED
                //               site. A ranking among the sites actually
                //               proposed, which is what the peak-picking path
                //               reports, so the two tables mean the same
                //               thing. Its top entry is 1.000 by construction
                //               and therefore says nothing on its own.
                //
                // `height` is kept as an alias of mapDensity so existing
                // callers keep working; new code should read mapDensity and
                // say which of the two it is showing.
                mapDensity,
                height: mapDensity,
                contacts: near.slice(0, 6),
                // The shortest contact of all, which is the single number that
                // says whether this site is physically possible.
                shortestContact: near.length ? near[0].distance : null
            };
        });

        // `relative` is each site's height against the tallest FITTED site,
        // matching what the peak-picking path reports, so the two routes'
        // tables mean the same thing.
        const topH = formattedSites.reduce(
            (m, s) => Number.isFinite(s.mapDensity) ? Math.max(m, s.mapDensity) : m, 0);
        for (const s of formattedSites) {
            s.relative = (topH > 0 && Number.isFinite(s.mapDensity)) ? s.mapDensity / topH : null;
        }

        return {
            sites: formattedSites,
            // Carried out of runWyckoffSearch so the message handler can label
            // the result. The search returns its per-assignment global bests
            // whether or not it ran to term, so a stopped run still has sites.
            stopped: !!(out && out.stopped),
            assignment: best.assignment,
            candidates: out.candidates,
            refinement,
            // The search's own numbers, kept separate from the refinement's.
            // searchCC is the FULL-RESOLUTION correlation from the quench, not
            // the ramped figure the progress readout shows.
            searchSites,
            searchCC: Number.isFinite(best.cc) ? best.cc : null,
            searchScore: Number.isFinite(best.score) ? best.score : null,
            z: Z
        };
    }
    
    return { sites: [] };
}
// ===========================================================================
//  MESSAGE HANDLING
// ===========================================================================

// Set by 'cancel-wyckoff', read by runWyckoffSearch through shouldStop().
// A stop is co-operative: the swarm keeps a per-assignment global best, which
// IS a structure, and killing the worker would throw it away.
let wyckoffStopRequested = false;

/**
 * One post, one message.
 *
 * An earlier draft of this file also emitted the old 'cf-*' names so that an
 * un-updated host would keep working. That is a trap: powder5.html listens for
 * BOTH, so every error arrived twice, ran done() twice and raised two toasts.
 * The host is updated in the same change as this file, so the names are simply
 * the new ones.
 */
function postBoth(msg) { postMessage(msg); }

/** Progress and commentary from the search, to the console and the host. */
function wySay(m) {
    console.log('[Wyckoff]', m);
    postMessage({ type: 'wy-log', message: String(m) });
}

/**
 * Is this failure one that a different Z would fix?
 *
 * Only these. A GPU error, a missing shader or a bad cell would fail the same
 * way for every Z, and retrying would just take three times as long to say so.
 */
function wyRetryableZError(e) {
    const m = String((e && e.message) || e || '');
    return /no wyckoff assignment|not a whole number|atoms per cell|assignment matches|multiplicit|workgroup budget/i.test(m);
}

globalThis.onmessage = async (e) => {
    const job = e.data || {};

    // Must not fall through into the switch and must not yield: this arrives
    // WHILE a search is running.
    if (job.type === 'cancel-wyckoff') { wyckoffStopRequested = true; return; }

    if (job.type !== 'wyckoff-search') return;

    try {
        const symops = job.symops || [];
        if (!symops.length) {
            postBoth({ type: 'wy-error', message: 'No symmetry operators for the chosen setting.' });
            return;
        }
        if (!job.targetComposition) {
            postBoth({ type: 'wy-error', message: 'A target composition is required for a Wyckoff search.' });
            return;
        }
        if (!job.wyckoffTable || !job.wyckoffTable.length) {
            postBoth({ type: 'wy-error', message: 'This setting carries no Wyckoff table.' });
            return;
        }

        wyckoffStopRequested = false;
        const res = await runWyckoffFromIntensities(job.cell, symops, job);

        if (res.error) { postBoth({ type: 'wy-error', message: res.error }); return; }
        if (!res.sites || !res.sites.length) {
            postBoth({ type: 'wy-error', message: res.stopped
                ? 'Stopped before any assignment produced a finite score, so there is nothing to show.'
                : 'The Wyckoff search produced no sites.' });
            return;
        }

        postBoth({
            type: 'wy-structure',
            results: {
                sites: res.sites,
                wyckoffResult: res,
                refinement: res.refinement,
                cell: job.cell,
                nOps: symops.length,
                symops: symops.map(o => o.xyz).filter(Boolean),
                rawPeakCount: res.sites.length,
                // No map on this route, by design. Explicit nulls so the UI can
                // tell "not applicable" from "zero".
                map: null, gridSize: null,
                originShift: null, hand: null,
                symmetryCorrelation: null, originScore: null,
                source: 'wyckoff',
                z: res.z,
                stopped: !!res.stopped
            }
        });
    } catch (err) {
        postBoth({ type: 'wy-error', message: (err && err.message) || String(err) });
    }
};
