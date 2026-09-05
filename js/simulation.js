/* ===========================================================================
   SIMULATION  --  a calculated powder pattern from a structural model.
   ===========================================================================

   WHAT THIS IS

   Every other mode in powder5 starts from a measured pattern and works
   backwards: Le Bail and Pawley extract intensities from it, charge flipping
   and the Wyckoff search turn those intensities into a structure. This module
   runs the same physics FORWARDS -- atoms in, pattern out -- and needs no data
   file at all.

   WHAT IT DELIBERATELY DOES NOT DO

   It does not contain a profile function, a background, a Lorentz-polarisation
   factor, an hkl generator or a peak-area integral. Every one of those already
   exists in this program, is already used by the fit, and is called here
   unchanged:

     generateAndCacheHklIndices()  reflections + systematic absences   (powder5.html)
     updateHklPositions()          2-theta and d                        (crystal.js)
     reflectionHeightsFromFsq()    |F|^2 -> peak height                 (data_io.js)
     calculatePatternCPU()         heights -> profile                   (profile.js)
     calculateTotalBackground()    background                           (profile.js)

   That is the whole point. The header of profile.js records what happened the
   last time this program had two implementations of the pattern: "the preview
   drew one pattern and the fit minimised another, and every report would have
   been self-consistent about the wrong answer." A simulation that used its own
   profile would be a third such copy, and the failure would be worse than a
   drifting fit -- a simulated pattern that cannot be compared with a measured
   one is not a simulation of anything.

   So the ONLY new physics here is the structure factor:

       F(hkl) = sum_j  occ_j . f_j(s) . exp(-B_j s^2 / 4) . exp(2 pi i h.r_j)

   summed over every atom in the unit cell, |F|^2 handed to
   reflectionHeightsFromFsq(), and the rest of the chain is the existing one.

   THE TWO THINGS TO GET RIGHT

   1. s versus stol. The Gaussian scattering coefficients are defined against
      stol = sin(theta)/lambda = s/2. scatterers.js's makeFormFactor() takes s
      and halves it internally, and this module does the same for the
      Debye-Waller factor: exp(-B stol^2) == exp(-B s^2 / 4). Getting it wrong
      leaves f(0) untouched, so it survives a casual check, while making every
      atom fall off at four times the correct rate.

   2. One tolerance, used twice. The atom orbit and the Wyckoff letter are
      computed with the SAME distance tolerance (SIM_MERGE_TOL_A). They have to
      be: if the orbit merged two images that the Wyckoff test called distinct,
      the panel would report a multiplicity the pattern was not computed with.
      A site on a mirror scatters half as strongly as one that is merely near
      it, so the disagreement would be a factor of two in that atom's
      contribution with nothing on screen to show for it.
   =========================================================================== */

(function (root) {
'use strict';

const SIM_VERSION = '1.0.0';

/* Distance below which two symmetry images are the SAME atom, in angstrom.
   Also the distance below which a site counts as sitting on a Wyckoff
   position. See note 2 in the header: these must be one number.

   0.02 A is well below any real interatomic distance, and comfortably above
   the rounding in a coordinate a user types by hand -- 1/3 entered as 0.3333
   is 3e-5 fractional out, which is 3e-4 A in a 10 A cell. */
const SIM_MERGE_TOL_A = 0.02;

/* Step of the calculated 2-theta axis when no data file is loaded, in degrees.
   With data the axis IS the data's, so the simulation lands on the same points
   as the measurement and the two can be subtracted; the step box is ignored. */
const SIM_STEP_MIN = 1e-4, SIM_STEP_MAX = 0.1, SIM_STEP_DEFAULT = 0.01;

/* Hard ceiling on the number of axis points, whatever step is asked for.
   0.0001 deg across 5-120 is 1.15 million points, and the profile is the most
   expensive part of a render -- that request would freeze the tab for the best
   part of a minute on every keystroke. The fine end of the step range is meant
   for a narrow window, where it fits under this cap comfortably; asking for it
   across the whole pattern is answered with the finest step that does fit, and
   the status line says which. */
const SIM_MAX_POINTS = 20000;

/* Flat background added under a calculated pattern, in the same units as the
   peaks. Only meaningful with no data file open: with one, the background is
   the measured one from the spline, and adding a synthetic constant on top
   would be inventing counts on a real measurement. */
const SIM_BKG_MIN = 0, SIM_BKG_MAX = 100, SIM_BKG_DEFAULT = 10;

/* COUNTING NOISE, IN UNITS OF THE COUNTING STATISTICS THEMSELVES.
   0 is a clean calculated curve. 1 is the scatter a real detector would give
   for that many counts. Above 1 is a noisier measurement at the same intensity
   scale -- a shorter count, a worse instrument. See simNoisy(). */
const SIM_NOISE_MIN = 0, SIM_NOISE_MAX = 10, SIM_NOISE_DEFAULT = 0;

/* Above this many reflections the Ka tick strip stops being information and
   starts being a denial of service on the renderer. Chart.js builds and lays
   out one element per bar, so a wide 2-theta range on a large cell hands it
   tens of thousands of them on every redraw -- and at that density the ticks
   are a solid band that tells the eye nothing anyway. Past the cap only the
   strongest are drawn, and the status line says so. */
const SIM_MAX_MARKERS = 2000;

/** Trailing debounce for edits that arrive in bursts (typing, slider drags). */
const SIM_DEBOUNCE_MS = 200;

const DEG2RAD_SIM = Math.PI / 180;

// ---------------------------------------------------------------------------
//  STATE
// ---------------------------------------------------------------------------

/** [{ id, label, x, y, z, occ, biso }] -- label is a key into the scattering
 *  table ("Fe", "Fe3+"), not necessarily a bare element symbol. */
let simAtoms = [];
let simNextId = 1;

/** From loadScatteringTables(); null until the first activation resolves. */
let simTables = null;
let simTablesError = null;
let simTablesPending = null;

/** Cache: sg identity -> { ops, wyckoff, settingSymbol, error }. */
const simSymCache = new Map();

/** Re-entrancy guard. render() writes to the DOM and to the chart, and
 *  anything that reacts to either must not start a second render inside the
 *  first. See simRender() for how a trigger arriving mid-render is handled. */
let simRendering = false;

/* THE REFLECTION LIST IS CACHED SEPARATELY FROM THE ATOMS.
   Which reflections exist, and where they sit in 2-theta, depends on the cell,
   the wavelengths, the zero point, the space group and the 2-theta range. It
   does NOT depend on the atoms. Editing an atom used to regenerate the whole
   list, re-sort it, and hand Chart.js several thousand freshly allocated bar
   objects for the tick strip -- all to arrive at exactly the geometry that was
   already on screen. */
let simGeomKey = null;
let simGeomList = null;
/** Geometry key the Ka marker datasets were last built for. */
let simMarkersKey = null;
/** Reused {x,y} array for the curve, so a redraw does not allocate one point
 *  object per axis step and leave the previous four thousand to the collector. */
let simPts = null;
/** The same curve WITHOUT counting noise. Kept so the export can write the
 *  noisy pattern as the observation and the clean one as the calculation --
 *  which is exactly the pair a Rietveld program needs to be tested against. */
let simCleanY = null;
let simNoiseApplied = false;

/** Counters for SIM.stats(), so a pathological session can be diagnosed from
 *  the console instead of guessed at. reentryStack holds where a render was
 *  re-triggered from, which is the one thing the counters alone cannot say. */
const simCounters = { renders: 0, geomBuilds: 0, markerBuilds: 0, lastMs: 0, maxMs: 0,
                      reentries: 0, reentryStack: null };

/** Session log. logEvent() lives inside powder5.html's DOMContentLoaded
 *  closure, so it is only reachable through window. */
function simLog(msg) {
    if (typeof root.logEvent === 'function') root.logEvent(msg);
}

/** One line describing an atom, for the log. */
function simAtomLine(a, i) {
    const wy = (() => {
        try {
            const sym = currentSG ? simSymmetry(currentSG) : null;
            const p = getAllParams();
            if (!sym || sym.error || !(p.a > 0)) return '';
            return ' ' + simWyckoffFor(a.x, a.y, a.z, sym, simMetric(p)).label;
        } catch (e) { return ''; }
    })();
    return `${simSiteName(a, i)} ${a.label}${wy} ` +
           `x=${simFmt(a.x, 5)} y=${simFmt(a.y, 5)} z=${simFmt(a.z, 5)} ` +
           `occ=${simFmt(a.occ, 3)} B=${simFmt(a.biso, 3)}`;
}

/** Set once the tab panel has been built and wired. */
let simUiReady = false;

/** Previous value of simActive(), so entering the mode can be told from
 *  staying in it. See simApplyMode(). */
let simWasActive = false;
/** Whether the log has an unmatched "entered Simulation" line. */
let simWasLogged = false;

/* WHICH ROWS ARE EXPANDED, BY ATOM ID.
   Kept here rather than as a property on the atom, for two reasons: the atom
   objects are what gets written to the CIF and read by SIM.atoms(), and a
   display flag has no business travelling with them; and the list is rebuilt
   wholesale on every add and remove, so the state has to outlive the DOM that
   showed it. */
const simOpenRows = new Set();

// ---------------------------------------------------------------------------
//  MODE
// ---------------------------------------------------------------------------

/**
 * True when the Method selector is on Simulation.
 *
 * Read from the DOM rather than cached: the selector is the single source of
 * truth for the mode, and a cached copy is one more thing that can be stale
 * when something else changes the value (a preset load, a restored session).
 */
function simActive() {
    const sel = document.getElementById('algorithm-select');
    return !!sel && sel.value === 'sim';
}

/* The controls that only mean something when a pattern is being FITTED.
   Hidden, never disabled and never removed: everything in them keeps its
   value, so switching to Simulation and back leaves the background anchors,
   the charge-flipping settings and the Wyckoff composition exactly as they
   were. That is what "data should not be lost" has to mean here -- there is no
   separate store to save them into, the DOM is the store.

   Background is on this list even though the simulation still ADDS the
   background to the curve when a file is loaded. The tab is a spline editor
   driven by clicking on measured data: with nothing measured there is nothing
   to place anchors against, and anchors set during a fit stay in force and
   stay correct without the tab being open. */
const SIM_FIT_ONLY_LEFT_TABS  = ['background', 'charge-flipping', 'wyckoff'];
const SIM_FIT_ONLY_RIGHT_TABS = ['lebail', 'pawley', 'cf', 'wyckoff'];

/**
 * Show or hide the fit-only controls and the Atoms tab.
 *
 * Called on every change of the Method selector, and once at startup.
 */
function simApplyMode() {
    const on = simActive();

    // Drives the stylesheet rule that hides the per-parameter "Fit" flags.
    // A class on <body> rather than a JS sweep, because those flags are
    // recreated whenever the crystal system or the profile type changes; see
    // the comment on the rule itself in powder5.html.
    document.body.classList.toggle('sim-active', on);

    // --- left panel tab strip -------------------------------------------
    SIM_FIT_ONLY_LEFT_TABS.forEach(name => {
        const btn = document.querySelector(`.tab-buttons .tab-btn[data-tab="${name}"]`);
        if (btn) btn.style.display = on ? 'none' : '';
    });
    const atomsBtn = document.getElementById('tab-btn-atoms');
    if (atomsBtn) atomsBtn.style.display = on ? '' : 'none';

    // --- right panel tab strip ------------------------------------------
    SIM_FIT_ONLY_RIGHT_TABS.forEach(name => {
        const btn = document.querySelector(`.rp-tab[data-rptab="${name}"]`);
        if (btn) btn.style.display = on ? 'none' : '';
    });

    // --- the run bar, the iteration count and the fit statistics ---------
    // Rp / Rwp / chi^2 are goodness-of-fit measures against a measurement.
    // With no measurement they are not "blank", they are meaningless, and a
    // dash in a results card reads as "not run yet" rather than "does not
    // apply".
    const bar = document.getElementById('action-button-bar');
    if (bar) bar.style.display = on ? 'none' : '';
    const iterLabel = document.querySelector('label[for="iterations-slider"]');
    const iterTrack = document.getElementById('iterations-slider');
    const iterRow   = iterTrack ? iterTrack.closest('.slider-value-track') : null;
    if (iterLabel) iterLabel.style.display = on ? 'none' : '';
    if (iterRow)   iterRow.style.display   = on ? 'none' : '';
    const grid = document.querySelector('.bottom-actions .results-grid');
    if (grid) grid.style.display = on ? 'none' : '';

    // --- if the tab in front is one that just disappeared, step off it ---
    //
    // On ENTERING simulation, go to Atoms unconditionally: the tab has just
    // appeared, it is the only place the new mode can be driven from, and a
    // mode switch that changes nothing visible in the left panel reads as
    // having done nothing at all. On every later call -- and simApplyMode()
    // runs on each change of the selector -- leave the tab alone, or the panel
    // would jump under the user's hands.
    const entering = on && !simWasActive;
    simWasActive = on;

    if (on) {
        const activeLeft = document.querySelector('.tab-buttons .tab-btn.active');
        if (entering || (activeLeft && SIM_FIT_ONLY_LEFT_TABS.includes(activeLeft.dataset.tab))) {
            simSwitchLeftTab('atoms');
        }
        const activeRight = document.querySelector('.rp-tab.active');
        if (activeRight && SIM_FIT_ONLY_RIGHT_TABS.includes(activeRight.dataset.rptab)) {
            if (typeof rpSwitchTab === 'function') rpSwitchTab('plot');
        }
    } else {
        const activeLeft = document.querySelector('.tab-buttons .tab-btn.active');
        if (activeLeft && activeLeft.dataset.tab === 'atoms') simSwitchLeftTab('sample');
    }

    simEnableTthRange(on);

    if (entering) {
        simLog(`Method set to Simulation. ${simAtoms.length} atom site(s) in the model` +
               (simAtoms.length ? ': ' + simAtoms.map(simAtomLine).join(' | ') : '.'));
    } else if (!on && simWasLogged) {
        simLog('Method left Simulation.');
    }
    simWasLogged = on;

    if (on) {
        simEnsureTables();
        simRender();
    } else {
        // Hand the chart back to whichever path owns it outside simulation.
        // The cached list is dropped with it: the fit paths write to
        // masterHklList themselves, and a stale simulation list surviving into
        // a refinement would be a list the fit never generated.
        simInvalidateGeometry();
        simClearPattern();
        if (typeof hasExperimentalData === 'function' && hasExperimentalData()) {
            if (typeof updatePreviewPattern === 'function') updatePreviewPattern();
        } else if (typeof root.drawTheoreticalPreview === 'function') {
            root.drawTheoreticalPreview();
        }
    }
}

/** Click the left-panel tab button, so the app's own handler does the work. */
function simSwitchLeftTab(name) {
    const btn = document.querySelector(`.tab-buttons .tab-btn[data-tab="${name}"]`);
    if (btn) btn.click();
}

/**
 * The 2-theta sliders are disabled until a file is loaded, because outside
 * simulation they select a slice OF the data and there is nothing to slice.
 * A simulation has the opposite problem: the range is the only thing that says
 * how much pattern to compute, so it has to be settable with no file open.
 *
 * With data loaded the sliders are already live and already carry the data's
 * own limits -- this leaves them completely alone, which is what makes the
 * simulated pattern come out on the same 2-theta interval as the fit would.
 */
function simEnableTthRange(on) {
    if (typeof hasExperimentalData === 'function' && hasExperimentalData()) return;

    const els = [
        document.getElementById('tth-min-range'),
        document.getElementById('tth-min-slider'),
        document.getElementById('tth-max-range'),
        document.getElementById('tth-max-slider')
    ];
    if (els.some(e => !e)) return;

    if (!on) {
        // Restore the no-data state exactly: disabled, and no stale numbers
        // left behind to be read by anything that only checks .value.
        els.forEach(e => { e.disabled = true; });
        return;
    }

    const LO = 1, HI = 160;
    els.forEach(e => { e.min = LO; e.max = HI; e.step = 0.01; e.disabled = false; });

    // Seed once, from the theoretical preview's own range, so the first
    // simulation covers the same interval the stick preview did. Re-entering
    // simulation mode later keeps whatever the user set.
    const lo = document.getElementById('tth-min-range');
    const hi = document.getElementById('tth-max-range');
    if (!lo.value || !hi.value || !(parseFloat(hi.value) > parseFloat(lo.value))) {
        const d0 = (typeof THEO_TTH_MIN === 'number') ? THEO_TTH_MIN : 5;
        const d1 = (typeof THEO_TTH_MAX === 'number') ? THEO_TTH_MAX : 120;
        lo.value = d0; hi.value = d1;
        document.getElementById('tth-min-slider').value = d0;
        document.getElementById('tth-max-slider').value = d1;
    }
}

/** The current 2-theta interval, from the sliders. */
function simRange() {
    const lo = parseFloat((document.getElementById('tth-min-slider') || {}).value);
    const hi = parseFloat((document.getElementById('tth-max-slider') || {}).value);
    const d0 = (typeof THEO_TTH_MIN === 'number') ? THEO_TTH_MIN : 5;
    const d1 = (typeof THEO_TTH_MAX === 'number') ? THEO_TTH_MAX : 120;
    const a = Number.isFinite(lo) ? lo : d0;
    const b = Number.isFinite(hi) ? hi : d1;
    return (b > a) ? { lo: a, hi: b } : { lo: d0, hi: d1 };
}

// ---------------------------------------------------------------------------
//  SCATTERING TABLES
// ---------------------------------------------------------------------------

/**
 * Load scatters/index.json and the preferred X-ray table, once.
 *
 * Deliberately lazy: someone who never opens Simulation should not pay for a
 * 50 kB fetch, and someone who does should not have the app refuse to start
 * because the tables are missing.
 *
 * A failure here is reported and is FATAL to the simulation rather than being
 * papered over with f = Z. Point atoms are not a degraded pattern, they are a
 * different one -- no angular fall-off means every high-angle reflection comes
 * out several times too strong -- and a user comparing that against a
 * measurement would read the discrepancy as a wrong structure.
 */
function simEnsureTables() {
    if (simTables || simTablesPending) return simTablesPending;
    if (typeof loadScatteringTables !== 'function') {
        simTablesError = 'js/scatterers.js is not loaded.';
        simSetStatus(simTablesError, 'warn');
        return null;
    }
    simTablesPending = loadScatteringTables('scatters')
        .then(t => {
            simTables = t;
            simTablesError = null;
            // The editor can be open and waiting on this: simOpenAtomModal()
            // kicks the load off and draws a "loading" placeholder in the list,
            // so the list has to be redrawn when the tables land or the dialog
            // stays empty until the user types.
            simFilterElements();
            simRender();
        })
        .catch(e => {
            simTablesError = `Could not load the scattering factors: ${e.message || e}`;
            simSetStatus(simTablesError, 'warn');
        })
        .finally(() => { simTablesPending = null; });
    return simTablesPending;
}

/**
 * Every label in the X-ray table, ordered for the picker.
 *
 * BY ELEMENT FIRST, THEN CHARGE -- so "T" gives Ti, Ti2+, Ti3+, Ti4+, Tc, Te,
 * Tb, Tm, Ta, Tl, Tl1+, Th. An ion belongs next to the atom it came from,
 * because someone typing "T" is looking for an element and then deciding its
 * charge, not scanning two separate alphabets. Sorting all the neutral atoms
 * ahead of all the ions is the other obvious grouping and it separates Ti from
 * Ti4+ by eight rows, which is exactly the lookup this list exists to make
 * easy.
 *
 * Within an element the neutral atom leads and the ions follow by label, which
 * for the tabulated charges is the same as by magnitude.
 */
function simTableLabels() {
    if (!simTables || !simTables.xray || !simTables.xray.table) return [];
    const entries = Object.entries(simTables.xray.table);
    const isNeutral = (label) => /^[A-Za-z]+$/.test(label);
    entries.sort((p, q) => {
        const [la, ea] = p, [lb, eb] = q;
        const za = ea.Z || 0, zb = eb.Z || 0;
        if (za !== zb) return za - zb;             // element first
        const na = isNeutral(la), nb = isNeutral(lb);
        if (na !== nb) return na ? -1 : 1;         // its neutral atom leads
        return la.localeCompare(lb);               // then its ions
    });
    return entries.map(([label]) => label);
}

/**
 * f(label, s) for the loaded table, or null if the label is not in it.
 *
 * s = 1/d. evaluateGaussian() wants stol = s/2 -- see note 1 in the header.
 */
function simFormFactor(label, s) {
    if (!simTables || !simTables.xray) return null;
    const e = simTables.xray.table[label];
    if (!e) return null;
    if (typeof evaluateGaussian !== 'function') return null;
    return evaluateGaussian(e, s / 2);
}

// ---------------------------------------------------------------------------
//  SYMMETRY
// ---------------------------------------------------------------------------

/**
 * The full operator list of the current setting: every sym_op crossed with
 * every centring translation, deduplicated.
 *
 * THE CROSS PRODUCT IS NOT REDUNDANT, and it is not harmful either. The v8
 * JSON writes order_z operators into sym_ops -- the full set, centring
 * included -- but it ALSO writes centring_translations, and nothing in the
 * schema promises which convention a regenerated file will use. Multiplying
 * and then deduplicating is correct under both: if sym_ops already carries the
 * centring, every product is a duplicate of an operator already in the list
 * and the dedup removes it; if it does not, the products are exactly the
 * operators that were missing.
 *
 * Getting this wrong is not a subtle error. An F-centred cell has four times
 * the atoms of its primitive subgroup, and dropping the centring would put a
 * quarter of the scattering matter in the cell while leaving every centring
 * absence in place -- a pattern with the right peak POSITIONS and uniformly
 * wrong intensities, which is the hardest kind of error to see.
 */
function simSymmetry(sg) {
    if (!sg) return { error: 'No space group is selected.' };
    const key = `${sg.number}|${sg.hall || ''}|${sg.symbol || ''}`;
    if (simSymCache.has(key)) return simSymCache.get(key);

    let out;
    const get = root.getSymopsForSetting;
    if (typeof get !== 'function') {
        out = { error: 'getSymopsForSetting() is not available.' };
    } else {
        const r = get(sg);
        if (r.error) {
            out = { error: r.error };
        } else {
            const cen = (r.setting && Array.isArray(r.setting.centring_translations)
                         && r.setting.centring_translations.length)
                      ? r.setting.centring_translations
                      : [{ t: [0, 0, 0] }];

            const wrap = v => ((v % 1) + 1) % 1;
            const seen = new Set();
            const ops = [];
            for (const op of r.symops) {
                for (const c of cen) {
                    // Exact rationals where the generator gave them: 1/3
                    // arriving as 0.333333 does not close a rhombohedral cell.
                    const ct = (c.t_num && c.t_den)
                        ? c.t_num.map(n => Number(n) / Number(c.t_den))
                        : (c.t || [0, 0, 0]).map(Number);
                    const t = [wrap(op.t[0] + ct[0]),
                               wrap(op.t[1] + ct[1]),
                               wrap(op.t[2] + ct[2])];
                    const k = op.r.join(',') + '|' +
                              t.map(v => Math.round(v * 24)).join(',');
                    if (seen.has(k)) continue;
                    seen.add(k);
                    ops.push({ r: op.r, t, xyz: op.xyz || '' });
                }
            }
            out = {
                ops,
                wyckoff: (r.setting && Array.isArray(r.setting.wyckoff))
                         ? r.setting.wyckoff : [],
                settingSymbol: (r.setting && r.setting.symbol) || sg.symbol || '',
                orderZ: (r.setting && r.setting.order_z) || ops.length
            };
        }
    }
    simSymCache.set(key, out);
    return out;
}

/** The real-space metric tensor for the cell currently on screen. */
function simMetric(params) {
    return metricTensor({
        a: params.a, b: params.b > 0 ? params.b : params.a,
        c: params.c > 0 ? params.c : params.a,
        alpha: params.alpha, beta: params.beta, gamma: params.gamma
    });
}

/**
 * The orbit of one fractional point: every DISTINCT image under the group.
 *
 * Distinctness is decided in angstrom, with the minimum-image convention, and
 * that matters twice over. A point at x = 1e-9 and its image at
 * x = 0.999999999 are the same atom one lattice translation apart, so a plain
 * |a - b| comparison calls them different and doubles that site's scattering
 * -- and sites on or near a coordinate of zero are exactly the special
 * positions a structure is most likely to use. crExpand() in coord_refine.js
 * carries the same warning for the same reason.
 */
function simOrbit(x, y, z, ops, G) {
    const wrap = v => ((v % 1) + 1) % 1;
    const out = [];
    for (const op of ops) {
        const r = op.r, t = op.t;
        const q = [
            wrap(r[0] * x + r[1] * y + r[2] * z + t[0]),
            wrap(r[3] * x + r[4] * y + r[5] * z + t[1]),
            wrap(r[6] * x + r[7] * y + r[8] * z + t[2])
        ];
        let dup = false;
        for (const p of out) {
            if (fracDistance(G, q[0] - p[0], q[1] - p[1], q[2] - p[2]) < SIM_MERGE_TOL_A) {
                dup = true; break;
            }
        }
        if (!dup) out.push(q);
    }
    return out;
}

/** Projection onto a Wyckoff position: P.x + T, from the JSON's rationals. */
function simProject(w, x, y, z) {
    let P, T;
    if (Array.isArray(w.P_num) && w.P_num.length === 9 && w.P_den) {
        P = w.P_num.map(v => v / w.P_den);
        T = (w.T_num || [0, 0, 0]).map(v => v / (w.T_den || 1));
    } else {
        P = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        T = [0, 0, 0];
    }
    return [
        P[0] * x + P[1] * y + P[2] * z + T[0],
        P[3] * x + P[4] * y + P[5] * z + T[1],
        P[6] * x + P[7] * y + P[8] * z + T[2]
    ];
}

/**
 * Which Wyckoff position a coordinate sits on. Derived, never entered.
 *
 * THE MULTIPLICITY DECIDES, THE PROJECTION ONLY BREAKS TIES.
 *
 * The orbit size computed above is not an estimate of the multiplicity, it IS
 * the multiplicity -- it is the number of atoms this site will actually
 * contribute to the structure factor, counted from the same operators and with
 * the same tolerance. So the letter is looked up among the positions of that
 * multiplicity and nowhere else, which means the letter on screen can never
 * disagree with the pattern being drawn.
 *
 * Several positions can share a multiplicity (Pnma has 4a, 4b and 4c), and
 * those are separated by the projection: the special_op of the right one
 * leaves some image of the point where it already is. The test runs over every
 * image because special_op describes ONE representative of the position -- a
 * point on 4b at (0,0,1/2) may only match after the operator that carries it
 * there.
 */
function simWyckoffFor(x, y, z, sym, G) {
    if (!sym || sym.error || !sym.ops) return { label: '?', mult: 0 };
    const orbit = simOrbit(x, y, z, sym.ops, G);
    const mult = orbit.length;

    const cands = (sym.wyckoff || []).filter(w => Number(w.multiplicity) === mult);
    if (!cands.length) return { label: `${mult}?`, mult, siteSym: '' };
    if (cands.length === 1) {
        return { label: `${mult}${cands[0].letter}`, mult,
                 siteSym: cands[0].site_symmetry || '' };
    }

    let best = null, bestD = Infinity;
    for (const w of cands) {
        let d = Infinity;
        for (const q of orbit) {
            const p = simProject(w, q[0], q[1], q[2]);
            const dd = fracDistance(G, p[0] - q[0], p[1] - q[1], p[2] - q[2]);
            if (dd < d) d = dd;
        }
        if (d < bestD) { bestD = d; best = w; }
    }
    return best
        ? { label: `${mult}${best.letter}`, mult, siteSym: best.site_symmetry || '' }
        : { label: `${mult}?`, mult, siteSym: '' };
}

// ---------------------------------------------------------------------------
//  STRUCTURE FACTORS
// ---------------------------------------------------------------------------

/**
 * |F|^2 for every reflection in `list`.
 *
 * The occupancy and the Debye-Waller factor both belong to the SITE and depend
 * on the reflection only through s, so they are folded into the site's
 * scattering factor once per reflection rather than once per atom. For a
 * high-symmetry structure that is the difference between 192 exponentials and
 * one.
 *
 * @param {object[]} list    positioned reflections (h_orig/k_orig/l_orig, d)
 * @param {object}   params  getAllParams()
 * @param {object}   sym     simSymmetry()
 * @returns {{rows: object[], warnings: string[], nAtoms: number}}
 */
function simStructureFactors(list, params, sym) {
    const warnings = [];
    const G = simMetric(params);
    const sites = simAtoms.filter(simAtomIsUsable);

    // Expand each site to its orbit once; the same atoms serve every
    // reflection.
    const atoms = [];
    const perSite = [];
    sites.forEach((s, i) => {
        const orbit = simOrbit(s.x, s.y, s.z, sym.ops, G);
        perSite.push(orbit.length);
        for (const q of orbit) atoms.push({ site: i, x: q[0], y: q[1], z: q[2] });
    });

    const missing = new Set();
    const TWO_PI = 2 * Math.PI;
    const f = new Float64Array(sites.length);

    const rows = list.map(r => {
        const out = { h: r.h_orig, k: r.k_orig, l: r.l_orig, Fsq: 0 };
        if (!Number.isFinite(r.d) || !(r.d > 0)) return out;
        const s = 1 / r.d;

        for (let i = 0; i < sites.length; i++) {
            const st = sites[i];
            const fi = simFormFactor(st.label, s);
            if (fi === null) { missing.add(st.label); f[i] = 0; continue; }
            // exp(-B stol^2) with stol = s/2, i.e. exp(-B s^2 / 4).
            f[i] = st.occ * fi * Math.exp(-st.biso * s * s / 4);
        }

        let A = 0, B = 0;
        for (let j = 0; j < atoms.length; j++) {
            const a = atoms[j];
            const fj = f[a.site];
            if (fj === 0) continue;
            const ph = TWO_PI * (r.h_orig * a.x + r.k_orig * a.y + r.l_orig * a.z);
            A += fj * Math.cos(ph);
            B += fj * Math.sin(ph);
        }
        out.Fsq = A * A + B * B;
        return out;
    });

    if (missing.size) {
        warnings.push(`No scattering factor for ${[...missing].join(', ')} -- ` +
                      `${missing.size > 1 ? 'those atoms contribute' : 'that atom contributes'} nothing.`);
    }
    return { rows, warnings, nAtoms: atoms.length, perSite, sites };
}

// ---------------------------------------------------------------------------
//  THE PATTERN
// ---------------------------------------------------------------------------

/**
 * Hand the Ka tick datasets to the chart, capped and cached.
 *
 * Two costs are avoided here. The first is rebuilding the datasets at all when
 * the geometry has not moved -- an atom edit changes intensities, not
 * positions, so the ticks are already right. The second is the cap: past
 * SIM_MAX_MARKERS the strongest reflections are kept and the rest dropped,
 * because handing Chart.js one layout element per reflection is what turns a
 * large cell over a wide range from slow into unusable.
 *
 * @returns {number} how many reflections were dropped from the tick strip.
 */
function simBuildMarkers(params, list) {
    if (typeof refreshHklMarkers !== 'function') return 0;

    if (list.length <= SIM_MAX_MARKERS) {
        // Positions only: safe to skip when the geometry is unchanged.
        if (simMarkersKey !== simGeomKey) {
            refreshHklMarkers(params, list);
            simMarkersKey = simGeomKey;
            simCounters.markerBuilds++;
        }
        return 0;
    }

    // Capped: which reflections survive depends on intensity, so this has to
    // be redone when the atoms change. Bounded at SIM_MAX_MARKERS, so the cost
    // is flat however many reflections the cell actually has.
    const keep = list.slice()
        .sort((a, b) => (b.intensity || 0) - (a.intensity || 0))
        .slice(0, SIM_MAX_MARKERS)
        .sort((a, b) => a.tth - b.tth);
    refreshHklMarkers(params, keep);
    simMarkersKey = null;
    simCounters.markerBuilds++;
    return list.length - keep.length;
}

/** Give the Ka1 / Ka2 bars a height.
 *
 * refreshHklMarkers() writes only an x per reflection; the height is set by
 * whoever knows the y scale. In data mode that is rescalePlot(), which sizes
 * them in PIXELS so they stay the same size on screen at any zoom -- and that
 * is the behaviour "as usual" means, so this leaves it alone. With no data
 * rescalePlot() returns before it reaches them, so the heights are set here
 * as a fraction of the plotted maximum instead.
 *
 * @param {number} refMax  the pattern maximum the ticks are scaled against.
 */
function simSizeMarkers(refMax) {
    const bar = (label, frac) => {
        const ds = mainChart.data.datasets.find(d => d.label === label);
        if (!ds || !ds.data) return;
        const lo = -frac * refMax;
        // Written in place. Allocating a fresh two-element array per bar per
        // render was thousands of short-lived objects on every keystroke, and
        // the arrays are private to this dataset.
        for (let i = 0; i < ds.data.length; i++) {
            const pt = ds.data[i];
            if (Array.isArray(pt.y)) { pt.y[0] = lo; pt.y[1] = 0; }
            else pt.y = [lo, 0];
        }
    };
    bar('Ka1', 0.055);
    bar('Ka2', 0.035);
}

/**
 * The reflection positions and nothing else: no atoms, so no intensities.
 *
 * Deliberately the same picture drawTheoreticalPreview() draws outside
 * simulation mode -- same sticks, same axis treatment -- because it is
 * answering the same question. The alternative was a blank chart, which says
 * "the simulation is broken" rather than "the model is empty".
 */
function simMarkersOnly(haveData) {
    const y = mainChart.options.scales.y;
    if (haveData) {
        // The data is still on screen and still sets the scale; let the
        // existing machinery place the ticks against it.
        if (typeof updatePlotRange === 'function') updatePlotRange(true);
        if (typeof rescalePlot === 'function') rescalePlot(true);
    } else {
        const find = (l) => mainChart.data.datasets.find(d => d.label === l);
        ['Experimental', 'Calculated', 'Background', 'Difference', 'Difference Zero']
            .forEach(l => { const d = find(l); if (d) d.data = []; });
        const { lo, hi } = simRange();
        mainChart.options.scales.x.min = lo;
        mainChart.options.scales.x.max = hi;
        y.min = -0.15;
        y.max = 1.15;
        if (y.ticks) y.ticks.display = false;
        if (y.title) y.title.text = 'Reflections (positions only)';
        mainChart.options.globalYMax = 1.15;
        simSizeMarkers(1.0);
    }
    mainChart.update('none');
}

/** The chart dataset the simulation draws into. */
function simDataset() {
    if (!mainChart) return null;
    return mainChart.data.datasets.find(d => d.label === 'Simulation (Manual)') || null;
}

/**
 * The simulation reuses the existing 'Simulation (Manual)' dataset rather than
 * adding a tenth one. That dataset already means "a pattern computed from the
 * model rather than fitted to the data", which is exactly what this is, and a
 * second curve with the same meaning could only ever disagree with the first.
 *
 * It is filtered OUT of the legend by initializeChart(), because outside
 * simulation it is a transient preview. Here it is the whole answer, so the
 * legend filter is relaxed and the label says so.
 */
function simShowLegendEntry(on) {
    if (!mainChart) return;
    const legend = mainChart.options.plugins && mainChart.options.plugins.legend;
    if (!legend || !legend.labels) return;

    if (on) {
        legend.labels.filter = item => item.text !== 'Difference Zero'
                                    && item.text !== 'Spline Points';
        // The dataset KEY stays 'Simulation (Manual)' -- eight call sites
        // across three files find this curve by that string, and renaming it
        // to make one legend read better would be a rename with eight chances
        // to be missed. The legend text is a separate thing, so relabel that.
        legend.labels.generateLabels = (chart) => {
            const base = Chart.defaults.plugins.legend.labels.generateLabels(chart);
            base.forEach(it => {
                if (it.text === 'Simulation (Manual)') it.text = 'Simulated';
            });
            return base;
        };
    } else {
        legend.labels.filter = item => item.text !== 'Difference Zero'
                                    && item.text !== 'Simulation (Manual)'
                                    && item.text !== 'Spline Points';
        delete legend.labels.generateLabels;
    }
}

function simClearPattern() {
    const ds = simDataset();
    if (ds) ds.data = [];
    simShowLegendEntry(false);
    if (mainChart) mainChart.update('none');
}

/** Read a clamped number out of one of the Pattern boxes. */
function simReadBox(id, lo, hi, dflt) {
    const el = document.getElementById(id);
    const v = el ? parseFloat(el.value) : NaN;
    if (!Number.isFinite(v)) return dflt;
    return simClamp(v, lo, hi);
}

const simConstantBackground = () =>
    simReadBox('sim-background', SIM_BKG_MIN, SIM_BKG_MAX, SIM_BKG_DEFAULT);
const simNoiseLevel = () =>
    simReadBox('sim-noise', SIM_NOISE_MIN, SIM_NOISE_MAX, SIM_NOISE_DEFAULT);

/**
 * A small deterministic PRNG (mulberry32).
 *
 * DETERMINISTIC ON PURPOSE, and for two separate reasons.
 *
 * On screen: the pattern is recomputed on every keystroke, and noise drawn
 * fresh each time would make the whole curve shimmer while an occupancy is
 * being typed -- the change the user is trying to see would be buried in the
 * change they are not.
 *
 * In the file: a synthetic dataset whose noise cannot be reproduced is a poor
 * test case. With a fixed seed the same model and the same settings give the
 * same file, so a disagreement between two programs reading it is a
 * disagreement about the programs and not about which random numbers each run
 * happened to draw.
 */
function simRng(seed) {
    let a = seed >>> 0;
    return function () {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

/**
 * Add counting noise to one point.
 *
 * WHY sqrt(I) AND NOT A FIXED SIGMA. A diffractometer counts photons, and
 * photon counting is a Poisson process: the variance of a measured intensity
 * IS the intensity, so the scatter goes as sqrt(I). That single fact is what
 * makes a simulated pattern behave like a real one -- strong peaks noisy in
 * absolute terms and clean in relative terms, weak peaks the other way round,
 * and a background whose scatter sets the detection limit. Gaussian noise of
 * one fixed width reproduces none of that: it would swamp the weak reflections
 * and be invisible on the strong ones, and any program weighting by 1/sigma^2
 * would be given a wrong weight at every point.
 *
 * The Gaussian approximation to Poisson is used rather than a true Poisson
 * draw, because that is what admits the `level` multiplier: a true Poisson
 * variate has no free width, its scatter is fixed by the count, so the only
 * way to ask for "twice the noise at this intensity" is to widen a Gaussian.
 * The approximation is good above roughly 10 counts and the default background
 * sits there deliberately.
 *
 * Negative results are clipped to zero. Counts cannot be negative, some
 * readers reject them, and at the default settings zero is more than three
 * sigma below the background, so the clip essentially never fires -- but when
 * it does, it biases that point upward, which is the honest trade against
 * writing a physically impossible number.
 */
function simNoisy(clean, level, rnd) {
    const sigma = level * Math.sqrt(Math.max(clean, 1));
    // Box-Muller.
    const u1 = Math.max(rnd(), 1e-12), u2 = rnd();
    const g = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    const v = clean + sigma * g;
    return v > 0 ? v : 0;
}

/** The step from the Atoms tab, clamped to the physical range. */
function simRequestedStep() {
    const el = document.getElementById('sim-step');
    const v = el ? parseFloat(el.value) : NaN;
    if (!Number.isFinite(v)) return SIM_STEP_DEFAULT;
    return simClamp(v, SIM_STEP_MIN, SIM_STEP_MAX);
}

/**
 * A synthetic ascending axis over [lo, hi] at the requested step.
 *
 * @returns {{axis: Float64Array, step: number, capped: boolean}}
 */
function simSyntheticAxis(lo, hi) {
    const want = simRequestedStep();
    let n = Math.floor((hi - lo) / want) + 1;
    let capped = false;
    if (n > SIM_MAX_POINTS) { n = SIM_MAX_POINTS; capped = true; }
    if (n < 2) n = 2;
    const step = (hi - lo) / (n - 1);
    // Only report a cap that actually coarsened the step. Asking for 0.0001
    // across exactly 2 degrees needs 20001 points and gets 20000, which is a
    // step of 0.00010001 -- clipping one point off the end is not something
    // worth a warning about, and a warning that fires when nothing went wrong
    // teaches the user to ignore the next one.
    if (capped && step <= want * 1.001) capped = false;
    const ax = new Float64Array(n);
    for (let i = 0; i < n; i++) ax[i] = lo + i * step;
    return { axis: ax, step, capped };
}

/**
 * Compute and draw the simulated pattern. The single entry point: every
 * control that can change the pattern ends up here.
 *
 * A RE-ENTRANT CALL IS IGNORED, NOT QUEUED.
 *
 * JavaScript is single-threaded, so a call arriving while simRendering is true
 * can only have come from inside simRenderInner's own call stack. It is never
 * a new request from the user -- it is this render's own side effects coming
 * back round. There is therefore nothing to catch up on, and the earlier
 * version's second pass simply recomputed the identical pattern and doubled
 * the cost of every render.
 *
 * There is one such path and it is legitimate. In data mode simRenderInner
 * calls updateWorkingData() to slice the pattern to the 2-theta range, and
 * that calls updateBackgroundForPreview(), which ends by calling
 * updatePreviewPattern(), which is routed back here. The nesting is correct
 * behaviour -- the background really does need recomputing, and it happens
 * before this render reads it -- so it is counted rather than warned about.
 * The stack of the first re-entry is kept for SIM.stats() in case a path
 * appears that is NOT this one.
 */
function simRender() {
    if (!simActive()) return;

    if (simRendering) {
        simCounters.reentries++;
        if (!simCounters.reentryStack) {
            simCounters.reentryStack = (new Error('SIM re-entry').stack || '').split('\n')
                .slice(1, 8).map(l => l.trim()).join(' | ');
            // Once per session, at debug level: the known path is benign, and a
            // warning that fires on every render for something working as
            // designed is how real warnings get ignored.
            if (typeof console !== 'undefined' && console.debug) {
                console.debug('SIM: render re-entered from within itself (expected via ' +
                              'updateWorkingData -> updateBackgroundForPreview -> ' +
                              'updatePreviewPattern). Ignored. Path:\n' +
                              simCounters.reentryStack);
            }
        }
        return;
    }

    simRendering = true;
    const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
    try {
        simRenderInner();
    } catch (e) {
        console.error('Simulation failed:', e);
        simSetStatus(`Simulation failed: ${e.message || e}`, 'warn');
    } finally {
        simRendering = false;
        const dt = (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0;
        simCounters.renders++;
        simCounters.lastMs = Math.round(dt);
        if (dt > simCounters.maxMs) simCounters.maxMs = Math.round(dt);
    }
}

/**
 * The reflections and their 2-theta positions, cached.
 *
 * The key is everything the geometry depends on and nothing else. The atoms
 * are deliberately absent from it: they set the intensities, which are applied
 * to this list afterwards, and rebuilding the list when an atom moves was pure
 * waste -- the expensive parts are the generator, the sort, and handing the
 * marker datasets to Chart.js.
 */
function simGeometry(sg, params, maxTth) {
    const key = [sg.number, sg.hall || '', sg.symbol || '', sg.system || '',
                 params.a, params.b, params.c, params.alpha, params.beta, params.gamma,
                 params.lambda, params.lambda2, params.ratio, params.zeroShift,
                 maxTth].join('|');
    if (key === simGeomKey && simGeomList) return simGeomList;

    const raw = generateAndCacheHklIndices(sg, maxTth, params);
    let list = raw.map(r => ({ ...r, hkl_list: r.hkl_list ? r.hkl_list.slice() : r.hkl_list }));
    updateHklPositions(list, params, sg.system);
    list = list.filter(p => Number.isFinite(p.tth) && p.tth <= maxTth)
               .sort((a, b) => a.tth - b.tth);

    simGeomKey = key;
    simGeomList = list;
    simMarkersKey = null;            // the tick strip belongs to the old geometry
    simCounters.geomBuilds++;
    return list;
}

/**
 * Restore the axes after a zoom, for the no-data case.
 *
 * The chart's own reset-view handler restores x from the sliders and then
 * calls updatePlotRange() and rescalePlot(), both of which return immediately
 * without a workingData -- so under Simulation with no file open the y axis
 * was simply left where the zoom had put it. This re-applies the extent the
 * last render computed, without recomputing the pattern: nothing about the
 * model changed, only the view.
 */
function simResetView() {
    if (!mainChart || !simPts || !simPts.length) return;
    let max = 0;
    for (const p of simPts) if (p.y > max) max = p.y;
    if (!(max > 0)) max = 1;
    const { lo, hi } = simRange();
    mainChart.options.scales.x.min = lo;
    mainChart.options.scales.x.max = hi;
    mainChart.options.scales.y.min = -0.10 * max;
    mainChart.options.scales.y.max = 1.10 * max;
    mainChart.options.globalYMax = 1.10 * max;
    simSizeMarkers(max);
    mainChart.update('none');
}

/** Drop the cached geometry: the cell, the group or the mode changed. */
function simInvalidateGeometry() {
    simGeomKey = null;
    simGeomList = null;
    simMarkersKey = null;
}

function simRenderInner() {
    const sg = currentSG;
    if (!sg) { simSetStatus('Select a space group first.', 'warn'); return; }

    if (!mainChart) {
        if (typeof initializeChart !== 'function') return;
        initializeChart();
    }
    // With no file loaded the placeholder covers the chart; the simulation is
    // the reason to look at it.
    if (controls && controls.placeholder && controls.resultsContainer) {
        controls.placeholder.style.display = 'none';
        controls.resultsContainer.style.display = 'flex';
    }

    const params = getAllParams();
    if (!(params.lambda > 0) || !(params.a > 0)) {
        simSetStatus('Set a wavelength and a lattice parameter first.', 'warn');
        return;
    }

    if (!simTables) {
        simEnsureTables();
        simSetStatus(simTablesError || 'Loading scattering factors...',
                     simTablesError ? 'warn' : '');
        return;
    }

    const sym = simSymmetry(sg);
    if (sym.error) { simSetStatus(sym.error, 'warn'); return; }

    // --- the axis --------------------------------------------------------
    const haveData = (typeof hasExperimentalData === 'function') && hasExperimentalData();
    let axis, stepNote = '';
    if (haveData) {
        if (typeof updateWorkingData === 'function') updateWorkingData();
        axis = workingData.tth;
        if (!axis || !axis.length) { simSetStatus('No points in the 2θ range.', 'warn'); return; }
    } else {
        const { lo, hi } = simRange();
        const built = simSyntheticAxis(lo, hi);
        axis = built.axis;
        stepNote = built.capped
            ? `A step of ${simRequestedStep()}° over ${(hi - lo).toFixed(2)}° would need ` +
              `more than ${SIM_MAX_POINTS} points, so ${built.step.toFixed(5)}° was used. ` +
              `Narrow the 2θ range to get the finer step.`
            : '';
    }
    const maxTth = axis[axis.length - 1] + 2.0;

    // --- reflections (the app's own generator, absences included) --------
    // Cached on the geometry, not on the atoms: see simGeometry().
    const list = simGeometry(sg, params, maxTth);

    // --- |F|^2, then heights through the app's own inverse ---------------
    const sf = simStructureFactors(list, params, sym);
    const heights = reflectionHeightsFromFsq(sf.rows, params, list, { scale: 1 });
    let anyIntensity = false;
    list.forEach(p => {
        const h = heights.get(p.hkl_list[0]);
        p.intensity = Number.isFinite(h) ? h : 0;
        if (p.intensity !== 0) anyIntensity = true;
    });

    // Downstream consumers -- the tooltip, the marker refresh, the reflection
    // report -- read these, so they describe the simulation while it is on
    // screen rather than a stale fit.
    masterHklList = list;
    lastGeneratedHklList = list;
    const capped = simBuildMarkers(params, list);

    // --- profile ---------------------------------------------------------
    const ds = simDataset();
    if (!ds) { simSetStatus('The chart has no simulation dataset.', 'warn'); return; }

    if (!anyIntensity) {
        // NO ATOMS YET, OR NOTHING THAT SCATTERS.
        //
        // This is the state the Atoms tab opens in, and it used to return here
        // without ever giving the Ka markers a height -- refreshHklMarkers()
        // writes an x and leaves y to whoever sizes the axes, and the only two
        // things that do (rescalePlot, and the branch in simApplyAxes below)
        // were both downstream of this return. So the reflection positions
        // were computed, were correct, were in the dataset, and drew as
        // nothing.
        //
        // Falling back to a stick plot is also the right ANSWER here: with no
        // atoms the positions are the whole of what the model predicts, which
        // is exactly what drawTheoreticalPreview() shows outside this mode.
        ds.data = [];
        simShowLegendEntry(false);
        simMarkersOnly(haveData);
        simSetStatus(simAtoms.length
            ? 'No reflection has any intensity. Check the elements and occupancies.'
            : 'Add at least one atom to simulate a pattern. The sticks are the reflection positions for this cell.',
            simAtoms.length ? 'warn' : '');
        simSetSummary(sf, list, sym);
        return;
    }

    const unscaled = calculatePatternCPU(axis, list, params);
    const bg = haveData
        ? calculateTotalBackground(axis, params, backgroundAnchors)
        : new Float64Array(axis.length);

    const scale = simScaleFactor(unscaled, bg, haveData);

    // The flat background and the noise are synthetic, so they apply only when
    // there is no measurement. With a file open the background under the curve
    // is the real one, and adding invented counts to it -- or scattering a
    // calculated curve that is being compared against data -- would make the
    // comparison meaningless.
    const flat  = haveData ? 0 : simConstantBackground();
    const level = haveData ? 0 : simNoiseLevel();
    const rnd   = level > 0 ? simRng(0x5EED) : null;

    // Reused across renders. Four thousand fresh {x,y} objects per keystroke
    // is the kind of garbage that shows up as the whole tab stuttering rather
    // than as any one slow function.
    if (!simPts || simPts.length !== axis.length) {
        simPts = new Array(axis.length);
        for (let i = 0; i < axis.length; i++) simPts[i] = { x: 0, y: 0 };
    }
    const pts = simPts;
    if (!simCleanY || simCleanY.length !== axis.length) simCleanY = new Float64Array(axis.length);
    for (let i = 0; i < axis.length; i++) {
        // Noise is applied to peaks PLUS background, because the detector
        // counts both -- a background of 10 carries its own scatter of about
        // 3, and that is what sets how weak a peak can be and still be seen.
        const clean = unscaled[i] * scale + bg[i] + flat;
        simCleanY[i] = clean;
        pts[i].x = axis[i];
        pts[i].y = rnd ? simNoisy(clean, level, rnd) : clean;
    }
    simNoiseApplied = level > 0;
    ds.data = pts;
    simShowLegendEntry(true);

    // The difference curve is a fit residual and has no meaning here unless
    // there is a measurement to take it against.
    const find = (l) => mainChart.data.datasets.find(d => d.label === l);
    if (haveData) {
        const diff = new Float64Array(axis.length);
        for (let i = 0; i < axis.length; i++) diff[i] = workingData.intensity[i] - pts[i].y;
        workingData.lastRawDifference = diff;
    } else {
        ['Experimental', 'Calculated', 'Background', 'Difference', 'Difference Zero']
            .forEach(l => { const d = find(l); if (d) d.data = []; });
    }
    const calc = find('Calculated'); if (calc) calc.data = [];

    simApplyAxes(haveData, pts);
    const notes = [];
    if (stepNote) notes.push(stepNote);
    if (capped) notes.push(`Tick marks show the ${SIM_MAX_MARKERS} strongest reflections; ` +
                           `${capped} weaker ones are omitted from the strip. The pattern uses all of them.`);
    simSetStatus(notes.join(' '), notes.length ? 'warn' : '');
    simSetSummary(sf, list, sym);
}

/**
 * The overall scale. Two different questions, so two different answers.
 *
 * With data on screen the useful scale is the one that puts the simulation
 * where the measurement is. With none there is nothing to scale to and the
 * absolute value is arbitrary (it absorbs the beam intensity, the sample
 * volume and the counting time, none of which exist), so the pattern is
 * normalised to a maximum the user sets.
 *
 * WHY THE TALLEST PEAK AND NOT LEAST SQUARES.
 *
 * A least-squares scale minimises the residual over every point, so it is
 * pulled by the whole pattern at once: a model missing a phase, or with one
 * badly wrong temperature factor, comes out uniformly short because the fit
 * splits the error across all the peaks. That is the correct thing to do when
 * the scale is a refined parameter and the model is nearly right. It is the
 * wrong thing here, where the scale is a viewing convenience and the model is
 * usually a first guess: what the eye wants is the tallest calculated peak
 * sitting on the tallest measured one, so the RELATIVE heights of everything
 * else can be read off against it. Anchoring on the maximum does exactly that,
 * and it does not move when a peak elsewhere is wrong.
 *
 * The background is subtracted from the measurement first, because the
 * simulated peaks are added on top of it further down -- comparing a net
 * calculated height against a gross measured one would scale the model down by
 * whatever the background happens to be under that peak.
 */
function simScaleFactor(unscaled, bg, haveData) {
    const el = document.getElementById('sim-scale-max');
    const auto = document.getElementById('sim-scale-auto');

    let peak = 0;
    for (let i = 0; i < unscaled.length; i++) if (unscaled[i] > peak) peak = unscaled[i];
    if (!(peak > 0)) return 1.0;

    if (haveData && auto && auto.checked) {
        let obs = 0;
        for (let i = 0; i < unscaled.length; i++) {
            const net = workingData.intensity[i] - bg[i];
            if (net > obs) obs = net;
        }
        return (obs > 0) ? obs / peak : 1.0;
    }

    const target = el ? parseFloat(el.value) : NaN;
    return (Number.isFinite(target) && target > 0 ? target : 10000) / peak;
}

/** Axis limits and labels. */
function simApplyAxes(haveData, pts) {
    const y = mainChart.options.scales.y;

    // drawTheoreticalPreview() leaves the y axis unlabelled and unticked,
    // because a stick pattern has no intensity scale. This one does.
    if (y.ticks) y.ticks.display = true;
    if (y.title) y.title.text = 'Intensity (a.u.)';

    if (haveData) {
        // rescalePlot() sizes the negative strip, draws the difference curve
        // from workingData.lastRawDifference and turns the Ka markers into
        // floating bars. All of that is data-mode machinery and it is already
        // correct; the simulation only had to fill in the difference.
        if (typeof updatePlotRange === 'function') updatePlotRange(true);
        if (typeof rescalePlot === 'function') rescalePlot(true);
    } else {
        let max = 0;
        for (const p of pts) if (p.y > max) max = p.y;
        if (!(max > 0)) max = 1;
        const { lo, hi } = simRange();
        mainChart.options.scales.x.min = lo;
        mainChart.options.scales.x.max = hi;
        y.min = -0.10 * max;
        y.max = 1.10 * max;
        mainChart.options.globalYMax = 1.10 * max;

        // With no data, rescalePlot() returns before it reaches the markers
        // (it needs workingData), so the ticks are sized here instead.
        simSizeMarkers(max);
    }
    mainChart.update('none');
}

// ---------------------------------------------------------------------------
//  THE ATOMS PANEL
//
//  The list in the tab is READ ONLY. Every add and every edit goes through the
//  modal below.
//
//  The first version put live inputs in the list, which read well as a spec and
//  badly in use: a site is six values that only mean anything together, and
//  editing them in place meant the pattern was recomputed against half-entered
//  states -- an atom on its way from (0,0,0) to (0.5,0.5,0.5) passes through
//  (0.5,0,0), which is a different Wyckoff position with a different
//  multiplicity, so the letter and the curve both flickered through answers
//  that were never asked for. A modal makes the commit explicit: nothing is
//  recomputed until Add or Save.
// ---------------------------------------------------------------------------

function simAtomIsUsable(a) {
    return !!a && !!a.label
        && Number.isFinite(a.x) && Number.isFinite(a.y) && Number.isFinite(a.z)
        && Number.isFinite(a.occ) && a.occ > 0
        && Number.isFinite(a.biso);
}

function simSetStatus(text, kind) {
    const el = document.getElementById('sim-status');
    if (!el) return;
    el.textContent = text || '';
    el.className = kind === 'warn' ? 'cf-warn' : 'cf-hint';
}

/** The line under the list: what is actually in the cell. */
function simSetSummary(sf, list, sym) {
    const el = document.getElementById('sim-summary');
    if (!el) return;
    if (!sf || !sf.sites || !sf.sites.length) { el.textContent = ''; return; }
    const parts = [`${sf.nAtoms} atom${sf.nAtoms === 1 ? '' : 's'} in the cell`,
                   `${list.length} reflection${list.length === 1 ? '' : 's'}`,
                   `${sym.ops.length} symmetry operator${sym.ops.length === 1 ? '' : 's'}`];
    el.textContent = parts.join(' · ') + '.';
    el.className = 'cf-hint';
    if (sf.warnings && sf.warnings.length) {
        el.textContent += ' ' + sf.warnings.join(' ');
        el.className = 'cf-warn';
    }
}

function simRemoveAtom(id) {
    simOpenRows.delete(id);
    const i = simAtoms.findIndex(a => a.id === id);
    if (i >= 0) simLog(`Simulation: removed atom ${simAtomLine(simAtoms[i], i)}`);
    simAtoms = simAtoms.filter(a => a.id !== id);
    simBuildRows();
    simRender();
}

/* PHYSICAL RANGES, ENFORCED.
   Occupancy is a fraction of a site: outside [0,1] it is not a small error but
   a different quantity. B is bounded above at 10 A^2 because past roughly that
   an atom's form factor has fallen so far by mid-angle that it stops
   contributing anything the pattern can show -- values above it are almost
   always a typo or a U written into a B box (U = B/8pi^2, so a plausible
   U of 0.006 typed as B is fine, but a B of 79 is a U mistaken for a B). */
const SIM_LIMITS = { occ: [0, 1], biso: [0, 10] };

const simClamp = (v, lo, hi) => (v < lo ? lo : v > hi ? hi : v);

const simFmt = (v, n) => (Number.isFinite(v) ? v.toFixed(n) : '—');

/**
 * Rebuild the atom list.
 *
 * EVERY FIELD IS LIVE. x, y, z, occupancy, B and the element are edited in
 * place and the pattern follows, which is the point: watching an intensity
 * change as an atom moves is how a structure gets understood, and it is not
 * something a modal that has to be opened and dismissed per edit can give.
 *
 * The earlier version made these read-only and put editing behind a dialog, to
 * stop the pattern being recomputed against half-typed coordinates. The
 * debounce below does that job without taking the feedback away.
 *
 * Full rebuild rather than a diff: the list is at most a few dozen entries and
 * it is only rebuilt when a row is added or removed. An edit updates the model
 * and the Wyckoff cell in place, so the box being typed into is never replaced
 * under the cursor.
 */
function simBuildRows() {
    const host = document.getElementById('sim-atom-rows');
    if (!host) return;
    host.innerHTML = '';

    if (!simAtoms.length) {
        const p = document.createElement('div');
        p.className = 'sg-explain-empty';
        p.textContent = 'No atoms yet. "Add atom" opens the picker, or load a CIF.';
        host.appendChild(p);
        simSetSummary(null);
        return;
    }

    const [oLo, oHi] = SIM_LIMITS.occ, [bLo, bHi] = SIM_LIMITS.biso;

    simAtoms.forEach((a, i) => {
        const card = document.createElement('div');
        const open = simOpenRows.has(a.id);
        card.className = 'sim-atom-card' + (open ? ' open' : '');
        card.dataset.id = a.id;
        card.innerHTML = `
            <div class="sim-atom-bar" role="button" tabindex="0"
                 aria-expanded="${open ? 'true' : 'false'}"
                 title="Click to ${open ? 'collapse' : 'expand'} this atom">
                <span class="sim-chev" aria-hidden="true">&rsaquo;</span>
                <span class="sim-atom-name">${simEsc(simSiteName(a, i))}</span>
                <span class="sim-atom-wy" title="Derived from the coordinates and the space group — the multiplicity the structure factor sums over."></span>
                <span class="sim-atom-sum"></span>
                <button type="button" class="sim-del" title="Remove this atom"
                        aria-label="Remove ${simEsc(simSiteName(a, i))}">&times;</button>
            </div>

            <div class="sim-atom-body">
                <div class="sim-atom-line">
                    <label class="sim-cap">Element</label>
                    <label class="sim-cap">Occupancy</label>
                    <label class="sim-cap">B<sub>iso</sub> (&Aring;&sup2;)</label>
                    <select class="control-select sim-el" data-field="label" data-value="${simEsc(a.label)}"></select>
                    <input type="number" class="control-input sim-num" data-field="occ"
                           step="0.01" min="${oLo}" max="${oHi}" value="${a.occ}">
                    <input type="number" class="control-input sim-num" data-field="biso"
                           step="0.05" min="${bLo}" max="${bHi}" value="${a.biso}">
                </div>
                <div class="sim-atom-line">
                    <label class="sim-cap">x</label>
                    <label class="sim-cap">y</label>
                    <label class="sim-cap">z</label>
                    <input type="number" class="control-input sim-num" data-field="x" step="0.001" value="${a.x}">
                    <input type="number" class="control-input sim-num" data-field="y" step="0.001" value="${a.y}">
                    <input type="number" class="control-input sim-num" data-field="z" step="0.001" value="${a.z}">
                </div>
                <div class="sim-atom-sp"></div>
            </div>`;
        host.appendChild(card);
    });

    simSyncExpandAll();
    simPopulateElementSelects();
    simRefreshWyckoff();
}

/** Open or close one card. */
function simSetRowOpen(card, open) {
    const id = parseInt(card.dataset.id, 10);
    if (open) simOpenRows.add(id); else simOpenRows.delete(id);
    card.classList.toggle('open', open);
    const bar = card.querySelector('.sim-atom-bar');
    if (bar) {
        bar.setAttribute('aria-expanded', open ? 'true' : 'false');
        bar.title = `Click to ${open ? 'collapse' : 'expand'} this atom`;
    }
    simSyncExpandAll();
}

/** Keep the Expand-all / Collapse-all button describing what it will do. */
function simSyncExpandAll() {
    const btn = document.getElementById('sim-expand-all');
    if (!btn) return;
    btn.style.display = simAtoms.length > 1 ? '' : 'none';
    const allOpen = simAtoms.length > 0 && simAtoms.every(a => simOpenRows.has(a.id));
    btn.textContent = allOpen ? 'Collapse all' : 'Expand all';
    btn.dataset.action = allOpen ? 'collapse' : 'expand';
}

/**
 * The site's display name: "O1", "Ti1".
 *
 * A CIF's own _atom_site_label is kept when there is one, because it is how
 * the structure is referred to in whatever the file came from and renaming it
 * breaks that link. Otherwise it is built from the element, with the charge
 * stripped -- "Ti4+1" reads as an oxidation state applied to a serial number.
 */
function simSiteName(a, i) {
    if (a.siteLabel) return a.siteLabel;
    const el = String(a.label || 'X').replace(/[^A-Za-z]/g, '') || 'X';
    return `${el}${i + 1}`;
}

/** Minimal escaping: labels come from the scattering table and from CIF files,
 *  and both end up in innerHTML. */
function simEsc(s) {
    return String(s).replace(/[&<>"]/g, c =>
        ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

/** Fill every inline element <select> from the loaded table. */
function simPopulateElementSelects() {
    const labels = simTableLabels();
    if (!labels.length) return;
    document.querySelectorAll('#sim-atom-rows select.sim-el').forEach(sel => {
        const want = sel.dataset.value || sel.value;
        sel.innerHTML = '';
        // A CIF can name a scatterer the table does not have. Dropping it from
        // the list would silently rewrite the structure to whatever sorts
        // first, so it is added as its own option and left selected.
        const opts = labels.includes(want) || !want ? labels : [want].concat(labels);
        const frag = document.createDocumentFragment();
        for (const l of opts) {
            const o = document.createElement('option');
            o.value = l; o.textContent = l;
            frag.appendChild(o);
        }
        sel.appendChild(frag);
        if (want) sel.value = want;
    });
}

/** The Wyckoff label for every atom, in list order. */
function simWyckoffLabels() {
    const sym = currentSG ? simSymmetry(currentSG) : null;
    let params = null;
    try { params = getAllParams(); } catch (e) { /* not ready yet */ }
    if (!sym || sym.error || !params || !(params.a > 0)) return simAtoms.map(() => '—');
    const G = simMetric(params);
    return simAtoms.map(a => {
        if (!Number.isFinite(a.x) || !Number.isFinite(a.y) || !Number.isFinite(a.z)) return '—';
        return simWyckoffFor(a.x, a.y, a.z, sym, G).label;
    });
}

/**
 * Recompute the Wyckoff cells in place.
 *
 * In place, not by rebuilding the rows: this runs on every keystroke in a
 * coordinate box, and replacing the row would take the focus and the caret
 * with it.
 *
 * The letter depends on the cell and the space group as well as the
 * coordinates, so it also goes stale when either of those changes -- which is
 * why handleSpaceGroupApplied() calls this.
 */
function simRefreshWyckoff() {
    const sym = currentSG ? simSymmetry(currentSG) : null;
    let params = null;
    try { params = getAllParams(); } catch (e) { /* not ready yet */ }
    const ok = sym && !sym.error && params && params.a > 0;
    const G = ok ? simMetric(params) : null;

    document.querySelectorAll('#sim-atom-rows .sim-atom-card').forEach((card, i) => {
        const a = simAtoms[i];
        if (!a) return;
        const wyEl = card.querySelector('.sim-atom-wy');
        const spEl = card.querySelector('.sim-atom-sp');
        if (wyEl) {
            if (!ok || ![a.x, a.y, a.z].every(Number.isFinite)) {
                wyEl.textContent = '—';
                if (spEl) spEl.textContent = '';
            } else {
                const w = simWyckoffFor(a.x, a.y, a.z, sym, G);
                wyEl.textContent = w.label;
                if (spEl) {
                    spEl.textContent = w.siteSym
                        ? `Site symmetry ${w.siteSym} · ${w.mult} atom${w.mult === 1 ? '' : 's'} in the cell`
                        : `${w.mult} atom${w.mult === 1 ? '' : 's'} in the cell`;
                }
            }
        }
        simUpdateSummary(card, a);
    });
}

/**
 * The one-line value summary shown on the collapsed row.
 *
 * WITHOUT IT THE COLLAPSED LIST IS A LIST OF NAMES, which is the wrong trade:
 * the reason to fold a twenty-atom structure down is to see all of it at once,
 * and a column of labels shows less than the expanded form did. The
 * coordinates are what gets compared down a list, so they stay visible;
 * occupancy and B appear only when they are not the ordinary values, so that
 * anything unusual stands out instead of being lost in a row of "1.000 0.500"
 * repeated twenty times.
 */
function simUpdateSummary(card, a) {
    const el = card.querySelector('.sim-atom-sum');
    if (!el) return;
    const bits = [`${simFmt(a.x, 4)} ${simFmt(a.y, 4)} ${simFmt(a.z, 4)}`];
    if (Number.isFinite(a.occ) && Math.abs(a.occ - 1) > 1e-9) bits.push(`occ ${simFmt(a.occ, 3)}`);
    if (Number.isFinite(a.biso) && Math.abs(a.biso - 0.5) > 1e-9) bits.push(`B ${simFmt(a.biso, 2)}`);
    el.textContent = bits.join('  ');
    el.title = `${a.label}  x=${simFmt(a.x, 5)} y=${simFmt(a.y, 5)} z=${simFmt(a.z, 5)} ` +
               `occ=${simFmt(a.occ, 3)} B=${simFmt(a.biso, 3)}`;
}

/**
 * Read one edited field back into the model.
 *
 * OUT-OF-RANGE VALUES ARE CLAMPED FOR THE CALCULATION IMMEDIATELY, AND THE BOX
 * IS REWRITTEN ONLY ON COMMIT.
 *
 * Rewriting the box on every keystroke fights the typist: clearing it to type
 * "0.75" passes through the empty string, and typing "10" into a field capped
 * at 10 passes through "1". So the value the pattern uses is always inside the
 * physical range, while the text stays whatever is being typed until focus
 * leaves or Enter is pressed, at which point the box is corrected to match
 * what was actually used. The two can never disagree once the user stops
 * typing, and they can never disagree about the PATTERN at any point.
 */
function simOnRowInput(e) {
    const el = e.target;
    const card = el.closest && el.closest('.sim-atom-card');
    if (!card) return;
    const a = simAtoms.find(x => String(x.id) === card.dataset.id);
    if (!a) return;
    const field = el.dataset.field;
    if (!field) return;

    if (field === 'label') {
        a.label = el.value;
        el.dataset.value = el.value;
        // An auto-generated name is built from the element, so it has to follow
        // when the element changes. A name that came from a CIF does not.
        const nameEl = card.querySelector('.sim-atom-name');
        if (nameEl && !a.siteLabel) {
            nameEl.textContent = simSiteName(a, simAtoms.indexOf(a));
        }
    } else {
        let v = parseFloat(el.value);
        const lim = SIM_LIMITS[field];
        if (Number.isFinite(v) && lim) {
            const c = simClamp(v, lim[0], lim[1]);
            // Only on commit, so typing is not interrupted. See above.
            if (c !== v && e.type === 'change') el.value = String(c);
            v = c;
        }
        a[field] = v;
        el.classList.toggle('sim-bad', !Number.isFinite(v));
    }

    simUpdateSummary(card, a);
    simRefreshWyckoff();

    // ON COMMIT ONLY. 'input' fires per keystroke, and a log with a line for
    // every character typed into a coordinate box is not a record of what the
    // user did, it is a record of them typing -- and it would bury the entries
    // that matter. 'change' fires on blur or Enter, which is where the value
    // actually settled.
    if (e.type === 'change') {
        simLog(`Simulation: edited atom ${simAtomLine(a, simAtoms.indexOf(a))}`);
    }
    simRenderDebounced();
}

// ---------------------------------------------------------------------------
//  THE ATOM EDITOR (modal)
// ---------------------------------------------------------------------------

/* The modal is now ADD-ONLY. Every existing atom is edited in place in the
   list, so the only thing the dialog is still needed for is choosing a
   scatterer out of the ~200 in the table, which no inline control does well. */
/** The scatterer label currently chosen in the picker. */
let simPickLabel = '';

function simModalEl(id) { return document.getElementById(id); }

/**
 * Open the editor. `id` null adds, otherwise edits that atom.
 */
function simOpenAtomModal() {
    const overlay = simModalEl('sim-atom-modal-overlay');
    if (!overlay) return;
    if (!simTables) { simEnsureTables(); }

    simPickLabel = '';
    simModalEl('sim-el-search').value = '';
    simModalEl('sim-f-x').value    = 0;
    simModalEl('sim-f-y').value    = 0;
    simModalEl('sim-f-z').value    = 0;
    simModalEl('sim-f-occ').value  = 1;
    simModalEl('sim-f-biso').value = 0.5;

    simFilterElements();
    simUpdateModalWyckoff();
    simModalError('');

    overlay.classList.add('open');
    setTimeout(() => {
        const f = simModalEl('sim-el-search');
        if (f) { f.focus(); f.select(); }
    }, 0);
}

function simCloseAtomModal() {
    const overlay = simModalEl('sim-atom-modal-overlay');
    if (overlay) overlay.classList.remove('open');
}

function simModalError(msg) {
    const el = simModalEl('sim-atom-error');
    if (!el) return;
    el.textContent = msg || '';
    el.style.display = msg ? 'block' : 'none';
}

/**
 * Filter the scatterer list against what has been typed.
 *
 * Prefix match on the label, so "T" gives Ti, Ti2+, Ti3+, Ti4+, Ta, Tb, Tc,
 * Te, Th, Tl, Tm -- every scatterer whose name starts that way, neutral atoms
 * before their own ions and the whole thing in atomic-number order, which is
 * the order a chemist scanning for an element expects. A substring match was
 * the alternative and is worse: typing "O" would surface Co, Mo, Ho and Po
 * above oxygen.
 */
function simFilterElements() {
    const box = simModalEl('sim-el-list');
    const q = (simModalEl('sim-el-search').value || '').trim().toLowerCase();
    if (!box) return;

    if (!simTables) {
        box.innerHTML = `<div class="sg-list-empty">${simEsc(simTablesError || 'Loading scattering factors…')}</div>`;
        return;
    }

    const table = simTables.xray.table;
    const all = simTableLabels();
    const hits = q ? all.filter(l => l.toLowerCase().startsWith(q)) : all;

    if (!hits.length) {
        box.innerHTML = `<div class="sg-list-empty">Nothing in the table starts with “${simEsc(q)}”.</div>`;
        return;
    }

    // An exact typed label selects itself, so a user who knows what they want
    // can type "Ti4+" and press Add without touching the list.
    const exact = hits.find(l => l.toLowerCase() === q);
    if (exact) simPickLabel = exact;
    else if (!hits.includes(simPickLabel)) simPickLabel = hits[0];

    box.innerHTML = hits.slice(0, 400).map(l => {
        const Z = (table[l] && table[l].Z) || '';
        const on = (l === simPickLabel) ? ' selected' : '';
        return `<div class="sg-list-item sim-el-item${on}" data-label="${simEsc(l)}">` +
               `<span class="num">${Z}</span><span>${simEsc(l)}</span>` +
               `<span class="sys">${simIsIon(l) ? 'ion' : 'atom'}</span></div>`;
    }).join('');

    const sel = box.querySelector('.sim-el-item.selected');
    if (sel && sel.scrollIntoView) sel.scrollIntoView({ block: 'nearest' });
}

function simIsIon(label) { return !/^[A-Za-z]+$/.test(label); }

/**
 * Show, live, which Wyckoff position the coordinates in the dialog land on.
 *
 * This is the reason the field is derived rather than typed: the user finds out
 * what the symmetry does to their coordinates BEFORE committing, instead of
 * asserting a position and discovering the mismatch in the pattern.
 */
function simUpdateModalWyckoff() {
    const out = simModalEl('sim-atom-wyck');
    if (!out) return;
    const v = simReadModalFields();
    const sym = currentSG ? simSymmetry(currentSG) : null;
    let params = null;
    try { params = getAllParams(); } catch (e) { /* not ready */ }

    if (!sym || sym.error || !params || !(params.a > 0)) { out.textContent = ''; return; }
    if (![v.x, v.y, v.z].every(Number.isFinite)) { out.textContent = ''; return; }

    const w = simWyckoffFor(v.x, v.y, v.z, sym, simMetric(params));
    out.innerHTML = `Wyckoff <b>${simEsc(w.label)}</b>` +
        (w.siteSym ? ` · site symmetry ${simEsc(w.siteSym)}` : '') +
        ` · ${w.mult} atom${w.mult === 1 ? '' : 's'} in the cell`;
}

function simReadModalFields() {
    const num = (id) => parseFloat(simModalEl(id).value);
    return {
        label: simPickLabel,
        x: num('sim-f-x'), y: num('sim-f-y'), z: num('sim-f-z'),
        occ: num('sim-f-occ'), biso: num('sim-f-biso')
    };
}

/** Validate and commit. */
function simCommitAtom() {
    const v = simReadModalFields();

    if (!v.label) { simModalError('Choose a scatterer from the list.'); return; }
    if (simTables && !simTables.xray.table[v.label]) {
        simModalError(`“${v.label}” is not in the scattering table.`); return;
    }
    for (const [k, name] of [['x', 'x'], ['y', 'y'], ['z', 'z']]) {
        if (!Number.isFinite(v[k])) { simModalError(`${name} must be a number.`); return; }
    }
    // The same limits the inline boxes enforce, so an atom cannot enter the
    // list through the dialog in a state the list would not allow.
    const [oLo, oHi] = SIM_LIMITS.occ, [bLo, bHi] = SIM_LIMITS.biso;
    if (!Number.isFinite(v.occ) || v.occ < oLo || v.occ > oHi) {
        simModalError(`Occupancy must be between ${oLo} and ${oHi}.`); return;
    }
    if (!Number.isFinite(v.biso) || v.biso < bLo || v.biso > bHi) {
        simModalError(`B must be between ${bLo} and ${bHi} A^2.`); return;
    }

    simAtoms.push({ id: simNextId++, siteLabel: '', ...v });
    // Just created, so almost certainly about to be adjusted: leave it open.
    simOpenRows.add(simAtoms[simAtoms.length - 1].id);
    simLog(`Simulation: added atom ${simAtomLine(simAtoms[simAtoms.length - 1], simAtoms.length - 1)}`);
    simCloseAtomModal();
    simBuildRows();
    simRender();
}

// ---------------------------------------------------------------------------
//  CIF IN
//
//  A structure reader, not a pattern reader. data_io.js already has
//  parsePdCifFile() for powder profiles; this reads the _atom_site_ loop and
//  the cell out of an ordinary crystal-structure CIF, which is a different
//  half of the same format.
// ---------------------------------------------------------------------------

/** Strip a CIF standard uncertainty: "0.1234(5)" is a value, not a syntax error. */
function simCifNum(tok) {
    if (tok == null) return NaN;
    const s = String(tok).trim().replace(/\(\d+\)\s*$/, '');
    if (s === '' || s === '?' || s === '.') return NaN;
    const v = parseFloat(s);
    return Number.isFinite(v) ? v : NaN;
}

/** Split a CIF data line into tokens, respecting quotes. */
function simCifTokens(line) {
    const out = [];
    const re = /'([^']*)'|"([^"]*)"|(\S+)/g;
    let m;
    while ((m = re.exec(line)) !== null) {
        out.push(m[1] !== undefined ? m[1] : m[2] !== undefined ? m[2] : m[3]);
    }
    return out;
}

/**
 * Parse a crystal-structure CIF.
 *
 * Only the first data block is read. Multi-block CIFs are common (a deposition
 * with several phases), and silently merging their atom lists would produce a
 * structure that exists in none of them; taking the first and saying so is
 * recoverable, and the user can split the file.
 *
 * @returns {{atoms: object[], cell: object|null, sgName: string|null,
 *            sgNumber: number|null, blocks: number, blockName: string}}
 */
function simParseStructureCif(text) {
    const rawLines = String(text).split(/\r?\n/);

    // Fold semicolon text fields into single tokens so they cannot be mistaken
    // for tags or loop rows. Their content is never needed here, but a
    // multi-line publication title containing the word "loop_" would otherwise
    // derail the parser.
    const lines = [];
    let inText = false;
    for (const raw of rawLines) {
        if (raw.startsWith(';')) { inText = !inText; lines.push(inText ? "';'" : ''); continue; }
        lines.push(inText ? '' : raw);
    }

    const tags = {};
    const atoms = [];
    let blocks = 0, blockName = '';
    let stop = false;

    for (let i = 0; i < lines.length && !stop; i++) {
        const line = lines[i];
        const t = line.trim();
        if (!t || t.startsWith('#')) continue;

        if (/^data_/i.test(t)) {
            blocks++;
            if (blocks === 1) blockName = t.slice(5);
            else { stop = true; break; }   // first block only
            continue;
        }

        if (/^loop_$/i.test(t)) {
            // Collect the header, then the rows.
            const headers = [];
            let j = i + 1;
            for (; j < lines.length; j++) {
                const h = lines[j].trim();
                if (!h || h.startsWith('#')) continue;
                if (h.startsWith('_')) { headers.push(h.split(/\s+/)[0].toLowerCase()); continue; }
                break;
            }
            const rows = [];
            for (; j < lines.length; j++) {
                const r = lines[j].trim();
                if (!r || r.startsWith('#')) continue;
                if (r.startsWith('_') || /^loop_$/i.test(r) || /^data_/i.test(r)) break;
                const toks = simCifTokens(r);
                if (toks.length) rows.push(toks);
            }
            i = j - 1;

            const col = (name) => headers.indexOf(name);
            const cx = col('_atom_site_fract_x');
            if (cx >= 0) {
                const cy = col('_atom_site_fract_y'), cz = col('_atom_site_fract_z');
                const cl = col('_atom_site_label'), ct = col('_atom_site_type_symbol');
                const co = col('_atom_site_occupancy');
                const cb = col('_atom_site_b_iso_or_equiv');
                const cu = col('_atom_site_u_iso_or_equiv');

                // Rows can wrap across lines in badly written files; only take
                // rows that actually have enough tokens rather than silently
                // reading a coordinate out of the next atom.
                const need = Math.max(cx, cy, cz, cl, ct, co, cb, cu) + 1;
                for (const r of rows) {
                    if (r.length < need) continue;
                    const rec = {
                        siteLabel: cl >= 0 ? r[cl] : '',
                        type: ct >= 0 ? r[ct] : '',
                        x: simCifNum(r[cx]),
                        y: cy >= 0 ? simCifNum(r[cy]) : NaN,
                        z: cz >= 0 ? simCifNum(r[cz]) : NaN,
                        occ: co >= 0 ? simCifNum(r[co]) : 1,
                        biso: cb >= 0 ? simCifNum(r[cb]) : NaN,
                        uiso: cu >= 0 ? simCifNum(r[cu]) : NaN
                    };
                    if ([rec.x, rec.y, rec.z].every(Number.isFinite)) atoms.push(rec);
                }
            }
            continue;
        }

        if (t.startsWith('_')) {
            const toks = simCifTokens(t);
            if (toks.length >= 2) tags[toks[0].toLowerCase()] = toks.slice(1).join(' ');
            continue;
        }
    }

    const g = (...names) => {
        for (const n of names) if (tags[n] !== undefined) return tags[n];
        return undefined;
    };
    const cellA = simCifNum(g('_cell_length_a'));
    const cell = Number.isFinite(cellA) ? {
        a: cellA,
        b: simCifNum(g('_cell_length_b')),
        c: simCifNum(g('_cell_length_c')),
        alpha: simCifNum(g('_cell_angle_alpha')),
        beta: simCifNum(g('_cell_angle_beta')),
        gamma: simCifNum(g('_cell_angle_gamma'))
    } : null;

    const sgName = g('_space_group_name_h-m_alt', '_symmetry_space_group_name_h-m',
                     '_space_group_name_hall', '_symmetry_space_group_name_hall') || null;
    const sgNum = simCifNum(g('_space_group_it_number', '_symmetry_int_tables_number'));

    return { atoms, cell, sgName, sgNumber: Number.isFinite(sgNum) ? sgNum : null,
             blocks, blockName };
}

/**
 * Map a CIF site onto a label in the scattering table.
 *
 * _atom_site_type_symbol is already in the table's spelling for ions ("Ti4+",
 * "O2-"), so it is tried verbatim first. Failing that the charge is dropped,
 * because a table without Ti4+ still scatters perfectly well as Ti and a site
 * that contributes nothing is a worse answer than one with a slightly wrong
 * form factor. The site LABEL is the last resort: "Ca1" carries the element
 * only by convention, and that convention is broken often enough ("M1", "OW")
 * that guessing from it silently would be worse than reporting the site as
 * unidentified.
 */
function simResolveScatterer(rec) {
    const table = (simTables && simTables.xray && simTables.xray.table) || {};
    const tryLabel = (s) => (s && table[s]) ? s : null;

    const type = (rec.type || '').trim();
    if (type) {
        const exact = tryLabel(type);
        if (exact) return { label: exact, note: null };
        const bare = type.replace(/[0-9+\-.]/g, '');
        const b = tryLabel(bare) || tryLabel(bare.charAt(0).toUpperCase() + bare.slice(1).toLowerCase());
        if (b) return { label: b, note: `${type} is not tabulated; used ${b}` };
    }

    const lab = (rec.siteLabel || '').trim();
    const m = lab.match(/^([A-Za-z]{1,2})/);
    if (m) {
        // Two letters then one. Labels like "OW" (water oxygen) and "CA"
        // (a carbon, not calcium) are both common, and "Ow" is not an element
        // while "O" is -- so the longer reading is tried first and the shorter
        // one catches what it misses.
        const cap = (s) => s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
        for (const cand of [cap(m[1]), cap(m[1].charAt(0))]) {
            const c = tryLabel(cand);
            if (c) return { label: c, note: `no type symbol for ${lab}; read as ${c}` };
        }
    }
    return { label: null, note: `could not identify the scatterer for ${lab || type || 'an atom'}` };
}

/** Read a CIF into the atom list. */
function simLoadCifText(text, fileName) {
    if (!simTables) {
        simSetStatus('The scattering factors are still loading. Try again in a moment.', 'warn');
        simEnsureTables();
        return;
    }

    let parsed;
    try { parsed = simParseStructureCif(text); }
    catch (e) { simSetStatus(`Could not read the CIF: ${e.message || e}`, 'warn'); return; }

    if (!parsed.atoms.length) {
        simSetStatus('No _atom_site_fract_x loop in that file, so there is no structure to load.', 'warn');
        return;
    }

    const notes = [];
    const next = [];
    for (const rec of parsed.atoms) {
        const r = simResolveScatterer(rec);
        if (r.note) notes.push(r.note);
        if (!r.label) continue;

        // B, not U. Files carry one or the other; B = 8 pi^2 U.
        let biso = rec.biso;
        if (!Number.isFinite(biso) && Number.isFinite(rec.uiso)) biso = 8 * Math.PI * Math.PI * rec.uiso;
        if (!Number.isFinite(biso)) biso = 0.5;

        let occ = Number.isFinite(rec.occ) ? rec.occ : 1;
        if (!(occ > 0)) occ = 1;
        if (occ > 1) occ = 1;

        next.push({ id: simNextId++, label: r.label, siteLabel: rec.siteLabel || '',
                    x: rec.x, y: rec.y, z: rec.z, occ, biso: Math.max(0, biso) });
    }

    if (!next.length) {
        simSetStatus('None of the sites in that file could be matched to a scatterer. ' +
                     notes.join('; '), 'warn');
        return;
    }

    // THE CELL AND THE SPACE GROUP COME WITH THE ATOMS.
    //
    // Fractional coordinates are fractions OF A PARTICULAR CELL under A
    // PARTICULAR set of operators; they mean nothing apart from those two. So
    // the whole triple is applied together. The space group goes first,
    // because it decides which lattice boxes exist -- a cubic group leaves
    // only `a`, and writing b, c or the angles before the system is set puts
    // them into inputs that are about to be destroyed.
    const said = [];
    let sgApplied = false;

    if ((parsed.sgNumber || parsed.sgName) && root.SG_ENGINE &&
        typeof root.handleSpaceGroupApplied === 'function') {
        try {
            const want = String(parsed.sgName || '').replace(/\s+/g, '');
            let rec = null;
            if (parsed.sgNumber) {
                const settings = SG_ENGINE.listSettings(parsed.sgNumber) || [];
                // A NUMBER DOES NOT NAME A SETTING. Space group 14 has six of
                // them, with different operators and different absences from
                // the same cell, so the H-M string picks which one and the
                // number only narrows the search. Falling straight to
                // settings[0] would silently load P21/c for a file that said
                // P21/n -- a structure that still produces a plausible pattern.
                rec = settings.find(r => String(r.symbol || '').replace(/\s+/g, '') === want)
                   || settings.find(r => String(r.standard_symbol || '').replace(/\s+/g, '') === want)
                   || settings.find(r => String(r.hall || '').replace(/\s+/g, '') === want)
                   || settings[0] || null;
                if (rec && want && String(rec.symbol || '').replace(/\s+/g, '') !== want) {
                    said.push(`the file's symbol "${parsed.sgName}" did not match a known setting of No. ${parsed.sgNumber}; used ${rec.symbol}`);
                }
            } else if (want) {
                // No number: search every group for a setting with this symbol.
                for (let n = 1; n <= 230 && !rec; n++) {
                    const s = SG_ENGINE.listSettings(n) || [];
                    rec = s.find(r => String(r.symbol || '').replace(/\s+/g, '') === want)
                       || s.find(r => String(r.standard_symbol || '').replace(/\s+/g, '') === want)
                       || null;
                }
            }
            if (rec) {
                root.handleSpaceGroupApplied(rec);
                sgApplied = true;
                said.push(`space group ${rec.symbol} (No. ${rec.number})`);
            } else {
                said.push(`space group "${parsed.sgName || parsed.sgNumber}" was not found in the database and was left unchanged`);
            }
        } catch (e) {
            said.push('the space group could not be applied: ' + (e.message || e));
        }
    }

    if (parsed.cell) {
        let n = 0;
        document.querySelectorAll('#lattice-parameters-container input[type="number"]')
            .forEach(input => {
                if (!input.id.startsWith('lattice-param-')) return;
                const name = input.id.replace('lattice-param-', '');
                if (Number.isFinite(parsed.cell[name])) { input.value = parsed.cell[name].toFixed(5); n++; }
            });
        if (n) said.push(`${n} lattice parameter${n === 1 ? '' : 's'}`);
        else if (!sgApplied) said.push('the cell could not be applied: no crystal system is active');
    }

    // REPLACES the list rather than appending. A CIF is one whole structure,
    // and merging it into whatever was already there would build a phase that
    // exists in neither file.
    simAtoms = next;
    simNextId = simAtoms.length + 1;
    // A loaded structure arrives collapsed. Expanding fifty sites at once
    // would bury the thing the user actually wants to see, which is that the
    // structure loaded and how many sites it has.
    simOpenRows.clear();

    // The cell and the group both changed, so every cached quantity keyed on
    // them is stale: the operator list, the reflection list, and the hkl
    // generator's own cache.
    simSymCache.clear();
    simInvalidateGeometry();
    if (typeof root.invalidateHklCache === 'function') root.invalidateHklCache();

    simBuildRows();
    simRender();

    const extra = [];
    if (parsed.blocks > 1) {
        extra.push(`the file has ${parsed.blocks} data blocks; only the first (${parsed.blockName}) was read`);
    }
    if (notes.length) extra.push(notes.join('; '));

    const head = `Loaded ${next.length} site${next.length === 1 ? '' : 's'} from ${fileName || 'the CIF'}` +
                 (said.length ? `, with ${said.join(', ')}` : '') + '.';
    if (extra.length) simSetStatus(`${head} Note: ${extra.join('. ')}.`, 'warn');
    else simSetStatus(head, '');

    simLog(`Simulation: loaded ${next.length} atom site(s) from ${fileName || 'a CIF'}` +
           (said.length ? `, with ${said.join(', ')}` : '') + '.');
    simAtoms.forEach((a, i) => simLog(`Simulation:   ${simAtomLine(a, i)}`));
}

// ---------------------------------------------------------------------------
//  pdCIF OUT
// ---------------------------------------------------------------------------

/**
 * Write the simulated pattern and its structure as a pdCIF.
 *
 * This is what the report button produces on the Plot tab while the Method is
 * Simulation. A PDF of a calculated curve is a picture of something the
 * program can state exactly, and the numbers are the point of a simulation --
 * they are what gets compared against a measurement, fed to another program,
 * or checked by hand. Every other mode still prints its PDF; see
 * printActiveTab().
 */
function simExportPdCif() {
    const say = (m, kind) => (typeof showToast === 'function')
        ? showToast(m, kind || 'error') : console.warn(m);

    if (!simGeomList || !simPts || !simPts.length) {
        say('There is no simulated pattern yet. Add at least one atom.');
        return;
    }
    if (typeof exportPattern !== 'function' || typeof downloadTextFile !== 'function') {
        say('The export writer is not available.');
        return;
    }

    const params = getAllParams();
    const haveData = (typeof hasExperimentalData === 'function') && hasExperimentalData();
    const wy = simWyckoffLabels();
    const sym = currentSG ? simSymmetry(currentSG) : null;
    const G = (sym && !sym.error) ? simMetric(params) : null;

    const tth = simPts.map(p => p.x);

    // WITH NOISE ON, THE PATTERN IS A SYNTHETIC OBSERVATION.
    //
    // That is what noise means: a curve with counting statistics on it is not
    // a calculation, it is what a detector would have recorded. Writing it as
    // the observed column, with the noise-free curve alongside as the
    // calculated one, is what makes the file useful for its stated purpose --
    // another program can fit the observation and be checked against the model
    // that produced it. Without noise the pattern is exactly a calculation and
    // is written as one. writePdCIF() switches the angle tag to match, and
    // labels the block so nobody mistakes the counts for a real measurement.
    let obs, calc, synthetic = false;
    if (haveData) {
        obs = Array.from(workingData.intensity);
        calc = simPts.map(p => p.y);
    } else if (simNoiseApplied && simCleanY) {
        obs = simPts.map(p => p.y);
        calc = Array.from(simCleanY);
        synthetic = true;
    } else {
        obs = null;
        calc = simPts.map(p => p.y);
    }
    const bkgDs = mainChart && mainChart.data.datasets.find(d => d.label === 'Background');
    const bkg = (haveData && bkgDs && bkgDs.data.length === tth.length)
        ? bkgDs.data.map(p => p.y) : null;

    const bundle = {
        sourceName: 'simulation',
        tth, obs, calc, bkg, syntheticObs: synthetic,
        params,
        spaceGroup: currentSG || null,
        hklList: simGeomList,
        backgroundAnchors: null,
        stats: null,
        esds: null,
        scaleFactor: 1,
        atoms: simAtoms.map((a, i) => {
            const w = (G && sym) ? simWyckoffFor(a.x, a.y, a.z, sym, G) : null;
            return {
                label: a.label,
                siteLabel: a.siteLabel || `${String(a.label).replace(/[^A-Za-z]/g, '')}${i + 1}`,
                x: a.x, y: a.y, z: a.z, occ: a.occ, biso: a.biso,
                wyckoff: wy[i] || '?',
                multiplicity: w ? w.mult : NaN
            };
        }),
        polarisationUi: {
            mode: (document.getElementById('polarization-mode') || {}).value || null,
            monoTth: parseFloat((document.getElementById('mono-2theta') || {}).value),
            fraction: parseFloat((document.getElementById('pol-fraction') || {}).value)
        }
    };

    try {
        const { text, ext } = exportPattern('pdcif', bundle);
        const stamp = new Date().toISOString().slice(0, 19).replace(/[:T]/g, '-');
        const name = `simulation-${stamp}.${ext}`;
        downloadTextFile(text, name);
        say(`Saved ${name}`, 'success');
        simLog(`Simulation: wrote ${name} (${tth.length} points, ` +
               `${simGeomList.length} reflections, ${simAtoms.length} sites)`);
    } catch (e) {
        say('Could not write the pdCIF: ' + (e.message || e));
    }
}

// ---------------------------------------------------------------------------
//  WIRING
// ---------------------------------------------------------------------------

let simDebounceTimer = null;
function simRenderDebounced() {
    if (simDebounceTimer !== null) clearTimeout(simDebounceTimer);
    simDebounceTimer = setTimeout(() => { simDebounceTimer = null; simRender(); },
                                  SIM_DEBOUNCE_MS);
}

function simInit() {
    if (simUiReady) return;
    simUiReady = true;

    // --- the read-only list: Edit and Remove ----------------------------
    const host = document.getElementById('sim-atom-rows');
    if (host) {
        // Delegated, so rows added later are wired automatically.
        // 'input' gives the live update; 'change' additionally commits the
        // clamp back into the box. See simOnRowInput().
        host.addEventListener('input', simOnRowInput);
        host.addEventListener('change', simOnRowInput);
        host.addEventListener('click', (e) => {
            const card = e.target.closest('.sim-atom-card');
            if (!card) return;
            // Remove is checked FIRST: it sits inside the header, so a plain
            // "did the click land in the header" test would also collapse the
            // row on its way to deleting it.
            if (e.target.closest('.sim-del')) {
                simRemoveAtom(parseInt(card.dataset.id, 10));
                return;
            }
            if (e.target.closest('.sim-atom-bar')) {
                simSetRowOpen(card, !card.classList.contains('open'));
            }
        });
        // The header is a div, so it gets the keyboard behaviour a button
        // would have had. Space and Enter both toggle; Space is prevented from
        // scrolling the panel.
        host.addEventListener('keydown', (e) => {
            if (e.key !== 'Enter' && e.key !== ' ') return;
            const bar = e.target.closest && e.target.closest('.sim-atom-bar');
            if (!bar) return;
            e.preventDefault();
            const card = bar.closest('.sim-atom-card');
            if (card) simSetRowOpen(card, !card.classList.contains('open'));
        });
    }

    const add = document.getElementById('sim-add-atom');
    if (add) add.addEventListener('click', () => simOpenAtomModal());

    const expand = document.getElementById('sim-expand-all');
    if (expand) expand.addEventListener('click', () => {
        const wantOpen = expand.dataset.action !== 'collapse';
        simOpenRows.clear();
        if (wantOpen) simAtoms.forEach(a => simOpenRows.add(a.id));
        document.querySelectorAll('#sim-atom-rows .sim-atom-card')
                .forEach(card => simSetRowOpen(card, wantOpen));
        simSyncExpandAll();
    });

    // --- Load CIF --------------------------------------------------------
    const loadBtn = document.getElementById('sim-load-cif');
    const fileIn = document.getElementById('sim-cif-file');
    if (loadBtn && fileIn) {
        loadBtn.addEventListener('click', () => fileIn.click());
        fileIn.addEventListener('change', () => {
            const f = fileIn.files && fileIn.files[0];
            if (!f) return;
            const fr = new FileReader();
            fr.onload = () => simLoadCifText(String(fr.result || ''), f.name);
            fr.onerror = () => simSetStatus(`Could not read ${f.name}.`, 'warn');
            fr.readAsText(f);
            // Cleared so that re-picking the SAME file fires 'change' again --
            // otherwise a user who edits a CIF and reloads it gets silence.
            fileIn.value = '';
        });
    }

    // --- the modal -------------------------------------------------------
    const ok = document.getElementById('sim-atom-ok');
    if (ok) ok.addEventListener('click', simCommitAtom);
    const cancel = document.getElementById('sim-atom-cancel');
    if (cancel) cancel.addEventListener('click', simCloseAtomModal);

    const overlay = document.getElementById('sim-atom-modal-overlay');
    if (overlay) overlay.addEventListener('click', (e) => {
        if (e.target === overlay) simCloseAtomModal();   // click outside = cancel
    });

    const search = document.getElementById('sim-el-search');
    if (search) {
        search.addEventListener('input', () => { simFilterElements(); simModalError(''); });
        // Enter in the search box commits, so adding an atom can be done
        // without reaching for the mouse: type, tab through the numbers, Enter.
        search.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') { e.preventDefault(); simCommitAtom(); }
        });
    }

    const list = document.getElementById('sim-el-list');
    if (list) list.addEventListener('click', (e) => {
        const item = e.target.closest('.sim-el-item');
        if (!item) return;
        simPickLabel = item.dataset.label;
        document.getElementById('sim-el-search').value = simPickLabel;
        simFilterElements();
        simModalError('');
    });

    const mOcc = document.getElementById('sim-f-occ');
    const mB = document.getElementById('sim-f-biso');
    if (mOcc) { mOcc.min = SIM_LIMITS.occ[0]; mOcc.max = SIM_LIMITS.occ[1]; }
    if (mB) { mB.min = SIM_LIMITS.biso[0]; mB.max = SIM_LIMITS.biso[1]; }

    ['sim-f-x', 'sim-f-y', 'sim-f-z'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener('input', simUpdateModalWyckoff);
    });
    ['sim-f-x', 'sim-f-y', 'sim-f-z', 'sim-f-occ', 'sim-f-biso'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') { e.preventDefault(); simCommitAtom(); }
        });
    });

    document.addEventListener('keydown', (e) => {
        if (e.key !== 'Escape') return;
        const ov = document.getElementById('sim-atom-modal-overlay');
        if (ov && ov.classList.contains('open')) simCloseAtomModal();
    });

    // --- the pattern controls --------------------------------------------
    ['sim-step', 'sim-background', 'sim-noise', 'sim-scale-max', 'sim-scale-auto'].forEach(id => {
        const el = document.getElementById(id);
        if (!el) return;
        el.addEventListener('input', simRenderDebounced);
        el.addEventListener('change', simRender);
    });

    // The 2-theta sliders. Outside simulation their handler slices the data;
    // with no data that handler returns early, so the simulation needs its own.
    ['tth-min-range', 'tth-min-slider', 'tth-max-range', 'tth-max-slider'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener('input', () => { if (simActive()) simRenderDebounced(); });
    });

    const sel = document.getElementById('algorithm-select');
    if (sel) sel.addEventListener('change', simApplyMode);

    simBuildRows();
    simApplyMode();
}

// ---------------------------------------------------------------------------
//  PUBLIC SURFACE
//
//  Small on purpose. powder5.html calls isActive() to decide whether the fit's
//  preview paths should stand aside, and render() when something it owns
//  changes (the cell, the space group, the wavelength, the profile). Everything
//  else in this file is reached through the DOM.
// ---------------------------------------------------------------------------
root.SIM = {
    VERSION: SIM_VERSION,
    init: simInit,
    isActive: simActive,
    render: simRender,
    renderDebounced: simRenderDebounced,
    applyMode: simApplyMode,
    refreshWyckoff: simRefreshWyckoff,
    /** The report button on the Plot tab routes here while simulating. */
    exportPdCif: simExportPdCif,
    resetView: simResetView,
    loadCifText: simLoadCifText,
    /** For the space-group change path: the operators depend on the setting. */
    invalidateSymmetry: () => { simSymCache.clear(); simInvalidateGeometry(); },
    /** Console diagnostic: SIM.stats() in devtools. */
    stats: () => ({ ...simCounters,
                    atoms: simAtoms.length,
                    reflections: simGeomList ? simGeomList.length : 0,
                    cachedGeometry: !!simGeomList }),
    /** Read-only view, for a report or a CIF export later. */
    atoms: () => simAtoms.map(a => ({ ...a }))
};

})(typeof window !== 'undefined' ? window : globalThis);
