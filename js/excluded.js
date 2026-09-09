// excluded.js
// version 1, september 2026
//
// User-defined EXCLUDED 2-theta REGIONS: intervals the refinement must not
// see. A sample holder line, a detector artefact, a known impurity peak, the
// beam stop shadow at very low angle -- anything present in the measurement
// that the model cannot account for and that would otherwise be charged to
// the structure.
//
// ---------------------------------------------------------------------------
//  HOW EXCLUSION IS APPLIED, AND WHY THAT WAY
// ---------------------------------------------------------------------------
//  The points inside an excluded region are REMOVED from workingData before
//  anything else touches it, rather than being kept and given zero weight.
//
//  Zero weighting is the obvious alternative and it is not enough here. It
//  would correct chi-square, Rwp and the normal equations, because those are
//  all weighted sums -- but Rp is an UNWEIGHTED sum of |y_obs - y_calc| over
//  the array, and the degrees of freedom are N - P with N the array length.
//  A zero-weighted point still contributes its whole residual to Rp and still
//  counts as an observation in the dof. The fit would be right and the numbers
//  reported for it would be wrong, which is the worse of the two failures:
//  nothing on screen would say so.
//
//  Removing the points instead makes every downstream sum correct with no
//  change at all in refinement_worker.js. The arrays stay ascending in
//  2-theta, so the binary searches over the axis keep working; nothing in the
//  fit path assumes a uniform step. The one thing that DID assume it,
//  setMinProfileFwhmFromAxis, now takes the median step rather than the mean
//  so a gap cannot inflate the resolvable-FWHM floor -- see constants.js.
//
//  Because the arrays are compacted, working index no longer equals
//  full-axis index minus startIndex. updateWorkingData publishes an
//  `indexMap` alongside them for the one consumer that needs to map back
//  (the background trace on the chart).
//
// ---------------------------------------------------------------------------
//  WHAT IS NOT EXCLUDED
// ---------------------------------------------------------------------------
//  The excluded points are still DRAWN. Hiding them would leave the user
//  looking at a pattern that is not the one they loaded, with no way to see
//  what the exclusion actually caught; the region is shaded on the chart
//  instead, so what was removed stays visible as removed.
//
//  Reflections whose position falls in an excluded region are NOT dropped from
//  the hkl list. A reflection is a property of the cell and the space group,
//  not of the data, and its neighbours' tails still reach into the kept range.
//  It simply has no data of its own to be determined by, which the degrees of
//  freedom already account for.
// ---------------------------------------------------------------------------

const EXCLUDED = (() => {
    'use strict';

    /**
     * @typedef {{min: number, max: number}} ExcludedRegion
     *   Inclusive interval in degrees 2-theta, always min < max.
     */

    /** @type {ExcludedRegion[]} Sorted by min, never overlapping. */
    let regions = [];

    /** Master switch. Off keeps the list intact but applies none of it, so a
     *  user can check what a region is worth without retyping it. */
    let enabled = true;

    /** Called after any change, so the page can re-slice and redraw. */
    let onChange = null;

    // -----------------------------------------------------------------------
    //  Model
    // -----------------------------------------------------------------------

    /**
     * Sorts and MERGES overlapping or touching intervals.
     *
     * Merging matters for more than tidiness: `count` below reports how many
     * points a region removes, and two overlapping regions would otherwise
     * report the shared points twice and read as excluding more of the pattern
     * than they do.
     *
     * @param {ExcludedRegion[]} list
     * @returns {ExcludedRegion[]} A new normalised list.
     */
    function normalise(list) {
        const clean = list
            .filter(r => r && Number.isFinite(r.min) && Number.isFinite(r.max) && r.max > r.min)
            .map(r => ({ min: r.min, max: r.max }))
            .sort((a, b) => a.min - b.min);

        const out = [];
        for (const r of clean) {
            const last = out[out.length - 1];
            if (last && r.min <= last.max) last.max = Math.max(last.max, r.max);
            else out.push(r);
        }
        return out;
    }

    /** @returns {ExcludedRegion[]} A copy; the caller cannot mutate the model. */
    function list() {
        return regions.map(r => ({ min: r.min, max: r.max }));
    }

    /** @returns {boolean} Whether exclusion is being applied at all. */
    function isEnabled() { return enabled && regions.length > 0; }

    /**
     * @param {number} tth Degrees 2-theta.
     * @returns {boolean} True if this angle lies in an active excluded region.
     */
    function contains(tth) {
        if (!enabled || !Number.isFinite(tth)) return false;
        // Linear: the list is a handful of entries, and a binary search here
        // would be a second thing to keep sorted correctly for no measurable
        // gain.
        for (let i = 0; i < regions.length; i++) {
            if (tth >= regions[i].min && tth <= regions[i].max) return true;
        }
        return false;
    }

    /**
     * The kept indices of an axis, in order.
     *
     * Returns null when nothing would be removed. That is the signal for the
     * caller to skip compaction entirely and keep its original typed arrays,
     * which is the common case and the one worth not paying for.
     *
     * @param {ArrayLike<number>} tthAxis Ascending 2-theta axis.
     * @returns {Int32Array|null} Indices to keep, or null for "keep everything".
     */
    function keepIndices(tthAxis) {
        if (!enabled || regions.length === 0 || !tthAxis || !tthAxis.length) return null;

        const n = tthAxis.length;
        const keep = new Int32Array(n);
        let k = 0;
        for (let i = 0; i < n; i++) {
            if (!contains(tthAxis[i])) keep[k++] = i;
        }
        if (k === n) return null;      // nothing landed in a region
        return keep.subarray(0, k);
    }

    /**
     * How many points of an axis the current regions remove.
     * @param {ArrayLike<number>} tthAxis
     * @returns {number}
     */
    function count(tthAxis) {
        if (!tthAxis || !tthAxis.length) return 0;
        const keep = keepIndices(tthAxis);
        return keep ? tthAxis.length - keep.length : 0;
    }

    /**
     * Adds a region. Silently normalises a reversed pair, because a user who
     * typed the bounds the other way round meant the interval between them.
     *
     * @param {number} min
     * @param {number} max
     * @returns {boolean} False if the pair does not describe an interval.
     */
    function add(min, max) {
        const a = Number(min), b = Number(max);
        if (!Number.isFinite(a) || !Number.isFinite(b) || a === b) return false;
        regions = normalise(regions.concat([{ min: Math.min(a, b), max: Math.max(a, b) }]));
        changed();
        return true;
    }

    /** @param {number} i Index into the list as rendered. */
    function remove(i) {
        if (i < 0 || i >= regions.length) return;
        regions.splice(i, 1);
        changed();
    }

    function clear() {
        if (!regions.length) return;
        regions = [];
        changed();
    }

    /** @param {boolean} on */
    function setEnabled(on) {
        const next = !!on;
        if (next === enabled) return;
        enabled = next;
        changed();
    }

    /**
     * Replaces the whole list at once, for a loader or a preset.
     * @param {ExcludedRegion[]} list
     */
    function setAll(newList) {
        regions = normalise(Array.isArray(newList) ? newList : []);
        changed();
    }

    /** @param {Function} fn Called after every model change. */
    function setOnChange(fn) { onChange = (typeof fn === 'function') ? fn : null; }

    /**
     * A frozen record of what was excluded, to travel WITH a set of results.
     *
     * WHY A SNAPSHOT AND NOT A LIVE READ. Reports are regenerated long after
     * the run that produced them -- from the history selector, from a PDF
     * button pressed an hour later -- and the region list is editable the
     * whole time. A report that read the CURRENT list would caption an old
     * refinement with regions that were never applied to it, and it would do
     * so in the one place a reader has no way to check: the section that
     * exists precisely to say which data the numbers came from. Wrong
     * provenance is worse than absent provenance, because it is believed.
     *
     * `enabled` is recorded separately from the list because the master switch
     * can be off with regions still defined, and "three regions, none applied"
     * and "no regions" are different facts about a run.
     *
     * @param {ArrayLike<number>} [axis] Full 2-theta axis, to record how many
     *        points the regions actually removed.
     * @returns {{regions: ExcludedRegion[], enabled: boolean,
     *            removed: number, axisPoints: number}}
     */
    function snapshot(axis) {
        const src = axis || ((typeof fullExperimentalData !== 'undefined' &&
                              fullExperimentalData && fullExperimentalData.tth) || null);
        return {
            regions: list(),
            enabled: enabled,
            removed: (src && enabled) ? count(src) : 0,
            axisPoints: src ? src.length : 0
        };
    }

    function changed() {
        render();
        if (onChange) {
            try { onChange(); } catch (e) { console.warn('excluded onChange failed:', e); }
        }
    }

    // -----------------------------------------------------------------------
    //  UI
    // -----------------------------------------------------------------------

    // Guarded, not a bare document.getElementById. render() is reached from
    // every model change, so without this the model could not be exercised at
    // all outside a document -- which contradicts the note above about this
    // file being readable from a realm that has no DOM.
    function el(id) {
        return (typeof document !== 'undefined') ? document.getElementById(id) : null;
    }

    /**
     * Redraws the region list and the summary line.
     *
     * The summary reports the count of REMOVED POINTS, not the width of the
     * regions in degrees. Degrees say how much of the axis is gone; points say
     * how many observations the fit lost, and that is the number that moves
     * the degrees of freedom and every R factor derived from them.
     */
    function render() {
        const listEl = el('excluded-list');
        if (!listEl) return;   // the tab is not in the document

        if (!regions.length) {
            listEl.innerHTML =
                '<p class="control-label" style="text-align: center; padding: 4px 8px 8px 0;">' +
                'No excluded regions.</p>';
        } else {
            // Editable, like the spline anchors: a region read off the chart by
            // eye almost always wants nudging afterwards, and retyping it as a
            // delete-and-re-add loses the reference point you were adjusting
            // against.
            listEl.innerHTML = regions.map((r, i) => `
                <div class="excluded-row" data-idx="${i}">
                    <label class="excluded-cap">From</label>
                    <input type="number" class="control-input excluded-min-input"
                           data-idx="${i}" value="${r.min.toFixed(3)}" step="0.01" title="Start of the region, degrees 2-theta">
                    <label class="excluded-cap">To</label>
                    <input type="number" class="control-input excluded-max-input"
                           data-idx="${i}" value="${r.max.toFixed(3)}" step="0.01" title="End of the region, degrees 2-theta">
                    <button type="button" class="excluded-del" data-del="${i}"
                            title="Remove this region">&times;</button>
                </div>`).join('');
        }

        const chk = el('excluded-enable');
        if (chk) chk.checked = enabled;

        const sum = el('excluded-summary');
        if (sum) {
            if (!regions.length) {
                sum.textContent = 'Nothing excluded; the fit uses the whole 2\u03b8 range.';
            } else if (!enabled) {
                sum.textContent = `${regions.length} region(s) defined but switched off.`;
            } else {
                // Against the full pattern, not the current slice: the slice
                // moves with the 2-theta sliders and the number would change
                // for a reason that has nothing to do with the regions.
                const axis = (typeof fullExperimentalData !== 'undefined' &&
                              fullExperimentalData && fullExperimentalData.tth) || null;
                const pts = axis ? count(axis) : 0;
                sum.textContent = axis && axis.length
                    ? `${regions.length} region(s), removing ${pts} of ${axis.length} points from the fit.`
                    : `${regions.length} region(s). Load a pattern to see how many points they remove.`;
            }
        }

        if (typeof mainChart !== 'undefined' && mainChart) chartUpdateSafe(mainChart);
    }

    /**
     * Wires the controls. Safe to call more than once: the listeners are
     * attached to elements that exist for the life of the page, and a second
     * call would double them, so it latches.
     */
    let wired = false;
    function init() {
        if (wired) return;
        const listEl = el('excluded-list');
        if (!listEl) return;
        wired = true;

        const addBtn = el('excluded-add-btn');
        const minEl = el('excluded-min');
        const maxEl = el('excluded-max');

        const doAdd = () => {
            const a = parseFloat(minEl && minEl.value);
            const b = parseFloat(maxEl && maxEl.value);
            if (!add(a, b)) {
                if (typeof showToast === 'function') {
                    showToast('Give a From and a To angle, and make them different.', 'error');
                }
                return;
            }
            // Cleared so the next region is typed from scratch rather than
            // edited out of the last one, which is how a duplicate gets added
            // by accident.
            if (minEl) minEl.value = '';
            if (maxEl) maxEl.value = '';
            if (minEl) minEl.focus();
        };

        if (addBtn) addBtn.addEventListener('click', doAdd);
        [minEl, maxEl].forEach(input => {
            if (!input) return;
            input.addEventListener('keydown', ev => {
                if (ev.key === 'Enter') { ev.preventDefault(); doAdd(); }
            });
        });

        const clearBtn = el('excluded-clear-btn');
        if (clearBtn) clearBtn.addEventListener('click', () => clear());

        const chk = el('excluded-enable');
        if (chk) chk.addEventListener('change', () => setEnabled(chk.checked));

        // Delegated, because the rows are rebuilt on every render.
        listEl.addEventListener('click', ev => {
            const btn = ev.target.closest('[data-del]');
            if (!btn) return;
            remove(parseInt(btn.dataset.del, 10));
        });

        // Edits commit on `change`, not on `input`: a half-typed "3" on its way
        // to "37.2" is a valid number, and committing it would re-slice the
        // pattern and re-sort the list under the cursor between keystrokes.
        listEl.addEventListener('change', ev => {
            const input = ev.target;
            const isMin = input.classList.contains('excluded-min-input');
            const isMax = input.classList.contains('excluded-max-input');
            if (!isMin && !isMax) return;

            const i = parseInt(input.dataset.idx, 10);
            if (!(i >= 0 && i < regions.length)) return;

            const v = parseFloat(input.value);
            if (!Number.isFinite(v)) { render(); return; }   // put the old value back

            // Edited through add(), so a bound dragged past its partner or
            // into a neighbouring region goes through the same normalise and
            // merge as any other region rather than leaving the list in a
            // state the model says cannot happen.
            const others = regions.filter((_, k) => k !== i);
            const edited = isMin ? { min: v, max: regions[i].max }
                                 : { min: regions[i].min, max: v };
            if (!(Math.abs(edited.max - edited.min) > 0)) { render(); return; }
            regions = normalise(others.concat([edited]));
            changed();
        });

        render();
    }

    // -----------------------------------------------------------------------
    //  Ctrl+Shift-drag on the chart
    // -----------------------------------------------------------------------
    //  THE THREE GESTURES ARE NOW DISJOINT AT THE MODIFIER, with no ordering
    //  between them and no shared state:
    //
    //      left drag                        rectangle zoom   (charting.js)
    //      Shift / Alt + left drag          pan              (charting.js)
    //      middle drag                      pan              (charting.js)
    //      Ctrl + click, no other modifier  background point (powder5.html)
    //      Ctrl + Shift + left drag         excluded region  (here)
    //
    //  Each handler tests for its own combination and bails on every other
    //  one, so no gesture is ever half-owned by two of them.
    //
    //  THIS FILE USED TO OWN PLAIN Ctrl-DRAG, and that was a mistake worth
    //  recording. Ctrl-click was already the background-point gesture, so a
    //  plain Ctrl-CLICK ran this module's mousedown, mousemove and mouseup --
    //  each calling chart.render() -- interleaved with the click handler that
    //  mutates the chart's data and calls update(). Two code paths driving one
    //  chart through a single click is not a thing to make safe with a
    //  suppression flag; the fix is for the gestures not to overlap. With
    //  Shift added, nothing here runs on a Ctrl-click at all, and the flag
    //  that used to tell the click handler to stand down is gone with it.
    // -----------------------------------------------------------------------

    /** Live drag state, in DATA units. Null when not dragging. */
    let drag = null;

    function chartOf() {
        return (typeof mainChart !== 'undefined' && mainChart) ? mainChart : null;
    }

    /** @returns {number|null} 2-theta under the pointer, or null if outside the plot. */
    function tthAtEvent(chart, ev) {
        const xs = chart.scales.x, ys = chart.scales.y;
        if (!xs || !ys) return null;
        const rect = chart.canvas.getBoundingClientRect();
        const px = ev.clientX - rect.left;
        const py = ev.clientY - rect.top;
        if (py < ys.top - 4 || py > ys.bottom + 4) return null;
        // Clamped rather than rejected: a drag that runs off the side of the
        // plot means "to the edge", which is what the user is showing you.
        const clamped = Math.min(Math.max(px, xs.left), xs.right);
        const v = xs.getValueForPixel(clamped);
        return Number.isFinite(v) ? v : null;
    }

    function wireChartDrag(chart) {
        // Takes the chart as an ARGUMENT rather than reading mainChart.
        // afterInit fires from inside `new Chart(...)`, so the assignment
        // `mainChart = new Chart(...)` has not happened yet and mainChart is
        // still undefined at this point -- reading it here would fail the guard
        // below every time and the drag would never be wired at all. The
        // handlers themselves resolve mainChart lazily, by which time it is set.
        if (!chart || !chart.canvas || chart.canvas.dataset.excludedDrag === '1') return;
        chart.canvas.dataset.excludedDrag = '1';
        const canvas = chart.canvas;

        canvas.addEventListener('mousedown', ev => {
            // Ctrl (or Cmd) AND Shift, and nothing else. Alt is the pan
            // modifier and is left alone.
            if (ev.button !== 0) return;
            if (!(ev.ctrlKey || ev.metaKey) || !ev.shiftKey || ev.altKey) return;
            const c = chartOf();
            if (!c) return;
            const t = tthAtEvent(c, ev);
            if (t === null) return;
            // Stops the browser turning the drag into a text selection, which
            // makes the whole page flash blue as the pointer moves.
            ev.preventDefault();
            drag = { from: t, to: t };
            c.render();
        });

        // On window, not the canvas: a drag that leaves the plot still has to
        // track and still has to finish, and a mouseup outside a canvas that
        // only listens to itself leaves the band stuck on screen.
        window.addEventListener('mousemove', ev => {
            if (!drag) return;
            const c = chartOf();
            if (!c) { drag = null; return; }
            const t = tthAtEvent(c, ev);
            if (t !== null) { drag.to = t; c.render(); }
        });

        window.addEventListener('mouseup', ev => {
            if (!drag) return;
            const c = chartOf();
            const d = drag;
            drag = null;
            if (!c) return;
            const lo = Math.min(d.from, d.to), hi = Math.max(d.from, d.to);

            // A shift-CLICK is a zero-width drag and means nothing here. The
            // threshold is in data units and scaled to the visible span, so it
            // is the same gesture whether the chart shows 5 degrees or 150.
            const xs = c.scales.x;
            const span = Math.abs((xs.max ?? 180) - (xs.min ?? 0)) || 180;
            if (hi - lo < span * 0.002) { c.render(); return; }

            add(lo, hi);
            if (typeof showToast === 'function') {
                showToast(`Excluded ${lo.toFixed(2)}\u2013${hi.toFixed(2)}\u00b0.`, 'info');
            }
        });

        // Escape abandons a drag in progress without adding anything.
        window.addEventListener('keydown', ev => {
            if (ev.key !== 'Escape' || !drag) return;
            drag = null;
            const c = chartOf();
            if (c) c.render();
        });
    }

    // -----------------------------------------------------------------------
    //  Chart shading
    // -----------------------------------------------------------------------
    //  Drawn UNDER the datasets, so an excluded peak is still legible through
    //  the shading. A region the user cannot see is a region they will forget
    //  is there, and then wonder why a peak is not being fitted.
    // -----------------------------------------------------------------------
    let drawFailed = false;
    const excludedRegionsPlugin = {
        id: 'excludedRegions',
        beforeDatasetsDraw(chart) {
            // WRAPPED. A throw from a plugin hook unwinds out of Chart.js's own
            // draw loop and leaves the chart half-updated; every later hit test
            // then fails inside Chart.js, once per mouse move, which reads as
            // an avalanche of errors from a file that did nothing wrong.
            // Reported once so a real fault here is still visible.
            try { drawRegions(chart); }
            catch (e) {
                if (!drawFailed) { drawFailed = true; console.warn('excludedRegions draw failed:', e); }
            }
        }
    };

    function drawRegions(chart) {
            if (!drag && (!enabled || !regions.length)) return;
            const xs = chart.scales.x, ys = chart.scales.y;
            if (!xs || !ys) return;

            const ctx = chart.ctx;
            ctx.save();
            ctx.beginPath();
            ctx.rect(xs.left, ys.top, xs.right - xs.left, ys.bottom - ys.top);
            ctx.clip();

            const band = (lo, hi, fill, stroke) => {
                const x0 = xs.getPixelForValue(lo);
                const x1 = xs.getPixelForValue(hi);
                if (!Number.isFinite(x0) || !Number.isFinite(x1)) return;
                const left = Math.min(x0, x1), width = Math.abs(x1 - x0);
                if (left > xs.right || left + width < xs.left) return;

                ctx.fillStyle = fill;
                ctx.fillRect(left, ys.top, Math.max(width, 1), ys.bottom - ys.top);

                ctx.strokeStyle = stroke;
                ctx.lineWidth = 1;
                ctx.setLineDash([4, 3]);
                ctx.beginPath();
                ctx.moveTo(left, ys.top);         ctx.lineTo(left, ys.bottom);
                ctx.moveTo(left + width, ys.top); ctx.lineTo(left + width, ys.bottom);
                ctx.stroke();
            };

            if (enabled) {
                for (const r of regions) {
                    band(r.min, r.max, 'rgba(239, 68, 68, 0.10)', 'rgba(239, 68, 68, 0.45)');
                }
            }
            // The drag in progress, drawn stronger than a committed region so
            // it is obvious which one is still under the pointer.
            if (drag) {
                band(Math.min(drag.from, drag.to), Math.max(drag.from, drag.to),
                     'rgba(239, 68, 68, 0.22)', 'rgba(239, 68, 68, 0.85)');
            }
            ctx.restore();
    }

    // The chart is destroyed and rebuilt on every file load, so this fires
    // again each time. The canvas element itself survives, so the dataset flag
    // inside keeps the listeners from stacking up.
    excludedRegionsPlugin.afterInit = function (chart) { wireChartDrag(chart); };

    // Registered here rather than in charting.js so that everything about
    // excluded regions is in one file. Guarded because this module is also
    // readable from a worker realm, where there is no Chart.
    if (typeof Chart !== 'undefined' && Chart.register) {
        Chart.register(excludedRegionsPlugin);
    }

    if (typeof document !== 'undefined') {
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', init);
        } else {
            init();
        }
    }

    return {
        list, setAll, add, remove, clear,
        isEnabled, setEnabled, contains, keepIndices, count,
        setOnChange, render, init, snapshot
    };
})();

if (typeof window !== 'undefined') window.EXCLUDED = EXCLUDED;
