/*
 * plot_styles.js -- colour, line style, width and point markers for the series
 * of the main plot, edited from its legend.
 *
 *   click on a legend entry           show / hide the series (unchanged)
 *   right-click on a legend entry     open the style editor
 *   long-press on a legend entry      open the style editor (touch screens)
 *
 * Choices are kept in localStorage under 'powder5-plot-styles', keyed by the
 * dataset LABEL, and re-applied by a Chart.js plugin at the start of every
 * update. That is the one place they can be applied reliably: initializeChart()
 * destroys the chart and rebuilds it from its defaults whenever a data file is
 * loaded, and a style written once onto the old datasets would vanish with
 * them. It also means the PDF report, which copies this canvas, gets the same
 * styles with no extra work.
 *
 * Right-click in the PLOT AREA still resets the zoom (charting.js). Only a
 * right-click that lands on a legend entry is taken here, and it is stopped
 * in the capture phase so the reset-zoom handler never sees it.
 *
 * Load after lib/chart.umd.min.js; placed after js/charting.js in the page.
 */
(function () {
    'use strict';
    if (typeof Chart === 'undefined' || typeof Chart.register !== 'function') return;

    const CANVAS_ID = 'main-chart';
    const STORE_KEY = 'powder5-plot-styles';
    const LONG_PRESS_MS = 550;
    const HINT = 'Click to show or hide. Right-click to change the style.';

    // The only dataset properties this module ever writes.
    const PROPS = ['borderColor', 'backgroundColor', 'borderWidth', 'borderDash',
                   'borderCapStyle', 'showLine', 'pointStyle', 'pointRadius',
                   'pointHoverRadius', 'pointBackgroundColor', 'pointBorderColor',
                   'barThickness'];
    const ABSENT = Symbol('absent');

    const LINE_OPTS = [['solid', 'Solid'], ['dashed', 'Dashed'], ['dotted', 'Dotted'],
                       ['dashdot', 'Dash-dot'], ['none', 'No line']];
    const MARK_OPTS = [['none', 'No points'], ['circle', 'Circles'], ['rect', 'Squares'],
                       ['triangle', 'Triangles'], ['cross', 'Crosses']];
    const DASHES = new Set(LINE_OPTS.map(o => o[0]));
    const MARKERS = new Set(MARK_OPTS.map(o => o[0]));

    // ---------------------------------------------------------------------
    //  Storage
    // ---------------------------------------------------------------------
    let styles = {};
    try {
        const s = JSON.parse(localStorage.getItem(STORE_KEY) || '{}');
        if (s && typeof s === 'object' && !Array.isArray(s)) styles = s;
    } catch (e) { /* private mode or corrupt entry: start clean */ }

    let saveTimer = 0;
    function save() {
        clearTimeout(saveTimer);
        saveTimer = setTimeout(() => {
            try { localStorage.setItem(STORE_KEY, JSON.stringify(styles)); } catch (e) { /* quota */ }
        }, 200);
    }

    // ---------------------------------------------------------------------
    //  Applying a style
    // ---------------------------------------------------------------------
    //  Each dataset's own values are photographed the first time the plugin
    //  sees it -- inside `new Chart(...)`, before anything has been changed --
    //  so a style can always be taken back off exactly. Datasets that were
    //  never styled are never touched.
    const originals = new WeakMap();
    const styled = new WeakSet();

    const copy = v => (Array.isArray(v) ? v.slice() : v);

    function snapshot(ds) {
        if (originals.has(ds)) return;
        const s = {};
        for (const p of PROPS) s[p] = Object.prototype.hasOwnProperty.call(ds, p) ? copy(ds[p]) : ABSENT;
        originals.set(ds, s);
    }

    function restore(ds) {
        const s = originals.get(ds);
        if (!s) return;
        for (const p of PROPS) {
            if (s[p] === ABSENT) delete ds[p];
            else ds[p] = copy(s[p]);
        }
    }

    const isHex = c => typeof c === 'string' && /^#[0-9a-f]{6}$/i.test(c);
    function num(v, lo, hi) {
        const x = Number(v);
        return (v === null || v === undefined || !Number.isFinite(x)) ? null : Math.min(hi, Math.max(lo, x));
    }

    function dashPattern(kind, w) {
        switch (kind) {
            case 'dashed':  return [4 + 3 * w, 3 + 2 * w];
            case 'dotted':  return [Math.max(1, w), 2 + 2 * w];
            case 'dashdot': return [4 + 3 * w, 2 + 2 * w, Math.max(1, w), 2 + 2 * w];
            default:        return [];
        }
    }

    function applyStyle(ds) {
        const o = styles[ds.label];
        if (!o || typeof o !== 'object') {
            if (styled.has(ds)) { restore(ds); styled.delete(ds); }
            return;
        }
        restore(ds);                       // start from the defaults every time: idempotent
        styled.add(ds);
        const bar = ds.type === 'bar';

        if (isHex(o.color)) {
            const a = num(o.alpha, 0, 1);
            const c = rgba(o.color, a === null ? 1 : a);
            if (bar) ds.backgroundColor = c;
            else ds.borderColor = ds.pointBackgroundColor = ds.pointBorderColor = c;
        }
        const w = num(o.width, 0.1, 12);
        if (bar) {
            if (w !== null) ds.barThickness = w;
            return;
        }
        if (w !== null) ds.borderWidth = w;
        if (DASHES.has(o.dash)) {
            if (o.dash === 'none') {
                ds.showLine = false;
            } else {
                ds.showLine = true;
                ds.borderDash = dashPattern(o.dash, typeof ds.borderWidth === 'number' ? ds.borderWidth : 1);
                ds.borderCapStyle = 'butt';
            }
        }
        if (MARKERS.has(o.marker)) {
            if (o.marker === 'none') {
                ds.pointRadius = 0;
            } else {
                const r = num(o.size, 0.5, 12);
                ds.pointStyle = o.marker;
                ds.pointRadius = r === null ? 2 : r;
                ds.pointHoverRadius = ds.pointRadius + 1.5;
            }
        }
    }

    /** The dataset's current look, in the editor's terms. */
    function readStyle(ds) {
        const bar = ds.type === 'bar';
        const col = parseColour(bar ? ds.backgroundColor : ds.borderColor) || { hex: '#888888', alpha: 1 };
        const s = { bar, color: col.hex, alpha: col.alpha };
        if (bar) {
            s.width = typeof ds.barThickness === 'number' ? ds.barThickness : 1;
            return s;
        }
        s.width = typeof ds.borderWidth === 'number' ? ds.borderWidth : 3;      // Chart.js default
        s.dash = dashKind(ds);
        const r = typeof ds.pointRadius === 'number' ? ds.pointRadius : 3;       // Chart.js default
        s.marker = r > 0 ? (MARKERS.has(ds.pointStyle) ? ds.pointStyle : 'circle') : 'none';
        s.size = r > 0 ? r : 2;
        return s;
    }

    function dashKind(ds) {
        if (ds.showLine === false) return 'none';
        const d = ds.borderDash;
        if (!Array.isArray(d) || d.length < 2 || d.every(v => !v)) return 'solid';
        if (d.length >= 4) return 'dashdot';
        const w = typeof ds.borderWidth === 'number' ? ds.borderWidth : 1;
        return d[0] <= Math.max(2, 1.5 * w) ? 'dotted' : 'dashed';
    }

    // Colours arrive as rgba() strings, hex, or names. The 2D context is the
    // one parser that accepts everything a Chart.js colour can be.
    const cx = document.createElement('canvas').getContext('2d');
    function parseColour(c) {
        if (typeof c !== 'string' || !cx) return null;
        cx.fillStyle = '#000000';
        cx.fillStyle = c;
        const n = String(cx.fillStyle);
        if (n[0] === '#') return { hex: n.toLowerCase(), alpha: 1 };
        const m = n.match(/rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*(?:,\s*([\d.]+)\s*)?\)/);
        if (!m) return null;
        const hex = '#' + [m[1], m[2], m[3]].map(v => (+v).toString(16).padStart(2, '0')).join('');
        return { hex, alpha: m[4] === undefined ? 1 : +m[4] };
    }
    function rgba(hex, a) {
        const n = parseInt(hex.slice(1), 16);
        return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${Math.round(a * 1000) / 1000})`;
    }

    // ---------------------------------------------------------------------
    //  The plugin
    // ---------------------------------------------------------------------
    //  Identified by canvas id, not by `chart === mainChart`: the first update
    //  runs inside `new Chart(...)`, before that assignment has happened.
    const isMain = chart => !!(chart && chart.canvas && chart.canvas.id === CANVAS_ID);

    Chart.register({
        id: 'powder5PlotStyles',
        afterInit(chart) {
            if (isMain(chart)) wireCanvas(chart.canvas);
        },
        beforeUpdate(chart) {
            if (!isMain(chart)) return;
            // WRAPPED. A throw here would unwind out of Chart.js's update and
            // leave the chart half-built; a style is never worth that.
            try {
                for (const ds of chart.data.datasets) {
                    if (!ds) continue;
                    snapshot(ds);
                    applyStyle(ds);
                }
            } catch (e) {
                console.error('plot_styles:', e);
            }
        }
    });

    let raf = 0;
    function redraw() {
        if (raf) return;
        raf = requestAnimationFrame(() => {
            raf = 0;
            const chart = currentChart();
            if (!chart) return;
            try {
                if (typeof window.chartUpdateSafe === 'function') window.chartUpdateSafe(chart);
                else chart.update();
            } catch (e) {
                console.error('plot_styles:', e);
            }
        });
    }

    function currentChart() {
        const canvas = document.getElementById(CANVAS_ID);
        return canvas ? (Chart.getChart(canvas) || null) : null;
    }

    // ---------------------------------------------------------------------
    //  Finding the legend entry under the pointer
    // ---------------------------------------------------------------------
    function legendHit(chart, clientX, clientY) {
        const lg = chart && chart.legend;
        if (!lg || !Array.isArray(lg.legendHitBoxes) || !Array.isArray(lg.legendItems)) return null;
        const r = chart.canvas.getBoundingClientRect();
        const x = clientX - r.left, y = clientY - r.top;
        for (let i = 0; i < lg.legendHitBoxes.length; i++) {
            const b = lg.legendHitBoxes[i], it = lg.legendItems[i];
            if (!b || !it || typeof it.datasetIndex !== 'number') continue;
            if (x >= b.left && x <= b.left + b.width && y >= b.top && y <= b.top + b.height) {
                const ds = chart.data.datasets[it.datasetIndex];
                if (!ds) return null;
                return {
                    key: ds.label,
                    text: it.text,
                    anchor: { left: r.left + b.left, top: r.top + b.top, bottom: r.top + b.top + b.height }
                };
            }
        }
        return null;
    }

    function mainChartOf(ev) {
        const t = ev.target;
        return (t && t.id === CANVAS_ID) ? (Chart.getChart(t) || null) : null;
    }

    // Hover: pointer cursor and a tooltip, so the right-click can be found.
    function wireCanvas(canvas) {
        if (canvas.dataset.plotStyles === '1') return;
        canvas.dataset.plotStyles = '1';
        let over = false;
        const leave = () => {
            if (!over) return;
            over = false;
            if (canvas.getAttribute('title') === HINT) canvas.removeAttribute('title');
            if (canvas.style.cursor === 'pointer') canvas.style.cursor = '';
        };
        canvas.addEventListener('mousemove', e => {
            const hit = legendHit(Chart.getChart(canvas), e.clientX, e.clientY);
            if (hit && !over) {
                over = true;
                canvas.setAttribute('title', HINT);
                if (!canvas.style.cursor) canvas.style.cursor = 'pointer';
            } else if (!hit) {
                leave();
            }
        });
        canvas.addEventListener('mouseleave', leave);
    }

    // Right-click on an entry. Capture phase on document, so it runs before
    // the canvas's own contextmenu listener (reset zoom) and can stop it.
    document.addEventListener('contextmenu', e => {
        const chart = mainChartOf(e);
        if (!chart) return;
        const hit = legendHit(chart, e.clientX, e.clientY);
        if (!hit) return;
        e.preventDefault();
        e.stopPropagation();
        cancelLongPress();
        openEditor(hit);
    }, true);

    // Long-press for touch and pen. The click that may follow the release is
    // swallowed, or it would also toggle the series off.
    let lp = null, swallowClickUntil = 0;
    function cancelLongPress() {
        if (lp) { clearTimeout(lp.timer); lp = null; }
    }
    document.addEventListener('pointerdown', e => {
        if (e.pointerType === 'mouse') return;
        const chart = mainChartOf(e);
        const hit = chart && legendHit(chart, e.clientX, e.clientY);
        if (!hit) return;
        cancelLongPress();
        lp = {
            id: e.pointerId, x: e.clientX, y: e.clientY,
            timer: setTimeout(() => {
                lp = null;
                swallowClickUntil = performance.now() + 800;
                openEditor(hit);
            }, LONG_PRESS_MS)
        };
    }, true);
    document.addEventListener('pointermove', e => {
        if (lp && e.pointerId === lp.id && Math.hypot(e.clientX - lp.x, e.clientY - lp.y) > 10) cancelLongPress();
    }, true);
    for (const type of ['pointerup', 'pointercancel']) {
        document.addEventListener(type, e => { if (lp && e.pointerId === lp.id) cancelLongPress(); }, true);
    }
    document.addEventListener('click', e => {
        if (swallowClickUntil && performance.now() < swallowClickUntil && mainChartOf(e)) {
            e.stopPropagation();
            e.preventDefault();
        }
        swallowClickUntil = 0;
    }, true);

    // ---------------------------------------------------------------------
    //  The editor
    // ---------------------------------------------------------------------
    injectCss();
    let ed = null;          // { key, text, state, el }

    const svgLine = dash => (dash === 'none')
        ? '<span class="ps-word">None</span>'
        : `<svg viewBox="0 0 28 10" width="28" height="10" aria-hidden="true"><line x1="2" y1="5" x2="26" y2="5"
             stroke="currentColor" stroke-width="2" stroke-dasharray="${({ solid: '', dashed: '6 3', dotted: '2 3', dashdot: '6 3 2 3' })[dash]}"/></svg>`;
    const svgMark = m => ({
        none:     '<span class="ps-word">None</span>',
        circle:   '<svg viewBox="0 0 14 14" width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="3.6" fill="currentColor"/></svg>',
        rect:     '<svg viewBox="0 0 14 14" width="14" height="14" aria-hidden="true"><rect x="3.5" y="3.5" width="7" height="7" fill="currentColor"/></svg>',
        triangle: '<svg viewBox="0 0 14 14" width="14" height="14" aria-hidden="true"><path d="M7 2.8 11.4 11H2.6z" fill="currentColor"/></svg>',
        cross:    '<svg viewBox="0 0 14 14" width="14" height="14" aria-hidden="true"><path d="M7 2.5v9M2.5 7h9" stroke="currentColor" stroke-width="1.7"/></svg>'
    })[m];

    const seg = (name, label, opts, icon) =>
        `<div class="ps-seg" role="radiogroup" aria-label="${label}" data-seg="${name}">` +
        opts.map(([v, t]) => `<button type="button" role="radio" data-v="${v}" title="${t}" aria-label="${t}">${icon(v)}</button>`).join('') +
        '</div>';

    function openEditor(hit) {
        const chart = currentChart();
        const ds = chart && chart.data.datasets.find(d => d && d.label === hit.key);
        if (!ds) return;
        closeEditor();

        const el = document.createElement('div');
        el.className = 'ps-pop';
        el.setAttribute('role', 'dialog');
        el.setAttribute('aria-label', `Style of ${hit.text}`);
        el.tabIndex = -1;
        el.innerHTML = `
          <div class="ps-head">
            <svg class="ps-sample" viewBox="0 0 44 16" width="44" height="16" aria-hidden="true"></svg>
            <span class="ps-name"></span>
            <button type="button" class="ps-x" aria-label="Close">&times;</button>
          </div>
          <div class="ps-row">
            <span class="ps-lab">Colour</span>
            <div class="ps-ctl ps-colour">
              <input type="color" class="ps-color" aria-label="Colour">
              <input type="range" class="ps-alpha" min="10" max="100" step="5" aria-label="Opacity">
            </div>
            <output class="ps-val ps-alpha-v"></output>
          </div>
          <div class="ps-row ps-only-line">
            <span class="ps-lab">Line</span>
            <div class="ps-ctl ps-span">${seg('dash', 'Line style', LINE_OPTS, svgLine)}</div>
          </div>
          <div class="ps-row">
            <span class="ps-lab">Width</span>
            <div class="ps-ctl"><input type="range" class="ps-width" aria-label="Width"></div>
            <output class="ps-val ps-width-v"></output>
          </div>
          <div class="ps-row ps-only-line">
            <span class="ps-lab">Points</span>
            <div class="ps-ctl ps-span">${seg('marker', 'Point markers', MARK_OPTS, svgMark)}</div>
          </div>
          <div class="ps-row ps-only-line ps-size-row">
            <span class="ps-lab">Size</span>
            <div class="ps-ctl"><input type="range" class="ps-size" min="0.5" max="6" step="0.5" aria-label="Point size"></div>
            <output class="ps-val ps-size-v"></output>
          </div>
          <div class="ps-foot">
            <label class="ps-show"><input type="checkbox" class="ps-vis"> Show on plot</label>
            <button type="button" class="ps-reset" title="Go back to the default style for this series">Reset</button>
          </div>`;
        document.body.appendChild(el);
        el.querySelector('.ps-name').textContent = hit.text;

        ed = { key: hit.key, text: hit.text, state: readStyle(ds), el };
        el.classList.toggle('ps-bar', ed.state.bar);
        const w = el.querySelector('.ps-width');
        if (ed.state.bar) { w.min = '1'; w.max = '6'; w.step = '1'; }
        else              { w.min = '0.5'; w.max = '5'; w.step = '0.5'; }

        bindEditor(el);
        syncEditor();

        el.style.visibility = 'hidden';
        place(el, hit.anchor);
        el.style.visibility = '';
        el.focus({ preventScroll: true });
    }

    function bindEditor(el) {
        const q = s => el.querySelector(s);
        q('.ps-x').addEventListener('click', closeEditor);
        q('.ps-color').addEventListener('input', e => change({ color: e.target.value, alpha: ed.state.alpha }));
        q('.ps-alpha').addEventListener('input', e => change({ color: ed.state.color, alpha: +e.target.value / 100 }));
        q('.ps-width').addEventListener('input', e => change({ width: +e.target.value }));
        q('.ps-size').addEventListener('input', e => {
            change({ marker: ed.state.marker === 'none' ? 'circle' : ed.state.marker, size: +e.target.value });
        });
        el.querySelectorAll('.ps-seg').forEach(g => {
            g.addEventListener('click', e => {
                const b = e.target.closest('button[data-v]');
                if (!b) return;
                if (g.dataset.seg === 'dash') change({ dash: b.dataset.v });
                else change({ marker: b.dataset.v, size: ed.state.size });
            });
            // Arrow keys move along the group, as a radio group should.
            g.addEventListener('keydown', e => {
                if (!['ArrowLeft', 'ArrowRight'].includes(e.key)) return;
                const bs = [...g.querySelectorAll('button')];
                const i = bs.indexOf(document.activeElement);
                if (i < 0) return;
                e.preventDefault();
                const n = bs[(i + (e.key === 'ArrowRight' ? 1 : bs.length - 1)) % bs.length];
                n.focus();
                n.click();
            });
        });
        q('.ps-vis').addEventListener('change', e => {
            const chart = currentChart();
            const i = chart ? chart.data.datasets.findIndex(d => d && d.label === ed.key) : -1;
            if (i < 0) return;
            chart.setDatasetVisibility(i, e.target.checked);
            redraw();
        });
        q('.ps-reset').addEventListener('click', () => {
            delete styles[ed.key];
            save();
            const chart = currentChart();
            if (chart) {
                try {
                    if (typeof window.chartUpdateSafe === 'function') window.chartUpdateSafe(chart);
                    else chart.update();
                } catch (err) { console.error('plot_styles:', err); }
                const ds = chart.data.datasets.find(d => d && d.label === ed.key);
                if (ds) ed.state = readStyle(ds);
            }
            syncEditor();
        });
    }

    function change(patch) {
        if (!ed) return;
        styles[ed.key] = Object.assign({}, styles[ed.key], patch);
        Object.assign(ed.state, patch);
        save();
        syncEditor();
        redraw();
    }

    function syncEditor() {
        if (!ed) return;
        const el = ed.el, s = ed.state, q = sel => el.querySelector(sel);
        q('.ps-color').value = s.color;
        q('.ps-alpha').value = String(Math.round(s.alpha * 100));
        q('.ps-alpha-v').textContent = `${Math.round(s.alpha * 100)}%`;
        q('.ps-width').value = String(s.width);
        q('.ps-width-v').textContent = `${+Number(s.width).toFixed(1)} px`;
        if (!s.bar) {
            setSeg(q('[data-seg="dash"]'), s.dash);
            setSeg(q('[data-seg="marker"]'), s.marker);
            const off = s.marker === 'none';
            q('.ps-size').value = String(s.size);
            q('.ps-size').disabled = off;
            q('.ps-size-row').classList.toggle('ps-off', off);
            q('.ps-size-v').textContent = off ? '' : `${+Number(s.size).toFixed(1)} px`;
        }
        const chart = currentChart();
        const i = chart ? chart.data.datasets.findIndex(d => d && d.label === ed.key) : -1;
        q('.ps-vis').checked = i >= 0 ? chart.isDatasetVisible(i) : true;
        q('.ps-sample').innerHTML = sampleSvg(s);
    }

    function setSeg(g, v) {
        g.querySelectorAll('button').forEach(b => {
            const on = b.dataset.v === v;
            b.setAttribute('aria-checked', on ? 'true' : 'false');
            b.tabIndex = on ? 0 : -1;
        });
    }

    // The header swatch: the series as it will be drawn.
    function sampleSvg(s) {
        const c = rgba(s.color, s.alpha);
        if (s.bar) return `<rect x="${22 - s.width / 2}" y="1" width="${s.width}" height="14" fill="${c}"/>`;
        let out = '';
        if (s.dash !== 'none') {
            const d = dashPattern(s.dash, s.width).join(' ');
            out += `<line x1="1" y1="8" x2="43" y2="8" stroke="${c}" stroke-width="${s.width}" stroke-dasharray="${d}"/>`;
        }
        if (s.marker !== 'none') {
            const r = Math.min(6, s.size);
            const shape = {
                circle:   `<circle cx="22" cy="8" r="${r}" fill="${c}"/>`,
                rect:     `<rect x="${22 - r}" y="${8 - r}" width="${2 * r}" height="${2 * r}" fill="${c}"/>`,
                triangle: `<path d="M22 ${8 - r}L${22 + r} ${8 + r}H${22 - r}z" fill="${c}"/>`,
                cross:    `<path d="M22 ${8 - r}v${2 * r}M${22 - r} 8h${2 * r}" stroke="${c}" stroke-width="1.5"/>`
            }[s.marker] || '';
            out += shape;
        }
        return out;
    }

    function place(el, a) {
        const m = 8, w = el.offsetWidth, h = el.offsetHeight;
        let left = a.left, top = a.bottom + 6;
        if (left + w > window.innerWidth - m) left = window.innerWidth - m - w;
        if (left < m) left = m;
        if (top + h > window.innerHeight - m) top = Math.max(m, a.top - 6 - h);
        el.style.left = `${Math.round(left)}px`;
        el.style.top = `${Math.round(top)}px`;
    }

    function closeEditor() {
        if (!ed) return;
        ed.el.remove();
        ed = null;
    }

    document.addEventListener('pointerdown', e => {
        if (ed && !ed.el.contains(e.target)) closeEditor();
    }, true);
    document.addEventListener('keydown', e => {
        if (ed && e.key === 'Escape') closeEditor();
    });
    window.addEventListener('resize', closeEditor);

    // ---------------------------------------------------------------------
    //  Look. Values are var()s from style.css, as in the page's own style
    //  blocks, so the editor follows the light/dark switch by itself.
    // ---------------------------------------------------------------------
    function injectCss() {
        if (document.getElementById('plot-styles-css')) return;
        const st = document.createElement('style');
        st.id = 'plot-styles-css';
        st.textContent = `
.ps-pop {
    position: fixed; z-index: 10000; width: 284px; box-sizing: border-box;
    padding: 10px 12px 10px;
    background: var(--bg-panel); color: var(--fg);
    border: 1px solid var(--rule); border-radius: var(--r-2);
    box-shadow: 0 14px 34px -10px rgba(0, 0, 0, 0.45);
    font-size: 12.5px; line-height: 1.3;
}
.ps-pop:focus { outline: none; }
.ps-head { display: flex; align-items: center; gap: 9px; margin-bottom: 8px; }
.ps-sample { flex: none; }
.ps-name { flex: 1; font-size: 13px; font-weight: 600; }
.ps-x {
    width: 22px; height: 22px; padding: 0; border: 0; border-radius: var(--r-1);
    background: transparent; color: var(--fg-3); font-size: 17px; line-height: 1; cursor: pointer;
}
.ps-x:hover { background: var(--bg-hover); color: var(--fg); }
.ps-row {
    display: grid; grid-template-columns: 50px 1fr 42px;
    align-items: center; gap: 8px; min-height: 30px;
}
.ps-lab { color: var(--fg-3); }
.ps-ctl { display: flex; align-items: center; gap: 8px; min-width: 0; }
.ps-span { grid-column: 2 / 4; }
.ps-val {
    color: var(--fg-3); font-size: 11.5px; text-align: right;
    font-variant-numeric: tabular-nums; white-space: nowrap;
}
.ps-pop input[type=range] { flex: 1; width: 100%; min-width: 0; margin: 0; accent-color: var(--accent); }
.ps-pop input[type=color] {
    flex: none; width: 30px; height: 22px; padding: 1px;
    border: 1px solid var(--rule); border-radius: var(--r-1);
    background: var(--bg-inset); cursor: pointer;
}
.ps-seg {
    display: flex; flex: 1; gap: 2px; padding: 2px;
    background: var(--bg-inset); border-radius: var(--r-1);
}
.ps-seg button {
    flex: 1; height: 22px; padding: 0; border: 0; border-radius: var(--r-1);
    display: grid; place-items: center;
    background: transparent; color: var(--fg-3); cursor: pointer;
}
.ps-seg button:hover { color: var(--fg); }
.ps-seg button[aria-checked="true"] {
    background: var(--accent-tint); color: var(--accent);
}
.ps-word { font-size: 10.5px; }
.ps-pop button:focus-visible, .ps-pop input:focus-visible {
    outline: 2px solid var(--accent); outline-offset: 1px;
}
.ps-off { opacity: 0.45; }
.ps-bar .ps-only-line { display: none; }
.ps-foot {
    display: flex; align-items: center; justify-content: space-between;
    margin-top: 8px; padding-top: 8px; border-top: 1px solid var(--rule-soft);
}
.ps-show { display: flex; align-items: center; gap: 6px; cursor: pointer; }
.ps-show input { margin: 0; accent-color: var(--accent); }
.ps-reset {
    padding: 3px 10px; border: 1px solid var(--rule); border-radius: var(--r-1);
    background: transparent; color: var(--fg); font: inherit; cursor: pointer;
}
.ps-reset:hover { background: var(--bg-hover); }
`;
        document.head.appendChild(st);
    }

    // Console access, e.g. PlotStyles.resetAll() to go back to every default.
    window.PlotStyles = {
        get: () => JSON.parse(JSON.stringify(styles)),
        reset(label) { delete styles[label]; save(); redraw(); },
        resetAll() { styles = {}; save(); closeEditor(); redraw(); }
    };
})();
