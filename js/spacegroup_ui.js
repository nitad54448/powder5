        // ===================================================================
        // Space-Group modal: filter chips + search + explanation panel.
        // All settings (≈530 rows) are pre-materialised into a flat
        // array `ALL_SETTINGS` once; every re-filter walks that array.
        // With an AND'd filter chain this is cheap (< 2 ms on any machine).
        // ===================================================================
        const SGM = (() => {
            //   Build the pooled list once.  
            //  Listed high → low symmetry so the most-used systems appear first.
            const SYSTEMS = ['cubic','hexagonal','trigonal','tetragonal','orthorhombic','monoclinic','triclinic'];
            const CENTERINGS = ['P','A','B','C','I','F','R'];

 
            function buildAllSettings() {
                const out = [];
                for (let n = 1; n <= 230; n++) {
                    const entries = SG_ENGINE.listSettings(n);
                    entries.forEach(rec => {
                        // rec already has: number, symbol, standard_symbol,
                        // setting_description, hall, system, laue_class,
                        // point_group, centering, centrosymmetric,
                        // reflection_conditions, display_label (optional),
                        // alternate_halls (optional).
                        //
                        // The HKL generator only supports hexagonal-axes
                        // indexing for R-centred groups (the d-spacing
                        // formula and asymmetric-unit loop are written
                        // that way). Drop the rhomb-axes duplicates so
                        // the user can't pick a setting the engine
                        // wouldn't handle correctly.
                        const lbl = rec.display_label || '';
                        if (/rhomb\. axes/.test(lbl)) return;

                        // Strip the now-redundant "[hex. axes]" tag from
                        // the remaining label, so an R-centred entry
                        // just reads e.g. "R-3c" like every other group.
                        if (/hex\. axes/.test(lbl)) {
                            rec = { ...rec,
                                    display_label: lbl.replace(/\s*\[hex\. axes\]\s*/, '').trim() };
                        }
                        out.push(rec);
                    });
                }
                return out;
            }

            let ALL_SETTINGS = null;

            //   Modal DOM handles (resolved lazily).

            let nodes = null;
            function grabNodes() {
                if (nodes) return nodes;
                nodes = {
                    overlay: document.getElementById('sg-modal-overlay'),
                    systemRow: document.getElementById('sg-chip-row-system'),
                    centRow:   document.getElementById('sg-chip-row-centering'),
                    axesLabel: document.getElementById('sg-chip-row-axes-label'),
                    axesRow:   document.getElementById('sg-chip-row-axes'),
                    numInput:  document.getElementById('sg-search-number'),
                    txtInput:  document.getElementById('sg-search-text'),
                    listWrap:  document.getElementById('sg-list-wrap'),
                    explain:   document.getElementById('sg-explain'),
                    btnCancel: document.getElementById('sg-btn-cancel'),
                    btnReset:  document.getElementById('sg-btn-reset'),
                    btnApply:  document.getElementById('sg-btn-apply')
                };
                return nodes;
            }

            //   Filter state. Reset() brings everything back to this.  
            const state = {
                system: null,       // one of SYSTEMS or null
                centering: null,    // one of CENTERINGS or null
                axis: null,         // 'hex' | 'rhomb' | null
                numText: '',
                txtText: '',
                selected: null      // the chosen setting record, or null
            };

            function resetState() {
                state.system = null;
                state.centering = null;
                state.axis = null;
                state.numText = '';
                state.txtText = '';
                state.selected = null;
            }

            function mkChip(label, active, disabled, onClick) {
                const chip = document.createElement('button');
                chip.type = 'button';
                chip.className = 'sg-chip' + (active ? ' active' : '') + (disabled ? ' disabled' : '');
                chip.textContent = label;
                chip.disabled = !!disabled;
                if (!disabled) chip.addEventListener('click', onClick);
                return chip;
            }

            // Which centerings are possible, given the current filter set
            // (excluding `centering` itself — we derive its enabled options
            // from what the OTHER filters allow). Same idea for axes.
            function availableCenteringsGiven(exclSystem) {
                const out = new Set();
                ALL_SETTINGS.forEach(rec => {
                    if (exclSystem && rec.system !== exclSystem) return;
                    out.add(rec.centering);
                });
                return out;
            }
            function axisChoiceApplies() {
                // Show the axis row only when the selection could be an
                // R-centred trigonal group (that's the only case where we
                // surface both hex- and rhomb-axis indexings).
                // Criterion: any row with an axis tag survives the current
                // non-axis filters.
                return ALL_SETTINGS.some(rec => {
                    if (!rec.axis) return false;
                    if (state.system && rec.system !== state.system) return false;
                    if (state.centering && rec.centering !== state.centering) return false;
                    return true;
                });
            }

            function renderChips() {
                const n = grabNodes();

                // System chips. Every system is always enabled because
                // toggling one simply reshapes the rest of the UI.
                n.systemRow.innerHTML = '';
                SYSTEMS.forEach(sys => {
                    n.systemRow.appendChild(mkChip(
                        sys.charAt(0).toUpperCase() + sys.slice(1),
                        state.system === sys,
                        false,
                        () => {
                            state.system = state.system === sys ? null : sys;
                            // Clear centering/axis if they no longer fit.
                            const avail = availableCenteringsGiven(state.system);
                            if (state.centering && !avail.has(state.centering)) state.centering = null;
                            if (!axisChoiceApplies()) state.axis = null;
                            renderAll();
                        }
                    ));
                });

                // Centering chips. Enabled set depends on current system.
                n.centRow.innerHTML = '';
                const avail = availableCenteringsGiven(state.system);
                CENTERINGS.forEach(c => {
                    const enabled = avail.has(c);
                    n.centRow.appendChild(mkChip(
                        c,
                        state.centering === c,
                        !enabled,
                        () => {
                            state.centering = state.centering === c ? null : c;
                            if (!axisChoiceApplies()) state.axis = null;
                            renderAll();
                        }
                    ));
                });

                // Axis-choice chips (shown only when relevant).
                const axisShown = axisChoiceApplies();
                n.axesLabel.style.display = axisShown ? '' : 'none';
                n.axesRow.style.display   = axisShown ? '' : 'none';
                n.axesRow.innerHTML = '';
                if (axisShown) {
                    [['hex','Hex. axes'],['rhomb','Rhomb. axes']].forEach(([key, label]) => {
                        n.axesRow.appendChild(mkChip(
                            label,
                            state.axis === key,
                            false,
                            () => {
                                state.axis = state.axis === key ? null : key;
                                renderAll();
                            }
                        ));
                    });
                }
            }

            function applyFilters() {
                const numQ = (state.numText || '').trim();
                const txtQ = (state.txtText || '').trim().toLowerCase();
                return ALL_SETTINGS.filter(rec => {
                    if (state.system    && rec.system     !== state.system)    return false;
                    if (state.centering && rec.centering  !== state.centering) return false;
                    if (state.axis) {
                        if (rec.axis !== state.axis) return false;
                    }
                    if (numQ) {
                        // Exact integer match, or prefix match on the number.
                        if (/^\d+$/.test(numQ)) {
                            if (String(rec.number) !== numQ) return false;
                        } else {
                            return false;
                        }
                    }
                    if (txtQ) {
                        const hay = (
                            rec.symbol + ' ' +
                            rec.standard_symbol + ' ' +
                            (rec.hall || '') + ' ' +
                            (rec.display_label || '')
                        ).toLowerCase();
                        if (hay.indexOf(txtQ) === -1) return false;
                    }
                    return true;
                });
            }

            function renderList(list) {
                const n = grabNodes();
                if (list.length === 0) {
                    n.listWrap.innerHTML = '<div class="sg-list-empty">No space groups match these filters.</div>';
                    return;
                }
                // Make sure the selection (if any) is still in the filtered list;
                // if not, drop it.


if (state.selected) {
    const keep = list.some(r =>
        r.number === state.selected.number
        && r.hall === state.selected.hall
    );
    if (!keep) state.selected = null;
}


                const frag = document.createDocumentFragment();
                list.slice(0, 600).forEach(rec => {   // cap — 530 max anyway
                    const row = document.createElement('div');
                    row.className = 'sg-list-item';
                    
                    const isSel = state.selected
    && rec.number === state.selected.number
    && rec.hall === state.selected.hall;

                    if (isSel) row.classList.add('selected');
                    row.innerHTML =
                        `<span class="num">#${rec.number}</span>`
                        + `<span>${escapeHtml(rec.display_label || rec.symbol)}</span>`
                        + `<span class="sys">${rec.system}</span>`;
                    row.addEventListener('click', () => {
                        state.selected = rec;
                        renderAll();
                    });
                    frag.appendChild(row);
                });
                n.listWrap.innerHTML = '';
                n.listWrap.appendChild(frag);
            }

            function escapeHtml(s) {
                return String(s)
                    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
            }

            function renderExplain(list) {
                const n = grabNodes();
                // Precedence:
                //   1. Empty list  -> "No match" message in the panel.
                //   2. Has list, nothing selected -> hint.
                //   3. Has selection -> full details.
                if (list.length === 0) {
                    n.explain.innerHTML =
                        '<div class="sg-explain-warning">'
                        + '<strong>No space groups match these filters.</strong><br>'
                        + 'Adjust the selections or clear a search field. '
                        + 'Apply stays disabled until a row is selected.'
                        + '</div>';
                    return;
                }
                if (!state.selected) {
                    n.explain.innerHTML =
                        '<div class="sg-explain-empty">'
                        + `<strong>${list.length}</strong> space-group setting`
                        + (list.length === 1 ? '' : 's')
                        + ' match the current filters. Pick one from the list to see details.'
                        + '</div>';
                    return;
                }
                const rec = state.selected;
                const condHtml = renderConditions(rec.reflection_conditions);
                n.explain.innerHTML =
                    '<div class="sg-explain">'
                    + `<h3>${escapeHtml(rec.display_label || rec.symbol)}</h3>`
                    + '<dl class="sg-explain-grid">'
                    + row('Number',    '#' + rec.number)
                    + row('Symbol',    rec.symbol)
                    + row('Setting',   rec.setting_description || '(default)')
                    + row('System',    rec.system)
                    + row('Laue class', rec.laue_class)
                    + row('Point group', rec.point_group)
                    + row('Centering', rec.centering)
                    + row('Hall',      rec.hall || '—')
                    + row('Centrosymmetric', rec.centrosymmetric ? 'yes' : 'no')
                    + '</dl>'
                    + '<div style="margin-top:12px;">'
                    +   '<div style="font-size:11px;text-transform:uppercase;letter-spacing:0.04em;color:var(--system-secondary-label);margin-bottom:4px;">Reflection conditions</div>'
                    +   condHtml
                    + '</div>'
                    + '</div>';
                function row(k, v) {
                    return `<dt>${escapeHtml(k)}</dt><dd>${escapeHtml(v)}</dd>`;
                }
            }

            function renderConditions(conds) {
                const keys = conds ? Object.keys(conds) : [];
                if (keys.length === 0) {
                    return '<div style="font-size:12px;color:var(--system-secondary-label);">No reflection conditions (every hkl is allowed).</div>';
                }
                let html = '<table class="sg-explain-cond-table"><thead><tr><th>Family</th><th>Rule</th></tr></thead><tbody>';
                keys.forEach(fam => {
                    (conds[fam] || []).forEach(rule => {
                        html += `<tr><td>${escapeHtml(fam)}</td><td>${escapeHtml(rule)}</td></tr>`;
                    });
                });
                html += '</tbody></table>';
                return html;
            }

            function renderAll() {
                renderChips();
                const list = applyFilters();
                renderList(list);
                renderExplain(list);
                const n = grabNodes();
                n.btnApply.disabled = !state.selected;
            }

           //   Public API  
            function open(prefill) {
                if (!ALL_SETTINGS) ALL_SETTINGS = buildAllSettings();
                const n = grabNodes();
                // Either preload from the current SG or start blank.

if (prefill) {
    state.system = null;
    state.centering = null;
    state.axis = null;
    state.numText = '';
    state.txtText = '';



  state.selected = ALL_SETTINGS.find(r =>
    r.number === prefill.number
    && r.hall === prefill.hall
) || null;


} else {
    resetState();
}


                n.numInput.value = state.numText;
                n.txtInput.value = state.txtText;

                
                renderAll();
n.overlay.classList.add('open');
setTimeout(() => {
    n.txtInput.focus();
    // Scroll the list to the selected item
    if (state.selected) {
        const selectedRow = n.listWrap.querySelector('.sg-list-item.selected');
        if (selectedRow) {
            selectedRow.scrollIntoView({ block: 'center', behavior: 'instant' });
        }
    }
}, 40);


            }
            function close() {
                grabNodes().overlay.classList.remove('open');
            }

            //   Event wiring (idempotent).  
            function wire(onApply) {
                const n = grabNodes();
                n.numInput.addEventListener('input', (e) => {
                    state.numText = e.target.value;
                    renderAll();
                });
                n.txtInput.addEventListener('input', (e) => {
                    state.txtText = e.target.value;
                    renderAll();
                });
                n.btnCancel.addEventListener('click', close);
                n.btnReset.addEventListener('click', () => {
                    resetState();
                    n.numInput.value = '';
                    n.txtInput.value = '';
                    renderAll();
                });
                n.btnApply.addEventListener('click', () => {
                    if (!state.selected) return;
                    onApply(state.selected);
                    close();
                });
                // Close on Escape; backdrop click.
                n.overlay.addEventListener('click', (e) => {
                    if (e.target === n.overlay) close();
                });
                document.addEventListener('keydown', (e) => {
                    if (e.key === 'Escape' && n.overlay.classList.contains('open')) close();
                });
            }

            return { open, close, wire };
        })();



        // =====================================================================
        //  Space-group probability test
        //  Markvardsen, David, Johnson & Shankland (2001),
        //  Acta Cryst. A57, 47-54:
        //  "A Probabilistic Approach to Space-Group Determination from
        //   Powder Diffraction Data."
        //
        //  IMPORTANT LIMITATION
        //  --------------------
        //  A reflection can only be tested if it was actually a parameter of
        //  the Pawley fit. If the Pawley refinement was run in a group that
        //  already imposes reflection conditions, the reflections it forbids
        //  were never refined and carry no information. For a complete test,
        //  run the Pawley refinement in the condition-free group of the same
        //  Laue class (P-1, P2/m, Pmmm, P4/mmm, P6/mmm, Pm-3m, ...). The test
        //  detects this case and says so.
        // =====================================================================
        const SGTEST = (() => {

            //   ---------- small dense linear algebra (symmetric, PD) ----------
            function chol(A) {
                const n = A.length;
                const L = [];
                for (let i = 0; i < n; i++) L.push(new Float64Array(n));
                for (let i = 0; i < n; i++) {
                    for (let j = 0; j <= i; j++) {
                        let s = A[i][j];
                        for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
                        if (i === j) {
                            if (!(s > 0) || !isFinite(s)) return null;   // not PD
                            L[i][i] = Math.sqrt(s);
                        } else {
                            L[i][j] = s / L[j][j];
                        }
                    }
                }
                return L;
            }
            function cholLogDet(L) {
                let s = 0;
                for (let i = 0; i < L.length; i++) s += 2 * Math.log(L[i][i]);
                return s;
            }
            function cholSolve(L, b) {
                const n = L.length;
                const y = new Float64Array(n), x = new Float64Array(n);
                for (let i = 0; i < n; i++) {
                    let s = b[i];
                    for (let k = 0; k < i; k++) s -= L[i][k] * y[k];
                    y[i] = s / L[i][i];
                }
                for (let i = n - 1; i >= 0; i--) {
                    let s = y[i];
                    for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k];
                    x[i] = s / L[i][i];
                }
                return x;
            }
            function cholInverse(L) {
                const n = L.length;
                const inv = [];
                for (let i = 0; i < n; i++) inv.push(new Float64Array(n));
                const e = new Float64Array(n);
                for (let c = 0; c < n; c++) {
                    e.fill(0); e[c] = 1;
                    const col = cholSolve(L, e);
                    for (let r = 0; r < n; r++) inv[r][c] = col[r];
                }
                for (let i = 0; i < n; i++)
                    for (let j = i + 1; j < n; j++) {
                        const v = 0.5 * (inv[i][j] + inv[j][i]);
                        inv[i][j] = v; inv[j][i] = v;
                    }
                return inv;
            }

            //   ---------- space-group helpers ----------
            // A group is "condition-free" if it forbids nothing: those are the
            // groups a Pawley fit must be run in for the test to see every
            // reflection.
            const PROBE_HKL = (() => {
                const out = [];
                for (let h = -3; h <= 3; h++)
                    for (let k = -3; k <= 3; k++)
                        for (let l = -3; l <= 3; l++)
                            if (h || k || l) out.push([h, k, l]);
                return out;
            })();

            function isConditionFree(rec) {
                if (!rec) return false;
                try {
                    for (const [h, k, l] of PROBE_HKL)
                        if (!SG_ENGINE.isReflectionAllowed(h, k, l, rec)) return false;
                } catch { return false; }
                return true;
            }

            // Every setting sharing the Laue class of the refinement. The Laue
            // class is fixed by the diffraction symmetry the user has already
            // committed to (it sets the reflection orbits and the
            // multiplicities), so it is not itself under test here -- only the
            // absences are. SG_ENGINE.listSettings already drops the
            // rhombohedral-axes settings the HKL generator cannot index.
            // Monoclinic settings off unique axis b are no longer excluded:
            // the generator reads the unique axis from the setting and the
            // metric follows it, so they are scored on their own terms.
            function candidateSettings(laue, system) {
                const out = [];
                for (let num = 1; num <= 230; num++) {
                    let entries;
                    try { entries = SG_ENGINE.listSettings(num) || []; } catch { entries = []; }
                    for (const rec of entries) {
                        if (!rec || rec.laue_class !== laue) continue;
                        out.push(rec);
                    }
                }
                return out;
            }

            function label(rec) {
                // listSettings() returns the setting symbol, which is the most
                // informative label here: Pnma and Pbnm are the same group but
                // predict different absences.
                return String(rec.symbol || rec.standard_symbol || `#${rec.number}`)
                    .replace(/\s+/g, ' ').trim();
            }

            //   ---------- gather the Pawley result ----------
            function collect() {
                if (!fitResults)
                    return { error: 'No refinement results are available yet. Run a Pawley refinement first.' };
                if (fitResults.refinementMode !== 'pawley')
                    return { error: 'The last refinement was a Le Bail fit. Le Bail intensities are not least-squares parameters and have no normal matrix, so the test cannot be applied. Run a Pawley refinement.' };
                if (!fitResults.JtJ || !fitResults.parameterInfo)
                    return { error: 'The test needs the normal-equations matrix. The Levenberg-Marquardt refiner builds one, and a Parallel Tempering run gets one from its LM polish; neither was available here. Re-run the Pawley fit.' };
                if (!workerWorkingData || !workerWorkingData.intensity || !workerWorkingData.intensity.length)
                    return { error: 'The fitted data slice is no longer in memory. Re-run the Pawley refinement.' };

                const { hklList, JtJ, parameterInfo, ss_res, params } = fitResults;
                const P = parameterInfo.length;
                const N = workerWorkingData.intensity.length;
                const dof = N - P;
                if (!(dof > 0)) return { error: `Not enough data points (${N}) for the ${P} refined parameters; no error estimates exist.` };
                if (!(ss_res >= 0) || !isFinite(ss_res)) return { error: 'The refinement returned an invalid sum of squared residuals.' };
                if (!Array.isArray(JtJ) || JtJ.length !== P || !JtJ[0] || JtJ[0].length !== P)
                    return { error: 'The normal-equations matrix has unexpected dimensions.' };
                const s2 = ss_res / dof;

                const nameToParam = new Map();
                parameterInfo.forEach((p, i) => { if (p.isIntensity) nameToParam.set(p.name, i); });
                if (nameToParam.size === 0)
                    return { error: 'The refinement contains no refined intensity parameters.' };

                const tthMin = workerWorkingData.tth[0];
                const tthMax = workerWorkingData.tth[workerWorkingData.tth.length - 1];

                // Convert each refined peak HEIGHT into an integrated intensity
                // per symmetry-equivalent reflection. The profile area is taken
                // analytically so that it stays valid for reflections that
                // refined to (near) zero -- which are precisely the ones the
                // test cares about, and for which the "area / height" route
                // used by the report would collapse to 0.
                // The test compares |F|^2 against a Wilson distribution, so the
                // refined HEIGHT has to be turned into something proportional
                // to |F|^2: multiply by the profile area, then divide by m
                // AND by Lp. Dividing by m alone (the old behaviour) left the
                // whole Lorentz-polarisation fall-off in the "structure
                // factors", which at Cu Ka spans more than an order of
                // magnitude across a normal scan -- the low-angle reflections
                // then look systematically strong and the extinction test is
                // biased against groups whose conditions kill them.
                const sgtPol = polarizationFromParams(params);
                const scaleOf = (hkl) => {
                    let w = 1;
                    try {
                        w = getPseudoVoigtArea(hkl.tth, hkl, params);
                        const l1 = params.lambda, l2 = params.lambda2, r21 = params.ratio;
                        if (r21 > 1e-6 && l2 > 1e-6 && Math.abs(l1 - l2) > 1e-6) {
                            const st2 = Math.sin(hkl.tth * Math.PI / 360) * (l2 / l1);
                            if (Math.abs(st2) < 1) {
                                const tth2 = 2 * Math.asin(st2) * 180 / Math.PI;
                                w += r21 * getPseudoVoigtArea(tth2, hkl, params);
                            }
                        }
                    } catch { w = 1; }
                    if (!isFinite(w) || w <= 0) w = 1;
                    const lp = cfLorentzPolarization(hkl.tth, sgtPol);
                    return w / (Math.max(1, hkl.multiplicity || 1) * lp);
                };

                const refl = [];
                (hklList || []).forEach(hkl => {
                    if (!hkl || typeof hkl.tth !== 'number') return;
                    if (!(hkl.tth >= tthMin && hkl.tth <= tthMax)) return;
                    const pi = nameToParam.get(`I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`);
                    if (pi === undefined) return;
                    const f = scaleOf(hkl);
                    refl.push({
                        h: hkl.h_orig, k: hkl.k_orig, l: hkl.l_orig,
                        tth: hkl.tth, pi, f,
                        I: Math.max(0, hkl.intensity || 0) * f
                    });
                });

                if (refl.length < 4)
                    return { error: `Only ${refl.length} usable refined reflections -- far too few to say anything about the space group.` };

                // Low-angle reflections carry nearly all of the extinction
                // information and are the best determined, so if the list is
                // enormous the test keeps the low-angle end rather than
                // choking on an n^3 that buys nothing.
                const MAX_REFL = 220;
                let truncatedTo = 0;
                refl.sort((a, b) => a.tth - b.tth);
                if (refl.length > MAX_REFL) { truncatedTo = MAX_REFL; refl.length = MAX_REFL; }
                const n = refl.length;

                //   Marginal precision of the intensities.
                // Schur complement over every other refined parameter, so the
                // lattice, profile and background uncertainties are folded in
                // rather than treated as exactly known.
                const inI = new Uint8Array(P);
                const Iidx = refl.map(r => r.pi);
                Iidx.forEach(i => { inI[i] = 1; });
                const Oidx = [];
                for (let i = 0; i < P; i++) if (!inI[i]) Oidx.push(i);
                const no = Oidx.length;

                const Q = [];
                for (let i = 0; i < n; i++) {
                    const row = new Float64Array(n);
                    const src = JtJ[Iidx[i]];
                    for (let j = 0; j < n; j++) row[j] = src[Iidx[j]];
                    Q.push(row);
                }

                let marginalised = false;
                if (no > 0) {
                    const Aoo = [];
                    for (let i = 0; i < no; i++) {
                        const row = new Float64Array(no);
                        const src = JtJ[Oidx[i]];
                        for (let j = 0; j < no; j++) row[j] = src[Oidx[j]];
                        Aoo.push(row);
                    }
                    const Loo = chol(Aoo);
                    if (Loo) {
                        // X[i] = Aoo^-1 * (column i of Q_OI)
                        const X = [];
                        for (let i = 0; i < n; i++) {
                            const b = new Float64Array(no);
                            for (let o = 0; o < no; o++) b[o] = JtJ[Oidx[o]][Iidx[i]];
                            X.push(cholSolve(Loo, b));
                        }
                        for (let i = 0; i < n; i++) {
                            const rowI = JtJ[Iidx[i]];
                            for (let j = i; j < n; j++) {
                                let s = 0;
                                for (let o = 0; o < no; o++) s += rowI[Oidx[o]] * X[j][o];
                                Q[i][j] -= s;
                                if (j !== i) Q[j][i] = Q[i][j];
                            }
                        }
                        marginalised = true;
                    }
                    // If Aoo is singular the raw sub-block is used instead; the
                    // resulting errors are optimistic but the test still runs.
                }

                // Rescale from refined heights to the integrated intensities in
                // `refl`, and divide out the variance scale: I = D h means
                // Q_I = D^-1 Q_h D^-1, and the fit's variance estimate is s2.
                for (let i = 0; i < n; i++)
                    for (let j = 0; j < n; j++) {
                        const v = Q[i][j] / (s2 * refl[i].f * refl[j].f);
                        Q[i][j] = isFinite(v) ? v : 0;
                    }
                for (let i = 0; i < n; i++)
                    for (let j = i + 1; j < n; j++) {
                        const v = 0.5 * (Q[i][j] + Q[j][i]);
                        Q[i][j] = v; Q[j][i] = v;
                    }

                // Diagnostic sigmas. The proper marginal variance needs Q^-1,
                // which does not exist when two reflections coincide exactly in
                // 2-theta. Fall back to the conditional variance 1/Q_jj in that
                // case and say so, rather than refusing to run.
                const sigma2 = new Float64Array(n);
                let sigmaExact = false;
                const LQ = chol(Q);
                if (LQ) {
                    const Qinv = cholInverse(LQ);
                    let good = true;
                    for (let i = 0; i < n; i++) {
                        if (!(Qinv[i][i] > 0) || !isFinite(Qinv[i][i])) { good = false; break; }
                        sigma2[i] = Qinv[i][i];
                    }
                    sigmaExact = good;
                }
                if (!sigmaExact) {
                    for (let i = 0; i < n; i++) {
                        const q = Q[i][i];
                        sigma2[i] = (q > 0 && isFinite(q)) ? 1 / q : Infinity;
                    }
                }

                return { refl, Q, sigma2, sigmaExact, marginalised, n, s2, truncatedTo };
            }

            //   ---------- the actual test ----------
            function run() {
                const sgFit = sgUsedForFit || currentSG;
                if (!sgFit) return { error: 'No space group is set.' };

                const data = collect();
                if (data.error) return data;
                const { refl, Q, sigma2, n } = data;

                // c = Q * Ihat does not depend on the candidate, so it is
                // computed once.
                const Ivec = Float64Array.from(refl, r => r.I);
                const cvec = new Float64Array(n);
                for (let i = 0; i < n; i++) {
                    let s = 0;
                    for (let j = 0; j < n; j++) s += Q[i][j] * Ivec[j];
                    cvec[i] = s;
                }

                const sig = new Float64Array(n);
                for (let i = 0; i < n; i++) sig[i] = Math.sqrt(Math.max(1e-300, sigma2[i]));

                // tau^2: the Wilson-like intensity scale of the reflections a
                // candidate says are PRESENT, from a sliding window in 2-theta
                // over that candidate's free set only (refl is sorted by
                // 2-theta). Estimating it once over all reflections instead is
                // a trap: in a C-centred or rhombohedral cell most reflections
                // are genuinely absent, the window fills up with zeros, tau
                // collapses, and the Occam term stops rewarding the candidate
                // that correctly explains the most absences -- which is exactly
                // the candidate one is trying to find. The window mean of
                // Ihat^2 estimates <I^2> + <var>, so the measurement variance
                // is subtracted back off; floored at var(I) so a region where
                // nothing is measurable cannot masquerade as evidence.
                function tauForFreeSet(free) {
                    const nf = free.length;
                    const tau2 = new Float64Array(nf);
                    if (nf === 0) return tau2;
                    let globalMS = 0;
                    for (let q = 0; q < nf; q++) {
                        const I = refl[free[q]].I;
                        globalMS += I * I;
                    }
                    globalMS = Math.max(1e-12, globalMS / nf);
                    const win = Math.max(10, Math.min(nf, Math.round(nf / 5)));
                    for (let i = 0; i < nf; i++) {
                        const lo = Math.max(0, Math.min(nf - win, i - (win >> 1)));
                        const hi = Math.min(nf, lo + win);
                        let s = 0, v = 0, c = 0;
                        for (let q = lo; q < hi; q++) {
                            const I = refl[free[q]].I;
                            s += I * I;
                            const sv = sigma2[free[q]];
                            if (isFinite(sv)) v += sv;
                            c++;
                        }
                        let t = c ? (s - v) / c : 0;
                        const fl = isFinite(sigma2[free[i]]) ? sigma2[free[i]] : 0;
                        if (!(t > fl)) t = fl;
                        if (!(t > 0) || !isFinite(t)) t = globalMS;
                        tau2[i] = t;
                    }
                    return tau2;
                }

                // Evidence for one candidate, given the indices it declares
                // absent. B = Q_ff + diag(1/tau_f^2) is positive definite even
                // when Q is singular, so exact reflection coincidences are
                // harmless here.
                function evidence(ext) {
                    const isExt = new Uint8Array(n);
                    for (const a of ext) isExt[a] = 1;
                    const free = [];
                    for (let i = 0; i < n; i++) if (!isExt[i]) free.push(i);
                    const nf = free.length;
                    if (nf === 0) return 0;          // everything absent: no free parameters

                    const tau2 = tauForFreeSet(free);
                    const B = [];
                    for (let a = 0; a < nf; a++) {
                        const row = new Float64Array(nf);
                        const qa = Q[free[a]];
                        for (let b = 0; b < nf; b++) row[b] = qa[free[b]];
                        row[a] += 1 / tau2[a];
                        B.push(row);
                    }
                    const LB = chol(B);
                    if (!LB) return null;
                    const bvec = Float64Array.from(free, i => cvec[i]);
                    const y = cholSolve(LB, bvec);
                    let quad = 0;
                    for (let a = 0; a < nf; a++) quad += bvec[a] * y[a];
                    let lnTau = 0;
                    for (let a = 0; a < nf; a++) lnTau += Math.log(tau2[a]);
                    return -0.5 * (lnTau + cholLogDet(LB) - quad);
                }

                const cands = candidateSettings(sgFit.laue_class, sgFit.system);
                if (!cands.length) return { error: `No space-group settings found for Laue class ${sgFit.laue_class}.` };

                const openRec = cands.find(isConditionFree) || null;
                const pawleyIsOpen = isConditionFree(sgFit);

                // Reflections the Pawley model never contained, so never tested.
                let missingSet = null;
                if (!pawleyIsOpen && openRec) {
                    try {
                        const have = new Set(refl.map(r => `${r.h},${r.k},${r.l}`));
                        const loTth = refl[0].tth, hiTth = refl[n - 1].tth;
                        const full = generateAndCacheHklIndices(
                            openRec, hiTth, fitResults.params) || [];
                        // Give them 2-theta values so the comparison is limited
                        // to the range that was actually fitted.
                        updateHklPositions(full, fitResults.params, openRec.system);
                        missingSet = full
                            .filter(r => typeof r.tth === 'number'
                                      && r.tth >= loTth && r.tth <= hiTth
                                      && !have.has(`${r.h_orig},${r.k_orig},${r.l_orig}`))
                            .map(r => [r.h_orig, r.k_orig, r.l_orig]);
                    } catch { missingSet = null; }
                }

                const groups = new Map();   // extinction signature -> row
                for (const rec of cands) {
                    let sigStr = '';
                    const ext = [];
                    try {
                        for (let i = 0; i < n; i++) {
                            if (SG_ENGINE.isReflectionAllowed(refl[i].h, refl[i].k, refl[i].l, rec)) {
                                sigStr += '1';
                            } else {
                                sigStr += '0'; ext.push(i);
                            }
                        }
                    } catch { continue; }

                    let row = groups.get(sigStr);
                    if (!row) {
                        row = { ext, recs: [], untested: 0 };
                        if (missingSet) {
                            let u = 0;
                            for (const [h, k, l] of missingSet)
                                if (!SG_ENGINE.isReflectionAllowed(h, k, l, rec)) u++;
                            row.untested = u;
                        }
                        groups.set(sigStr, row);
                    }
                    row.recs.push(rec);
                }

                const rows = [];
                for (const row of groups.values()) {
                    const e = row.ext, m = e.length;
                    const lnL = evidence(e);
                    if (lnL === null || !isFinite(lnL)) continue;

                    // Plain diagnostics on the would-be-absent set, using only
                    // the diagonal errors. Not what drives the ranking, but it
                    // is the number a crystallographer looks at first.
                    let meanIsig = 0, maxIsig = 0;
                    for (const i of e) {
                        const z = refl[i].I / sig[i];
                        meanIsig += z;
                        if (Math.abs(z) > Math.abs(maxIsig)) maxIsig = z;
                    }
                    if (m) meanIsig /= m;

                    rows.push({ recs: row.recs, m, untested: row.untested, lnL, meanIsig, maxIsig });
                }
                if (!rows.length) return { error: 'No candidate space group could be evaluated.' };

                // Flat prior over distinct extinction signatures -- powder data
                // cannot distinguish settings that predict the same absences.
                let maxLn = -Infinity;
                rows.forEach(r => { if (r.lnL > maxLn) maxLn = r.lnL; });
                let norm = 0;
                rows.forEach(r => { r.w = Math.exp(r.lnL - maxLn); norm += r.w; });
                rows.forEach(r => { r.p = r.w / norm; });
                rows.sort((a, b) => b.p - a.p);

                return {
                    rows, n,
                    laue: sgFit.laue_class,
                    system: sgFit.system,
                    sgFit, openRec, pawleyIsOpen,
                    missingCount: missingSet ? missingSet.length : 0,
                    sigmaExact: data.sigmaExact,
                    marginalised: data.marginalised,
                    truncatedTo: data.truncatedTo
                };
            }

            //   ---------- presentation ----------
            let lastResult = null;
            let lastRecIndex = [];

            function esc(s) {
                return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
            }

            function renderError(msg) {
                return `<div class="sg-explain-warning">${esc(msg)}</div>`;
            }

            function renderResult(res) {
                lastRecIndex = [];
                const parts = [];

                parts.push(
                    `<div class="sgt-intro">Probability of each distinct set of systematic absences, given the `
                    + `${res.n} refined Pawley intensities and their full correlation structure. `
                    + `Method: Markvardsen, David, Johnson &amp; Shankland, <i>Acta Cryst.</i> (2001) <b>A57</b>, 47&ndash;54.</div>`
                );

                if (!res.pawleyIsOpen) {
                    parts.push(
                        `<div class="sg-explain-warning" style="margin:10px 0;">`
                        + `<b>Incomplete test.</b> The Pawley refinement was run in `
                        + `<b>${esc(label(res.sgFit))}</b>, which already forbids `
                        + `${res.missingCount ? res.missingCount + ' ' : ''}reflections in this range. `
                        + `Those reflections were never refined, so nothing is known about them and any `
                        + `candidate that relies on them is scored on partial evidence (marked below). `
                        + (res.openRec
                            ? `For a proper test, re-run the Pawley fit in <b>${esc(label(res.openRec))}</b> `
                            + `&mdash; same Laue class, no reflection conditions &mdash; then come back here.`
                            : ``)
                        + `</div>`
                    );
                }

                const notes = [];
                if (res.truncatedTo)
                    notes.push(`Restricted to the ${res.truncatedTo} lowest-angle reflections for tractability.`);
                if (!res.sigmaExact)
                    notes.push('Some reflections coincide exactly in 2\u03B8, so the I/\u03C3 columns are conditional (optimistic) estimates. The probabilities are unaffected.');
                if (!res.marginalised)
                    notes.push('The non-intensity parameters could not be marginalised out; errors are slightly optimistic.');
                if (notes.length)
                    parts.push(`<div class="sgt-note" style="margin-top:8px;">${notes.map(esc).join(' ')}</div>`);

                parts.push(`<table class="sgt-table"><thead><tr>`
                    + `<th style="width:64px;">P(SG&nbsp;|&nbsp;data)</th>`
                    + `<th>Space group(s)</th>`
                    + `<th style="width:52px;text-align:right;">absent</th>`
                    + `<th style="width:66px;text-align:right;">&lang;I/&sigma;&rang;</th>`
                    + `<th style="width:72px;text-align:right;">worst&nbsp;I/&sigma;</th>`
                    + `</tr></thead><tbody>`);

                const shown = res.rows.filter((r, i) => i < 12 || r.p > 0.001);
                shown.forEach((r, rank) => {
                    const pct = r.p >= 0.0005 ? (100 * r.p).toFixed(1) + '%' : '&lt;0.1%';
                    const bar = Math.max(1, Math.round(100 * r.p));
                    const chips = r.recs.slice(0, 14).map(rec => {
                        lastRecIndex.push(rec);
                        return `<button type="button" class="sgt-chip" data-rec="${lastRecIndex.length - 1}" `
                            + `title="Select this space group">${esc(label(rec))}</button>`;
                    }).join('');
                    const more = r.recs.length > 14 ? `<span class="sgt-more">+${r.recs.length - 14} more</span>` : '';
                    const flag = r.untested
                        ? `<span class="sgt-flag" title="${r.untested} of this group's conditions involve reflections that were not in the Pawley model">partial</span>`
                        : '';
                    parts.push(`<tr${rank === 0 ? ' class="sgt-top"' : ''}>`
                        + `<td><div class="sgt-bar"><span style="width:${bar}%"></span></div>`
                        + `<div class="sgt-pct">${pct}</div></td>`
                        + `<td>${chips}${more}${flag}</td>`
                        + `<td class="sgt-num">${r.m}</td>`
                        + `<td class="sgt-num">${r.m ? r.meanIsig.toFixed(2) : '&ndash;'}</td>`
                        + `<td class="sgt-num">${r.m ? r.maxIsig.toFixed(1) : '&ndash;'}</td>`
                        + `</tr>`);
                });
                parts.push(`</tbody></table>`);

                parts.push(`<div class="sgt-note" style="margin-top:12px;">`
                    + `Rows are distinct <i>extinction symbols</i>: powder data cannot separate space groups that `
                    + `predict identical absences, so those are listed together and share one probability. `
                    + `&lang;I/&sigma;&rang; and worst I/&sigma; use diagonal errors only and are diagnostics; the `
                    + `probability itself uses the full correlation structure, which is what handles overlapped `
                    + `reflections correctly. Click any symbol to adopt that setting.`
                    + `</div>`);

                return parts.join('');
            }

            function asText(res) {
                const L = [];
                L.push('Space-Group Probability Test');
                L.push('Markvardsen, David, Johnson & Shankland (2001), Acta Cryst. A57, 47-54');
                L.push('');
                L.push(`Laue class          : ${res.laue}`);
                L.push(`Pawley group used   : ${label(res.sgFit)}${res.pawleyIsOpen ? '' : '  (imposes conditions - test is incomplete)'}`);
                if (!res.pawleyIsOpen && res.openRec)
                    L.push(`Recommended re-run  : ${label(res.openRec)}  (condition-free, same Laue class)`);
                L.push(`Reflections tested  : ${res.n}`);
                L.push('');
                L.push('  P(SG|data)  absent  <I/sig>  worst   space group(s)');
                L.push('  ' + '-'.repeat(72));
                res.rows.slice(0, 25).forEach(r => {
                    L.push('  '
                        + (100 * r.p).toFixed(1).padStart(8) + '%'
                        + String(r.m).padStart(8)
                        + (r.m ? r.meanIsig.toFixed(2) : '-').padStart(9)
                        + (r.m ? r.maxIsig.toFixed(1) : '-').padStart(8)
                        + '   ' + r.recs.map(label).join(', ')
                        + (r.untested ? '   [partial]' : ''));
                });
                return L.join('\n');
            }

            //   ---------- modal plumbing ----------
            function nodes() {
                return {
                    overlay: document.getElementById('sgt-modal-overlay'),
                    body: document.getElementById('sgt-body'),
                    copy: document.getElementById('sgt-btn-copy'),
                    close: document.getElementById('sgt-btn-close')
                };
            }

            function open() {
                const n = nodes();
                if (!n.overlay) return;
                lastResult = null;
                n.body.innerHTML = `<div class="sgt-note">Evaluating candidate space groups&hellip;</div>`;
                n.copy.disabled = true;
                n.overlay.classList.add('open');

                // Let the modal paint before the (synchronous) linear algebra.
                setTimeout(() => {
                    let res;
                    try { res = run(); }
                    catch (err) {
                        console.error('Space-group test failed:', err);
                        res = { error: `The test failed: ${err && err.message ? err.message : err}` };
                    }
                    if (res.error) { n.body.innerHTML = renderError(res.error); return; }
                    lastResult = res;
                    n.body.innerHTML = renderResult(res);
                    n.copy.disabled = false;
                    n.body.querySelectorAll('.sgt-chip').forEach(btn => {
                        btn.addEventListener('click', () => {
                            const rec = lastRecIndex[parseInt(btn.dataset.rec, 10)];
                            if (!rec) return;
                            close();
                            handleSpaceGroupApplied(rec);
                            showToast(`Space group set to ${label(rec)}.`, 'success');
                        });
                    });
                }, 30);
            }

            function close() {
                const n = nodes();
                if (n.overlay) n.overlay.classList.remove('open');
            }

            function wire() {
                const n = nodes();
                if (!n.overlay) return;
                n.close.addEventListener('click', close);
                n.copy.addEventListener('click', () => {
                    if (!lastResult) return;
                    const txt = asText(lastResult);
                    if (navigator.clipboard && navigator.clipboard.writeText) {
                        navigator.clipboard.writeText(txt)
                            .then(() => showToast('Result copied to the clipboard.', 'success'))
                            .catch(() => showToast('Could not access the clipboard.', 'error'));
                    } else {
                        showToast('Clipboard not available in this browser.', 'error');
                    }
                });
                n.overlay.addEventListener('click', (e) => { if (e.target === n.overlay) close(); });
                document.addEventListener('keydown', (e) => {
                    if (e.key === 'Escape' && n.overlay.classList.contains('open')) close();
                });
            }

            return { open, close, wire };
        })();

