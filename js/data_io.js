// data_io.js
// Writing patterns out, and reading pdCIF back in.
//
// ---------------------------------------------------------------------------
//  WHICH FORMATS ARE WRITTEN, AND WHY NOT THE OTHERS
// ---------------------------------------------------------------------------
//  parsers.js READS nine families. Most of them cannot sensibly be written:
//
//    .brml, .rasx   proprietary ZIP+XML containers. A synthetic one is of no
//                   use to anybody -- nobody feeds a fabricated Bruker file
//                   back to Bruker software -- and producing one means
//                   inventing instrument metadata that was never measured.
//    .rd, .sd       Philips binary.
//    .xrdml         a large required schema; almost every mandatory field
//                   would have to be fabricated.
//    .ras           writable in principle, but the header is the point of the
//                   format and we do not have it.
//
//  What is written is the set that round-trips honestly: the two-column
//  interchange formats, the two legacy text formats this program can read back
//  itself, and pdCIF -- which is the only one of the group actually DESIGNED
//  to carry a powder pattern together with its calculated profile, its
//  wavelength, its cell and its reflections. Everything else throws that away
//  and keeps two columns.
// ---------------------------------------------------------------------------

/**
 * @typedef {object} ExportBundle
 * @property {number[]|Float64Array} tth        Observed 2-theta.
 * @property {number[]|Float64Array} obs        Observed intensity.
 * @property {number[]|Float64Array} [calc]     Calculated total, if fitted.
 * @property {number[]|Float64Array} [bkg]      Background, if available.
 * @property {object} [params]                  Profile/cell parameters.
 * @property {object} [spaceGroup]              { symbol, number }.
 * @property {object[]} [hklList]               Reflections with intensities.
 * @property {object} [stats]                   Figures of merit.
 * @property {object} [esds]                    Parameter esds, keyed as the
 *                                              report keys them ("I_(h,k,l)").
 *                                              Without it the reflection esds
 *                                              in a Pawley pdCIF are "?".
 * @property {number} [scaleFactor]             Intensity scale of the fit.
 * @property {object} [polarisationUi]           { mode, monoTth, fraction } as
 *                                              set in the control panel.
 * @property {string} [sourceName]              Original data file name.
 */

/** Formats offered in the save dialog. `needsCalc` entries are hidden until a
 *  fit exists, because they would otherwise write a column of zeros and call
 *  it a calculated pattern. */
const DATA_EXPORT_FORMATS = [
    { id: 'xy',      label: 'Two-column (2\u03B8, intensity)',      ext: 'xy',  needsCalc: false },
    { id: 'xycalc',  label: 'Multi-column (obs, calc, bkg, diff)',  ext: 'xy',  needsCalc: true  },
    { id: 'uxd',     label: 'Bruker UXD',                            ext: 'uxd', needsCalc: false },
    { id: 'xra',     label: 'GSAS (BANK / XRA)',                     ext: 'xra', needsCalc: false },
    { id: 'pdcif',   label: 'pdCIF (profile + reflections)',         ext: 'cif', needsCalc: true  }
];

/**
 * The formats that can be written from the data currently in hand.
 * @param {boolean} haveCalc
 * @returns {object[]}
 */
function availableExportFormats(haveCalc) {
    return DATA_EXPORT_FORMATS.filter(f => haveCalc || !f.needsCalc);
}

/** Fixed-decimal helper that never emits "NaN" or "undefined" into a data file. */
function ioNum(v, dp) {
    return Number.isFinite(v) ? v.toFixed(dp) : (0).toFixed(dp);
}

/**
 * Uniform-step detection. UXD and GSAS both describe the axis by a start and a
 * step rather than listing it, so a non-uniform axis has to be refused rather
 * than silently written as if it were uniform -- which would misplace every
 * point after the first irregularity.
 *
 * @param {number[]|Float64Array} tth
 * @returns {{uniform: boolean, start: number, step: number, n: number}}
 */
function ioAxisStep(tth) {
    const n = tth.length;
    if (n < 2) return { uniform: false, start: n ? tth[0] : 0, step: 0, n };
    const step = (tth[n - 1] - tth[0]) / (n - 1);
    let uniform = step > 0;
    // A tenth of a step is generous; real instruments write far tighter axes,
    // and anything looser is a merged or trimmed pattern rather than a scan.
    const tol = Math.abs(step) * 0.1;
    for (let i = 1; i < n && uniform; i++) {
        if (Math.abs((tth[i] - tth[i - 1]) - step) > tol) uniform = false;
    }
    return { uniform, start: tth[0], step, n };
}

// ---------------------------------------------------------------------------
//  WRITERS
// ---------------------------------------------------------------------------

/** Two columns, space separated. The universal interchange. */
function writeXY(b) {
    const out = [`# ${b.sourceName || 'pattern'} exported by powder5`,
                 '# 2theta  intensity'];
    for (let i = 0; i < b.tth.length; i++) {
        out.push(`${ioNum(b.tth[i], 5)}  ${ioNum(b.obs[i], 4)}`);
    }
    return out.join('\n') + '\n';
}

/**
 * Observed, calculated, background and difference.
 *
 * The difference is written as obs - calc, the same sign convention the
 * difference curve uses on screen, so a plot made from this file looks like
 * the plot it came from.
 */
function writeXYCalc(b) {
    const out = [`# ${b.sourceName || 'pattern'} exported by powder5`,
                 '# 2theta  y_obs  y_calc  y_bkg  y_obs-y_calc'];
    for (let i = 0; i < b.tth.length; i++) {
        const o = b.obs[i], c = b.calc ? b.calc[i] : NaN, g = b.bkg ? b.bkg[i] : NaN;
        out.push([ioNum(b.tth[i], 5), ioNum(o, 4), ioNum(c, 4), ioNum(g, 4),
                  ioNum(Number.isFinite(c) ? o - c : NaN, 4)].join('  '));
    }
    return out.join('\n') + '\n';
}

/**
 * Bruker UXD. Counts only, with the axis described by _START and _STEPSIZE --
 * which is why a non-uniform axis is refused rather than written.
 */
function writeUXD(b) {
    const ax = ioAxisStep(b.tth);
    if (!ax.uniform) {
        throw new Error('UXD describes the 2-theta axis by a start and a fixed step, ' +
                        'and this pattern is not evenly spaced. Use the two-column ' +
                        'format instead, which lists every 2-theta explicitly.');
    }
    const wl = (b.params && b.params.lambda) || 1.5406;
    const out = [
        '; UXD written by powder5',
        `; source: ${b.sourceName || 'unknown'}`,
        '_FILEVERSION=1',
        '_SAMPLE=' + (b.sourceName || 'pattern'),
        '_WL1=' + wl.toFixed(6),
        '_DRIVE=THETA_2THETA',
        '_STEPTIME=1.0',
        `_STEPSIZE=${ax.step.toFixed(6)}`,
        `_START=${ax.start.toFixed(6)}`,
        `_STEPCOUNT=${ax.n}`,
        '_COUNTS'
    ];
    // Ten values per line, which is what the format's own writers produce.
    for (let i = 0; i < b.obs.length; i += 10) {
        out.push(Array.from(b.obs.slice(i, i + 10), v => ioNum(v, 2)).join(' '));
    }
    return out.join('\n') + '\n';
}

/**
 * GSAS constant-step powder data, ESD flavour written under the .xra name.
 *
 * BANK ... RALF/CONST records are fixed-width by specification: ten fields of
 * eight characters per line. Writing them as free-form text produces a file
 * that looks right and that GSAS cannot read.
 */
function writeXRA(b) {
    const ax = ioAxisStep(b.tth);
    if (!ax.uniform) {
        throw new Error('The GSAS constant-step format describes the axis by a start ' +
                        'and a step, and this pattern is not evenly spaced. Use the ' +
                        'two-column format instead.');
    }
    const n = ax.n;
    const perLine = 10;
    const nRec = Math.ceil(n / perLine);
    // GSAS counts 2-theta in centidegrees.
    const start100 = Math.round(ax.start * 100);
    const step100 = Math.round(ax.step * 100);
    const out = [
        (b.sourceName || 'powder5 export').slice(0, 80),
        `BANK 1 ${n} ${nRec} CONST ${start100} ${step100} 0 0 STD`
    ];
    for (let i = 0; i < n; i += perLine) {
        let line = '';
        for (let k = i; k < Math.min(i + perLine, n); k++) {
            line += String(Math.max(0, Math.round(b.obs[k]))).padStart(8);
        }
        out.push(line);
    }
    return out.join('\n') + '\n';
}

/**
 * pdCIF.
 *
 * The only format here that can carry the whole result: measured profile,
 * calculated profile, background, wavelength, cell, space group and the
 * reflection list with its standard uncertainties. Everything else in this
 * module keeps two columns and discards the rest.
 *
 * Standard uncertainties are written in the CIF parenthesised form where they
 * exist. A value printed without one asserts that it is exact, so a missing
 * sigma is left as a bare number rather than being invented as zero.
 */
function writePdCIF(b) {
    const p = b.params || {};
    const sg = b.spaceGroup || {};
    const esc = (v) => {
        const s = String(v == null ? '' : v);
        return (s === '') ? '?' : (/[\s'"]/.test(s) ? `'${s.replace(/'/g, "\\'")}'` : s);
    };
    const L = [
        '#',
        '# Powder diffraction data written by powder5',
        `# ${new Date().toISOString()}`,
        '#',
        `data_${(b.sourceName || 'powder5').replace(/[^A-Za-z0-9_]/g, '_')}`,
        '',
        `_audit_creation_date            ${new Date().toISOString().slice(0, 10)}`,
        "_audit_creation_method          'powder5'",
        `_pd_block_id                    ${esc(b.sourceName || 'powder5')}`,
        ''
    ];

    if (Number.isFinite(p.a)) {
        L.push('_cell_length_a                  ' + ioNum(p.a, 5),
               '_cell_length_b                  ' + ioNum(p.b, 5),
               '_cell_length_c                  ' + ioNum(p.c, 5),
               '_cell_angle_alpha               ' + ioNum(p.alpha, 4),
               '_cell_angle_beta                ' + ioNum(p.beta, 4),
               '_cell_angle_gamma               ' + ioNum(p.gamma, 4),
               '');
    }
    if (sg.symbol || sg.number) {
        L.push('_space_group_name_H-M_alt       ' + esc(sg.symbol || ''),
               '_space_group_IT_number          ' + (sg.number || '?'),
               '');
    }
    if (Number.isFinite(p.zeroShift)) {
        // THE TWO ZERO TAGS CARRY OPPOSITE SIGNS, DELIBERATELY.
        //
        // powder5 places a calculated Bragg peak at (2theta_calc + zeroShift)
        // on the observed axis and never modifies the data, so its own
        // convention is
        //
        //     zeroShift = 2theta_observed - 2theta_calculated
        //
        // i.e. the instrument zero error: positive means the diffractometer
        // reads high.
        //
        // The pdCIF dictionary defines its offset the other way round,
        //
        //     2theta_calibrated = 2theta_measured + 2theta_offset
        //
        // and "calibrated" is the true angle, so
        //
        //     _pd_calib_2theta_offset = 2theta_true - 2theta_measured
        //                             = -zeroShift
        //
        // Writing powder5's value into the standard tag unchanged would export
        // a zero point that every other program applies backwards -- a peak
        // displaced by twice the offset, in a file that looks perfectly valid.
        // The standard tag therefore gets the negated value, and the private
        // tag keeps powder5's own so a round trip through this program is
        // exact regardless.
        const cifOffset = -p.zeroShift;
        L.push('# Zero point.',
               '#   _pd_calib_2theta_offset follows the dictionary:',
               '#       2theta_calibrated = 2theta_measured + 2theta_offset',
               '#   _powder5_2theta_zero_shift is powder5\'s own convention:',
               '#       zeroShift = 2theta_observed - 2theta_calculated  (= -offset)',
               '_pd_calib_2theta_offset         ' + ioNum(cifOffset, 5),
               '_powder5_2theta_zero_shift      ' + ioNum(p.zeroShift, 5),
               '');
    }
    if (Number.isFinite(p.lambda)) {
        // THE LOOP FORM, ALWAYS -- because the reflection loop below writes
        // _pd_refln_wavelength_id and that has to refer to something. With the
        // scalar tag there was no wavelength list for the id to index, which
        // made the file quietly self-inconsistent.
        //
        // It is also the only way to state a Ka1/Ka2 pair. The dictionary's
        // _wt is a relative weight, so the Ka2 row carries the ratio directly
        // and no private tag is needed: ratio = wt(2) / wt(1). Neither was
        // written before, so a doublet measurement exported as if it were
        // monochromatic -- and read back as one.
        const dbl = (p.ratio > 1e-6 && p.lambda2 > 1e-6 &&
                     Math.abs(p.lambda - p.lambda2) > 1e-6);
        L.push('loop_',
               '    _diffrn_radiation_wavelength_id',
               '    _diffrn_radiation_wavelength',
               '    _diffrn_radiation_wavelength_wt',
               '  1  ' + p.lambda.toFixed(6) + '  1.0');
        if (dbl) L.push('  2  ' + p.lambda2.toFixed(6) + '  ' + ioNum(p.ratio, 5));
        L.push('', "_diffrn_radiation_probe         x-ray", '');
    }

    // ---- the profile model, in private tags -------------------------------
    //
    // The pdCIF dictionary has _pd_proc_ls_profile_function as free text and
    // nothing at all for Caglioti terms, Lorentzian broadening, asymmetry or
    // Stephens coefficients. They are therefore written under a _powder5_
    // prefix, which is what a private extension prefix is for and what GSAS and
    // TOPAS both do with the same problem.
    //
    // Without them a powder5 pdCIF could not reproduce its own fit: reading one
    // back restored the pattern and lost the model that described it.
    {
        // Everything already carried by a standard tag above, plus the
        // non-numeric bookkeeping that would not survive a toFixed().
        const STANDARD = new Set(['a', 'b', 'c', 'alpha', 'beta', 'gamma',
                                  'lambda', 'lambda2', 'ratio', 'zeroShift',
                                  'profileType', 'system', 'polarization',
                                  '__invalidFields']);
        const keys = Object.keys(p).filter(
            k => !STANDARD.has(k) && !k.startsWith('__') && Number.isFinite(p[k])).sort();
        if (p.profileType || keys.length) {
            L.push('# The profile model. Private tags: the dictionary has no items for these.');
            if (p.profileType) {
                L.push('_powder5_profile_type           ' + esc(p.profileType));
            }
            // The polarisation model, which decides Lp and therefore every
            // |F|^2 in the reflection loop below. Taken from the three controls
            // verbatim rather than from params.polarization, whose internal key
            // names are not something this writer should have to know.
            const pol = b.polarisationUi || {};
            if (pol.mode) L.push('_powder5_pol_mode               ' + esc(pol.mode));
            if (Number.isFinite(pol.monoTth)) {
                L.push('_powder5_pol_mono_2theta        ' + ioNum(pol.monoTth, 4));
            }
            if (Number.isFinite(pol.fraction)) {
                L.push('_powder5_pol_fraction           ' + ioNum(pol.fraction, 4));
            }
            for (const k of keys) {
                L.push(('_powder5_param_' + k).padEnd(32) + ioNum(p[k], 6));
            }
            L.push('');
        }
    }
    if (b.stats) {
        if (Number.isFinite(b.stats.rwp))  L.push('_pd_proc_ls_prof_wR_factor      ' + (b.stats.rwp / 100).toFixed(5));
        if (Number.isFinite(b.stats.r_p))  L.push('_pd_proc_ls_prof_R_factor       ' + (b.stats.r_p / 100).toFixed(5));
        if (Number.isFinite(b.stats.chi2)) L.push('_refine_ls_goodness_of_fit_all  ' + b.stats.chi2.toFixed(4));
        L.push('');
    }

    // ---- the profile ------------------------------------------------------
    const haveCalc = !!(b.calc && b.calc.length === b.tth.length);
    const haveBkg = !!(b.bkg && b.bkg.length === b.tth.length);
    // _pd_meas_2theta_scan, NOT _pd_proc_2theta_corrected.
    //
    // These are the angles as measured. powder5 applies the zero-point
    // correction to the CALCULATED peak positions, not to the data, so the
    // observed axis has had nothing done to it. Labelling it "corrected" would
    // assert that the zero had already been taken out -- and a reader that
    // believed the tag would then apply its own refined zero on top, shifting
    // the pattern twice. The zero is written separately below as the refinable
    // quantity it is.
    // The scan extent, stated as scalars as well as listed as a column.
    //
    // parsePdCifFile() BELOW already looks for these three tags -- they are how
    // it reconstructs an axis from a file whose profile loop has no angle
    // column -- and this writer never emitted them, so a powder5 pdCIF was
    // missing the summary of its own range. They are also the fastest way for
    // a reader (human or otherwise) to see what was actually saved, which is
    // exactly what was in doubt when the exported range did not follow the
    // sliders. _inc is written only for a uniform axis; asserting a step for an
    // irregular one would misplace every point after the first gap.
    if (b.tth.length) {
        const ax = ioAxisStep(b.tth);
        L.push('_pd_meas_2theta_range_min       ' + ioNum(b.tth[0], 5),
               '_pd_meas_2theta_range_max       ' + ioNum(b.tth[b.tth.length - 1], 5));
        if (ax.uniform) L.push('_pd_meas_2theta_range_inc       ' + ioNum(ax.step, 6));
        L.push('_pd_meas_number_of_points       ' + b.tth.length, '');
    }

    L.push('loop_',
           '    _pd_proc_point_id',
           '    _pd_meas_2theta_scan',
           '    _pd_meas_intensity_total');
    if (haveCalc) L.push('    _pd_calc_intensity_total');
    if (haveBkg)  L.push('    _pd_proc_intensity_bkg_calc');
    // A MISSING VALUE IS "?", NOT ZERO.
    //
    // ioNum() turns a non-finite value into 0.000, which in a calculated
    // intensity column is not a missing value at all -- it is the assertion
    // that the model predicts no counts there. Over a region outside the fitted
    // range that reads as a calculated pattern which collapses to the axis, and
    // any program computing a residual from this file gets the observed counts
    // back as the difference. CIF has a null token; this uses it.
    const cifNum = (v, dp) => Number.isFinite(v) ? v.toFixed(dp) : '?';
    for (let i = 0; i < b.tth.length; i++) {
        let row = `  ${i + 1}  ${ioNum(b.tth[i], 5)}  ${cifNum(b.obs[i], 3)}`;
        if (haveCalc) row += `  ${cifNum(b.calc[i], 3)}`;
        if (haveBkg)  row += `  ${cifNum(b.bkg[i], 3)}`;
        L.push(row);
    }
    L.push('');

    // ---- the reflections --------------------------------------------------
    //
    // _refln_F_squared_meas MEANS |F|^2, so it gets |F|^2.
    //
    // This loop used to write ioNum(r.intensity, 3) -- the raw refined peak
    // height -- into that tag, and '?' into the esd for every reflection,
    // because r.sigma is a field nothing ever sets. The result disagreed with
    // the same program's own report by 163x to 1562x depending on the
    // reflection. reflectionStructureFactors() is now the single source of the
    // conversion and the report reads it too, so the two cannot drift apart
    // again.
    const refl = (b.hklList || []).filter(
        r => r && Number.isFinite(r.h_orig) && Number.isFinite(r.intensity));
    if (refl.length) {
        // The esd of the fitted HEIGHT. Pawley heights are least-squares
        // parameters and carry a covariance entry; Le Bail heights are not, but
        // the extraction propagates counting statistics into hkl.I_sigma.
        // reflectionStructureFactors() converts height -> area -> |F|^2.
        const esds = b.esds || null;
        const sigmaArea = (hkl, widthFactor) => {
            if (esds) {
                const key = `I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`;
                const sh = esds[key];
                if (Number.isFinite(sh)) return sh * widthFactor;
            }
            if (Number.isFinite(hkl.I_sigma)) {
                return hkl.I_sigma * (Number.isFinite(b.scaleFactor) ? b.scaleFactor : 1);
            }
            return NaN;
        };
        const sf = reflectionStructureFactors(refl, b.params, {
            scale: Number.isFinite(b.scaleFactor) ? b.scaleFactor : 1,
            sigmaArea
        });

        L.push('# |F|^2 = I(area) / ( m * Lp * (1 + I2/I1) ), on an arbitrary scale.',
               '#   m  is the powder multiplicity, Lp the Lorentz-polarisation factor at the',
               '#   corrected angle, and the last term removes the Ka2 contribution to the',
               '#   fitted area (it is 1 when the Ka2/Ka1 ratio is 0).',
               '# A NEGATIVE value is a real measurement, not an error: the refined intensity',
               '#   of a weak reflection can fall below zero when the background runs above',
               '#   the data. It is written with its esd so that it can be treated properly',
               '#   rather than silently replaced by nothing.',
               'loop_',
               '    _refln_index_h',
               '    _refln_index_k',
               '    _refln_index_l',
               '    _refln_d_spacing',
               '    _pd_refln_wavelength_id',
               '    _refln_F_squared_meas',
               '    _refln_F_squared_sigma');
        for (const r of sf) {
            // A reflection whose |F|^2 could not be formed at all (no Lp, no
            // profile area) is '?', which is different from a negative value
            // and has to stay different.
            L.push(['  ' + r.h, r.k, r.l,
                    ioNum(r.d, 5), 1,
                    Number.isFinite(r.Fsq) ? r.Fsq.toFixed(4) : '?',
                    Number.isFinite(r.Fsq_sigma) ? r.Fsq_sigma.toFixed(4) : '?'].join('  '));
        }
        L.push('');
    }
    return L.join('\n') + '\n';
}

/**
 * Observed structure factors from refined peak intensities.
 *
 * THE ONLY PLACE THIS ARITHMETIC EXISTS. It used to live inline in
 * generateReportContent() while writePdCIF() wrote the raw refined parameter
 * into _refln_F_squared_meas, and the two disagreed by a factor that ranged
 * from 163 to 1562 across a single pattern. Four separate steps were missing
 * from the CIF, three of them reflection-dependent, so no scale factor could
 * have reconciled them:
 *
 *   1. hkl.intensity is a peak HEIGHT, not an integrated area. The refinement
 *      parameterises the height; integratedPeakArea() supplies the profile
 *      integral that converts one into the other, and that integral grows with
 *      2theta as the peak broadens.
 *   2. m, the powder multiplicity, varies from 1 to 8 between neighbouring
 *      reflections for reasons that are pure lattice geometry.
 *   3. Lp, the Lorentz-polarisation factor, spans a factor of four over a
 *      twenty-degree window and an order of magnitude over a full pattern.
 *   4. (1 + I2/I1), because the fitted area is the sum over the Ka1/Ka2 pair.
 *      The only constant of the four, and 1 exactly when ratio is 0.
 *
 *   |F|^2 = area / ( m * Lp * (1 + I2/I1) )
 *
 * The overall scale is arbitrary and cancels out of everything downstream, so
 * `scale` only matters for reproducing a particular report's numbers.
 *
 * NEGATIVE VALUES ARE RETURNED, NOT CLAMPED. A weak reflection whose refined
 * intensity came out below zero is a real measurement -- the background ran
 * above the data there -- and CIF accepts a negative _refln_F_squared_meas.
 * Clamping to zero would assert a measurement nobody made, and dropping the
 * row would bias the set towards positive noise. The caller decides what to do
 * with it; the report has no square root to print, a CIF reader with the esd
 * beside it can apply a French-Wilson correction.
 *
 * @param {object[]} hklList
 * @param {object} params        refined profile/cell parameters
 * @param {object} [opts]
 * @param {number} [opts.scale]         intensity scale factor (default 1)
 * @param {object} [opts.polarisation]  model from polarizationFromParams()
 * @param {function} [opts.sigmaArea]   (hkl, widthFactor) -> sigma(area), or NaN
 * @returns {object[]} one entry per reflection, in input order
 */
function reflectionStructureFactors(hklList, params, opts = {}) {
    const p = params || {};
    const scale = Number.isFinite(opts.scale) ? opts.scale : 1;
    const pol = opts.polarisation
        || (typeof polarizationFromParams === 'function' ? polarizationFromParams(p) : null);

    // 1 unless there really is a second wavelength with a non-zero weight.
    // The same guard as everywhere else in powder5, so ratio = 0 disables it
    // exactly rather than approximately.
    const dbl = (p.ratio > 1e-6 && p.lambda2 > 1e-6 &&
                 Math.abs(p.lambda - p.lambda2) > 1e-6);
    const doubletSum = dbl ? (1 + p.ratio) : 1;

    return (hklList || []).map(hkl => {
        const out = {
            hkl, h: hkl.h_orig, k: hkl.k_orig, l: hkl.l_orig, d: hkl.d,
            tthCorr: NaN, m: 1, lp: NaN, area: NaN,
            Fsq: NaN, Fsq_sigma: NaN, doubletSum
        };
        if (!hkl || !Number.isFinite(hkl.tth)) return out;

        // Lp is evaluated at the CORRECTED angle, where the peak actually
        // sits. Using the ideal Bragg angle would put Lp out of step with the
        // position reported beside it.
        let tthCorr = hkl.tth + (p.zeroShift || 0);
        try {
            if (typeof calculatePeakShift === 'function') tthCorr += calculatePeakShift(hkl.tth, p);
        } catch { /* a shift model that cannot evaluate here leaves the angle alone */ }
        out.tthCorr = tthCorr;

        out.m = (Number.isFinite(hkl.multiplicity) && hkl.multiplicity > 0) ? hkl.multiplicity : 1;
        out.lp = (typeof lorentzPolarization === 'function') ? lorentzPolarization(tthCorr, pol) : NaN;
        out.area = (typeof integratedPeakArea === 'function')
            ? integratedPeakArea(hkl, p, scale) : NaN;

        if (!Number.isFinite(out.lp) || !(out.lp > 0) || !Number.isFinite(out.area)) return out;

        const denom = out.m * out.lp * doubletSum;
        out.Fsq = out.area / denom;

        if (typeof opts.sigmaArea === 'function') {
            // Area = Height * widthFactor, so sigma(Area) = sigma(Height) * widthFactor.
            const h = hkl.intensity || 0;
            const widthFactor = Math.abs(h) > 1e-9 ? (out.area / h) : 0;
            const sa = opts.sigmaArea(hkl, widthFactor);
            // Linear in the same denominator: no square root, no sign question.
            if (Number.isFinite(sa)) out.Fsq_sigma = Math.abs(sa) / denom;
        }
        return out;
    });
}

/**
 * Render a bundle in the requested format.
 * @param {string} formatId
 * @param {ExportBundle} bundle
 * @returns {{text: string, ext: string}}
 */
function exportPattern(formatId, bundle) {
    const spec = DATA_EXPORT_FORMATS.find(f => f.id === formatId);
    if (!spec) throw new Error(`Unknown export format "${formatId}".`);
    const writers = { xy: writeXY, xycalc: writeXYCalc, uxd: writeUXD, xra: writeXRA, pdcif: writePdCIF };
    return { text: writers[formatId](bundle), ext: spec.ext };
}

// ---------------------------------------------------------------------------
//  pdCIF READER
// ---------------------------------------------------------------------------

/**
 * Read a powder profile out of a pdCIF.
 *
 * Deliberately tolerant about which tags carry the angle and the intensity:
 * different producers write _pd_proc_2theta_corrected, _pd_meas_2theta_scan or
 * _pd_meas_2theta_range_* and a file that uses the range form has no angle
 * column at all -- the axis is implied by min, max and the point count. A
 * reader that insists on one spelling rejects most real files.
 *
 * Standard uncertainties in parentheses are stripped: "1234(12)" is a value of
 * 1234, not a parse error.
 *
 * Everything the writer states is read back, not just the profile: the cell,
 * the space group, both wavelengths with the Ka2 ratio, the zero point and the
 * profile model. A file that describes a refinement and a reader that returns
 * only two columns of numbers meant a powder5 pdCIF could not restore the fit
 * it came from -- and the loss was silent, because the tags were all there.
 *
 * @param {string} content
 * @returns {{tth: number[], intensity: number[], wavelength: number|null,
 *            wavelength2: number|null, ratio21: number|null,
 *            zeroShift: number|null, cell: object|null,
 *            spaceGroup: object|null, profileType: string|null,
 *            profileParams: object, polarisation: object|null}}
 */
function parsePdCifFile(content) {
    const lines = content.split(/\r?\n/);
    const num = (tok) => {
        if (tok == null) return NaN;
        const t = String(tok).replace(/\(\d+\)$/, '');   // strip the s.u.
        if (t === '?' || t === '.') return NaN;
        return parseFloat(t);
    };

    // ---- scalars we care about -------------------------------------------
    let wavelength = null, rangeMin = NaN, rangeMax = NaN, rangeInc = NaN;
    let zeroShift = null, cifOffset = null;
    let wavelength2 = null, ratio21 = null;
    const cell = {};
    const sgIn = { symbol: null, number: null };
    let profileType = null;
    const profileParams = {};
    const polarisation = { mode: null, monoTth: NaN, fraction: NaN };
    // Quoted CIF strings arrive with their delimiters attached.
    const unq = (v) => v.replace(/^['"]|['"]$/g, '').trim();
    const CELL_TAG = {
        '_cell_length_a': 'a', '_cell_length_b': 'b', '_cell_length_c': 'c',
        '_cell_angle_alpha': 'alpha', '_cell_angle_beta': 'beta', '_cell_angle_gamma': 'gamma'
    };
    for (const raw of lines) {
        const s = raw.trim();
        const m = s.match(/^(_[\w.\-\[\]]+)\s+(.*)$/);
        if (!m) continue;
        const tag = m[1].toLowerCase(), val = m[2].trim();
        if (CELL_TAG[tag] && Number.isFinite(num(val))) cell[CELL_TAG[tag]] = num(val);
        // Both the modern and the deprecated spellings: files in the wild use
        // either, and rejecting the old one would fail on most of them.
        else if ((tag === '_space_group_name_h-m_alt' || tag === '_symmetry_space_group_name_h-m')
                 && unq(val) && unq(val) !== '?') sgIn.symbol = unq(val);
        else if ((tag === '_space_group_it_number' || tag === '_symmetry_int_tables_number')
                 && Number.isFinite(num(val))) sgIn.number = num(val);
        else if (tag === '_powder5_profile_type') profileType = unq(val) || null;
        else if (tag === '_powder5_pol_mode') polarisation.mode = unq(val) || null;
        else if (tag === '_powder5_pol_mono_2theta') polarisation.monoTth = num(val);
        else if (tag === '_powder5_pol_fraction') polarisation.fraction = num(val);
        else if (tag.startsWith('_powder5_param_') && Number.isFinite(num(val))) {
            profileParams[m[1].slice('_powder5_param_'.length)] = num(val);
        }
        else if (tag === '_diffrn_radiation_wavelength' && Number.isFinite(num(val))) wavelength = num(val);
        else if (tag === '_pd_meas_2theta_range_min') rangeMin = num(val);
        else if (tag === '_pd_meas_2theta_range_max') rangeMax = num(val);
        else if (tag === '_pd_meas_2theta_range_inc') rangeInc = num(val);
        // The private tag wins where present -- it is already in powder5's
        // convention. The standard tag is NEGATED on the way in, because the
        // dictionary defines
        //     2theta_calibrated = 2theta_measured + 2theta_offset
        // and powder5 wants 2theta_observed - 2theta_calculated, which is the
        // same quantity with the opposite sign.
        else if (tag === '_powder5_2theta_zero_shift' && Number.isFinite(num(val))) zeroShift = num(val);
        else if (tag === '_pd_calib_2theta_offset' && Number.isFinite(num(val)) &&
                 cifOffset === null) cifOffset = num(val);
    }

    // ---- the wavelength loop, if the file uses one ------------------------
    //
    // Scanned separately from the profile loop because it comes first and is
    // tiny. A file written by this program always uses the loop form (the
    // reflection loop refers to a wavelength id), and the loop is the only
    // place a Ka2 component can be stated at all.
    {
        let j = 0;
        while (j < lines.length) {
            if (lines[j].trim().toLowerCase() !== 'loop_') { j++; continue; }
            j++;
            const tags = [];
            while (j < lines.length && lines[j].trim().startsWith('_')) {
                tags.push(lines[j].trim().toLowerCase()); j++;
            }
            const iW = tags.indexOf('_diffrn_radiation_wavelength');
            if (iW < 0) continue;
            const iId = tags.indexOf('_diffrn_radiation_wavelength_id');
            const iWt = tags.indexOf('_diffrn_radiation_wavelength_wt');
            const rows = [];
            while (j < lines.length) {
                const t = lines[j].trim();
                if (t === '' || t.startsWith('_') || t.startsWith('#') ||
                    t.toLowerCase() === 'loop_' || t.toLowerCase().startsWith('data_')) break;
                const tok = t.split(/\s+/);
                if (tok.length >= tags.length) {
                    rows.push({ id: iId >= 0 ? String(tok[iId]) : String(rows.length + 1),
                                lam: num(tok[iW]),
                                wt: iWt >= 0 ? num(tok[iWt]) : 1 });
                }
                j++;
            }
            const good = rows.filter(r => Number.isFinite(r.lam));
            if (!good.length) continue;
            // Sorted by id so "1" is the primary line whatever order it was
            // written in; the weights are relative to that one.
            good.sort((x, y) => String(x.id).localeCompare(String(y.id), undefined, { numeric: true }));
            wavelength = good[0].lam;
            if (good.length > 1) {
                wavelength2 = good[1].lam;
                const w1 = Number.isFinite(good[0].wt) && good[0].wt !== 0 ? good[0].wt : 1;
                const w2 = Number.isFinite(good[1].wt) ? good[1].wt : NaN;
                if (Number.isFinite(w2)) ratio21 = w2 / w1;
            }
            break;
        }
    }

    // ---- find the profile loop -------------------------------------------
    const ANGLE = ['_pd_proc_2theta_corrected', '_pd_meas_2theta_scan',
                   '_pd_meas_2theta_corrected'];
    const INTEN = ['_pd_meas_intensity_total', '_pd_meas_counts_total',
                   '_pd_proc_intensity_total', '_pd_proc_intensity_net'];

    const tth = [], intensity = [];
    let i = 0;
    while (i < lines.length) {
        if (lines[i].trim().toLowerCase() !== 'loop_') { i++; continue; }
        i++;
        const tags = [];
        while (i < lines.length && lines[i].trim().startsWith('_')) {
            tags.push(lines[i].trim().toLowerCase()); i++;
        }
        const iA = tags.findIndex(t => ANGLE.includes(t));
        const iI = tags.findIndex(t => INTEN.includes(t));
        if (iI < 0) continue;      // not the profile loop; skip its body

        while (i < lines.length) {
            const s = lines[i].trim();
            if (s === '' || s.startsWith('_') || s.startsWith('#') ||
                s.toLowerCase() === 'loop_' || s.toLowerCase().startsWith('data_')) break;
            const tok = s.split(/\s+/);
            if (tok.length >= tags.length) {
                const y = num(tok[iI]);
                if (Number.isFinite(y)) {
                    intensity.push(y);
                    tth.push(iA >= 0 ? num(tok[iA]) : NaN);
                }
            }
            i++;
        }
        if (intensity.length) break;
    }

    if (!intensity.length) {
        throw new Error('No powder profile loop found. The file needs a loop_ containing ' +
                        '_pd_meas_intensity_total (or _pd_proc_intensity_total).');
    }

    // The range form: no angle column, the axis is implied.
    if (!Number.isFinite(tth[0])) {
        let step = rangeInc;
        if (!Number.isFinite(step) && Number.isFinite(rangeMin) && Number.isFinite(rangeMax)) {
            step = (rangeMax - rangeMin) / Math.max(1, intensity.length - 1);
        }
        if (!Number.isFinite(rangeMin) || !Number.isFinite(step)) {
            throw new Error('The profile loop has no 2-theta column, and the file does not ' +
                            'give _pd_meas_2theta_range_min / _inc to reconstruct the axis.');
        }
        for (let k = 0; k < intensity.length; k++) tth[k] = rangeMin + k * step;
    }
    // The private tag has priority; otherwise convert the dictionary's offset.
    if (zeroShift === null && cifOffset !== null) zeroShift = -cifOffset;

    return {
        tth, intensity, wavelength, zeroShift,
        wavelength2, ratio21,
        // Only handed back when the cell is actually complete enough to use.
        // A partial cell applied to the controls would leave a mix of this
        // file's edges and the previous sample's.
        cell: Number.isFinite(cell.a) ? cell : null,
        spaceGroup: (sgIn.symbol || Number.isFinite(sgIn.number)) ? sgIn : null,
        profileType,
        profileParams,
        polarisation: (polarisation.mode || Number.isFinite(polarisation.monoTth) ||
                       Number.isFinite(polarisation.fraction)) ? polarisation : null
    };
}

// ---------------------------------------------------------------------------
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        DATA_EXPORT_FORMATS, availableExportFormats, exportPattern, parsePdCifFile,
        writeXY, writeXYCalc, writeUXD, writeXRA, writePdCIF, ioAxisStep
    };
}
