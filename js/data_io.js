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
        L.push('_diffrn_radiation_wavelength    ' + p.lambda.toFixed(6),
               "_diffrn_radiation_probe         x-ray",
               '');
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
    L.push('loop_',
           '    _pd_proc_point_id',
           '    _pd_meas_2theta_scan',
           '    _pd_meas_intensity_total');
    if (haveCalc) L.push('    _pd_calc_intensity_total');
    if (haveBkg)  L.push('    _pd_proc_intensity_bkg_calc');
    for (let i = 0; i < b.tth.length; i++) {
        let row = `  ${i + 1}  ${ioNum(b.tth[i], 5)}  ${ioNum(b.obs[i], 3)}`;
        if (haveCalc) row += `  ${ioNum(b.calc[i], 3)}`;
        if (haveBkg)  row += `  ${ioNum(b.bkg[i], 3)}`;
        L.push(row);
    }
    L.push('');

    // ---- the reflections --------------------------------------------------
    const refl = (b.hklList || []).filter(
        r => r && Number.isFinite(r.h_orig) && Number.isFinite(r.intensity));
    if (refl.length) {
        L.push('loop_',
               '    _refln_index_h',
               '    _refln_index_k',
               '    _refln_index_l',
               '    _refln_d_spacing',
               '    _pd_refln_wavelength_id',
               '    _refln_F_squared_meas',
               '    _refln_F_squared_sigma');
        for (const r of refl) {
            L.push(['  ' + r.h_orig, r.k_orig, r.l_orig,
                    ioNum(r.d, 5), 1,
                    ioNum(r.intensity, 3),
                    Number.isFinite(r.sigma) ? ioNum(r.sigma, 3) : '?'].join('  '));
        }
        L.push('');
    }
    return L.join('\n') + '\n';
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
 * @param {string} content
 * @returns {{tth: number[], intensity: number[], wavelength: number|null}}
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
    for (const raw of lines) {
        const s = raw.trim();
        const m = s.match(/^(_[\w.\-\[\]]+)\s+(.*)$/);
        if (!m) continue;
        const tag = m[1].toLowerCase(), val = m[2].trim();
        if (tag === '_diffrn_radiation_wavelength' && Number.isFinite(num(val))) wavelength = num(val);
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

    return { tth, intensity, wavelength, zeroShift };
}

// ---------------------------------------------------------------------------
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        DATA_EXPORT_FORMATS, availableExportFormats, exportPattern, parsePdCifFile,
        writeXY, writeXYCalc, writeUXD, writeXRA, writePdCIF, ioAxisStep
    };
}
