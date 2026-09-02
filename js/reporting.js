/**
 * How many significant figures the standard uncertainty is quoted to.
 *
 * Crystallographic practice is one figure, except when the leading digit is 1:
 * rounding 0.00012 to 0.0001 misstates the uncertainty by 17%, and the reader
 * cannot tell 0.00010 from 0.00019. With a leading digit of 9 the same
 * rounding is worth at most 5%, so one figure carries its own precision.
 *
 * Set to 3 if you prefer the "1 or 2 leads to two figures" variant.
 * @type {number}
 */
const ESD_TWO_FIGURE_THRESHOLD = 2;

/** Decimal places beyond which a value is noise regardless of its su. */
const ESD_MAX_DECIMALS = 12;

/**
 * Round to a given number of decimal places, negative places included.
 * @param {number} v
 * @param {number} dp Negative rounds to tens, hundreds, ...
 * @returns {number}
 */
function roundTo(v, dp) {
    const f = Math.pow(10, dp);
    return Math.round(v * f) / f;
}

/**
 * The decimal place a standard uncertainty determines.
 *
 * THE EXPONENT COMES FROM toExponential(), NOT Math.log10. log10(0.001) can
 * evaluate to -2.9999999999999996, and floor() then returns -3 for one input
 * and -2 for its neighbour, which would move the decimal point of the reported
 * value by a factor of ten for a hairsbreadth change in the su.
 *
 * @param {number} esd
 * @returns {{dp:number, esd:number, digits:number}|null} null if the su is not
 *          a positive finite number.
 */
function esdDecimals(esd) {
    if (!isFinite(esd) || esd <= 0) return null;
    let [mStr, eStr] = esd.toExponential().split('e');
    let e = parseInt(eStr, 10);
    let sig = (parseFloat(mStr) < ESD_TWO_FIGURE_THRESHOLD) ? 2 : 1;
    let dp = sig - 1 - e;
    let r = roundTo(esd, dp);

    // Rounding can promote: 0.096 to one figure is 0.1, whose leading digit is
    // now 1 and which therefore wants two figures after all. One re-check
    // converges, because the second rounding cannot promote again.
    const parts = r.toExponential().split('e');
    if (parseInt(parts[1], 10) !== e) {
        e = parseInt(parts[1], 10);
        sig = (parseFloat(parts[0]) < ESD_TWO_FIGURE_THRESHOLD) ? 2 : 1;
        dp = sig - 1 - e;
        r = roundTo(esd, dp);
    }
    if (dp > ESD_MAX_DECIMALS) return null;
    return { dp, esd: r, digits: Math.round(r * Math.pow(10, dp)) };
}

/**
 * A value and its uncertainty in the crystallographic compact form:
 * 5.5139(3), where the parenthesis is the su in units of the last decimal.
 *
 * WHY NOT JUST PRINT MORE DIGITS. A cell edge shown as 5.439771388246059 with
 * its su in a separate column asks the reader to do this rounding themselves,
 * and most will copy the full string into a paper. Everything past the digit
 * the su reaches is the binary representation of a double, not a measurement.
 *
 * Falls back to `value +/- su` when the su is 10 or more, where the compact
 * form has no last decimal to refer to and would be ambiguous.
 *
 * @param {number} value
 * @param {number} esd
 * @returns {string|null} null when no su is available; the caller then
 *          formats the value on its own.
 */
function formatWithEsd(value, esd) {
    if (value === null || value === undefined || !isFinite(value)) return null;
    const d = esdDecimals(esd);
    if (!d) return null;
    if (d.dp < 0) return `${roundTo(value, d.dp)} +/- ${d.esd}`;
    return `${value.toFixed(d.dp)}(${d.digits})`;
}

/**
 * The data file the reports refer to.
 *
 * `controls.fileName.textContent` is the on-screen label, which reads
 * "No data file loaded." when there is none; `dataset.file` is the actual
 * name and is absent in that case. Preferring the dataset means a report
 * either names a real file or says N/A, rather than embedding a sentence
 * meant for a status line.
 *
 * The first of the three report headers dereferenced controls.fileName
 * without a guard, so a report generated before the controls were wired
 * would have thrown here rather than anywhere informative.
 *
 * @returns {string}
 */
function reportDataFileName() {
    if (typeof controls === 'undefined' || !controls || !controls.fileName) return 'N/A';
    const el = controls.fileName;
    const name = (el.dataset && el.dataset.file) || '';
    if (name) return name;
    const shown = (el.textContent || '').trim();
    return (shown && !/^no data/i.test(shown)) ? shown : 'N/A';
}

/**
 * A number at a sensible number of digits, with trailing zeros trimmed.
 *
 * `${p.a}` printed a cell edge as 5.439771388246059 -- sixteen significant
 * figures for a quantity whose uncertainty is in the fourth decimal. Beyond
 * about the fifth decimal those digits are the binary representation of the
 * double, not a measurement, and printing them invites someone to quote them.
 *
 * Trailing zeros go too, so a right angle reads 90 rather than 90.0000.
 *
 * @param {number} v
 * @param {number} dp Maximum decimal places.
 * @returns {string}
 */
function fmtNum(v, dp) {
    if (v === null || v === undefined || !isFinite(v)) return '-';
    let out = Number(v).toFixed(dp);
    if (out.indexOf('.') >= 0) out = out.replace(/0+$/, '').replace(/\.$/, '');
    return out;
}

/**
 * Reduce report text to pure ASCII, for jsPDF's built-in fonts.
 *
 * WHY THIS HAS TO BE EXHAUSTIVE, NOT A LIST OF THE CHARACTERS SOMEONE HIT.
 * The previous version named eleven characters and missed the rest. A line
 * containing an unlisted one -- alpha and gamma in the refined-cell block --
 * came out as
 *
 *     ±  =  9 0 d e g,   B e t a  =  9 0 d e g,   ³  =  9 0 d e g
 *
 * two failures at once: alpha (U+03B1) and gamma (U+03B3) rendered as the
 * WinAnsi characters at 0xB1 and 0xB3, and jsPDF switched that whole line to a
 * two-byte encoding, whose high bytes printed as the gaps between every
 * letter. Beta survived only because it happened to be on the list.
 *
 * So the rule is inverted: everything is mapped, and whatever is left after
 * mapping is forced into ASCII rather than passed through to be guessed at.
 * The catch-all is what makes this safe against report text nobody has written
 * yet -- a data file named "echantillon-2theta.raw", a mineral name with an
 * umlaut, a future label with a subscript.
 *
 * @param {string} text
 * @returns {string} ASCII only, 0x20 to 0x7E plus newline.
 */
function pdfAscii(text) {
    let out = String(text);

    // Multi-character forms first, before their components are replaced.
    out = out.replace(/χ²/g, 'Chi^2').replace(/Å⁻¹/g, '1/A');

    // Degree sign. `90°` needs the space that `deg` does not carry; `(°)` as a
    // column unit does not, or the header reads "2theta ( deg)".
    out = out.replace(/([\d)])\s*°/g, '$1 deg').replace(/°/g, 'deg');

    const MAP = {
        // Greek, lower case
        'α':'alpha','β':'beta','γ':'gamma','δ':'delta','ε':'epsilon','ζ':'zeta',
        'η':'eta','θ':'theta','ι':'iota','κ':'kappa','λ':'lambda','μ':'mu',
        'ν':'nu','ξ':'xi','ο':'o','π':'pi','ρ':'rho','ς':'s','σ':'sigma',
        'τ':'tau','υ':'upsilon','φ':'phi','χ':'chi','ψ':'psi','ω':'omega',
        // Greek, upper case
        'Α':'A','Β':'B','Γ':'Gamma','Δ':'Delta','Ε':'E','Ζ':'Z','Η':'H',
        'Θ':'Theta','Ι':'I','Κ':'K','Λ':'Lambda','Μ':'M','Ν':'N','Ξ':'Xi',
        'Ο':'O','Π':'Pi','Ρ':'P','Σ':'Sigma','Τ':'T','Υ':'Y','Φ':'Phi',
        'Χ':'Chi','Ψ':'Psi','Ω':'Omega',
        // Letters and units
        'Å':'A','å':'a','µ':'u',
        // Superscripts and subscripts
        '²':'^2','³':'^3','¹':'^1','⁰':'^0','⁴':'^4','⁻':'^-',
        '₀':'0','₁':'1','₂':'2','₃':'3','₄':'4',
        // Maths
        '×':'x','·':'.','±':'+/-','−':'-','÷':'/','∞':'inf','√':'sqrt',
        '≈':'~','≠':'!=','≤':'<=','≥':'>=','∑':'sum','∫':'int','∂':'d',
        '→':'->','←':'<-','↔':'<->','∙':'.',
        // Punctuation and marks
        '–':'-','—':'-','‘':"'",'’':"'",'“':'"','”':'"','…':'...',
        '«':'<<','»':'>>','•':'*','·':'.','†':'+','‡':'++','′':"'",'″':'"',
        '✓':'y','✗':'x','✕':'x','⚠':'!','○':'o','◉':'*','●':'*','■':'#','□':'-'
    };
    out = out.replace(/[^\x00-\x7F]/g, ch => (ch in MAP) ? MAP[ch] : ch);

    // Anything still outside ASCII: strip the accent if it has one (so a file
    // named "echantillon" keeps its letters), and only then give up.
    out = out.replace(/[^\x00-\x7F]/g, (ch) => {
        const bare = ch.normalize('NFD').replace(/[\u0300-\u036f]/g, '');
        return /^[\x00-\x7F]+$/.test(bare) ? bare : '?';
    });

    // Control characters other than newline would end the text object early.
    return out.replace(/[\x00-\x08\x0b-\x1f\x7f]/g, ' ');
}

/**
 * A PDF from monospaced text, optionally followed by the main chart.
 *
 * The body of generatePdfReport() was a text-layout loop and a canvas grab
 * with no dependence on what produced the text, so it is now shared. The
 * theoretical-pattern export (no data loaded, no fit) needs exactly this and
 * nothing else; duplicating forty lines of jsPDF layout to get it would have
 * meant two places to keep in step over page breaks and the WinAnsi
 * transliteration.
 *
 * @param {string} text        Report body. Lines starting '---' are headed.
 * @param {string} filename    Without the .pdf extension.
 * @param {boolean} [withChart=true] Append the main chart on its own page.
 */
function renderTextPdf(text, filename, withChart = true) {
    if (!window.jspdf || !window.jspdf.jsPDF) {
        throw new Error('jspdf.umd.min.js was not loaded. It belongs at lib/jspdf.umd.min.js, next to the other vendor bundles.');
    }
    const { jsPDF } = window.jspdf;
    const doc = new jsPDF({ orientation: 'p', unit: 'mm', format: 'a4' });
    const margin = 15;
    const contentWidth = doc.internal.pageSize.getWidth() - 2 * margin;
    let yPosition = 20;

    const pdfText = pdfAscii(text);

    doc.setFont('Courier');
    for (const line of pdfText.split('\n')) {
        if (yPosition > doc.internal.pageSize.getHeight() - 20) { doc.addPage(); yPosition = 20; }
        const isHeader = line.startsWith('---');
        const isTitle = line.includes('Refinement Report') || line.includes('Theoretical HKL');
        let fontSize = 9, fontStyle = 'normal';
        if (isTitle) { fontSize = 14; fontStyle = 'bold'; yPosition += 6; }
        else if (isHeader) { fontSize = 10; fontStyle = 'bold'; yPosition += 4; }
        doc.setFontSize(fontSize);
        doc.setFont(undefined, fontStyle);

        // A Courier line wider than the text width is silently CLIPPED by
        // jsPDF -- columns would just vanish off the right edge. Shrink the
        // offending line instead; monospace keeps the row aligned.
        try {
            const w = doc.getTextWidth(line);
            if (w > contentWidth && w > 0) doc.setFontSize(Math.max(4.5, fontSize * (contentWidth / w)));
        } catch { /* getTextWidth unavailable: keep the full size */ }

        doc.text(line, margin, yPosition);
        doc.setFontSize(fontSize);
        yPosition += fontSize * 0.4;
    }

    if (withChart && typeof controls !== 'undefined' && controls.mainChartCanvas) {
        const chartCanvas = controls.mainChartCanvas;
        doc.addPage();
        // Composite onto white first: the chart canvas is transparent, and a
        // transparent PNG in a PDF renders as black in several viewers.
        const tmp = document.createElement('canvas');
        tmp.width = chartCanvas.width;
        tmp.height = chartCanvas.height;
        const ctx = tmp.getContext('2d');
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, tmp.width, tmp.height);
        ctx.drawImage(chartCanvas, 0, 0);
        const img = tmp.toDataURL('image/png', 1.0);
        const h = Math.min((tmp.height * contentWidth) / tmp.width, 250);
        doc.addImage(img, 'PNG', margin, 20, contentWidth, h);
    }

    // Date AND time. Two reports of the same tab on the same day collided on
    // one name, and the browser silently appended "(1)" -- which tells you
    // there are two files but not which is which.
    const st = new Date();
    const stamp = st.getFullYear()
        + String(st.getMonth() + 1).padStart(2, '0')
        + String(st.getDate()).padStart(2, '0')
        + '-' + String(st.getHours()).padStart(2, '0')
        + String(st.getMinutes()).padStart(2, '0')
        + String(st.getSeconds()).padStart(2, '0');
    doc.save(`${filename}-${stamp}.pdf`);
}

// Reusable PDF generator. `results` is the snapshot to report; `modeName` is
// 'Pawley' or 'LeBail'. Refactored out of the old button handler so the
// per-run history tabs can export any stored run, not just the live fit.
/**
 * @param {object} results   The snapshot to report.
 * @param {string} modeName  'Pawley', 'LeBail', 'Wyckoff', 'ChargeFlipping'.
 * @param {boolean} [withChart=true] Append the powder pattern on its own page.
 *        FALSE for charge flipping and Wyckoff: those report a density map and
 *        a set of atomic positions, and the diffraction pattern behind them is
 *        the same picture the Plot tab already exports. A page of it at the
 *        back of a structure report is noise.
 */
async function generatePdfReport(results, modeName, withChart = true) {
    window.__activeReportResults = results;
    try {
        renderTextPdf(generateReportContent('summary'), `${modeName}-Report`, withChart);
    } finally {
        window.__activeReportResults = null;
    }
}
/**
 * The covariance of a fit, factored ONCE and memoised on the results object.
 *
 * covarianceFromJtJ used to run twice for every report: here, and again inside
 * intensityOverlapClusters on the same matrix. At 2500 reflections that is
 * ~460 ms each, on the main thread, in the same task as the history snapshot.
 *
 * The cached value holds Float64Array rows, so strip `__cov` before any
 * JSON.stringify of the results (refinement_controller.js already does).
 *
 * @param {object} results
 * @returns {{cov: Array<Float64Array>, undetermined: Set<number>}|null}
 */
/**
 * One factual line about the parallel-tempering budget, or '' for a non-PT run.
 *
 * PT proposes one scalar per replica per iteration, so a parameter receives
 * only maxIter * PT_NUM_REPLICAS / n_params proposals. Stating it costs
 * nothing and is the one number that explains a PT run that went nowhere.
 *
 * @param {object} results
 * @returns {string}
 */
function ptIterationsNote(results) {
    const pt = results && results.ptIterations;
    if (!pt || !pt.used) return '';
    return `${pt.used} iterations, ${pt.proposalsPerParam.toFixed(0)} proposals/parameter `
         + `over ${pt.nParams} parameters`;
}

/**
 * A square numeric matrix of the given size, with rows that may be plain
 * arrays OR typed arrays.
 *
 * The three call sites used to test `Array.isArray(JtJ[0])`, which is false
 * for a Float64Array. That conflates "is this the right shape" with "was it
 * built with the right constructor", and the two are not the same question --
 * covarianceFromJtJ only ever indexes it.
 *
 * @param {*} m
 * @param {number} n
 * @returns {boolean}
 */
function isSquareMatrix(m, n) {
    if (!m || m.length !== n) return false;
    const r0 = m[0];
    return !!r0 && typeof r0 !== 'string' && r0.length === n;
}

function reportCovariance(results) {
    if (!results || !results.JtJ) return null;
    if (results.__cov === undefined) {
        try { results.__cov = covarianceFromJtJ(results.JtJ) || null; }
        catch (err) { console.error('covarianceFromJtJ failed:', err); results.__cov = null; }
    }
    return results.__cov;
}

function generateReportContent(format = 'summary', resultsArg = null) {
    // resultsArg lets the per-run history export a specific snapshot;
    // callers that pass nothing keep using the live global fitResults.
    const fitResults = resultsArg || window.__activeReportResults || getLiveFitResults();
    
    // -------------------------------------------------------------
    // Intercept Charge Flipping results
    // -------------------------------------------------------------
    if (fitResults && fitResults.gridSize && fitResults.bestR !== undefined) {
        return generateChargeFlippingReport(fitResults);
    }

    // A Wyckoff solution has neither a gridSize nor a bestR -- it never built a
    // map -- so it fell past the charge-flipping branch and into the Pawley
    // guard below, which correctly reported that a Pawley fit was missing. The
    // PDF came out as a single line saying the fit results were incomplete.
    if (fitResults && Array.isArray(fitResults.sites) &&
        (fitResults.source === 'wyckoff' || fitResults.sites.some(x => x && x.wyckoff))) {
        return generateWyckoffReport(fitResults);
    }

    // Ensure fit results and necessary data are available for refinement runs
    if (!fitResults || !fitResults.params || !fitResults.stats || !fitResults.fitFlags || !fitResults.hklList) {
        return "Error: Fit results are incomplete or not available.";
    }
     if (fitResults.JtJ && (!workerWorkingData || !workerWorkingData.intensity || workerWorkingData.intensity.length === 0)) {
          console.warn("generateReportContent: workingData slice not available for ESD calculation.");
     }

    const now = new Date();
    const { params: finalParams, stats, fitFlags, hklList, algorithm, refinementMode, JtJ, parameterInfo, ss_res } = fitResults;
    const mainScaleFactor = stats.scaleFactor || 1.0;

    const formatLine = (cols, widths) => {
        return cols.map((col, i) => {
             const colStr = (col === null || col === undefined) ? '-' : String(col);
             // Pad numbers to the left within their column width for alignment
             if (typeof col === 'number' && i > 0 && widths[i] > 0) {
                 return colStr.padStart(widths[i]);
             }
             return colStr.padEnd(widths[i]);
        }).join(' ');
    };

    const modeName = refinementMode === 'pawley' ? 'Pawley' : 'Le Bail';

    const header = [
        `${modeName} Refinement Report`,
        '========================================================================',
        '',
        APP_VERSION_TEXT,
        `Report Generated: ${now.toLocaleString()}`,
        `Data File: ${reportDataFileName()}`,
        ''
    ];

    const algorithmNames = {
        lm: 'Levenberg-Marquardt',
        pt: 'Parallel Tempering'
    };
    const algorithmName = algorithmNames[algorithm] || 'Unknown';

    const profileNameMap = {
        "simple_pvoigt": "Simple pVoigt",
        "tch_aniso": "TCH (Size/Strain/Aniso)",
        "split_pvoigt": "Split pVoigt (Asymmetric)"
    };
    const profileType = finalParams.profileType || "simple_pvoigt";
    const profileName = profileNameMap[profileType] || "Unknown Profile";

    const selectedSg = currentSG;
    // FIX: report the SELECTED setting's symbol (e.g. "Pcnm"), not the
    // general/standard symbol (e.g. "Pnma"). Fall back to standard if
    // no setting-specific symbol is present. Note the setting name too.
   
    const spaceGroupName = selectedSg
        ? `${selectedSg.number} – ${selectedSg.symbol || selectedSg.standard_symbol}`
          + (selectedSg.setting_description && selectedSg.setting_description !== 'standard'
               ? ` (setting: ${selectedSg.setting_description})` : '')
        : 'N/A';

    // Calculate Rexp and safely grab Rp just like in the UI
    const rexp = stats.r_exp ?? stats.rexp ?? ((stats.rwp && stats.chi2) ? (stats.rwp / Math.sqrt(stats.chi2)) : undefined);
    const rpVal = stats.r_p ?? stats.rp;

    const statsSection = [
        '--- Refinement Statistics ---',
        `Rp (%):      ${rpVal !== undefined ? rpVal.toFixed(3) : 'N/A'}`,
        `Rwp (%):     ${stats.rwp !== undefined ? stats.rwp.toFixed(3) : 'N/A'}`,
        `Rexp (%):    ${rexp !== undefined ? rexp.toFixed(3) : 'N/A'}`,
        `χ² (GOF):    ${stats.chi2 !== undefined ? stats.chi2.toFixed(3) : 'N/A'}`,
        `Reflections:  ${hklList ? hklList.length : 'N/A'}`,
        `Algorithm:    ${algorithmName}`,
        // Factual, not advisory: the iteration count run and what it bought
        // per parameter.
        ...(ptIterationsNote(fitResults) ? [`PT search:    ${ptIterationsNote(fitResults)}`] : []),
        `Refinement:   ${modeName}`,
        `Profile:      ${profileName} (${profileType})`,
        `Space Group:  ${spaceGroupName}`,
        // The background is interpolated through fixed anchors, never refined.
        // Stating it here rather than leaving it implicit: it is the reason no
        // sigma below carries background uncertainty, and the reason chi-square
        // is optimistic by roughly the anchor count.
        // The worker composes backgroundNote and distinguishes the fitted case.
        // The fallback is only reached when a result predates that field, so it
        // must not assert either way.
        `Background:   ${stats.backgroundNote || 'spline through user anchor points'}`,
        ''
    ];

    //   Lorentz-polarisation model.
    // Declared explicitly, and with the K that was actually used, so a
    // reader can reproduce every |Fo| in the reflection table below
    // from the printed intensities without guessing the convention.
    const reportPol = polarizationFromParams(finalParams);
    const lpSection = [
        '--- Lorentz-Polarisation ---',
        `Model:        ${polarizationDescription(reportPol)}`,
        ...polarizationFormulaLines(reportPol).map(s => '  ' + s),
        '  I(hkl) = s * m * Lp * |F(hkl)|^2',
        '(Lp does not enter the profile: in Le Bail and Pawley the intensity is a',
        ' free parameter and absorbs it. It is applied when converting an extracted',
        ' intensity into |F|, both in the table below and in the charge-flipping input.)',
        ''
    ];

    let whAnalysisSection = [];
    if (profileType === "tch_aniso" && hklList && hklList.length > 0 && workerWorkingData && workerWorkingData.tth) {
        try {
            const whResults = calculateWilliamsonHall(finalParams, hklList);
            if (whResults) {
                whAnalysisSection.push(
                    '--- Williamson-Hall Size/Strain Analysis (Approximate) ---',
                    `Apparent Crystallite Size (nm): ${whResults.size_nm !== undefined ? whResults.size_nm.toFixed(1) : 'N/A'}`,
                    `Apparent Microstrain (%):      ${whResults.strain_percent !== undefined ? whResults.strain_percent.toFixed(4) : 'N/A'}`,
                    `  from Y directly (nm):        ${whResults.size_direct_nm !== undefined ? whResults.size_direct_nm.toFixed(1) : 'N/A'}`,
                    `  from X,U directly (%):       ${whResults.strain_direct_percent !== undefined ? whResults.strain_direct_percent.toFixed(4) : 'N/A'}`,
                    `Linear Fit R²:                 ${whResults.r_squared !== undefined ? whResults.r_squared.toFixed(5) : 'N/A'}`,
                    '(Note: Assumes isotropic size/strain broadening via TCH U, X, Y)',
                    ''
                );
            }
         } catch (whError) {
              console.warn("Could not calculate Williamson-Hall:", whError);
              whAnalysisSection.push('--- Williamson-Hall Size/Strain Analysis ---', '(Calculation failed)', '');
         }
    }


    let esds = {};
    let esdWarning = null;
    // FIX: gate on the presence of the normal-equations matrix, not on
    // the requested algorithm. A Parallel Tempering run is now followed
    // by a short Levenberg-Marquardt polish, so it supplies a JtJ too.
    // The data slice THIS result was fitted against, frozen by rpSnapshotFit,
    // in preference to the live one. Selecting an earlier run from the history
    // dropdown after changing the 2-theta range used to recompute its ESDs
    // against the new range -- N would be wrong, and with it every sigma.
    const esdData = (fitResults && fitResults.workingData) || workerWorkingData;

    if (JtJ && parameterInfo && ss_res !== undefined && esdData && esdData.intensity && esdData.intensity.length > 0) {
         const P = parameterInfo.length;
         const N = esdData.intensity.length;
         // FIX: use the EFFECTIVE degrees of freedom reported by the
         // worker, which include the one free intensity per reflection
         // that a Le Bail decomposition also extracts from this data.
         // N - parameterInfo.length counted only the profile parameters,
         // so chi-square was too small and every ESD scaled from it was
         // optimistic -- by a large factor on a pattern with many
         // reflections.
         const reportedDof = (fitResults && typeof fitResults.dof === 'number' && fitResults.dof > 0)
                             ? fitResults.dof : null;
         if (P === 0) {
              esdWarning = "No parameters were refined.";
         } else if (N > P) {
              const degreesOfFreedom = reportedDof !== null ? reportedDof : (N - P);

                if (ss_res >= 0 && isFinite(ss_res)) {
                    const reduced_chi_sq = (degreesOfFreedom > 0) ? (ss_res / degreesOfFreedom) : 0;
                    let cov_matrix = null;

                    try {
                         // Verify Matrix Dimensions
                         // Shape, not constructor. Array.isArray() is false
                         // for a Float64Array, so a perfectly valid matrix
                         // with typed rows failed this test and was reported
                         // as malformed.
                         if (isSquareMatrix(JtJ, P)) {
                             
                             const covResult = reportCovariance(fitResults);
                             if (!covResult) {
                                 throw new Error("Normal-equations matrix is not positive definite.");
                             }
                             // PERF FIX: the old line was
                             //     covResult.cov.map(row => row.map(v => v * reduced_chi_sq))
                             // which built a second full n x n copy -- 219 ms
                             // and ~50 MB at 2500 reflections -- so that the
                             // loop below could read 2500 diagonal entries.
                             // Scale on access instead.
                             cov_matrix = covResult.cov;
                         } else {
                             const rows = Array.isArray(JtJ) ? JtJ.length : 'not an array';
                             const cols = (JtJ && JtJ[0] && JtJ[0].length !== undefined) ? JtJ[0].length : '?';
                             throw new Error(`Normal-equations matrix is ${rows}x${cols}, expected ${P}x${P}.`);
                         }
                    } catch (invError) {
                         console.error("Matrix inversion error:", invError);
                         cov_matrix = null;
                         esdWarning = "Matrix singularity prevented error calculation (severe peak overlap).";
                    }

                    // If covariance matrix was successfully calculated, extract diagonal elements
                    if (cov_matrix) {
                        const undeterminedNames = [];
                        parameterInfo.forEach((p_info, i) => {
                             if (cov_matrix[i] && cov_matrix[i][i] !== undefined) {
                                 const variance = cov_matrix[i][i] * reduced_chi_sq;
                                  if (variance >= 0 && isFinite(variance)) {
                                      const sigma_scaled = Math.sqrt(variance);
                                      const scale = p_info.scale || 1.0;
                                      esds[p_info.name] = sigma_scaled * scale;
                                  } else {
                                      // Cannot happen for a well-conditioned
                                      // J^T J; means this parameter is not
                                      // determined by the data. Say so
                                      // rather than printing "NaN".
                                      undeterminedNames.push(p_info.name);
                                      esds[p_info.name] = NaN;
                                  }
                             } else {
                                 undeterminedNames.push(p_info.name);
                                 esds[p_info.name] = NaN;
                             }
                        });
                        if (undeterminedNames.length > 0) {
                            esdWarning = `Not determined by the data (no ESD): ${undeterminedNames.join(', ')}. `
                                       + `Usually means these parameters are fully correlated with others in the fit.`;
                        }
                    }
               } else {
                    esdWarning = "Invalid sum of squared residuals (ss_res).";
               }


         } else {
              esdWarning = `N (${N}) <= P (${P}), cannot calculate reliable errors.`;
         }
    } else if (!JtJ || !parameterInfo || ss_res === undefined || !esdData || !esdData.intensity || esdData.intensity.length === 0) {
         // Name the piece that is actually absent. "Required data is missing"
         // followed by a list of four candidates tells the user nothing they
         // can act on, and when the history snapshot lost its normal-equations
         // matrix it took a code read to work out which of the four it was.
         const absent = [];
         if (!JtJ) absent.push('normal-equations matrix (J^T J)');
         if (!parameterInfo) absent.push('parameterInfo');
         if (ss_res === undefined) absent.push('ss_res');
         if (!esdData || !esdData.intensity || esdData.intensity.length === 0) absent.push('the fitted data slice');
         esdWarning = `Cannot calculate ESDs: ${absent.join(', ')} not present in this result.`;
    }

    // A Le Bail ESD is a LOWER BOUND, and saying so is not optional.
    // The extracted intensities are counted in the degrees of freedom
    // (above), but they are not in the normal-equations matrix, so its
    // inverse contains no profile-to-intensity correlation at all -- and
    // the cell parameters are strongly correlated with the intensities
    // they help partition. A Pawley fit refines the intensities as
    // proper least-squares parameters and does not have this problem.
    if (refinementMode === 'le-bail' && Object.keys(esds).length > 0) {
        const nI = (fitResults && fitResults.nIntensities) || 0;
        esdWarning = (esdWarning ? esdWarning + ' ' : '')
            + `Le Bail ESDs are lower bounds.`;
    }

    // ------------------------------------------------------------------
    //  SYMMETRY-SLAVED PARAMETERS
    //
    //  In a hexagonal cell b is not refined -- it is set equal to a, by
    //  enforceSymmetryConstraintsWorker, every time the parameters change.
    //  The table reported that as "Fitted: No" with a blank ESD, which reads
    //  as "this number was assumed" when in fact it is determined exactly as
    //  well as a is. It IS a: same parameter, same random variable, so
    //  sigma(b) = sigma(a) identically, not merely approximately.
    //
    //  A slaved parameter is therefore shown as "= a" with a's uncertainty.
    //  A parameter fixed by symmetry to a CONSTANT -- gamma = 120 in the same
    //  cell -- is genuinely not determined by the data and stays "No" with no
    //  ESD. The two cases look the same in fitFlags and are not the same
    //  thing.
    //
    //  Mirrors enforceSymmetryConstraintsWorker exactly; if the constraints
    //  there change, change these with them.
    //
    //  @param {string} system
    //  @returns {Object<string,string>} slave parameter -> master parameter
    // ------------------------------------------------------------------
    const symmetrySlaves = (system) => {
        switch (system) {
            case 'cubic':
                return { b: 'a', c: 'a',
                         S040: 'S400', S004: 'S400', S202: 'S220', S022: 'S220' };
            case 'tetragonal': case 'hexagonal': case 'trigonal': case 'rhombohedral':
                return { b: 'a', S040: 'S400', S022: 'S202' };
            default:
                return {};
        }
    };
    const slavedTo = symmetrySlaves(
        (finalParams && finalParams.system) ||
        (fitResults && fitResults.system) ||
        (typeof currentSG !== 'undefined' && currentSG ? currentSG.system : null));

    const paramWidths = [24, 18, 10, 18]; // Name, Value, Fitted, ESD
    const paramHeader = formatLine(['Parameter', 'Value', 'Fitted', 'ESD'], paramWidths);
    const paramLines = [];
    const paramGroups = {};

    const addParam = (group, name, value, flag, esd_key) => {
         // Check if value is defined and not null before adding
         if (value !== undefined && value !== null) {
              group.push({ name: name, value: value, flag: flag, esd_key: esd_key });
         }
    };

    paramGroups['Structural & Instrumental'] = [];
    addParam(paramGroups['Structural & Instrumental'], 'a (Å)', finalParams.a, fitFlags.a, 'a');
    addParam(paramGroups['Structural & Instrumental'], 'b (Å)', finalParams.b, fitFlags.b, 'b');
    addParam(paramGroups['Structural & Instrumental'], 'c (Å)', finalParams.c, fitFlags.c, 'c');
    addParam(paramGroups['Structural & Instrumental'], 'alpha (°)', finalParams.alpha, fitFlags.alpha, 'alpha');
    addParam(paramGroups['Structural & Instrumental'], 'beta (°)', finalParams.beta, fitFlags.beta, 'beta');
    addParam(paramGroups['Structural & Instrumental'], 'gamma (°)', finalParams.gamma, fitFlags.gamma, 'gamma');
    addParam(paramGroups['Structural & Instrumental'], 'Radiation 1 (Å)', finalParams.lambda, undefined, undefined);
    addParam(paramGroups['Structural & Instrumental'], 'Radiation 2 (Å)', finalParams.lambda2, undefined, undefined);
    addParam(paramGroups['Structural & Instrumental'], 'Ratio (I2/I1)', finalParams.ratio, undefined, undefined);
    addParam(paramGroups['Structural & Instrumental'], 'Polarisation K', reportPol.K, undefined, undefined);
    addParam(paramGroups['Structural & Instrumental'], 'Zero Shift (°)', finalParams.zeroShift, fitFlags.zeroShift, 'zeroShift');
    addParam(paramGroups['Structural & Instrumental'], '2theta Min (°)', workerWorkingData ? workerWorkingData.tth[0] : null, undefined, undefined);
    addParam(paramGroups['Structural & Instrumental'], '2theta Max (°)', workerWorkingData ? workerWorkingData.tth[workerWorkingData.tth.length - 1] : null, undefined, undefined);

    paramGroups['Profile Parameters'] = [];
    switch (profileType) {
        case "simple_pvoigt":
            addParam(paramGroups['Profile Parameters'], 'GU', finalParams.GU, fitFlags.GU, 'GU');
            addParam(paramGroups['Profile Parameters'], 'GV', finalParams.GV, fitFlags.GV, 'GV');
            addParam(paramGroups['Profile Parameters'], 'GW', finalParams.GW, fitFlags.GW, 'GW');
            addParam(paramGroups['Profile Parameters'], 'GP', finalParams.GP, fitFlags.GP, 'GP');
            addParam(paramGroups['Profile Parameters'], 'LX', finalParams.LX, fitFlags.LX, 'LX');
            addParam(paramGroups['Profile Parameters'], 'eta (Mixing)', finalParams.eta, fitFlags.eta, 'eta');
            addParam(paramGroups['Profile Parameters'], 'shft (Displ.)', finalParams.shft, fitFlags.shft, 'shft');
            addParam(paramGroups['Profile Parameters'], 'trns (Transp.)', finalParams.trns, fitFlags.trns, 'trns');
            break;
        case "tch_aniso":
            addParam(paramGroups['Profile Parameters'], 'U', finalParams.U, fitFlags.U, 'U');
            addParam(paramGroups['Profile Parameters'], 'V', finalParams.V, fitFlags.V, 'V');
            addParam(paramGroups['Profile Parameters'], 'W', finalParams.W, fitFlags.W, 'W');
            addParam(paramGroups['Profile Parameters'], 'X', finalParams.X, fitFlags.X, 'X');
            addParam(paramGroups['Profile Parameters'], 'Y', finalParams.Y, fitFlags.Y, 'Y');
            addParam(paramGroups['Profile Parameters'], 'S/L (Asymm)', finalParams.SL, fitFlags.SL, 'SL');
            addParam(paramGroups['Profile Parameters'], 'H/L (Asymm)', finalParams.HL, fitFlags.HL, 'HL');
            addParam(paramGroups['Profile Parameters'], 'S400', finalParams.S400, fitFlags.S400, 'S400');
            addParam(paramGroups['Profile Parameters'], 'S040', finalParams.S040, fitFlags.S040, 'S040');
            addParam(paramGroups['Profile Parameters'], 'S004', finalParams.S004, fitFlags.S004, 'S004');
            addParam(paramGroups['Profile Parameters'], 'S220', finalParams.S220, fitFlags.S220, 'S220');
            addParam(paramGroups['Profile Parameters'], 'S202', finalParams.S202, fitFlags.S202, 'S202');
            addParam(paramGroups['Profile Parameters'], 'S022', finalParams.S022, fitFlags.S022, 'S022');
            break;
        case "split_pvoigt":
            addParam(paramGroups['Profile Parameters'], 'GU-L', finalParams.GU_L, fitFlags.GU_L, 'GU_L');
            addParam(paramGroups['Profile Parameters'], 'GV-L', finalParams.GV_L, fitFlags.GV_L, 'GV_L');
            addParam(paramGroups['Profile Parameters'], 'GW-L', finalParams.GW_L, fitFlags.GW_L, 'GW_L');
            addParam(paramGroups['Profile Parameters'], 'LX-L', finalParams.LX_L, fitFlags.LX_L, 'LX_L');
            addParam(paramGroups['Profile Parameters'], 'GU-R', finalParams.GU_R, fitFlags.GU_R, 'GU_R');
            addParam(paramGroups['Profile Parameters'], 'GV-R', finalParams.GV_R, fitFlags.GV_R, 'GV_R');
            addParam(paramGroups['Profile Parameters'], 'GW-R', finalParams.GW_R, fitFlags.GW_R, 'GW_R');
            addParam(paramGroups['Profile Parameters'], 'LX-R', finalParams.LX_R, fitFlags.LX_R, 'LX_R');
            addParam(paramGroups['Profile Parameters'], 'eta (Mixing)', finalParams.eta_split, fitFlags.eta_split, 'eta_split');
            addParam(paramGroups['Profile Parameters'], 'shft (Displ.)', finalParams.shft_split, fitFlags.shft_split, 'shft_split');
            addParam(paramGroups['Profile Parameters'], 'trns (Transp.)', finalParams.trns_split, fitFlags.trns_split, 'trns_split');
            break;
    }

    // Removed Background Parameters section (Chebyshev/Hump)

    // ------------------------------------------------------------------
    //  THE QUOTABLE CELL
    //
    //  One line, in the form a paper would print. The table below carries the
    //  same numbers, but spread over six rows with the su in a separate
    //  column, which is not what anyone copies. Slaved edges are shown with
    //  the master's su, and angles fixed by symmetry are given without one.
    // ------------------------------------------------------------------
    const cellSummary = (() => {
        const part = (key, label, unit) => {
            const v = finalParams[key];
            if (v === null || v === undefined || !isFinite(v)) return null;
            const master = slavedTo[key];
            const esdKey = (master && fitFlags[master]) ? master
                         : (fitFlags[key] ? key : null);
            const compact = esdKey ? formatWithEsd(v, esds[esdKey]) : null;
            return `${label} = ${compact !== null ? compact : fmtNum(v, 5)}${unit}`;
        };
        const edges = ['a', 'b', 'c'].map(k => part(k, k, ' \u00c5')).filter(Boolean);
        const angs = [['alpha', '\u03b1'], ['beta', '\u03b2'], ['gamma', '\u03b3']]
            .map(([k, g]) => part(k, g, '\u00b0')).filter(Boolean);
        if (edges.length === 0) return [];
        const out = ['--- Refined Cell ---', edges.join(',  ')];
        if (angs.length) out.push(angs.join(',  '));
        if (typeof finalParams.zeroShift === 'number' && isFinite(finalParams.zeroShift)) {
            const z = fitFlags.zeroShift ? formatWithEsd(finalParams.zeroShift, esds.zeroShift) : null;
            out.push(`zero shift = ${z !== null ? z : fmtNum(finalParams.zeroShift, 5)}\u00b0`);
        }
        // Whether the background contributes to these depends on whether it was
        // refined, so say which -- the statistics block above states it too,
        // but this block is the one that gets copied out on its own.
        out.push('Uncertainties are in units of the last digit shown, and come from the',
                 'diagonal of the covariance matrix scaled by the reduced chi-squared.',
                 fitFlags.fitBackground
                    ? 'They include background uncertainty: the anchor heights were refined.'
                    : 'They exclude background uncertainty: the anchors were held fixed.', '');
        return out;
    })();
    cellSummary.forEach(l => paramLines.push(l));

    for (const groupName in paramGroups) {
        // Ensure group exists and has items before adding headers
        if (paramGroups[groupName] && paramGroups[groupName].length > 0) {
            paramLines.push(`--- ${groupName} ---`, paramHeader, '-'.repeat(paramHeader.length));
            paramGroups[groupName].forEach(p => {
                 // Check value is valid number before formatting
                 if (p.value !== null && !isNaN(p.value)) {
                    // A slaved parameter borrows its master's row: '= a' in
                    // place of Yes/No, and the master's ESD, which is its own
                    // ESD exactly. Only when the master was actually refined
                    // -- with `a` fixed, b is fixed too and has no ESD.
                    const master = slavedTo[p.esd_key];
                    const masterFitted = master ? !!fitFlags[master] : false;

                    const fitStr = (p.flag === undefined) ? ''
                                 : master ? `= ${master}`
                                 : (p.flag ? 'Yes' : 'No');

                    const esdKey = (master && masterFitted) ? master
                                 : (p.flag ? p.esd_key : null);
                    const esdValue = esdKey ? esds[esdKey] : undefined;

                    let esdStr = '';
                    if (esdKey) {
                        if (typeof esdValue === 'number' && isFinite(esdValue)) {
                             esdStr = `(${esdValue.toExponential(2)})`;
                        } else if (isNaN(esdValue)) {
                             esdStr = '(NaN)';
                        }
                    }

                    // VALUE IN THE COMPACT CRYSTALLOGRAPHIC FORM where an su
                    // exists: 5.5139(3), rounded to the digit the su reaches.
                    // This is the string that gets copied into a paper, so it
                    // should be the one that is defensible. The raw su stays
                    // in its own column for anyone who needs it arithmetically.
                    //
                    // Without an su, a plain decimal -- toExponential(6) wrote
                    // a right angle as 9.000000e+1. Exponential is kept only
                    // for magnitudes where a decimal would be unreadable.
                    let valStr;
                    if (typeof p.value !== 'number') {
                        valStr = String(p.value);
                    } else {
                        const compact = formatWithEsd(p.value, esdValue);
                        if (compact !== null) {
                            valStr = compact;
                        } else {
                            const m = Math.abs(p.value);
                            valStr = (m !== 0 && (m >= 1e6 || m < 1e-4))
                                ? p.value.toExponential(4)
                                : fmtNum(p.value, 6);
                        }
                    }
                    paramLines.push(formatLine([p.name, valStr, fitStr, esdStr], paramWidths));
                }
            });
            paramLines.push(''); // Add blank line after each group
        }
    }
    if (esdWarning) {
         paramLines.push(`NOTE regarding ESDs: ${esdWarning}`, ''); // 
    }


    //   ADDED: Background Points Section, version 115, 25 oct 2025  
    const backgroundPointsSection = [];
    if (backgroundAnchors && backgroundAnchors.length > 0) {
         const bgWidths = [15, 18, 18]; // 2theta, Intensity, ESD
         const isBgRefined = fitFlags && fitFlags.fitBackground;
         const bgHeader = formatLine(['2theta (°)', 'Intensity', 'ESD'], bgWidths);
         backgroundPointsSection.push('', '--- Background Spline Points ---', bgHeader, '-'.repeat(bgHeader.length));
         // Use a copy sorted by tth for the report
         const sortedAnchors = [...backgroundAnchors].sort((a, b) => a.tth - b.tth);
         sortedAnchors.forEach((point, i) => {
              let esdStr = '-';
              if (isBgRefined && esds[`bg_y_${i}`] !== undefined) {
                  const esdValue = esds[`bg_y_${i}`];
                  if (typeof esdValue === 'number' && isFinite(esdValue)) {
                      esdStr = `(${esdValue.toExponential(2)})`;
                  } else if (isNaN(esdValue)) {
                      esdStr = '(NaN)';
                  }
              }
              const yVal = (isBgRefined && finalParams[`bg_y_${i}`] !== undefined) ? finalParams[`bg_y_${i}`] : point.y;
              backgroundPointsSection.push(formatLine([
                  point.tth.toFixed(4),
                  yVal.toFixed(2),
                  esdStr
              ], bgWidths));
         });
         backgroundPointsSection.push(''); // Add blank line after section
    }
    

    const reflectionsSection = [];
    // PERF FIX: this block was built unconditionally, so the SUMMARY preview
    // that rpRenderRun draws after every fit walked every reflection, ran
    // calculateAllObservedIntensities and reflectionStructureFactors over the
    // whole list, and asked for the covariance -- all on the main thread, in
    // the same task as the history snapshot. The summary does not display any
    // of it. Anything that is actually read by a human or written to a file
    // still asks for 'full'.
    const wantReflections = (format !== 'summary');
    if (wantReflections && hklList && hklList.length > 0 && workerWorkingData && workerWorkingData.tth) {
        const isPawley = (refinementMode === 'pawley');

        // The doublet: the integrated area below is the sum over the
        // Ka1/Ka2 pair, so it is larger than the single-wavelength
        // intensity by (1 + I2/I1). |Fo| divides that back out, which
        // makes the printed structure factors directly comparable
        // with a monochromatic (or synchrotron) measurement instead
        // of carrying a source-dependent constant.
        const dbl = (finalParams.ratio > 1e-6 && finalParams.lambda2 > 1e-6 &&
                     Math.abs(finalParams.lambda - finalParams.lambda2) > 1e-6);
        const doubletSum = dbl ? (1 + finalParams.ratio) : 1;

        // THE SAME FUNCTION THE pdCIF WRITER USES, keyed by the first hkl
        // string so the loop below can look each reflection up.
        //
        // The arithmetic was inline here and absent from writePdCIF(), which
        // wrote the raw refined height into _refln_F_squared_meas instead --
        // 163x to 1562x wrong depending on the reflection. Having one
        // implementation is the only thing that keeps the report, the PDF and
        // the CIF telling the same story, since all three now come through
        // here.
        const sfByKey = new Map();
        if (typeof reflectionStructureFactors === 'function') {
            const sigmaArea = (hkl, widthFactor) => {
                if (isPawley) {
                    const sh = esds[`I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`];
                    return Number.isFinite(sh) ? sh * widthFactor : NaN;
                }
                // I_sigma_full, not I_sigma: the area this is paired with is
                // the untruncated one, so the uncertainty has to be too. See
                // the same note in data_io.js.
                if (Number.isFinite(hkl.I_sigma_full)) return hkl.I_sigma_full * mainScaleFactor;
                return Number.isFinite(hkl.I_sigma) ? hkl.I_sigma * mainScaleFactor : NaN;
            };
            for (const r of reflectionStructureFactors(hklList, finalParams,
                    { scale: mainScaleFactor, polarisation: reportPol, sigmaArea })) {
                if (r.hkl && r.hkl.hkl_list) sfByKey.set(r.hkl.hkl_list[0], r);
            }
        }

        // Le Bail has ONE intensity per reflection: the decomposed
        // observed intensity. Pawley additionally has the refined
        // parameter, which is worth comparing against the observed.
        //
        // Both now report the full set required to get from a peak to
        // a structure factor: the integrated intensity of the whole
        // reflection, its multiplicity, its Lorentz-polarisation
        // factor, and the resulting |Fo|.
        const reflCols = ['h,k,l', '2th_corr (°)', 'd (Å)', 'm', 'Lp',
                          'I_hkl (Area)', 'sigma(I)', '|Fo|'];
        const reflWidths = [13, 12, 9, 5, 12, 14, 12, 12];
        if (isPawley) { reflCols.push('I_obs (Area)'); reflWidths.push(14); }
        const reflHeader = formatLine(reflCols, reflWidths);

        reflectionsSection.push(
            '',
            '--- Reflections List (Integrated Intensities) ---',
            `Lp model : ${polarizationDescription(reportPol)}`,
            ...polarizationFormulaLines(reportPol).map(s => '  ' + s),
            `I_hkl    : integrated intensity of the reflection${dbl ? ' (Ka1 + Ka2)' : ''}, scale factor applied`,
            'm        : powder multiplicity of the reflection',
            `|Fo|     : sqrt( I_hkl / ( m * Lp${dbl ? ` * ${doubletSum.toFixed(4)}` : ''} ) )    [arbitrary scale]`,
            ...(isPawley
                ? ['I_obs    : same reflection re-integrated from the observed pattern (cross-check)']
                : []),
            '',
            reflHeader, '-'.repeat(reflHeader.length));

        let i_obs_map = new Map();
        if (isPawley) {
            try {
                 i_obs_map = calculateAllObservedIntensities(hklList, finalParams, workerWorkingData).area;
            } catch(iobsError) {
                 console.error("Error calculating observed intensities:", iobsError);
                 reflectionsSection.push("(Error calculating I_obs)");
            }
        }

        const fitted_tth_min = workerWorkingData.tth[0];
        const fitted_tth_max = workerWorkingData.tth[workerWorkingData.tth.length - 1];

        hklList.filter(hkl => hkl && Number.isFinite(hkl.tth) && hkl.tth >= fitted_tth_min && hkl.tth <= fitted_tth_max)
            .forEach(hkl => {
                 let tthCorr = hkl.tth;
                 if (finalParams.zeroShift) tthCorr += finalParams.zeroShift;
                 try {
                     tthCorr += calculatePeakShift(hkl.tth, finalParams);
                 } catch { /* ignore shift error */ }
                const i_calc_height = hkl.intensity || 0;

                // The integrated intensity is whatever the Le Bail
                // decomposition actually produced, using the integral of
                // the profile as it was really drawn (hkl.shapeArea,
                // accumulated on the data grid during extraction).
                // Reconstructing it here from the analytic pseudo-Voigt
                // area gave a different number, because the drawn peak
                // is truncated and pedestal-subtracted. Fall back to the
                // analytic route only for reflections the extraction
                // never visited (e.g. a Pawley run).
                const total_calc_area = integratedPeakArea(hkl, finalParams, mainScaleFactor);

                // 5. Calculate Effective Width Factor for ESD Scaling
                // The matrix gives ESD for Height. We need ESD for Area.
                // Relation: Area = Height * EffectiveWidth. Therefore, sigma(Area) = sigma(Height) * EffectiveWidth.
                // EffectiveWidth = TotalArea / Height.
                const effective_width_factor = (i_calc_height > 1e-9) ? (total_calc_area / i_calc_height) : 0;

                const i_obs_raw = i_obs_map.get(hkl.hkl_list[0]);
                const i_obs_known = typeof i_obs_raw === 'number' && isFinite(i_obs_raw);
                const i_obs = i_obs_known ? i_obs_raw : 0;
                const i_obs_str = i_obs_known ? i_obs_raw.toFixed(1) : 'n/a';
                let esd_I_str = '';

                if (isPawley) {
                     // Pawley intensities are least-squares parameters,
                     // so their ESD comes from the covariance matrix.
                     const hkl_name = `I_(${hkl.h_orig},${hkl.k_orig},${hkl.l_orig})`;
                     const esd_I_val_height = esds[hkl_name];

                     if (typeof esd_I_val_height === 'number' && isFinite(esd_I_val_height)) {
                         // Convert Height ESD to Area ESD
                         const esd_I_val_area = esd_I_val_height * effective_width_factor;
                         esd_I_str = esd_I_val_area.toFixed(1);
                     } else {
                         esd_I_str = '-';
                     }
                } else {
                     // FIX: this column used to be blank for every Le Bail
                     // run. Le Bail intensities are not least-squares
                     // parameters so they have no covariance entry, but
                     // they DO have a well-defined uncertainty: the
                     // counting statistics of the observed points,
                     // propagated through the partition fractions. That
                     // is computed during the extraction (peak.I_sigma)
                     // and is what other Rietveld packages report here.
                     // I_sigma_full so this column matches the I_hkl column
                     // beside it, which is built from the untruncated area.
                     const sig = Number.isFinite(hkl.I_sigma_full) ? hkl.I_sigma_full
                                                                   : hkl.I_sigma;
                     const scaled = (typeof sig === 'number') ? sig * mainScaleFactor : NaN;
                     esd_I_str = isFinite(scaled) ? scaled.toFixed(1) : '-';
                }

                //   Lorentz-polarisation and the structure factor.
                // Lp is evaluated at the CORRECTED 2theta -- the same
                // angle printed in the table -- because that is where
                // the peak actually sits; using the ideal Bragg angle
                // would put Lp and the reported position out of step.
                const lpVal = lorentzPolarization(tthCorr, reportPol);
                const mult  = (Number.isFinite(hkl.multiplicity) && hkl.multiplicity > 0)
                              ? hkl.multiplicity : 1;
                const lp_str = Number.isFinite(lpVal) ? lpVal.toFixed(4) : 'n/a';

                // |Fo| exists only where the intensity is positive. A
                // negative extracted area is a real measurement (the
                // background ran above the data there), so it is
                // printed as such in I_hkl and simply has no square
                // root -- inventing a zero would hide it.
                let fo_str = '-';
                const sfRow = sfByKey.get(hkl.hkl_list[0]);
                const fsq = sfRow ? sfRow.Fsq
                                  : (Number.isFinite(lpVal) && lpVal > 0
                                     ? total_calc_area / (mult * lpVal * doubletSum) : NaN);
                if (Number.isFinite(fsq) && fsq > 0) fo_str = Math.sqrt(fsq).toFixed(3);

                if (Math.abs(i_obs) > 0.01 || total_calc_area > 0.01) {
                    const row = [
                        hkl.hkl_list[0],
                        tthCorr.toFixed(4),
                        Number.isFinite(hkl.d) ? hkl.d.toFixed(4) : '-',
                        String(mult),
                        lp_str,
                        total_calc_area.toFixed(1),
                        esd_I_str,
                        fo_str
                    ];
                    if (isPawley) row.push(i_obs_str);
                    reflectionsSection.push(formatLine(row, reflWidths));
                }
            });

        // ------------------------------------------------------------------
        //  EXACT COINCIDENCES.
        //
        //  Distinct from the overlap clusters below. An overlap cluster is a
        //  set of reflections the data separates BADLY; this is a set it
        //  cannot separate AT ALL, because their profiles are identical --
        //  333 against 511, 700 against 522 and their relatives in a cubic F
        //  cell, and every accidental coincidence elsewhere.
        //
        //  The linear solve puts the group's whole intensity on one member
        //  and leaves the rest near zero. That split is arbitrary; the SUM is
        //  exact. Without this section the arbitrary number looks ordinary:
        //  on a test pattern 775 was printed as 2340 against a true 406,
        //  because it had absorbed 11-1-1, and nothing said so.
        // ------------------------------------------------------------------
        const degGroups = (fitResults && fitResults.degenerateGroups) || [];
        if (isPawley && degGroups.length) {
            reflectionsSection.push('',
                '--- Exactly Coincident Reflections (not separable at all) ---',
                'These reflections have identical calculated profiles, so only the SUM',
                'within each group is determined by the data. The split printed in the',
                'table above is arbitrary and carries no information.');
            const dgW = [34, 16, 16];
            const dgHeader = formatLine(['group', 'I_group (Area)', '2theta (°)'], dgW);
            reflectionsSection.push(dgHeader, '-'.repeat(dgHeader.length));
            for (const grp of degGroups) {
                const names = grp.map(j => (hklList[j] && hklList[j].hkl_list)
                                            ? hklList[j].hkl_list[0] : `#${j}`);
                let sum = 0, tth0 = NaN;
                for (const j of grp) {
                    const h = hklList[j];
                    if (!h) continue;
                    if (!Number.isFinite(tth0) && Number.isFinite(h.tth)) tth0 = h.tth;
                    try { sum += integratedPeakArea(h, finalParams, mainScaleFactor); }
                    catch { /* leave the sum as far as it got */ }
                }
                reflectionsSection.push(formatLine(
                    [names.join(' = '), sum.toFixed(1),
                     Number.isFinite(tth0) ? tth0.toFixed(4) : '-'], dgW));
            }
        }

        // ------------------------------------------------------------------
        //  OVERLAP CLUSTERS
        //
        //  A Pawley intensity with an ESD larger than itself is not a failed
        //  fit, and printing it without context invites exactly that reading.
        //  It means the reflection overlaps its neighbours so closely that the
        //  data cannot say how the intensity divides between them -- while
        //  saying perfectly clearly what the TOTAL is. This section prints the
        //  total, so the number the measurement actually determined is on the
        //  page next to the ones it does not.
        //
        //  Only clusters containing at least one unresolved member are listed;
        //  a well-separated reflection is its own cluster and its sigma in the
        //  table above is already the whole story.
        // ------------------------------------------------------------------
        if (isPawley && typeof intensityOverlapClusters === 'function') {
            const oc = intensityOverlapClusters(fitResults, finalParams, reportCovariance(fitResults));
            const bad = oc ? oc.clusters.filter(c => c.members.length > 1 && c.nUnresolved > 0) : [];
            if (bad.length) {
                const ocCols = ['2th range (°)', 'n', 'unresolved', 'I_cluster', 'sigma(I_cluster)', 'members'];
                const ocW = [16, 5, 12, 14, 18, 40];
                const ocHeader = formatLine(ocCols, ocW);
                reflectionsSection.push(
                    '',
                    '--- Overlap Clusters (reflections the data cannot separate) ---',
                    'Two reflections closer than their profile width have nearly parallel',
                    'derivatives: intensity can move from one to the other with almost no',
                    'change to the calculated pattern. Their individual intensities are then',
                    'poorly determined even where the fit is visually excellent, while their',
                    'SUM is well determined. sigma(I_cluster) uses the full covariance,',
                    'including the large negative off-diagonal terms -- it is NOT the',
                    'quadrature sum of the individual sigmas, and is usually far smaller.',
                    'Ka2 satellites are included when deciding what overlaps.',
                    '',
                    'Use I_cluster, not the individual intensities, for anything downstream.',
                    '',
                    ocHeader, '-'.repeat(ocHeader.length));

                bad.forEach(c => {
                    const names = c.members.map(m => `(${m.h},${m.k},${m.l})${m.separable ? '' : '*'}`).join(' ');
                    reflectionsSection.push(formatLine([
                        `${c.tthMin.toFixed(3)}-${c.tthMax.toFixed(3)}`,
                        String(c.members.length),
                        String(c.nUnresolved),
                        Number.isFinite(c.sumI) ? c.sumI.toFixed(1) : '-',
                        Number.isFinite(c.sigmaSum) ? c.sigmaSum.toFixed(1) : '-',
                        names
                    ], ocW));
                });
                reflectionsSection.push('', '* = sigma exceeds the intensity: not separately determined.');
            }
        }
    }

    const dataSection = [];
    if (format === 'full' && workerWorkingData && workerWorkingData.tth) {
         dataSection.push('', '--- Point-by-Point Intensity Data (Fitted Range Only) ---',
             formatLine(['2theta', 'I_obs', 'I_calc', 'Difference', 'Background'], [15, 18, 18, 18, 18]),
             '-'.repeat(15+18+18+18+18+4));

         try {
              const finalBkg = calculateTotalBackground(workerWorkingData.tth, finalParams, backgroundAnchors); // Pass anchors
              const finalNetPattern = calculatePattern(workerWorkingData.tth, hklList, finalParams);
              const finalCalcPattern = new Float64Array(finalNetPattern.length);
              for (let i = 0; i < finalCalcPattern.length; i++) {
                  finalCalcPattern[i] = finalNetPattern[i] * mainScaleFactor + finalBkg[i];
              }

              for (let i = 0; i < workerWorkingData.tth.length; i++) {
                  const i_obs = workerWorkingData.intensity[i];
                  const i_calc = finalCalcPattern[i];
                  const i_bkg = finalBkg[i];
                   if (isFinite(i_obs) && isFinite(i_calc) && isFinite(i_bkg)) {
                      dataSection.push(formatLine([
                          workerWorkingData.tth[i].toFixed(4),
                          i_obs.toFixed(2),
                          i_calc.toFixed(2),
                          (i_obs - i_calc).toFixed(2),
                          i_bkg.toFixed(2)
                      ], [15, 18, 18, 18, 18]));
                   }
              }
         } catch (dataError) {
               console.error("Error generating point data for report:", dataError);
               dataSection.push("(Error generating point data)");
         }
    }

    // Assemble the report sections
    return [
        ...header,
        ...statsSection,
        ...lpSection,
        ...whAnalysisSection,
        ...paramLines,
        ...backgroundPointsSection, 
        ...reflectionsSection,
        ...dataSection
    ].join('\n');
}

/**
 * Text report for a Wyckoff solution.
 *
 * Deliberately not the charge-flipping report with fields blanked out. The two
 * describe different things: a map has an origin, a hand and a symmetry
 * correlation; a Wyckoff solution has an assignment, free-parameter count and
 * an R factor against the intensities. Printing one through the other's
 * template is how a reader ends up trusting a heading that does not apply.
 *
 * @param {object} st  The structure, as posted by the wyckoff-search worker.
 * @returns {string}
 */
function generateWyckoffReport(st) {
    const sg = st.spaceGroup || {};
    const c = st.cell || {};
    const wy = st.wyckoffResult || {};
    const ref = st.refinement || wy.refinement || {};
    const num = (v, n) => Number.isFinite(v) ? v.toFixed(n) : '-';
    const pct = (v) => Number.isFinite(v) ? (v * 100).toFixed(2) + '%' : '-';
    // Its own copy: the formatLine in generateReportContent is a local const,
    // not a module function, so referencing it from here would have been a
    // ReferenceError the first time anyone pressed Generate PDF.
    const formatLine = (cols, widths) =>
        cols.map((c, i) => String(c ?? '').padEnd(widths[i])).join('').trimEnd();

    const lines = [
        'Wyckoff Search Report',
        '='.repeat(72),
        APP_VERSION_TEXT,
        `Report Generated: ${new Date().toLocaleString()}`,
        `Data File: ${reportDataFileName()}`,
        '',
        '--- Model ---',
        `Space group            : ${sg.symbol || '-'}${sg.number ? ` (No. ${sg.number})` : ''}`,
        `Operators applied      : ${st.nOps || '-'}`,
        `Assignment             : ${wy.assignment || '-'}`,
        `Formula units Z        : ${wy.z || '-'}`,
        '',
        '--- Unit Cell ---',
        `a = ${num(c.a, 4)} A    b = ${num(c.b, 4)} A    c = ${num(c.c, 4)} A`,
        `alpha = ${num(c.alpha, 3)}    beta = ${num(c.beta, 3)}    gamma = ${num(c.gamma, 3)}`,
        '',
        '--- Agreement ---'
    ];

    // The search's own figure of merit, alongside the refinement's. They
    // optimise different quantities -- a correlation over reflection groups
    // versus a weighted residual -- so reporting only one hides whichever step
    // actually did the work.
    if (Number.isFinite(wy.searchCC)) {
        lines.push(
            `Search CC (full res.)  : ${wy.searchCC.toFixed(4)}`,
            `   The figure shown while a search runs is computed on a`,
            `   resolution-ramped SUBSET of the reflections, so it reads`,
            `   higher than this and is not comparable with wR.`);
    }

    if (Number.isFinite(ref.R)) {
        lines.push(
            `wR(F^2)                : ${pct(ref.R)}` +
                (Number.isFinite(ref.Rstart) ? `  (from ${pct(ref.Rstart)} before refinement)` : ''),
            `Observations           : ${ref.nObs ?? '-'}`,
            `Free coordinates       : ${ref.nParams ?? '-'} (after the Wyckoff constraints)`,
            `Converged              : ${ref.converged === false ? 'NO' : 'yes'}`,
            `Lp removed             : ${ref.lpApplied ? 'yes' : 'no'}`);
    } else {
        lines.push(`wR(F^2)               : not available` +
                   (ref.skipped ? ` -- ${ref.skipped}` : ''));
    }

    lines.push('', '--- Atom Positions ---');
    const w = [4, 10, 12, 12, 12, 14];
    const head = formatLine(['#', 'Element', 'x', 'y', 'z', 'Wyckoff/Mult'], w);
    lines.push(head, '-'.repeat(head.length));
    (st.sites || []).forEach(x => lines.push(formatLine([
        x.rank, x.element || '?', num(x.x, 4), num(x.y, 4), num(x.z, 4),
        x.wyckoff || (x.multiplicity + (x.special ? '*' : ''))
    ], w)));

    if (Array.isArray(wy.searchSites) && wy.searchSites.length === (st.sites || []).length) {
        const movedAny = wy.searchSites.some((x, i) =>
            Math.abs(x.x - st.sites[i].x) > 5e-5 ||
            Math.abs(x.y - st.sites[i].y) > 5e-5 ||
            Math.abs(x.z - st.sites[i].z) > 5e-5);
        if (movedAny) {
            lines.push('', '--- Search Positions, before refinement ---');
            const h2 = formatLine(['#', 'Element', 'x', 'y', 'z', 'Wyckoff/Mult'], w);
            lines.push(h2, '-'.repeat(h2.length));
            wy.searchSites.forEach((x, i) => lines.push(formatLine([
                st.sites[i].rank, x.element || '?', num(x.x, 4), num(x.y, 4), num(x.z, 4),
                st.sites[i].wyckoff || ''
            ], w)));
        }
    }

    // The restraints the search ran under, spelled out.
    //
    // Without this the report states an R factor and a set of coordinates and
    // leaves out what the coordinates were REQUIRED to satisfy -- and a
    // structure produced under a coordination demand is a different claim from
    // one produced without it. describeDistanceConstraint lives in powder5.html
    // and is the same formatter the Results panel uses, so the report and the
    // screen cannot describe the same rule two different ways.
    const cons = Array.isArray(st.distanceConstraints) ? st.distanceConstraints : [];
    lines.push('', '--- Distance Constraints ---');
    if (cons.length && typeof describeDistanceConstraint === 'function') {
        cons.forEach(w => lines.push('  ' + describeDistanceConstraint(w)));
        lines.push('',
            '  A constraint naming a coordination number is enforced on that many',
            '  nearest partners of the named centre, and only outward FROM the first',
            '  element: "S O 4" asks each S for four O, not each O for four S.',
            '  Partners are counted per symmetry image AND per lattice translation, so',
            '  an atom bonded to two translates of one neighbour counts both.');
    } else {
        lines.push('  None. Only the minimum-contact floor was applied.');
    }

    lines.push('',
        'Coordinates were fitted to the Pawley intensities and tabulated scattering',
        'factors, weighted on 1/sigma^2 from the Pawley decomposition rather than on',
        'counting statistics: the uncertainty on one of these intensities comes from',
        'the decomposition, and for an overlapped reflection it is far larger than',
        'sqrt(I). wR is the residual actually minimised.');

    // The quality table goes LAST, after the positions it judges.
    wyQualitySection(st, formatLine).forEach(l => lines.push(l));

    return lines.join('\n');
}

/**
 * Coordination, contact statistics and bond-valence sums per site.
 *
 * WHY THIS BELONGS IN THE REPORT. The refined coordinates and the wR say the
 * fit converged; they say nothing about whether the result is chemistry. A
 * sulfur with a mean S-O of 1.13 A and a valence sum of 8.6 where 6 is
 * expected is a converged fit to a wrong structure, and the numbers that show
 * it are these. The panel has shown them all along -- the report did not, so a
 * PDF handed to someone else carried the positions without the one table that
 * says whether to believe them.
 *
 * Computed by WyckoffReport.analyse(), the same call the panel makes, so the
 * two cannot drift apart.
 *
 * @param {object} st Structure, as passed to generateWyckoffReport.
 * @param {(cols:Array, widths:number[])=>string} formatLine
 * @returns {string[]}
 */
function wyQualitySection(st, formatLine) {
    if (typeof WyckoffReport === 'undefined' || typeof WyckoffReport.analyse !== 'function') return [];
    let A;
    try { A = WyckoffReport.analyse(st, {}); }
    catch (e) { return ['', '--- Structure quality ---', `(not available: ${e.message || e})`]; }
    if (!A || A.error) {
        return ['', '--- Structure quality ---', '(' + ((A && A.error) || 'not available') + ')'];
    }
    if (!Array.isArray(A.report) || A.report.length === 0) return [];

    const w = [14, 24, 11, 10, 24, 12];
    const head = formatLine(['Site', 'CN', 'Mean (A)', 'Spread', 'Bond valence', 'Verdict'], w);
    const out = ['', '--- Structure quality ---', head, '-'.repeat(head.length)];
    const f = (v, n) => Number.isFinite(v) ? v.toFixed(n) : '-';

    for (const r of A.report) {
        // The composition is shown beside CN only when the shell is MIXED.
        // "2 (1xS + 1xPb)" is the whole point for an oxygen bridging two
        // cations; "4 (4xO)" would just be noise.
        const cn = (r.shellComp && r.shellComp.indexOf('+') >= 0)
            ? `${r.cn} (${r.shellComp.replace(/\u00d7/g, 'x')})`
            : String(r.cn);
        const bv = Number.isFinite(r.bvs)
            ? `${r.bvs.toFixed(2)} of ${r.expected} expected`
            : (r.missingParams && r.missingParams.length
                ? `no R0 for ${r.missingParams[0]}`
                : (r.cn ? '-' : 'no coordination shell'));
        out.push(formatLine([
            `${r.element} ${r.wyckoff}`, cn, f(r.mean, 3), f(r.spread, 3), bv, r.verdict
        ], w));
    }

    out.push('',
        `Bond-valence parameters are Brese & O'Keeffe (1991), b = 0.37 A. A sum within`,
        '15% of the formal valence is reported as plausible, within 30% as borderline.',
        `Both the sum and the coordination number depend on the ${A.cutoff} A cutoff:`,
        'distant contacts contribute little to the valence but do raise CN.',
        'Oxidation states are the tabulated default for each element, so the expected',
        'valence is wrong for a mixed-valence compound and the verdict with it.');
    return out;
}

function generateChargeFlippingReport(cfResult) {
    const now = new Date();
    const header = [
        `Charge Flipping Report`,
        '========================================================================',
        '',
        typeof APP_VERSION_TEXT !== 'undefined' ? APP_VERSION_TEXT : 'Powder 5',
        `Report Generated: ${now.toLocaleString()}`,
        `Data File: ${reportDataFileName()}`,
        ''
    ];

    // Try to recover full run context (opts, structure) from the history
    let cfOpts = {};
    let cfSt = null;
    if (typeof rpHistory !== 'undefined' && rpHistory.cf) {
        const entry = rpHistory.cf.find(e => e.result === cfResult || (e.results && e.results === cfResult));
        if (entry) {
            cfOpts = entry.opts || {};
            if (entry.structure) cfSt = entry.structure.results;
        }
    }
    // Fallback to live globals if not in history
    if (!cfSt && typeof cfStructure !== 'undefined' && cfStructure) {
        cfSt = cfStructure;
    }

    const prov = (cfResult && cfResult.intensityProvenance)
              || cfOpts.intensityProvenance
              || (typeof cfIntensityProvenance !== 'undefined' && cfIntensityProvenance)
              || {};

    const formatLine = (cols, widths) => {
        return cols.map((col, i) => {
             const colStr = (col === null || col === undefined) ? '-' : String(col);
             if (typeof col === 'number' && i > 0 && widths[i] > 0) {
                 return colStr.padStart(widths[i]);
             }
             return colStr.padEnd(widths[i]);
        }).join(' ');
    };

    const num = (v, d = 4) => (v === undefined || v === null || Number.isNaN(v)) ? '–' : Number(v).toFixed(d);

    const ref = cfResult.reflections || {};
    const c = cfResult.cell || {};
    
    const summarySection = [
        '--- Solution Summary ---',
        `Best R (|F| agreement) : ${num(cfResult.bestR, 4)}`,
        `Reached at iteration   : ${cfResult.bestIter} ${cfOpts.maxIter ? 'of ' + cfOpts.maxIter : ''}`,
        `Grid Size              : ${cfResult.gridSize}³`,
        `Cell Volume            : ${num(cfResult.volume, 2)} Å³`,
        `Unique reflections     : ${ref.unique || '–'}`,
        `Grid points filled     : ${ref.gridPoints || '–'}`,
        `Expansion symmetry     : ${ref.symmetrySource === 'symops' ? (ref.nSymops + ' operator(s), Laue ' + (ref.laueClass || '?')) : ('Laue ' + (ref.laueClass || 'P1'))}`,
        `Symmetry in loop       : ${!cfResult.symLambda ? 'none (solved in P1)' : ('lambda = ' + cfResult.symLambda.toFixed(2))}`,
        // Three sources, in order of reliability: the result object itself, the
        // run-history options, then the live global. The history lookup is by
        // object identity and CAN miss, and when it did this line rendered as
        // a dash -- which reads as "nothing was done to the intensities" when
        // in fact every weak reflection had been replaced by its posterior
        // mean. A silent dash is worse than a wrong value here.
        `Lp model               : ${prov.lpModel || cfOpts.lpModel || '–'}`,
        `French-Wilson          : ${prov.fwApplied === true
              ? `applied to ${prov.fwCorrected} weak reflection(s)`
                + (prov.fwNegatives ? `, ${prov.fwNegatives} of them negative` : '')
              : (prov.fwApplied === false
                    ? `NOT APPLIED (${prov.fwReason || 'reason not recorded'})`
                    : (cfOpts.intensityNote || 'not recorded'))}`,
        ...(prov.clampedToZero
            ? [`Clamped to zero        : ${prov.clampedToZero} non-positive intensity/intensities`]
            : []),
        `Compute backend        : ${cfResult.backend === 'gpu' ? 'WebGPU' : 'CPU'}`,
        ...(ref.multiplicityMismatch
            ? [`Multiplicity WARNING   : ${ref.multiplicityMismatch} reflection(s) disagree with the`,
               `                         space-group operators, e.g. ${ref.multiplicityExample}.`,
               '                         Every |F| in the map is scaled by sqrt of that ratio.']
            : []),
        ''
    ];

    const cellSection = [
        '--- Unit Cell ---',
        `a = ${num(c.a, 4)} Å    b = ${num(c.b, 4)} Å    c = ${num(c.c, 4)} Å`,
        `alpha = ${num(c.alpha, 3)}°  beta = ${num(c.beta, 3)}°  gamma = ${num(c.gamma, 3)}°`,
        ''
    ];

    let structureSection = [];
    let sitesSection = [];

    if (cfSt) {
        const sg = cfSt.spaceGroup || {};
        structureSection = [
            '--- Structure ---',
            `Space Group          : ${sg.symbol || sg.standard_symbol || 'Unknown'} (No. ${sg.number || '?'})`,
            `Operators applied    : ${cfSt.nOps || '–'}`,
            `Origin shift         : (${(cfSt.originShift || [0,0,0]).map(v => v.toFixed(4)).join(', ')})`,
            `Solution hand        : ${cfSt.hand > 0 ? 'as solved' : 'inverted'}`,
            `Symmetry correlation : ${num(cfSt.symmetryCorrelation, 3)}`,
            `Phase agreement      : ${num(cfSt.originScore, 3)}`,
            `Peaks -> unique sites: ${cfSt.rawPeakCount} -> ${(cfSt.sites || []).length}`,
            ''
        ];

        const wyckoff = cfSt.wyckoffResult;
        if (wyckoff && wyckoff.candidates && wyckoff.candidates.length > 0) {
            structureSection.push('--- Wyckoff Target Composition ---');
            structureSection.push(`Formulas per cell (Z): ${wyckoff.z || '?'}`);
            
            if (wyckoff.refinement) {
                const wr = wyckoff.refinement;
                if (wr.skipped) {
                    structureSection.push(`Pawley Refinement    : Skipped (${wr.skipped})`);
                } else {
                    structureSection.push(`Pawley Refinement    : Converged in ${wr.iterations} iterations`);
                    structureSection.push(`Initial R-factor     : ${num(wr.Rstart, 4)}`);
                    structureSection.push(`Final R-factor       : ${num(wr.R, 4)}`);
                }
            }
            structureSection.push('');
        }

        if (cfSt.sites && cfSt.sites.length > 0) {
            // 'Shortest Contact' removed: nothing on the charge-flipping path
            // populates s.shortestContact, so the column was a header over a
            // full column of en-dashes. Contacts are computed for the Wyckoff
            // path and reported there, where the numbers actually exist.
            const siteWidths = [4, 8, 10, 10, 10, 14, 10, 10];
            const siteHeader = formatLine(['#', 'Element', 'x', 'y', 'z', 'Wyckoff/Mult', 'Height', 'Charge'], siteWidths);

            sitesSection = [
                '--- Atom Positions ---',
                'Height : the map value at the peak, relative to the tallest peak.',
                'Charge : the map INTEGRATED over a 0.65 A sphere, local baseline removed,',
                '         relative to the largest integral. This is the column that tracks',
                '         atomic number; Height does not. Both are on an arbitrary scale --',
                '         F(000) is held at zero during charge flipping, so only RATIOS mean',
                '         anything. Ratios are still compressed for light atoms, and a site',
                '         with partial occupancy is indistinguishable from a lighter atom.',
                '         Bond distances discriminate light elements better than either.',
                '',
                siteHeader,
                '-'.repeat(siteHeader.length)
            ];

            cfSt.sites.forEach(s => {
                const wm = s.wyckoff ? s.wyckoff : (s.multiplicity + (s.special ? '*' : ''));
                sitesSection.push(formatLine([
                    s.rank,
                    s.element || 'Q',
                    num(s.x, 4),
                    num(s.y, 4),
                    num(s.z, 4),
                    wm,
                    num(s.relative || s.height || s.mapDensity, 3),
                    Number.isFinite(s.chargeRel) ? num(s.chargeRel, 3) : '-'
                ], siteWidths));
            });
            sitesSection.push('');
        }
    } else {
        structureSection = [
            '--- Structure ---',
            'Structure not built. Map only.',
            ''
        ];
        
    }

    // ------------------------------------------------------------------
    //  ELECTRON-DENSITY PEAKS
    //
    //  This used to appear ONLY when no structure had been built -- the raw
    //  peak list was the else-branch of "did the structure builder run". But
    //  the peaks are the primary observation and the structure is an
    //  interpretation of them: if the assignment is wrong, the peak list is
    //  the only way to see that from the report. It is now always present,
    //  with the nearest-neighbour distance the on-screen table already shows.
    // ------------------------------------------------------------------
    const peaksSection = cfPeakSection(cfResult, formatLine, num);

    return [
        ...header,
        ...summarySection,
        ...cellSection,
        ...structureSection,
        ...sitesSection,
        ...peaksSection
    ].join('\n');
}

/**
 * The electron-density peak table, with each peak's distance to its nearest
 * other peak.
 *
 * NEAREST PEAK is the number worth reading here. The origin and the hand of a
 * charge-flipping solution are arbitrary, so an absolute coordinate says
 * nothing on its own; the distance to the closest other peak says whether a
 * peak is a plausible separate atom or a ripple sitting beside a stronger one.
 *
 * @param {object} cfResult  { peaks, cell, ... }
 * @param {(cols:Array, widths:number[])=>string} formatLine
 * @param {(v:number, n:number)=>string} num
 * @returns {string[]}
 */
function cfPeakSection(cfResult, formatLine, num) {
    const peaks = cfResult && cfResult.peaks;
    if (!Array.isArray(peaks) || peaks.length === 0) return [];

    const out = [
        '--- Electron-density peaks ---',
        'Fractional coordinates in P1. The origin and the hand of the solution are',
        'arbitrary -- compare interatomic distances, not absolute positions.',
        '"Nearest peak" is the distance to the closest other peak in this list: it says',
        'how isolated each peak is, and is neither a resolution limit nor a bond length.',
        'A value far below any real bond means the peak is a ripple beside a stronger',
        'neighbour, not an atom of its own.',
        '',
        'Height is the tallest single voxel of the peak, relative to the tallest peak in',
        'the map; it is sensitive to grid placement and to ripples. Charge is the density',
        'integrated over a fixed sphere, which is what actually ranks atoms by electron',
        'count. On PbSO4 Height puts sulfur above lead and Charge does not.',
        ''
    ];

    // cfMetricTensor and cfFracDistance are top-level in powder5.html. Guard
    // anyway: this file is also loaded by tooling that has no page around it,
    // and a missing distance is better than a report that throws.
    const haveGeom = (typeof cfMetricTensor === 'function')
                  && (typeof cfFracDistance === 'function')
                  && cfResult.cell;
    const G = haveGeom ? cfMetricTensor(cfResult.cell) : null;

    const w = [5, 10, 10, 10, 10, 10, 16];
    const head = formatLine(['#', 'x', 'y', 'z', 'Height', 'Charge', 'Nearest peak (A)'], w);
    out.push(head, '-'.repeat(head.length));

    for (const pk of peaks) {
        let dmin = Infinity;
        if (G) {
            for (const q of peaks) {
                if (q === pk) continue;
                const d = cfFracDistance(G, pk.x - q.x, pk.y - q.y, pk.z - q.z);
                if (d < dmin) dmin = d;
            }
        }
        out.push(formatLine([
            pk.rank,
            num(pk.x, 4), num(pk.y, 4), num(pk.z, 4),
            num(pk.relative, 3),
            Number.isFinite(pk.chargeRel) ? num(pk.chargeRel, 3) : '-',
            Number.isFinite(dmin) ? dmin.toFixed(2) : '-'
        ], w));
    }
    out.push('');
    return out;
}