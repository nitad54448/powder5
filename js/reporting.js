// Reusable PDF generator. `results` is the snapshot to report; `modeName` is
// 'Pawley' or 'LeBail'. Refactored out of the old button handler so the
// per-run history tabs can export any stored run, not just the live fit.
async function generatePdfReport(results, modeName) {
    window.__activeReportResults = results;
    try {
        if (!window.jspdf || !window.jspdf.jsPDF) {
            throw new Error('jspdf.umd.min.js was not loaded. It belongs at lib/jspdf.umd.min.js, next to the other vendor bundles.');
        }
        
        // Removed the html2canvas check here since we no longer need it.

        const { jsPDF } = window.jspdf;
        const doc = new jsPDF({ orientation: 'p', unit: 'mm', format: 'a4' });
        const margin = 15;
        const contentWidth = doc.internal.pageSize.getWidth() - 2 * margin;
        let yPosition = 20;

        const summaryText = generateReportContent('summary');
        // jsPDF's built-in fonts are WinAnsi, so anything outside it renders
        // as a blank or a mojibake box. The reflection table added a few more
        // non-ASCII characters (Angstrom, theta, superscripts), so the
        // transliteration list is extended to match rather than leaving holes
        // in the column headers.
        const pdfText = summaryText
            .replace(/χ²/g, 'Chi^2')
            .replace(/°/g, 'deg')
            .replace(/β/g, 'Beta')
            .replace(/θ/g, 'th')
            .replace(/η/g, 'eta')
            .replace(/λ/g, 'lambda')
            .replace(/σ/g, 'sigma')
            .replace(/Å/g, 'A')
            .replace(/²/g, '^2')
            .replace(/·/g, '.')
            .replace(/[–—]/g, '-')
            .replace(/[‘’]/g, "'")
            .replace(/[“”]/g, '"');

        const lines = pdfText.split('\n');
        doc.setFont('Courier');

        for (const line of lines) {
            if (yPosition > doc.internal.pageSize.getHeight() - 20) {
                doc.addPage();
                yPosition = 20;
            }
            const isHeader = line.startsWith('---');
            const isTitle = line.includes('Refinement Report');
            let fontSize = 9;
            let fontStyle = 'normal';

            if (isTitle) {
                fontSize = 14;
                fontStyle = 'bold';
                yPosition += 6;
            } else if (isHeader) {
                fontSize = 10;
                fontStyle = 'bold';
                yPosition += 4;
            }

            doc.setFontSize(fontSize);
            doc.setFont(undefined, fontStyle);

            // The reflection table is wider than it used to be (it now carries
            // m, Lp and |Fo|), and a Courier line that overflows the text
            // width is silently clipped by jsPDF -- the columns would just
            // vanish off the right edge of the page. Shrink the offending
            // line instead; monospace means the whole row stays aligned.
            let drawSize = fontSize;
            try {
                const w = doc.getTextWidth(line);
                if (w > contentWidth && w > 0) {
                    drawSize = Math.max(4.5, fontSize * (contentWidth / w));
                    doc.setFontSize(drawSize);
                }
            } catch { /* getTextWidth unavailable: fall back to the full size */ }

            doc.text(line, margin, yPosition);
            doc.setFontSize(fontSize);
            yPosition += fontSize * 0.4;
        }

        doc.addPage();
        
        // --- Native Canvas Image Extraction ---
        const chartCanvas = controls.mainChartCanvas;
        
        // Create a temporary canvas to apply a white background
        const tempCanvas = document.createElement('canvas');
        tempCanvas.width = chartCanvas.width;
        tempCanvas.height = chartCanvas.height;
        const ctx = tempCanvas.getContext('2d');
        
        // Fill white, then draw the chart over it
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, tempCanvas.width, tempCanvas.height);
        ctx.drawImage(chartCanvas, 0, 0);

        // Extract image data natively
        const mainImgData = tempCanvas.toDataURL('image/png', 1.0);
        const mainImgHeight = Math.min((tempCanvas.height * contentWidth) / tempCanvas.width, 250); 
        
        doc.addImage(mainImgData, 'PNG', margin, 20, contentWidth, mainImgHeight);
        // --------------------------------------
        
        doc.save(`${modeName}-Report-${new Date().toISOString().slice(0, 10)}.pdf`);
    } finally {
        window.__activeReportResults = null;
    }
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
        `Data File: ${controls.fileName.textContent || 'N/A'}`,
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
        `Refinement:   ${modeName}`,
        `Profile:      ${profileName} (${profileType})`,
        `Space Group:  ${spaceGroupName}`,
        // The background is interpolated through fixed anchors, never refined.
        // Stating it here rather than leaving it implicit: it is the reason no
        // sigma below carries background uncertainty, and the reason chi-square
        // is optimistic by roughly the anchor count.
        `Background:   ${stats.backgroundNote || 'fixed spline through user anchor points, not refined'}`,
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
    if (JtJ && parameterInfo && ss_res !== undefined && workerWorkingData && workerWorkingData.intensity && workerWorkingData.intensity.length > 0) {
         const P = parameterInfo.length;
         const N = workerWorkingData.intensity.length;
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
                         if (Array.isArray(JtJ) && JtJ.length === P && Array.isArray(JtJ[0]) && JtJ[0].length === P) {
                             
                             const covResult = covarianceFromJtJ(JtJ);
                             if (!covResult) {
                                 throw new Error("Normal-equations matrix is not positive definite.");
                             }
                             cov_matrix = covResult.cov.map(row => row.map(v => v * reduced_chi_sq));
                         } else {
                             throw new Error("JtJ matrix has incorrect dimensions or format.");
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
                                 const variance = cov_matrix[i][i];
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
    } else if (!JtJ || !parameterInfo || ss_res === undefined || !workerWorkingData || !workerWorkingData.intensity || workerWorkingData.intensity.length === 0) {
         esdWarning = "Required data for ESD calculation is missing (normal-equations matrix, parameterInfo, ss_res, or the fitted data slice).";
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

    for (const groupName in paramGroups) {
        // Ensure group exists and has items before adding headers
        if (paramGroups[groupName] && paramGroups[groupName].length > 0) {
            paramLines.push(`--- ${groupName} ---`, paramHeader, '-'.repeat(paramHeader.length));
            paramGroups[groupName].forEach(p => {
                 // Check value is valid number before formatting
                 if (p.value !== null && !isNaN(p.value)) {
                    const valStr = (typeof p.value === 'number') ? p.value.toExponential(6) : String(p.value);
                    const fitStr = (p.flag === undefined) ? '' : (p.flag ? 'Yes' : 'No');
                    let esdStr = '';
                    if (p.flag && p.esd_key) {
                        const esdValue = esds[p.esd_key];
                        if (typeof esdValue === 'number' && isFinite(esdValue)) {
                             esdStr = `(${esdValue.toExponential(2)})`;
                        } else if (isNaN(esdValue)) {
                             esdStr = '(NaN)';
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
         const bgWidths = [15, 18]; // 2theta, Intensity
         const bgHeader = formatLine(['2theta (°)', 'Intensity'], bgWidths);
         backgroundPointsSection.push('', '--- Background Spline Points ---', bgHeader, '-'.repeat(bgHeader.length));
         // Use a copy sorted by tth for the report
         const sortedAnchors = [...backgroundAnchors].sort((a, b) => a.tth - b.tth);
         sortedAnchors.forEach(point => {
              backgroundPointsSection.push(formatLine([
                  point.tth.toFixed(4),
                  point.y.toFixed(2)
              ], bgWidths));
         });
         backgroundPointsSection.push(''); // Add blank line after section
    }
    

    const reflectionsSection = [];
    if (hklList && hklList.length > 0 && workerWorkingData && workerWorkingData.tth) {
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
                     const sig = hkl.I_sigma;
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
            const oc = intensityOverlapClusters(fitResults, finalParams);
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
        `Data File: ${(controls && controls.fileName && controls.fileName.textContent) || 'N/A'}`,
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
        'sqrt(I). wR is the residual actually minimised. No charge-flipping map was used, so there is no origin shift, no',
        'hand and no symmetry correlation to report: the Wyckoff operators place every',
        'atom on a crystallographic site by construction.');
    return lines.join('\n');
}

function generateChargeFlippingReport(cfResult) {
    const now = new Date();
    const header = [
        `Charge Flipping Report`,
        '========================================================================',
        '',
        typeof APP_VERSION_TEXT !== 'undefined' ? APP_VERSION_TEXT : 'Powder 5',
        `Report Generated: ${now.toLocaleString()}`,
        `Data File: ${(typeof controls !== 'undefined' && controls.fileName) ? controls.fileName.textContent : 'N/A'}`,
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
            const siteWidths = [4, 8, 10, 10, 10, 14, 10, 10, 18];
            const siteHeader = formatLine(['#', 'Element', 'x', 'y', 'z', 'Wyckoff/Mult', 'Height', 'Charge', 'Shortest Contact'], siteWidths);

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
                const contact = s.shortestContact !== undefined && s.shortestContact !== null 
                                ? s.shortestContact.toFixed(2) + ' Å' : '–';
                sitesSection.push(formatLine([
                    s.rank,
                    s.element || 'Q',
                    num(s.x, 4),
                    num(s.y, 4),
                    num(s.z, 4),
                    wm,
                    num(s.relative || s.height || s.mapDensity, 3),
                    Number.isFinite(s.chargeRel) ? num(s.chargeRel, 3) : '-',
                    contact
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
        
        if (cfResult.peaks && cfResult.peaks.length > 0) {
            const peakWidths = [6, 10, 10, 10, 10, 10];
            const peakHeader = formatLine(['Rank', 'x', 'y', 'z', 'Height', 'Charge'], peakWidths);
            sitesSection = [
                '--- Raw Peaks (P1, arbitrary origin) ---',
                peakHeader,
                '-'.repeat(peakHeader.length)
            ];
            cfResult.peaks.forEach(p => {
                sitesSection.push(formatLine([
                    p.rank, num(p.x, 4), num(p.y, 4), num(p.z, 4), num(p.relative, 3),
                    Number.isFinite(p.chargeRel) ? num(p.chargeRel, 3) : '-'
                ], peakWidths));
            });
            sitesSection.push('');
        }
    }

    return [
        ...header,
        ...summarySection,
        ...cellSection,
        ...structureSection,
        ...sitesSection
    ].join('\n');
}