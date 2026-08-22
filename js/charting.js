const verticalCursorLine = { id: 'verticalCursorLine', afterDraw: chart => { if (chart.tooltip?._active?.length) { let x = chart.tooltip._active[0].element.x, yAxis = chart.scales.y, ctx = chart.ctx; ctx.save(); ctx.beginPath(); ctx.moveTo(x, yAxis.top); ctx.lineTo(x, yAxis.bottom); ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(156, 163, 175, 0.7)'; ctx.setLineDash([4, 4]); ctx.stroke(); ctx.restore(); } } };
Chart.register(verticalCursorLine);
Chart.Tooltip.positioners.experimentalAnchor = function(items) { if (!items.length) return false; const experimentalItem = items.find(item => item.datasetIndex === 0) || items[0]; return { x: experimentalItem.element.x, y: experimentalItem.element.y }; };

const THEO_TTH_MIN = 5;
const THEO_TTH_MAX = 120;

function hasExperimentalData() {
    return !!(fullExperimentalData && fullExperimentalData.tth && fullExperimentalData.tth.length > 0);
}

function getPeakInfoAt(tth) {
    if (!lastGeneratedHklList || lastGeneratedHklList.length === 0 || !mainChart) return { peak: null, inRegion: false };
    let closestPeak = null, minDiff = Infinity; const currentParams = getAllParams(); const zeroShift = currentParams.zeroShift || 0;
    for (const hkl of lastGeneratedHklList) { const peakShift = calculatePeakShift(hkl.tth, currentParams); const peakPos = hkl.tth + zeroShift + peakShift; const diff = Math.abs(tth - peakPos); if (diff < minDiff) { minDiff = diff; closestPeak = hkl; } }
    if (closestPeak) { 
        const threshold = 0.2;
           if (minDiff < threshold) return { peak: closestPeak, inRegion: true }; }
    return { peak: null, inRegion: false };
}

function initializeChart() {
    if (mainChart) mainChart.destroy();

    const experimentalPoints = [];
    for (let i = 0; i < fullExperimentalData.tth.length; i++) {
        experimentalPoints.push({ x: fullExperimentalData.tth[i], y: fullExperimentalData.intensity[i] });
    }
    
    const yMax = fullExperimentalData.intensity.reduce((m, v) => (v > m ? v : m), 0) || 1000;
    mainChart = new Chart(controls.mainChartCanvas, {
        type: 'line',
        data: { 
            datasets: [
                { label: 'Experimental', data: experimentalPoints, borderColor: 'rgba(107, 114, 128, 0.7)', borderWidth: 0.5, pointRadius: 1, pointBackgroundColor: 'rgba(107, 114, 128, 0.7)', order: 1 },
                { label: 'Simulation (Manual)', data: [], borderColor: 'rgba(249, 115, 22, 0.8)', borderWidth: 2, pointRadius: 0, borderDash: [8, 4], order: 2 },
                { label: 'Calculated', data: [], borderColor: 'rgba(59, 130, 246, 1)', borderWidth: 2, pointRadius: 0, order: 3 },
                { label: 'Background', data: [], borderColor: 'rgba(2, 9, 206, 0.8)', borderWidth: 1.5, pointRadius: 0, borderDash: [5, 5], order: 4 },
                { label: 'Difference', data: [], borderColor: 'rgba(239, 68, 68, 0.8)', borderWidth: 1.5, pointRadius: 0, order: 5 },
                { label: 'Difference Zero', data: [], borderColor: 'rgba(156, 163, 175, 0.8)', borderWidth: 1, pointRadius: 0, borderDash: [2, 2], order: 6 },
                { type: 'bar', label: 'Ka1', data: [], backgroundColor: 'rgba(22, 163, 74, 0.9)', barThickness: 1, categoryPercentage: 1.0, barPercentage: 1.0, order: 0 },
                { type: 'bar', label: 'Ka2', data: [], backgroundColor: 'rgba(134, 239, 172, 0.6)', barThickness: 1, categoryPercentage: 1.0, barPercentage: 1.0, order: 0 },
                {
                    label: 'Spline Points', 
                    type: 'scatter',
                    data: [],
                    showLine: false,
                    hidden: true,
                    pointBackgroundColor: 'rgba(34, 197, 94, 1)',
                    pointBorderColor: 'rgba(255, 255, 255, 0.9)',
                    pointRadius: 5,
                    pointHoverRadius: 7,
                    pointBorderWidth: 1.5,
                    order: -1
                }
            ]
        },
        options: {
            responsive: true, maintainAspectRatio: false, animation: false,
            scales: {
                x: { type: 'linear', title: { display: true, text: '2θ (degrees)', font: { size: 14 }}},
                y: { type: 'linear', position: 'left', title: { display: true, text: 'Intensity (a.u.)', font: { size: 14 }},
                    min: -yMax * 0.3, max: Math.ceil(yMax * 1.1),
                    ticks: { callback: function(value, index, ticks) { return value >= 0 ? value.toFixed(1) : null; }}
                }
            },
            plugins: {
                legend: { labels: { filter: item => item.text !== 'Difference Zero' && item.text !== 'Simulation (Manual)' && item.text !== 'Spline Points'}},
                tooltip: {
                    enabled: true, mode: 'nearest', intersect: false, position: 'experimentalAnchor',
                    filter: function(tooltipItem) {
                        const tth = tooltipItem.parsed.x;
                        const peakInfo = getPeakInfoAt(tth);
                        return peakInfo.inRegion;
                    },
                    callbacks: {
                        title: function(tooltipItems) {
                            if (!tooltipItems.length) return '';
                            const tth = tooltipItems[0].parsed.x;
                            const peakInfo = getPeakInfoAt(tth);
                            if (peakInfo.inRegion) {
                                const closestPeak = peakInfo.peak;
                                return `d: ${closestPeak.d.toFixed(2)} Å, HKL: ${closestPeak.hkl_list[0]}`;
                            }
                            return '';
                        },
                        label: function(context) { return null; }
                    }
                }
            }
        }
    });
    mainChart.options.globalYMax = yMax;
}

function updatePlotRange(recalculateYMax = false) {
    if (!mainChart || workingData.tth.length === 0) return;
    const minTth = parseFloat(controls.tthMinSlider.value);
    const maxTth = parseFloat(controls.tthMaxSlider.value);

    mainChart.options.scales.x.min = minTth;
    mainChart.options.scales.x.max = maxTth;

    if (recalculateYMax) {
        let currentMaxY = workingData.intensity.reduce((m, v) => (v > m ? v : m), -Infinity);
        
        const calcDataset = mainChart.data.datasets.find(d => d.label === 'Calculated');
        if (calcDataset && calcDataset.data && calcDataset.data.length > 0) {
            const calcMax = calcDataset.data.reduce((m, pt) => (pt.y > m ? pt.y : m), -Infinity);
            if (calcMax > currentMaxY) currentMaxY = calcMax;
        }

        const simDataset = mainChart.data.datasets.find(d => d.label === 'Simulation (Manual)');
        if (simDataset && simDataset.data && simDataset.data.length > 0) {
            const simMax = simDataset.data.reduce((m, pt) => (pt.y > m ? pt.y : m), -Infinity);
            if (simMax > currentMaxY) currentMaxY = simMax;
        }

        let newYMax = currentMaxY * 1.1; 
        if (newYMax < 100) newYMax = 100; 

        mainChart.options.globalYMax = newYMax;
    }
    
    mainChart.update('none');
}

function rescalePlot(updateY = false) {
    const DIFF_PLOT_PADDING = 1.1; 
    if (!mainChart || !workingData.isValid) return;

    const findDataset = (label) => mainChart.data.datasets.find(d => d.label === label);
    const globalYMax = mainChart.options.globalYMax || 1000;

    let diffPlotHeightRatio = 0.15;
    const diffData = workingData.lastRawDifference || [];
    let constantNegSpaceHeight = globalYMax * diffPlotHeightRatio;

    if (diffData.length > 0) {
        const maxAbsDiff = diffData.reduce((m, v) => (Math.abs(v) > m ? Math.abs(v) : m), 0);
        const requiredNegSpace = maxAbsDiff * DIFF_PLOT_PADDING; 
        constantNegSpaceHeight = Math.max(constantNegSpaceHeight, requiredNegSpace);
        if (constantNegSpaceHeight > globalYMax * 0.6) {
            constantNegSpaceHeight = globalYMax * 0.6;
        }
    }

    if (updateY) {
        mainChart.options.scales.y.max = globalYMax;
        mainChart.options.scales.y.min = -constantNegSpaceHeight;
    }
    
    const formatMarkers = (label, pixelHeight) => {
        const ds = findDataset(label);
        if (ds && ds.data.length > 0) {
            const y_scale = mainChart.scales.y;
            const zeroVal = 0; 
            const zeroPx = y_scale.getPixelForValue(zeroVal);
            const bottomPx = zeroPx + pixelHeight; 
            const bottomVal = y_scale.getValueForPixel(bottomPx);
            ds.data.forEach(pt => pt.y = [bottomVal, zeroVal]);
        }
    };

    formatMarkers('Ka1', 20);      
    formatMarkers('Ka2', 12);  

    if (diffData.length > 0 && diffData.length === workingData.tth.length) {
        const constantDiffPlotSpaceTop = -(constantNegSpaceHeight * 0.25);
        const diffPlotZeroLine = -constantNegSpaceHeight + (constantNegSpaceHeight - constantDiffPlotSpaceTop) / 2;

        const differenceDataPoints = [];
        for (let i = 0; i < workingData.tth.length; i++) {
            const difference = diffData[i];
            differenceDataPoints.push({ x: workingData.tth[i], y: difference + diffPlotZeroLine });
        }
        findDataset('Difference').data = differenceDataPoints;

        findDataset('Difference Zero').data = [
            { x: workingData.tth[0], y: diffPlotZeroLine },
            { x: workingData.tth[workingData.tth.length - 1], y: diffPlotZeroLine }
        ];
    } else {
        findDataset('Difference').data = [];
        findDataset('Difference Zero').data = [];
    }
}

function updateChart(netPeakPattern_sliced, background_sliced, hklList, params, scaleFactor = 1.0) {
    if (!mainChart || !workingData.isValid) return;

    const currentYMin = mainChart.scales.y.min;
    const currentYMax = mainChart.scales.y.max;
    const currentXMin = mainChart.scales.x.min;
    const currentXMax = mainChart.scales.x.max;
    const oldGlobalYMax = mainChart.options.globalYMax || 1000;
    const wasZoomedOutY = (currentYMax >= oldGlobalYMax * 0.99);

    const findDataset = (label) => mainChart.data.datasets.find(d => d.label === label);
    findDataset('Simulation (Manual)').data = [];

    refreshHklMarkers(params);

    const len = workingData.intensity.length;
    const totalCalcPattern = new Float64Array(len);
    const diff = new Float64Array(len);
    for (let i = 0; i < len; i++) {
        totalCalcPattern[i] = netPeakPattern_sliced[i] * scaleFactor + background_sliced[i];
        diff[i] = workingData.intensity[i] - totalCalcPattern[i];
    }
    workingData.lastRawDifference = diff;

    const calculatedData = [];
    const backgroundData = [];
    for(let i = 0; i < workingData.tth.length; i++) {
        calculatedData.push({ x: workingData.tth[i], y: totalCalcPattern[i] });
        backgroundData.push({ x: workingData.tth[i], y: background_sliced[i] });
    }
    findDataset('Calculated').data = calculatedData;
    findDataset('Background').data = backgroundData;

    updatePlotRange(true);

    if (wasZoomedOutY) {
        rescalePlot(true);
    } else {
        mainChart.options.scales.y.min = currentYMin;
        mainChart.options.scales.y.max = currentYMax;
        rescalePlot(false); 
    }

    mainChart.options.scales.x.min = currentXMin;
    mainChart.options.scales.x.max = currentXMax;
    mainChart.update('none');
}

function updatePreviewPattern() {
    if (!mainChart || !workingData.isValid || isFitting) return;

    const currentXMin = mainChart.scales.x.min;
    const currentXMax = mainChart.scales.x.max;
    const currentYMin = mainChart.scales.y.min;
    const currentYMax = mainChart.scales.y.max;
    
    const oldGlobalYMax = mainChart.options.globalYMax || 1000;
    const wasZoomedOutY = (currentYMax >= oldGlobalYMax * 0.99);

    const findDataset = (label) => mainChart.data.datasets.find(d => d.label === label);
    const params = getAllParams();
    if (isNaN(params.lambda) || params.lambda <= 0 || isNaN(params.a) || params.a <= 0) return;
    
    const selectedSg = currentSG;
    if (!selectedSg) return;

    let hklList = JSON.parse(JSON.stringify(masterHklList));
    updateHklPositions(hklList, params, selectedSg.system);
    lastGeneratedHklList = hklList; 

    if (lastFitResultsCache && lastFitResultsCache.hklList) {
        const intensityMap = new Map(
            lastFitResultsCache.hklList.map(p => [p.hkl_list[0], p.intensity])
        );
        hklList.forEach(p => {
            p.intensity = intensityMap.get(p.hkl_list[0]) || 0;
        });
    } else {
        const FALLBACK_I = 1000;
        if (workingData && workingData.tth && workingData.tth.length > 0) {
            const tthArr = workingData.tth;
            const yArr   = workingData.intensity;
            hklList.forEach(p => {
                if (!p || !p.tth) { if (p) p.intensity = FALLBACK_I; return; }
                const peakPos = p.tth + params.zeroShift;
                const win = 0.5;
                let localMax = 0;
                let localBg  = Infinity;
                for (let i = 0; i < tthArr.length; i++) {
                    const t = tthArr[i];
                    if (t < peakPos - win) continue;
                    if (t > peakPos + win) break;
                    const v = yArr[i];
                    if (v > localMax) localMax = v;
                    if (v < localBg)  localBg  = v;
                }
                const above = (isFinite(localBg) && localMax > localBg) ? (localMax - localBg) : localMax;
                p.intensity = (above > 1) ? above : FALLBACK_I;
            });
        } else {
            hklList.forEach(p => p.intensity = FALLBACK_I);
        }
    }
    
    refreshHklMarkers(params, hklList);

    const background_sliced = calculateTotalBackground(workingData.tth, params, backgroundAnchors);
    const unscaledPeakPattern_sliced = calculatePattern(workingData.tth, hklList, params);
   
    const len = workingData.intensity.length;
    let sum_obs_calc = 0, sum_calc_sq = 0;
    for (let i = 0; i < len; i++) {
        const y_net = Math.max(0, workingData.intensity[i] - background_sliced[i]);
        const y_calc = unscaledPeakPattern_sliced[i];
        sum_obs_calc += y_net * y_calc;
        sum_calc_sq += y_calc * y_calc;
    }
    const scaleFactor = (sum_calc_sq > 1e-9) ? sum_obs_calc / sum_calc_sq : 1.0;

    const totalCalcPattern = new Float64Array(len);
    const diff = new Float64Array(len);
    for (let i = 0; i < len; i++) {
        totalCalcPattern[i] = unscaledPeakPattern_sliced[i] * scaleFactor + background_sliced[i];
        diff[i] = workingData.intensity[i] - totalCalcPattern[i];
    }
    workingData.lastRawDifference = diff;

    const simulationData = [];
    for (let i = 0; i < workingData.tth.length; i++) {
        simulationData.push({ x: workingData.tth[i], y: totalCalcPattern[i] });
    }
    findDataset('Simulation (Manual)').data = simulationData;
    findDataset('Calculated').data = [];

    updatePlotRange(true);

    if (wasZoomedOutY) {
        rescalePlot(true);
    } else {
        mainChart.options.scales.y.min = currentYMin;
        mainChart.options.scales.y.max = currentYMax;
        rescalePlot(false);
    }

    mainChart.options.scales.x.min = currentXMin;
    mainChart.options.scales.x.max = currentXMax;

    mainChart.update('none');
}

function redrawFitForNewRange() {
    if (!mainChart) return; 

    const currentYMin = mainChart.scales.y.min;
    const currentYMax = mainChart.scales.y.max;
    const oldGlobalYMax = mainChart.options.globalYMax || 1000;
    const wasZoomedOutY = (currentYMax >= oldGlobalYMax * 0.99);

    updateWorkingData();
    if (!fitResults) {
        updatePreviewPattern();
    } else {
        const params = fitResults.params;
        const scaleFactor = fitResults.stats.scaleFactor;
        
        const background_sliced = calculateTotalBackground(workingData.tth, params, backgroundAnchors);
        const unscaledPeakPattern_sliced = calculatePattern(workingData.tth, fitResults.hklList, params);
       
        const len = workingData.intensity.length;
        const totalCalcPattern = new Float64Array(len);
        const diff = new Float64Array(len);
        for (let i = 0; i < len; i++) {
            totalCalcPattern[i] = unscaledPeakPattern_sliced[i] * scaleFactor + background_sliced[i];
            diff[i] = workingData.intensity[i] - totalCalcPattern[i];
        }
        workingData.lastRawDifference = diff;

        const findDataset = (label) => mainChart.data.datasets.find(d => d.label === label);
        
        const calculatedData = [];
        const backgroundData = [];
        for(let i = 0; i < workingData.tth.length; i++) {
            calculatedData.push({ x: workingData.tth[i], y: totalCalcPattern[i] });
            backgroundData.push({ x: workingData.tth[i], y: background_sliced[i] });
        }

        findDataset('Calculated').data = calculatedData;
        findDataset('Background').data = backgroundData;
        findDataset('Simulation (Manual)').data = [];

        refreshHklMarkers(params);

        updatePlotRange(true);

        if (wasZoomedOutY) {
            rescalePlot(true);
        } else {
            mainChart.options.scales.y.min = currentYMin;
            mainChart.options.scales.y.max = currentYMax;
            rescalePlot(false);
        }
    }
    
    mainChart.options.scales.x.min = parseFloat(controls.tthMinSlider.value);
    mainChart.options.scales.x.max = parseFloat(controls.tthMaxSlider.value);
    mainChart.update('none');
}

function refreshHklMarkers(params, positioned) {
    if (!mainChart) return;
    const findDataset = (label) => mainChart.data.datasets.find(d => d.label === label);
    const ka1Dataset = findDataset('Ka1');
    const ka2Dataset = findDataset('Ka2');
    if (!ka1Dataset) return;

    const clear = () => {
        ka1Dataset.data = [];
        if (ka2Dataset) ka2Dataset.data = [];
    };
    if (!params) { clear(); return; }

    let list = positioned;
    if (!list) {
        if (!masterHklList || masterHklList.length === 0) { clear(); return; }
        list = masterHklList.map(h => ({
            h_orig: h.h_orig, k_orig: h.k_orig, l_orig: h.l_orig,
            hkl_list: h.hkl_list, multiplicity: h.multiplicity
        }));
        const system = (currentSG && currentSG.system) || params.system;
        try { updateHklPositions(list, params, system); }
        catch (err) { console.warn('refreshHklMarkers: could not position reflections', err); clear(); return; }
    }


    const dataLen = fullExperimentalData.tth.length;
    const dataMin = dataLen > 0 ? parseFloat(controls.tthMinSlider.value) : -Infinity;
    const dataMax = dataLen > 0 ? parseFloat(controls.tthMaxSlider.value) : Infinity;

    const ka1Markers = [];
    const ka2Markers = [];
    const deg2rad = Math.PI / 180;
    const zeroShift = params.zeroShift || 0;
    const doubletEnabled = params.ratio > 1e-6 && params.lambda2 > 1e-6
                           && Math.abs(params.lambda - params.lambda2) > 1e-6;
    const lambdaRatio = doubletEnabled ? (params.lambda2 / params.lambda) : 1.0;

    list.forEach(hkl => {
        if (!hkl || typeof hkl.tth !== 'number' || !isFinite(hkl.tth)) return;
        const label = (hkl.hkl_list && hkl.hkl_list[0]) ? hkl.hkl_list[0] : '';
        const pos1 = hkl.tth + zeroShift;
        if (pos1 >= dataMin && pos1 <= dataMax) {
            ka1Markers.push({ x: pos1, hkl: `${label} [m=${hkl.multiplicity}]` });
        }
        if (doubletEnabled) {
            const sinTheta2 = Math.sin(hkl.tth * deg2rad / 2.0) * lambdaRatio;
            if (Math.abs(sinTheta2) < 1) {
                const pos2 = (2 * Math.asin(sinTheta2) / deg2rad) + zeroShift;
                if (pos2 >= dataMin && pos2 <= dataMax) {
                    ka2Markers.push({ x: pos2, hkl: `${label} Ka2` });
                }
            }
        }
    });

    ka1Dataset.data = ka1Markers;
    if (ka2Dataset) ka2Dataset.data = ka2Markers;
}

// ===========================================================================
//  drawTheoreticalPreview() USED TO BE DUPLICATED HERE.
//
//  powder5.html defines its own copy inside the DOMContentLoaded callback, and
//  every caller of it lives in that same scope, so the inline copy shadowed
//  this one and this one never ran.
//
//  The two had already drifted. The copy here skipped reflections with
//  tth <= 0 as well as non-finite ones; the live copy skips only the
//  non-finite. That guard is therefore NOT carried over by this deletion --
//  removing dead code should not change what the program does, and adding the
//  stricter test to the live copy would. If a preview peak at or below zero
//  2-theta is worth suppressing (it can only arise from a zero shift large
//  enough to push a low-angle reflection negative), add `|| hkl.tth <= 0` to
//  the guard in powder5.html deliberately, as a change in its own right.
// ===========================================================================


function initCharting() {
    const RECT_ZOOM_MIN_PX = 8;
    const WHEEL_ZOOM_STEP  = 1.15;
    let rectZoomState = null;

    const zoomSelectionBox = document.createElement('div');
    zoomSelectionBox.id = 'zoom-selection-box';
    (controls.mainChartCanvas.parentElement || document.body).appendChild(zoomSelectionBox);

    function drawZoomSelectionBox() {
        const host = zoomSelectionBox.parentElement;
        if (!host || !rectZoomState) return;
        const hostRect   = host.getBoundingClientRect();
        const canvasRect = controls.mainChartCanvas.getBoundingClientRect();
        const dx = canvasRect.left - hostRect.left;
        const dy = canvasRect.top  - hostRect.top;

        const s = rectZoomState;
        zoomSelectionBox.style.left    = (dx + Math.min(s.x0, s.x1)) + 'px';
        zoomSelectionBox.style.top     = (dy + Math.min(s.y0, s.y1)) + 'px';
        zoomSelectionBox.style.width   = Math.abs(s.x1 - s.x0) + 'px';
        zoomSelectionBox.style.height  = Math.abs(s.y1 - s.y0) + 'px';
        zoomSelectionBox.style.display = 'block';
    }

    function cancelRectZoom() {
        rectZoomState = null;
        zoomSelectionBox.style.display = 'none';
    }

    let panState = null;

    function beginPan(clientX, clientY) {
        const chart = mainChart;
        if (!chart || !chart.chartArea || !chart.scales || !chart.scales.x || !chart.scales.y) return false;
        const xs = chart.scales.x, ys = chart.scales.y;
        if (![xs.min, xs.max, ys.min, ys.max].every(isFinite)) return false;

        panState = {
            clientX, clientY,
            xMin: xs.min, xMax: xs.max,
            yMin: ys.min, yMax: ys.max,
            xPerPx: (xs.max - xs.min) / Math.max(1, xs.right - xs.left),
            yPerPx: (ys.max - ys.min) / Math.max(1, ys.bottom - ys.top)
        };
        controls.mainChartCanvas.style.cursor = 'grabbing';
        return true;
    }

    function movePan(clientX, clientY) {
        const p = panState;
        if (!p || !mainChart) return;
        const dx = (clientX - p.clientX) * p.xPerPx;
        const dy = (clientY - p.clientY) * p.yPerPx;
        mainChart.options.scales.x.min = p.xMin - dx;
        mainChart.options.scales.x.max = p.xMax - dx;
        mainChart.options.scales.y.min = p.yMin + dy;
        mainChart.options.scales.y.max = p.yMax + dy;
        mainChart.update('none');
    }

    function endPan() {
        if (!panState) return;
        panState = null;
        controls.mainChartCanvas.style.cursor = '';
        rescalePlot(false);
    }

    controls.mainChartCanvas.addEventListener('mousedown', e => {
        const chart = mainChart;
        if (!chart || !chart.chartArea) return;
        const { left, right, top, bottom } = chart.chartArea;
        const x = e.offsetX, y = e.offsetY;
        if (x < left || x > right || y < top || y > bottom) return; 

        if (e.button === 1 || (e.button === 0 && (e.altKey || e.shiftKey))) {
            cancelRectZoom();
            if (beginPan(e.clientX, e.clientY)) e.preventDefault();
            return;
        }

        if (e.button !== 0) return; 
        if (e.ctrlKey || e.metaKey) return;
        rectZoomState = { x0: x, y0: y, x1: x, y1: y, moved: false };
        e.preventDefault();
    });

    window.addEventListener('mousemove', e => {
        if (panState) { movePan(e.clientX, e.clientY); return; }
        const s = rectZoomState;
        if (!s) return;
        const chart = mainChart;
        if (!chart || !chart.chartArea) { cancelRectZoom(); return; }

        const { left, right, top, bottom } = chart.chartArea;
        const canvasRect = controls.mainChartCanvas.getBoundingClientRect();
        const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
        s.x1 = clamp(e.clientX - canvasRect.left, left, right);
        s.y1 = clamp(e.clientY - canvasRect.top,  top,  bottom);

        if (!s.moved && (Math.abs(s.x1 - s.x0) > 3 || Math.abs(s.y1 - s.y0) > 3)) s.moved = true;
        if (s.moved) drawZoomSelectionBox();
    });

    window.addEventListener('mouseup', e => {
        if (panState) { endPan(); return; }
        const s = rectZoomState;
        if (!s) return;
        cancelRectZoom();
        if (e.button !== 0 || !s.moved) return;

        const chart = mainChart;
        if (!chart || !chart.scales || !chart.scales.x || !chart.scales.y) return;

        const w = Math.abs(s.x1 - s.x0), h = Math.abs(s.y1 - s.y0);
        if (w < RECT_ZOOM_MIN_PX && h < RECT_ZOOM_MIN_PX) return;

        const xs = chart.scales.x, ys = chart.scales.y;
        let xMin = xs.min, xMax = xs.max, yMin = ys.min, yMax = ys.max;

        if (w >= RECT_ZOOM_MIN_PX) {
            xMin = xs.getValueForPixel(Math.min(s.x0, s.x1));
            xMax = xs.getValueForPixel(Math.max(s.x0, s.x1));
        }
        if (h >= RECT_ZOOM_MIN_PX) {
            yMax = ys.getValueForPixel(Math.min(s.y0, s.y1)); 
            yMin = ys.getValueForPixel(Math.max(s.y0, s.y1));
        }
        if (!isFinite(xMin) || !isFinite(xMax) || xMax <= xMin) return;
        if (!isFinite(yMin) || !isFinite(yMax) || yMax <= yMin) return;

        chart.options.scales.x.min = xMin;
        chart.options.scales.x.max = xMax;
        chart.options.scales.y.min = yMin;
        chart.options.scales.y.max = yMax;
        chart.update('none');
        rescalePlot(false);
    });

    controls.mainChartCanvas.addEventListener('wheel', e => {
        const chart = mainChart;
        if (!chart || !chart.chartArea) return;

        const { left, right, top, bottom } = chart.chartArea;
        if (e.offsetX < left || e.offsetX > right ||
            e.offsetY < top  || e.offsetY > bottom) return;

        e.preventDefault();
        cancelRectZoom();

        const axisId = e.shiftKey ? 'y' : 'x';
        const axis   = chart.scales[axisId];
        if (!axis) return;

        const focalValue = axis.getValueForPixel(e.shiftKey ? e.offsetY : e.offsetX);
        const range      = axis.max - axis.min;
        if (!isFinite(focalValue) || !(range > 0)) return;

        const newRange = range * (e.deltaY < 0 ? 1 / WHEEL_ZOOM_STEP : WHEEL_ZOOM_STEP);
        const ratio    = (focalValue - axis.min) / range;
        const newMin   = focalValue - newRange * ratio;
        const newMax   = newMin + newRange;
        if (!isFinite(newMin) || !isFinite(newMax) || newMax <= newMin) return;

        chart.options.scales[axisId].min = newMin;
        chart.options.scales[axisId].max = newMax;
        chart.update('none');
        rescalePlot(false);
    }, { passive: false });

    let pinchState = null;

    function touchSeparation(touches) {
        return Math.hypot(touches[0].clientX - touches[1].clientX,
                          touches[0].clientY - touches[1].clientY);
    }

    controls.mainChartCanvas.addEventListener('touchstart', e => {
        const chart = mainChart;
        if (!chart || !chart.scales) return;
        cancelRectZoom();

        if (e.touches.length === 1) {
            pinchState = null;
            if (beginPan(e.touches[0].clientX, e.touches[0].clientY)) e.preventDefault();
        } else if (e.touches.length === 2) {
            endPan();
            const sep = touchSeparation(e.touches);
            if (!(sep > 0)) return;
            const rect = controls.mainChartCanvas.getBoundingClientRect();
            const xs = chart.scales.x, ys = chart.scales.y;
            if (![xs.min, xs.max, ys.min, ys.max].every(isFinite)) return;
            pinchState = {
                sep,
                xMin: xs.min, xMax: xs.max, yMin: ys.min, yMax: ys.max,
                xFocal: xs.getValueForPixel((e.touches[0].clientX + e.touches[1].clientX) / 2 - rect.left),
                yFocal: ys.getValueForPixel((e.touches[0].clientY + e.touches[1].clientY) / 2 - rect.top)
            };
            e.preventDefault();
        }
    }, { passive: false });

    controls.mainChartCanvas.addEventListener('touchmove', e => {
        if (pinchState && e.touches.length === 2) {
            const chart = mainChart;
            if (!chart) return;
            const sep = touchSeparation(e.touches);
            if (!(sep > 0)) return;
            e.preventDefault();

            const p = pinchState;
            const factor = p.sep / sep;               
            const applyAxis = (id, focal, lo, hi) => {
                if (!isFinite(focal)) return;
                const range = hi - lo;
                if (!(range > 0)) return;
                const newRange = range * factor;
                const newMin = focal - newRange * ((focal - lo) / range);
                const newMax = newMin + newRange;
                if (!isFinite(newMin) || !isFinite(newMax) || newMax <= newMin) return;
                chart.options.scales[id].min = newMin;
                chart.options.scales[id].max = newMax;
            };
            applyAxis('x', p.xFocal, p.xMin, p.xMax);
            applyAxis('y', p.yFocal, p.yMin, p.yMax);
            chart.update('none');
            return;
        }
        if (panState && e.touches.length === 1) {
            e.preventDefault();
            movePan(e.touches[0].clientX, e.touches[0].clientY);
        }
    }, { passive: false });

    function endTouchNav() {
        if (pinchState) { pinchState = null; rescalePlot(false); }
        endPan();
    }
    controls.mainChartCanvas.addEventListener('touchend', endTouchNav);
    controls.mainChartCanvas.addEventListener('touchcancel', endTouchNav);

    function abortChartNav() { cancelRectZoom(); endPan(); pinchState = null; }
    window.addEventListener('blur', abortChartNav);
    window.addEventListener('keydown', e => { if (e.key === 'Escape') abortChartNav(); });

    controls.mainChartCanvas.addEventListener('contextmenu', e => {
        e.preventDefault();
        abortChartNav();
        if (!mainChart) return;

        mainChart.options.scales.x.min = parseFloat(controls.tthMinSlider.value);
        mainChart.options.scales.x.max = parseFloat(controls.tthMaxSlider.value);

        updatePlotRange(true);
        rescalePlot(true);

mainChart.update('none');
    });
}