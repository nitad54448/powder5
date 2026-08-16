// MOVED to profile.js: createMonotonicCubicSplineInterpolator.
// It was the THIRD copy of the same Fritsch-Carlson spline (the others
// were in powder5.html and refinement_worker.js). profile.js is loaded
// before this file in powder5.html and by importScripts in the worker,
// so the name is already in scope here.

         /**
 * Automatically finds points for the background spline, using local minima
 * refined within a sub-window and smoothed by averaging, ensuring points
 * exist at the current slider min/max tth.
 */
function autoFindSplinePoints(numPoints) {
    const AVERAGE_WINDOW_HALF_WIDTH = 3; // Window size = 2 * half_width + 1 = 7 points
    const REFINE_WINDOW_FRACTION = 0.30; // +/- 30% of chunk width for refinement search

    if (!fullExperimentalData || fullExperimentalData.tth.length < (2 * AVERAGE_WINDOW_HALF_WIDTH + 1)) {
        showToast("Not enough data loaded for background averaging.", "error");
        return;
    }

    const minTth = parseFloat(controls.tthMinSlider.value);
    const maxTth = parseFloat(controls.tthMaxSlider.value);

    if (numPoints < 2) numPoints = 2;

    const internalPointsToFind = numPoints - 2;
    const newSplinePoints = [];
    const tth = fullExperimentalData.tth;
    const intensity = fullExperimentalData.intensity;
    const totalDataPoints = tth.length;

    //   Find Edge Points, problem: how to deal with the edge
    const firstPoint = findClosestExperimentalPoint(minTth);
    const lastPoint = findClosestExperimentalPoint(maxTth);

    if (!firstPoint || !lastPoint) {
        showToast("Could not find edge points for background.", "error");
        return;
    }
    newSplinePoints.push({ index: firstPoint.index, tth: minTth, y: firstPoint.y });
    newSplinePoints.push({ index: lastPoint.index, tth: maxTth, y: lastPoint.y });

    //   Find Internal Points (with Refined Minimum and Averaging)  
    if (internalPointsToFind > 0) {
        let startIndex = fullExperimentalData.tth.findIndex(t => t >= minTth);
        let endIndex = fullExperimentalData.tth.findIndex(t => t > maxTth);
        if (startIndex === -1) startIndex = 0;
        if (endIndex === -1) endIndex = totalDataPoints;

        const relevantLength = endIndex - startIndex;
        if (relevantLength > internalPointsToFind + (2 * AVERAGE_WINDOW_HALF_WIDTH)) {
            const chunkSize = Math.floor(relevantLength / internalPointsToFind);

            for (let i = 0; i < internalPointsToFind; i++) {
                const chunkStartIdx = startIndex + i * chunkSize;
                const chunkEndIdx = (i === internalPointsToFind - 1) ? endIndex - 1 : startIndex + (i + 1) * chunkSize;

                if (chunkStartIdx >= chunkEndIdx) continue;

                //   Step 1: Find initial minimum index in the chunk  
                let initialMinVal = Infinity;
                let initialMinIndex = -1;
                for (let j = chunkStartIdx; j < chunkEndIdx; j++) {
                    if (j === firstPoint.index || j === lastPoint.index) continue;
                    if (intensity[j] < initialMinVal) {
                        initialMinVal = intensity[j];
                        initialMinIndex = j;
                    }
                }

                if (initialMinIndex === -1) continue; // No valid minimum found in chunk

                //   Step 2: Define and search within the refined window  
                const chunkStartTth = tth[chunkStartIdx];
                const chunkEndTth = tth[Math.max(chunkStartIdx, chunkEndIdx - 1)]; // Ensure valid index
                const chunkWidthTth = chunkEndTth - chunkStartTth;
                const refineRadiusTth = REFINE_WINDOW_FRACTION * chunkWidthTth;
                const refineCenterTth = tth[initialMinIndex];
                const refineMinTth = refineCenterTth - refineRadiusTth;
                const refineMaxTth = refineCenterTth + refineRadiusTth;

                let finalMinVal = Infinity;
                let finalMinIndex = initialMinIndex; // Default to initial minimum

                // Iterate through the original chunk indices
                for (let j = chunkStartIdx; j < chunkEndIdx; j++) {
                     if (j === firstPoint.index || j === lastPoint.index) continue;
                    if (tth[j] >= refineMinTth && tth[j] <= refineMaxTth) {
                        if (intensity[j] < finalMinVal) {
                            finalMinVal = intensity[j];
                            finalMinIndex = j;
                        }
                    }
                }
            

                // Check if this final minimum point already exists 
                const exists = newSplinePoints.some(p => p.index === finalMinIndex);
                if (!exists) {
                    //   Step 3: Calculate Average Intensity 
                    let sumY = 0;
                    let countY = 0;
                    const windowStart = Math.max(0, finalMinIndex - AVERAGE_WINDOW_HALF_WIDTH);
                    const windowEnd = Math.min(totalDataPoints, finalMinIndex + AVERAGE_WINDOW_HALF_WIDTH + 1);

                    for (let k = windowStart; k < windowEnd; k++) {
                        sumY += intensity[k];
                        countY++;
                    }
                    const averageY = (countY > 0) ? sumY / countY : intensity[finalMinIndex];
                    //   End Step 3  

                    newSplinePoints.push({
                        index: finalMinIndex,     
                        tth: tth[finalMinIndex],  
                        y: averageY                
                    });
                }
            }
        } else {
             console.warn("Not enough data points in the selected range to find distinct internal background points with averaging.");
        }
    }

    // Overwrite the global list
    backgroundAnchors = newSplinePoints.sort((a, b) => a.tth - b.tth);

    renderSplinePointList();
    updateBackgroundForPreview();
    updateSplinePointsOnChart();
}


function updateSplinePointsOnChart() {
    if (!mainChart) return;

    const anchorDataset = mainChart.data.datasets.find(d => d.label === 'Spline Points'); 
    if (anchorDataset) {
        // Map the anchors to the {x, y} format required by Chart.js 
        anchorDataset.data = backgroundAnchors.map(anchor => ({
            x: anchor.tth,
            y: anchor.y
        }));
       mainChart.update('none'); // Redraw the chart without animation, too slwo and useless
    }
}
        // function to find the nearest data point to a click
        function findClosestExperimentalPoint(targetTth) {
            if (!fullExperimentalData || fullExperimentalData.tth.length === 0) {
                return null;
            }

            // Binary search to the neighbourhood, then compare the two
            // candidates -- the axis is sorted, so the old full linear scan
            // over every data point was unnecessary.
            const arr = fullExperimentalData.tth;
            const hi = lowerBound(arr, targetTth);
            let bestIdx = -1, bestDiff = Infinity;
            for (const cand of [hi - 1, hi]) {
                if (cand < 0 || cand >= arr.length) continue;
                const d = Math.abs(arr[cand] - targetTth);
                if (d < bestDiff) { bestDiff = d; bestIdx = cand; }
            }
            if (bestIdx === -1) return null;
            return { index: bestIdx, tth: arr[bestIdx], y: fullExperimentalData.intensity[bestIdx] };
        }


        function renderSplinePointList() {
    const listContainer = document.getElementById('spline-points-list');
    listContainer.innerHTML = '';

    if (backgroundAnchors.length === 0) {
        listContainer.innerHTML = '<p class="control-label" style="text-align: center; padding-right: 8px;">No spline points defined.</p>';
    } else {
        backgroundAnchors.forEach((anchor, index) => {
            const item = document.createElement('div');
            item.className = 'anchor-point-item';

            // Check if it's an edge point
            const isEdgePoint = (index === 0 || index === backgroundAnchors.length - 1);
            const tthDisabled = isEdgePoint ? 'disabled' : '';
            const removeDisabled = isEdgePoint ? 'disabled' : '';
            const removeStyle = isEdgePoint ? 'opacity: 0.3; cursor: not-allowed;' : ''; 

            item.innerHTML = `
                <input type="number" class="control-input spline-tth-input" data-index="${index}" value="${anchor.tth.toFixed(4)}" step="0.01" title="2-theta" ${tthDisabled}>
                <input type="number" class="control-input spline-y-input" data-index="${index}" value="${anchor.y.toFixed(1)}" step="1" title="Intensity">
                <button class="anchor-remove-btn" data-index="${index}" title="${isEdgePoint ? 'Cannot remove edge point' : 'Remove point'}" ${removeDisabled} style="${removeStyle}">&times;</button>
            `;
            listContainer.appendChild(item);
        });
    }
    updateSplinePointsOnChart();
}

function handleSplinePointListInteraction(e) {
    let listNeedsUpdate = false;
    if (e.target.classList.contains('anchor-remove-btn')) {
        const index = parseInt(e.target.dataset.index, 10);
        // Prevent deleting first or last point (index 0 and length-1)
        if (!isNaN(index) && index > 0 && index < backgroundAnchors.length - 1) {
            backgroundAnchors.splice(index, 1);
            listNeedsUpdate = true; 
        } else if (!isNaN(index) && (index === 0 || index === backgroundAnchors.length - 1)) {
             showToast("Cannot remove edge points.", "error"); // Inform user
        }
    }

    if (e.type === 'input' && e.target.classList.contains('spline-tth-input')) {
        const index = parseInt(e.target.dataset.index, 10);
        const value = parseFloat(e.target.value);
        // Prevent editing tth of first or last point
        if (!isNaN(index) && index > 0 && index < backgroundAnchors.length - 1) {
            if (!isNaN(value) && backgroundAnchors[index]) {
                backgroundAnchors[index].tth = value;
                backgroundAnchors.sort((a, b) => a.tth - b.tth);
                listNeedsUpdate = true;
            }
        } else if (!isNaN(index) && backgroundAnchors[index]) {
             const currentMinTth = parseFloat(controls.tthMinSlider.value);
             const currentMaxTth = parseFloat(controls.tthMaxSlider.value);
             e.target.value = (index === 0 ? currentMinTth : currentMaxTth).toFixed(4);
             showToast("Cannot edit edge point 2θ.", "error");
        }
    }

    if (e.type === 'input' && e.target.classList.contains('spline-y-input')) {
        const index = parseInt(e.target.dataset.index, 10);
        const value = parseFloat(e.target.value);
        if (!isNaN(index) && !isNaN(value) && backgroundAnchors[index]) {
            backgroundAnchors[index].y = value;
            // No need to touch tth here, it's either fixed (edge) or editable via the other input
            updateBackgroundForPreview(); // Recalculate and redraw the background curve
            updateSplinePointsOnChart(); 
        }
    }

    if (listNeedsUpdate) {
        renderSplinePointList(); // Rebuilds the HTML list elements
        updateBackgroundForPreview(); // Ensure background curve reflects potential Tth changes
    }
}