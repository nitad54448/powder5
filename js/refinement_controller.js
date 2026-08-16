// FIX: previously there was no timeout and no way to recover. If the
// worker wedged (or died without firing onerror) the UI sat on
// "Refining..." forever, and any onerror permanently nulled the worker
// so every later fit was refused. The worker is now built by a factory
// that can be called again, and a watchdog re-armed on every progress
// message terminates a worker that has gone quiet.
const FIT_WATCHDOG_MS = 120000;   // 2 min with no progress => assume wedged
let fitWatchdog = null;

function disarmFitWatchdog() {
    if (fitWatchdog !== null) { clearTimeout(fitWatchdog); fitWatchdog = null; }
}

function armFitWatchdog() {
    disarmFitWatchdog();
    fitWatchdog = setTimeout(() => {
        fitWatchdog = null;
        if (!isFitting) return;
        console.error("Refinement watchdog fired: no progress from the worker.");
        try { if (refinementWorker) refinementWorker.terminate(); } catch (err) { /* already gone */ }
        refinementWorker = null;
        window.setUIState(false);
        controls.progressBar.style.width = '0%';
        showToast("Refinement stopped responding and was cancelled.", "error");
        createRefinementWorker();   // rebuild so the next fit can still run
    }, FIT_WATCHDOG_MS);
}

function createRefinementWorker() {
    try {
        refinementWorker = new Worker('js/refinement_worker.js');
        console.log("Refinement worker created.");

        //   Worker Message Handler
        refinementWorker.onmessage = function(e) {
            const { type, value, message, results } = e.data;

            if (type === 'warning') {
                console.warn("Worker warning:", message);
                showToast(message, 'error');
                return;
            }

            if (type === 'progress') {
                armFitWatchdog();
                controls.progressBar.style.transition = 'width 0.1s linear';
                controls.progressBar.style.width = `${Math.min(100, value * 100)}%`;
            } else if (type === 'result') {
                disarmFitWatchdog();
                fitResults = results; 
                if (fitResults && fitResults.params) window.enforceSymmetryConstraints(fitResults.params);
                lastFitResultsCache = JSON.parse(JSON.stringify(fitResults)); 

                try {
                    const histMode = (fitResults.refinementMode === 'pawley') ? 'pawley' : 'lebail';
                    window.rpSnapshotFit(histMode, fitResults);
                } catch (histErr) { console.warn('history snapshot failed:', histErr); }

                const displayParams = fitResults.params;
                const displayHklList = fitResults.hklList;
                const displayScaleFactor = fitResults.stats.scaleFactor;

                if (!workerWorkingData) {
                     console.error("Main thread doesn't have workerWorkingData to display results!");
                     showToast("Error displaying results: Data mismatch.", "error");
                     window.setUIState(false);
                     return;
                }

                 try {
                     const finalBackgroundDisplay = calculateTotalBackground(workerWorkingData.tth, displayParams, backgroundAnchors);
                     const finalNetPatternDisplay = calculatePattern(workerWorkingData.tth, displayHklList, displayParams);
                     window.updateUI(
                         displayParams,
                         fitResults.stats,
                         finalNetPatternDisplay, 
                         finalBackgroundDisplay, 
                         displayScaleFactor,
                         displayHklList
                     );
                 } catch (displayError) {
                      console.error("Error updating UI with worker results:", displayError);
                      showToast(`Error displaying results: ${displayError.message}`, "error");
                 }

                window.setUIState(false); 
                controls.progressBar.style.transition = 'width 0.3s ease';
                controls.progressBar.style.width = '100%';
            } else if (type === 'error') {
                disarmFitWatchdog();
                console.error("Worker reported error:", message);
                if (window.showToast) window.showToast(`Refinement Error: ${message}`, 'error');
                window.setUIState(false); 
                controls.progressBar.style.width = '0%'; 
                if (message && message.indexOf('initialization failed') !== -1) {
                    refinementWorker = null;
                }
            }
        };

        refinementWorker.onerror = function(e) {
            console.error("Error initializing or communicating with worker:", e.message, e);
            if (window.showToast) window.showToast(`Worker Error: ${e.message}`, 'error');
            disarmFitWatchdog();
            window.setUIState(false);
            try { if (refinementWorker) refinementWorker.terminate(); } catch (err) { /* already gone */ }
            refinementWorker = null;
        };

    } catch (e) {
        console.error("Failed to create Web Worker:", e);
        if (window.showToast) {
    window.showToast("Error initializing worker...", "error");
} else {
    console.error("Error initializing worker...");
}
        refinementWorker = null;
    }
}

createRefinementWorker();

// master HKL, 1.0.3 et ensuite
// This function now uses the cache to get raw indices and then calculates positions.
function generateMasterHklList() {
    if (fullExperimentalData.tth.length === 0) {
        masterHklList = [];
        updatePreviewPattern();
        return;
    }

    const selectedSg = currentSG;
    if (!selectedSg) {
        masterHklList = [];
        console.error("Cannot generate HKL list: no space group resolved.");
        updatePreviewPattern();
        return;
    }

    const maxTth = fullExperimentalData.tth.reduce((m, v) => (v > m ? v : m), -Infinity) + 2.0; 
    const params = getAllParams();

    const rawHklIndices = generateAndCacheHklIndices(selectedSg, maxTth, params);
    let workingHklList = JSON.parse(JSON.stringify(rawHklIndices));
    updateHklPositions(workingHklList, params, selectedSg.system);

    masterHklList = workingHklList
        .filter(peak => peak.tth !== null && peak.tth <= maxTth)
        .sort((a, b) => a.tth - b.tth);

    masterHklList.forEach(peak => {
        peak.intensity = peak.intensity || 1000;
    });

    updatePreviewPattern();
}

function runFit(refinementMode) {
    if (isFitting) {
         console.warn("Fit already in progress.");
         return;
    }
     if (!refinementWorker) {
         showToast("Refinement worker is not available. Cannot start fit.", "error");
         return; 
     }

    if (!currentSG) {
        showToast('Please select a space group first ("Space Group »" button).', "error");
        return;
    }

    isFitting = true;
    window.setUIState(true);
    updateWorkingData(); 

    if (masterHklList.length === 0) {
         showToast("No reflections in range. Check lattice parameters and 2-theta range.", "error");
         isFitting = false; window.setUIState(false); return;
    }

    const currentParams = getAllParams();
    if (currentParams.__invalidFields && currentParams.__invalidFields.length > 0) {
        showToast(
            `Cannot start: these fields are empty or not a number - ${currentParams.__invalidFields.join(', ')}.`,
            "error"
        );
        isFitting = false; window.setUIState(false); return;
    }
    delete currentParams.__invalidFields;

    const fitFlags = window.getFitFlags();
    const selectedSg = currentSG;
    const selectedSgNumber = selectedSg.number;
    sgUsedForFit = selectedSg;

    const system = selectedSg.system; 
    const maxIterations = parseInt(controls.iterationsSlider.value);
    const algorithm = controls.algorithmSelect.value;

    if (!workingData.isValid || workingData.tth.length === 0) {
        showToast("No data in the selected 2-theta range. Aborting fit.", "error");
        isFitting = false; window.setUIState(false); return;
    }

     workerWorkingData = {
         tth: workingData.tth.slice(), 
         intensity: workingData.intensity.slice(),
         weights: workingData.weights.slice(),
         startIndex: workingData.startIndex
     };

     const minTth = parseFloat(controls.tthMinSlider.value);
     const maxTth = parseFloat(controls.tthMaxSlider.value);
     const fittedHklList = masterHklList.filter(
         peak => peak.tth !== null && peak.tth >= minTth && peak.tth <= maxTth
     );

     const workerPayload = {
         initialParams: currentParams,
         fitFlags: fitFlags,
         workingData: workerWorkingData, 
         masterHklList: JSON.parse(JSON.stringify(fittedHklList)), 
         selectedSgNumber: selectedSgNumber,
         selectedSgQuery:  selectedSg.hall || String(selectedSgNumber),
         system: system,
         maxIterations: maxIterations,
         algorithm: algorithm,
         refinementMode: refinementMode,
         backgroundAnchors: JSON.parse(JSON.stringify(backgroundAnchors)) 
     };

    try {
         refinementWorker.postMessage(workerPayload);
         controls.progressBar.style.width = '0%'; 
         controls.progressBar.style.transition = 'none'; 
    } catch (error) {
        console.error("Error sending message to worker:", error);
        showToast(`Error starting refinement: ${error.message}`, "error");
        isFitting = false;
        window.setUIState(false);
    }
}