// FIX: previously there was no timeout and no way to recover. If the
// worker wedged (or died without firing onerror) the UI sat on
// "Refining..." forever, and any onerror permanently nulled the worker
// so every later fit was refused. The worker is now built by a factory
// that can be called again, and a watchdog re-armed on every progress
// message terminates a worker that has gone quiet.
// ADAPTIVE. A fixed two minutes is simultaneously too short and too long: it
// kills a healthy fit on a large pattern where one iteration legitimately
// takes minutes, and it takes two minutes to notice a wedged worker on a small
// one. The budget is now a multiple of the slowest gap between progress
// messages actually observed in THIS run, clamped at both ends.
const FIT_WATCHDOG_MIN_MS = 30000;
const FIT_WATCHDOG_MAX_MS = 600000;
const FIT_WATCHDOG_FACTOR = 20;
let fitWatchdog = null;
let fitLastProgressAt = 0;
let fitObservedGapMs = 0;

function disarmFitWatchdog() {
    if (fitWatchdog !== null) { clearTimeout(fitWatchdog); fitWatchdog = null; }
}

/** Call before posting a job, so one run's timings do not carry into the next. */
function resetFitWatchdogTimings() {
    fitLastProgressAt = 0;
    fitObservedGapMs = 0;
}

function armFitWatchdog() {
    disarmFitWatchdog();
    const now = Date.now();
    if (fitLastProgressAt) {
        fitObservedGapMs = Math.max(fitObservedGapMs, now - fitLastProgressAt);
    }
    fitLastProgressAt = now;
    const budget = Math.min(FIT_WATCHDOG_MAX_MS,
                            Math.max(FIT_WATCHDOG_MIN_MS, fitObservedGapMs * FIT_WATCHDOG_FACTOR));
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
    }, budget);
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
                // And to the unified footer, so a refinement stays visible
                // after switching to a results or log tab.
                if (typeof AppStatus !== 'undefined') {
                    AppStatus.set('fit', value,
                        message || `${Math.round(Math.min(100, value * 100))}%`);
                }
            } else if (type === 'result') {
                disarmFitWatchdog();
                fitResults = results;
                if (fitResults && fitResults.params) window.enforceSymmetryConstraints(fitResults.params);

                // ONE clone, and J^T J shared BY REFERENCE rather than copied.
                //
                // This used to be a full JSON round-trip here AND another
                // inside rpSnapshotFit, both carrying the matrix: three copies
                // of it per fit, ~13 MB of JSON and ~300 ms at 2500
                // reflections, retained for every run for the life of the tab.
                //
                // But the matrix cannot simply be dropped. generateReportContent
                // computes every ESD from it, and rpRenderRun draws the report
                // from the HISTORY SNAPSHOT, not from the live object -- so a
                // snapshot without J^T J means a report whose sigmas are all
                // N/A. Attaching it by reference gives one copy instead of
                // three and keeps the report intact: a fit result is never
                // mutated once it has been handed over, so there is no
                // aliasing hazard in sharing it.
                //
                // __cov IS dropped: it is a memoised inverse holding
                // Float64Array rows, cheap to recompute and wrong to clone.
                const { JtJ: sharedJtJ, __cov: _cov, ...cloneable } = fitResults;
                lastFitResultsCache = (typeof structuredClone === 'function')
                    ? structuredClone(cloneable)
                    : JSON.parse(JSON.stringify(cloneable));
                if (sharedJtJ) lastFitResultsCache.JtJ = sharedJtJ;

                if (fitResults.fitFlags && fitResults.fitFlags.fitBackground) {
                    backgroundAnchors.forEach((anchor, i) => {
                        if (fitResults.params[`bg_y_${i}`] !== undefined) {
                            anchor.y = fitResults.params[`bg_y_${i}`];
                        }
                    });
                    if (window.renderSplinePointList) window.renderSplinePointList();
                }

                try {
                    const histMode = (fitResults.refinementMode === 'pawley') ? 'pawley' : 'lebail';
                    // The clone above, not the live object: rpSnapshotFit no
                    // longer needs to make its own.
                    window.rpSnapshotFit(histMode, lastFitResultsCache);
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
            // REBUILT, exactly as the watchdog path already does.
            //
            // The comment at the top of this file says the factory exists
            // because "any onerror permanently nulled the worker so every
            // later fit was refused". The watchdog was given the rebuild;
            // THIS path was not -- so the bug the header describes as fixed
            // survived on the commonest route into it. An uncaught throw
            // anywhere inside the worker fires onerror, and from then on
            // runFit answers "Refinement worker is not available" until the
            // page is reloaded: one transient fault, permanently.
            //
            // Guarded, so a worker that cannot be constructed at all does
            // not throw out of an error handler.
            try { createRefinementWorker(); }
            catch (rebuildErr) { console.error("Worker rebuild failed:", rebuildErr); }
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
window.resetFitWatchdogTimings = resetFitWatchdogTimings;

// ===========================================================================
//  generateMasterHklList() AND runFit() USED TO BE DUPLICATED HERE.
//
//  powder5.html defines its own copy of each, inside the DOMContentLoaded
//  callback near the top of its inline script. Everything from there down
//  shares one function scope, the Le Bail and Pawley button handlers live in
//  it, and a function declared in that scope shadows the global of the same
//  name for every caller inside it. These two were therefore never reached.
//
//  They had already drifted apart. The inline copies write the background
//  anchor heights into the parameter set before posting the job; the copies
//  here did not, so they were an older revision preserved by the shadowing
//  rather than a second opinion about anything.
//
//  This is not theoretical. The watchdog fix for the unguarded startup window
//  was applied to the copy here first and had no effect at all, because
//  nothing runs it; it now lives in the inline copy, where it works.
//
//  Removed rather than annotated: nothing outside the shadowing scope
//  referenced either name, neither was placed on `window`, and the whole
//  point of deleting them is that the next person to fix runFit cannot
//  accidentally fix the wrong one. This file now owns the worker lifecycle
//  and the watchdog, and nothing else.
// ===========================================================================
