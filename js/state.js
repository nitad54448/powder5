// --- GLOBAL STATE ---
let controls = {};
let profileParamCache = { "simple_pvoigt": {}, "tch_aniso": {}, "split_pvoigt": {} };
let currentProfile = "simple_pvoigt";
let isFitting = false;
let lastGeneratedHklList = [];
let backgroundAnchors = [];
let masterHklList = [];
let hklIndexCache = {};
let fullExperimentalData = { tth: new Float64Array(), intensity: new Float64Array() };
let calculatedBackground = new Float64Array();
let fitResults = null;
const getLiveFitResults = () => fitResults;
let sgUsedForFit = null;
let lastFitResultsCache = null;
let mainChart;
let refinementWorker = null;
let workerWorkingData = null;
let currentSG = null;

// ---------------------------------------------------------------------------
//  STRUCTURE-SOLUTION RUN STATE
//
//  MOVED HERE (v151) from powder5.html, where these sat several hundred lines
//  BELOW the DOMContentLoaded block that initialises the controls. `let` is
//  hoisted but stays in its temporal dead zone until the declaration is
//  evaluated, so calling anything that reads them during setup threw
//  "Cannot access 'wyckoffRunning' before initialization". The existing code
//  only escaped it because it read cfRunning inside deferred callbacks; the
//  Wyckoff hint needs its value immediately, to decide the initial button
//  state.
//
//  state.js is loaded before every other script, so a declaration here cannot
//  be in the dead zone for any of them.
// ---------------------------------------------------------------------------
let chargeFlipWorker = null;
// The Wyckoff search has its OWN worker (js/wyckoff_worker.js). It is not a
// stage of charge flipping and does not consume its map, so it does not share
// its thread either: either run can be cancelled, crashed or terminated
// without touching the other. Declared here for the same reason as the rest of
// this block -- powder5.html reads it during setup.
let wyckoffWorker = null;
let cfRunning = false;
let cfCancelled = false;
let wyckoffRunning = false;
let wyckoffCancelled = false;
let cfLastResult = null;      // { map, gridSize, peaks, cell, ... }
let cfStructure = null;       // the solved structure, from CF or a Wyckoff search
let cfMapInputs = null;       // the observation set both solution methods share
let cfIntensityProvenance = null;   // how those observations were prepared
// Set by runChargeFlipping so the Stop button can reach that run's finish
// path. `best` lives in the run's closure, so without this the button could
// only terminate the worker -- discarding every trial that had completed.
let cfStopHandler = null;
let workingData = { tth: [], intensity: [], weights: [], startIndex: 0, lastRawDifference: [], isValid: false };