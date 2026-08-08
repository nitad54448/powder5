// test/test_load_integration.js
// ---------------------------------------------------------------------------
//  The other suites test the modules in isolation, via require(). That is NOT
//  how they load in production: the browser concatenates them into one script
//  scope via <script src>, and the workers do the same via importScripts.
//
//  Those two paths have a failure mode require() cannot reproduce: a `const`
//  declared in two files is a SyntaxError at load, and a name referenced
//  before its defining file is evaluated is a ReferenceError -- both of which
//  would take the whole app down on first paint, with the modules themselves
//  testing perfectly clean.
//
//  So: concatenate them the way each host does, evaluate, and use the result.
//
//  Run: node test/test_load_integration.js
// ---------------------------------------------------------------------------

'use strict';

const fs = require('fs');
const path = require('path');
const vm = require('vm');
const ROOT = path.join(__dirname, '..');

let failures = 0, checks = 0;
const check = (cond, msg) => {
    checks++;
    if (!cond) { failures++; console.error('  FAIL ' + msg); }
};

const read = f => fs.readFileSync(path.join(ROOT, f), 'utf8');

/**
 * Loads a list of files into one fresh realm, concatenated, the way a browser
 * script scope or importScripts would. `module` is left undefined so the
 * node-export footers are skipped, exactly as in a browser.
 *
 * NOTE on how names must then be probed. A top-level `const` or `let` in a
 * classic script goes into the realm's global LEXICAL environment, not onto
 * the global object -- so `ctx.PV_G0` is undefined even though PV_G0 is
 * perfectly visible to every other script in the realm. Only `function` and
 * `var` become properties. Use resolves() below rather than property access,
 * or the test reports failures that do not exist.
 */
function loadTogether(files) {
    const src = files.map(f => `/* ==== ${f} ==== */\n${read(f)}`).join('\n');
    const ctx = vm.createContext({ console, Math, Float64Array, Float32Array,
                                   Map, Set, Object, Array, JSON, isFinite, NaN,
                                   Infinity, Number, String });
    // Evaluate as ONE script, so a duplicate const across files is a
    // SyntaxError here just as it would be in the browser.
    vm.runInContext(src + '\n;globalThis.__ok = true;', ctx, { filename: files.join('+') });
    return ctx;
}

/**
 * True if `name` resolves inside the realm -- works for const/let as well as
 * function/var. See the note above loadTogether.
 * @param {object} ctx
 * @param {string} name
 * @returns {boolean}
 */
function resolves(ctx, name) {
    return vm.runInContext(`typeof ${name} !== 'undefined'`, ctx);
}

// ---------------------------------------------------------------------------
console.log('1. browser load order: constants.js + crystal.js + profile.js\n');
// ---------------------------------------------------------------------------
let browser;
try {
    browser = loadTogether(['constants.js', 'crystal.js', 'profile.js']);
    console.log('   three modules evaluate as one script scope: no redeclaration');
} catch (e) {
    failures++;
    console.error('  FAIL the three modules cannot coexist in one scope: ' + e.message);
    console.log(`\n${checks} checks, ${failures} failures`);
    process.exit(1);
}

// The cross-module reference that the load order exists to satisfy:
// profile.js reads constants.js's MIN_PROFILE_FWHM_DEG and softPositive.
for (const n of ['prepareVoigt', 'evalVoigt', 'updateHklPositions', 'lowerBound',
                 'PV_G0', 'PV_GAUSS_AREA', 'MIN_PROFILE_FWHM_DEG',
                 'MONOCLINIC_ANGLE_FOR_AXIS', 'buildInvDsqEvaluator']) {
    check(resolves(browser, n), `${n} should resolve in the shared script scope`);
}
console.log('   functions and const bindings alike resolve across the three files');

// ---------------------------------------------------------------------------
console.log('\n2. profile.js really does read constants.js at run time\n');
// ---------------------------------------------------------------------------
{
    // Change the constant, and the profile floor must follow. This is the
    // cross-module link that a require()-based test cannot check, because
    // there the two modules would have separate bindings.
    const params = { profileType: 'simple_pvoigt', GU: 0, GV: 0, GW: 0, GP: 0, LX: 0, eta: 0.5 };

    vm.runInContext('setMinProfileFwhmFromAxis([10, 10.5, 11.0, 11.5])', browser);
    const coarse = browser.calculateProfileWidths(40, { d: 2 }, params, 'center');

    vm.runInContext('setMinProfileFwhmFromAxis([10, 10.005, 10.010, 10.015])', browser);
    const fine = browser.calculateProfileWidths(40, { d: 2 }, params, 'center');

    console.log(`   0.5 deg step -> floor gives gamma_G = ${coarse.gamma_G.toFixed(4)} deg`);
    console.log(`   0.005 deg step -> floor gives gamma_G = ${fine.gamma_G.toFixed(4)} deg`);
    check(coarse.gamma_G > fine.gamma_G * 10,
        'the width floor must track MIN_PROFILE_FWHM_DEG across the module boundary');
}

// ---------------------------------------------------------------------------
console.log('\n3. refinement_worker load order (importScripts equivalent)\n');
// ---------------------------------------------------------------------------
{
    // The worker's own body needs fetch/postMessage/importScripts, which are
    // not worth stubbing. What matters is that the three shared modules load
    // in the order the worker asks for, and that everything the worker's
    // profile and geometry code references resolves.
    const workerSrc = read('refinement_worker.js');
    const m = workerSrc.match(/importScripts\(([^)]*)\)/);
    check(!!m, 'refinement_worker.js calls importScripts');
    const files = m[1].split(',').map(s => s.trim().replace(/['"]/g, ''));
    console.log(`   importScripts(${files.join(', ')})`);
    const ctx = loadTogether(files);
    for (const n of ['prepareVoigt', 'evalVoigt', 'calculateProfileWidths',
                     'calculatePeakShift', 'updateHklPositions', 'lowerBound',
                     'LEBAIL_FLAT_SEED', 'LEBAIL_REVIVAL_FRACTION',
                     'PT_NUM_REPLICAS', 'PT_MAX_TEMP', 'PT_MIN_TEMP',
                     'PT_SWAP_INTERVAL', 'PT_COST_SCALE_FLOOR',
                     'PT_RESCALE_THRESHOLD', 'FD_BASE_FRACTION', 'JTJ_RANK_TOL',
                     'LM_TRUST_GROWTH', 'GSAS_GAUSSIAN_TO_DEG']) {
        check(resolves(ctx, n), `${n} must resolve after the worker's importScripts`);
    }
    console.log('   every shared name the worker uses resolves');
}

// ---------------------------------------------------------------------------
console.log('\n4. charge_flipping_worker load order\n');
// ---------------------------------------------------------------------------
{
    const src = read('charge_flipping_worker.js');
    const m = src.match(/importScripts\(([^)]*)\)/);
    check(!!m, 'charge_flipping_worker.js calls importScripts');
    const files = m[1].split(',').map(s => s.trim().replace(/['"]/g, ''));
    console.log(`   importScripts(${files.join(', ')})`);
    const ctx = loadTogether(files);
    for (const n of ['metricTensor', 'fracDistance', 'cellVolume']) {
        check(resolves(ctx, n), `${n} must resolve`);
    }
    // It must NOT need profile.js -- it does not load it.
    check(!resolves(ctx, 'prepareVoigt'),
        'the CF worker does not load profile.js, so it must not depend on it');
    const stripped = src.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*/g, '');
    for (const n of ['prepareVoigt', 'evalVoigt', 'calculateProfileWidths']) {
        check(!new RegExp(`\\b${n}\\s*\\(`).test(stripped),
            `CF worker calls ${n} but never loads profile.js`);
    }
    console.log('   CF worker resolves its geometry names and needs nothing from profile.js');
}

// ---------------------------------------------------------------------------
console.log('\n5. the full browser chain, in document order\n');
// ---------------------------------------------------------------------------
{
    const html = read('powder5.html');
    const order = [...html.matchAll(/<script src="([^"]+)"><\/script>/g)]
        .map(m => m[1])
        .filter(f => !f.startsWith('lib/'));
    console.log(`   ${order.join(' -> ')}`);
    // density3d.js is an IIFE needing `global`; give it one and THREE absent.
    const src = order.map(f => `/* ==== ${f} ==== */\n${read(f)}`).join('\n');
    const ctx = vm.createContext({
        console, Math, Float64Array, Float32Array, Uint8Array, Uint32Array,
        Int32Array, Map, Set, Object, Array, JSON, isFinite, NaN, Infinity,
        Number, String, document: undefined, window: undefined
    });
    ctx.globalThis = ctx;
    try {
        vm.runInContext(src, ctx, { filename: 'powder5 head' });
        console.log('   the whole head chain evaluates cleanly');
    } catch (e) {
        failures++;
        console.error('  FAIL head chain: ' + e.message);
    }
    check(ctx.Density3D && typeof ctx.Density3D === 'object',
        'density3d.js should still publish Density3D');
    check(resolves(ctx, 'cellBasis'),
        'cellBasis should come from crystal.js, in scope for density3d.js');
    if (ctx.Density3D && ctx.Density3D._internals) {
        const b = ctx.Density3D._internals.cellBasis({ a: 5, b: 6, c: 7, beta: 100 });
        check(Math.abs(b.cv[0] - 7 * Math.cos(100 * Math.PI / 180)) < 1e-12,
            'density3d gets the crystal.js cellBasis, and it computes correctly');
        console.log('   density3d.js picks up crystal.js cellBasis through the scope chain');
    }
}

console.log(`\n${checks} checks, ${failures} failures`);
process.exit(failures ? 1 : 0);
