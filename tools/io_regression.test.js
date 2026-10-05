#!/usr/bin/env node
'use strict';
// Powder 5 v179: import/export and worker-startup regressions, no dependencies.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const quiet = { log() {}, warn() {}, error() {} };
const ctx = vm.createContext({ console: quiet, TextDecoder });
for (const file of ['data_io.js', 'parsers.js']) {
    vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8'), ctx);
}
function call(expression, input) {
    ctx.input = input;
    return vm.runInContext(expression, ctx);
}
let passed = 0;
function test(label, fn) { fn(); passed++; console.log(`ok  ${label}`); }
function close(actual, expected) { assert.ok(Math.abs(actual - expected) < 1e-9, `${actual} != ${expected}`); }

test('sample XRA channel counts and angular ranges', () => {
    for (const [file, count, end] of [['PBSO4.XRA', 6001, 160], ['FAP.XRA', 5753, 130.04]]) {
        const r = call('detectAndParseFile("sample.xra", input)', fs.readFileSync(path.join(root, 'test_files', file), 'utf8'));
        assert.equal(r.intensity.length, count);
        close(r.tth.at(-1), end);
    }
});
test('fractional centidegrees survive XRA round trips', () => {
    for (const [start, step] of [[10, 0.025], [10.003, 0.005], [15, 0.02]]) {
        const b = { tth: Array.from({ length: 13 }, (_, i) => start + step * i), obs: Array.from({ length: 13 }, (_, i) => i + 100) };
        const exported = call('writeXRA(input)', b);
        assert.ok(exported.trimEnd().split('\n').slice(2).every(line => line.length === 80));
        const r = call('parseGsasXraFile(input)', exported);
        assert.equal(r.tth.length, b.tth.length);
        r.tth.forEach((x, i) => close(x, b.tth[i]));
        assert.deepEqual(Array.from(r.intensity), b.obs);
    }
});
test('adjacent full-width STD fields stay separate', () => {
    const r = call('parseGsasXraFile(input)', 'test\nBANK 1 2 1 CONST 1000 2.5 0 0 STD\n1234567887654321\n');
    assert.deepEqual(Array.from(r.intensity), [12345678, 87654321]);
});
test('ESD padding and adjacent fixed-width fields', () => {
    const row = [12345678, 11111111, 87654321, 22222222, 0, 0].map(v => String(v).padStart(8)).join('');
    const r = call('parseGsasEsdFile(input)', `test\nBANK 1 2 1 CONST 1000 2.5 0 0 ESD\n${row}\n`);
    assert.deepEqual(Array.from(r.intensity), [12345678, 87654321]);
});
test('FXYE keeps the supplied axis and declared channel count', () => {
    const r = call('parseGsasEsdFile(input)', 'test\nBANK 1 2 2 CONST 1000 2.5 0 0 FXYE\n1000 12 3\n1002.5 15 4\n1005 0 0\n');
    assert.deepEqual(Array.from(r.intensity), [12, 15]);
    close(r.tth[1], 10.025);
});
test('first GSAS bank is isolated; truncated scans are rejected', () => {
    const r = call('parseGsasXraFile(input)', 'test\nBANK 1 2 1 CONST 1000 2.5 0 0 STD\n12 15\nBANK 2 2 1 CONST 2000 5 0 0 STD\n90 91\n');
    assert.deepEqual(Array.from(r.intensity), [12, 15]);
    assert.throws(() => call('parseGsasXraFile(input)', 'test\nBANK 1 3 1 CONST 1000 2.5 0 0 STD\n12 15\n'), /declares 3 channels/);
});
test('oversized STD intensity fields cannot corrupt exports', () => {
    assert.throws(() => call('writeXRA(input)', { tth: [10, 11], obs: [1e8, 12] }), /eight-character/);
});
test('all pdCIF loops tolerate comments and blank lines', () => {
    const tables = [
        ['_pd_calib_xcoord.nominal_2theta\n_pd_calib_xcoord.actual_2theta', ['10 10.1', '20 20.1']],
        ['_diffrn_radiation_wavelength.id\n_diffrn_radiation_wavelength.value\n_diffrn_radiation_wavelength.wt', ['1 1.5406 1', '2 1.5444 0.5']],
        ['_refln.index_h\n_refln.index_k\n_refln.index_l\n_refln.f_squared_meas', ['1 0 0 100', '0 1 0 200']],
        ['_powder5_bkg_anchor_2theta\n_powder5_bkg_anchor_intensity', ['10 5', '20 7']],
        ['_pd_meas.2theta_scan\n_pd_meas.intensity_total', ['10 100', '11 110', '12 120', '13 130']]
    ];
    const text = 'data_test\n' + tables.map(([tags, rows]) => 'loop_\n# before tags\n\n' + tags.replaceAll('\n', '\n# between tags\n\n') + '\n\n# before rows\n' + rows.join('\n\n# between rows\n') + '\n\n').join('');
    const r = call('normalizeParsedData(parsePdCifFile(input))', text);
    assert.equal(r.tth.length, 4);
    assert.equal(r.hklList.length, 2);
    assert.equal(r.backgroundAnchors.length, 2);
    close(r.zeroShift, -0.1);
    close(r.wavelength2, 1.5444);
    close(r.ratio21, 0.5);
});
test('four CommonJS modules export their version without throwing', () => {
    for (const [file, version] of [['observations', 'OBSERVATIONS'], ['scatterers', 'SCATTERERS'], ['swarm_wyckoff', 'SWARM_WYCKOFF'], ['coord_refine', 'COORD_REFINE']]) {
        assert.equal(typeof require(path.join(root, 'js', `${file}.js`))[`SHARKO_${version}_VERSION`], 'string');
    }
});
function workerContext() {
    const workers = [];
    const c = vm.createContext({ console: quiet, window: {}, controls: {}, isFitting: false, refinementWorker: null,
        setTimeout() { return 1; }, clearTimeout() {}, Worker: function () { this.dead = false; this.terminate = () => { this.dead = true; }; workers.push(this); } });
    vm.runInContext(fs.readFileSync(path.join(root, 'js/refinement_controller.js'), 'utf8'), c);
    return { c, workers };
}
test('early worker failures clean up and stop after three rebuilds', () => {
    const { c, workers } = workerContext();
    for (let i = 0; i < 4; i++) vm.runInContext('refinementWorker.onerror({message: "early failure"})', c);
    assert.equal(workers.length, 4);
    assert.ok(workers.every(w => w.dead));
    assert.equal(vm.runInContext('refinementWorker', c), null);
});
test('initialization-error messages use the same bounded recovery', () => {
    const { c, workers } = workerContext();
    for (let i = 0; i < 4; i++) vm.runInContext('refinementWorker.onmessage({data: {type: "error", message: "Worker initialization failed: test"}})', c);
    assert.equal(workers.length, 4);
    assert.ok(workers.every(w => w.dead));
});
console.log(`${passed} passed, 0 failed`);
