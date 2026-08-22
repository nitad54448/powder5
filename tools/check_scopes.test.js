#!/usr/bin/env node
/* Regression test for tools/check_scopes.js.
 *
 * Re-introduces each of the four duplicate-name bugs that were actually found
 * in this codebase, into a scratch copy of the tree, and asserts the checker
 * reports it and exits non-zero. Then asserts the clean tree passes.
 */
'use strict';
const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');

const SRC = process.argv[2] || path.join(__dirname, '..');
const TMP = path.join(require('os').tmpdir(), 'check_scopes_regress');

/* A RECURSIVE copy, because the real tree keeps its modules in js/.
 *
 * The first version of this copied only the files sitting at the top level and
 * skipped directories outright. On a flat scratch directory that happens to be
 * the whole tree, so it looked correct. On the actual repository it copied
 * powder5.html and nothing else -- and then the checker, finding no scripts at
 * all, correctly reported no duplicates, and case 0 PASSED. A test that passes
 * because it tested nothing is worse than one that fails, which is why case 0
 * now also asserts that the scripts were found. */
const SKIP_DIRS = new Set(['.git', 'node_modules', '.vscode', '.idea']);

function copyTree(from, to) {
    fs.mkdirSync(to, { recursive: true });
    for (const entry of fs.readdirSync(from, { withFileTypes: true })) {
        if (entry.isDirectory()) {
            if (SKIP_DIRS.has(entry.name)) continue;
            copyTree(path.join(from, entry.name), path.join(to, entry.name));
        } else if (entry.isFile()) {
            fs.copyFileSync(path.join(from, entry.name), path.join(to, entry.name));
        }
    }
}

function freshCopy() {
    fs.rmSync(TMP, { recursive: true, force: true });
    copyTree(SRC, TMP);
    // The tool must exist in the copy even if SRC/tools was excluded somehow.
    const dst = path.join(TMP, 'tools');
    fs.mkdirSync(dst, { recursive: true });
    fs.copyFileSync(path.join(SRC, 'tools', 'check_scopes.js'),
                    path.join(dst, 'check_scopes.js'));
}

/* Locate a file by relative path OR by bare basename, so the same test reads
 * the same whether the tree is flat or has its modules under js/. */
function locate(rel) {
    const direct = path.join(TMP, rel);
    if (fs.existsSync(direct)) return direct;
    const base = path.basename(rel);
    const stack = [TMP];
    while (stack.length) {
        const dir = stack.pop();
        for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
            const full = path.join(dir, e.name);
            if (e.isDirectory()) { if (!SKIP_DIRS.has(e.name)) stack.push(full); }
            else if (e.name === base) return full;
        }
    }
    throw new Error(`test harness: ${rel} is not anywhere in the copied tree`);
}

function run() {
    try {
        const out = execFileSync('node', [path.join(TMP, 'tools', 'check_scopes.js'), TMP],
                                 { encoding: 'utf8' });
        return { code: 0, out };
    } catch (e) {
        return { code: e.status, out: (e.stdout || '') + (e.stderr || '') };
    }
}

function edit(file, fn) {
    const p = locate(file);
    fs.writeFileSync(p, fn(fs.readFileSync(p, 'utf8')));
}

let pass = 0, fail = 0;
function check(label, { expectCode, mustContain = [] }) {
    const r = run();
    const problems = [];
    if (r.code !== expectCode) problems.push(`exit ${r.code}, expected ${expectCode}`);
    for (const m of mustContain) if (!r.out.includes(m)) problems.push(`missing "${m}"`);
    if (problems.length) { fail++; console.log(`FAIL  ${label}`); problems.forEach(p => console.log('        ' + p)); }
    else { pass++; console.log(`ok    ${label}`); }
}

/* ---- 0. the clean tree must pass, AND must actually have looked ------- */
freshCopy();
{
    const r = run();
    const m = /page scripts : (\d+)/.exec(r.out);
    const n = m ? Number(m[1]) : 0;
    if (n < 5) {
        console.log('FAIL  the harness did not copy the tree properly');
        console.log(`        the checker saw ${n} page script(s); expected the whole set.`);
        console.log('        Every case below would pass vacuously, so stopping here.');
        console.log(r.out.split('\n').slice(0, 6).map(l => '        ' + l).join('\n'));
        process.exit(1);
    }
    console.log(`      (checker sees ${n} page scripts, ` +
                `${/closure names : (\d+)/.exec(r.out)?.[1] || '?'} closure names)`);
}
check('clean tree reports nothing and exits 0', { expectCode: 0, mustContain: ['OK  no duplicate'] });

/* ---- 1. cfErfc: closure copy shadowing a module copy, with drift ------- */
freshCopy();
edit('powder5.html', s => {
    const bad = `function cfErfc(x) {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const ans = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (-9.9 +
        t * (0.09678418 + t * (-1.18628806)))));
    return x >= 0 ? ans : 2 - ans;
}
`;
    const i = s.indexOf('// The five helpers that used to sit here');
    return s.slice(0, i) + bad + s.slice(i);
});
check('cfErfc duplicated into the closure, body differs', {
    expectCode: 1,
    mustContain: ['cfErfc', 'shadows a global', 'STATUS: the copies DIFFER']
});

/* ---- 2. runFit: module copy shadowed by the closure copy --------------- */
freshCopy();
edit('refinement_controller.js', s => s + `
function runFit(refinementMode) {
    // the stale copy that the buttons never reach
    return null;
}
`);
check('runFit re-added to refinement_controller.js', {
    expectCode: 1,
    mustContain: ['runFit', 'refinement_controller.js', 'THIS is what runs']
});

/* ---- 3. drawTheoreticalPreview: same, in charting.js ------------------- */
freshCopy();
edit('charting.js', s => s + `
function drawTheoreticalPreview() {
    if (hasExperimentalData()) return;
    return 'stale';
}
`);
check('drawTheoreticalPreview re-added to charting.js', {
    expectCode: 1,
    mustContain: ['drawTheoreticalPreview', 'charting.js']
});

/* ---- 4. two module files defining the same function -------------------- */
freshCopy();
edit('data_io.js',   s => s + '\nfunction ioSharedHelper(x) { return x + 1; }\n');
edit('reporting.js', s => s + '\nfunction ioSharedHelper(x) { return x + 2; }\n');
check('same function in two module files, bodies differ', {
    expectCode: 1,
    mustContain: ['ioSharedHelper', 'this one wins (loaded last)', 'STATUS: the copies DIFFER']
});

/* ---- 5. identical duplicate is still reported -------------------------- */
freshCopy();
edit('data_io.js',   s => s + '\nfunction ioSharedHelper(x) { return x + 1; }\n');
edit('reporting.js', s => s + '\nfunction ioSharedHelper(x) { return x + 1; }\n');
check('identical duplicate is reported as not-yet-wrong', {
    expectCode: 1,
    mustContain: ['ioSharedHelper', 'identical today']
});

/* ---- 6. duplicate inside a worker's importScripts set ------------------ */
freshCopy();
edit('crystal.js', s => s + '\nfunction swPackReflections(a) { return a; }\n');
check('duplicate across a worker and one of its imports', {
    expectCode: 1,
    mustContain: ['swPackReflections']
});

/* ---- 7. a nested local of the same name must NOT be reported ----------- */
freshCopy();
edit('data_io.js', s => s + `
function ioWrapper() {
    function cellVolume(c) { return 0; }   // a LOCAL, legitimately shadowing
    return cellVolume({});
}
`);
check('a nested local shadowing a global is not reported', { expectCode: 0 });

/* ---- 8. names inside an IIFE must NOT be reported ---------------------- */
freshCopy();
edit('data_io.js', s => s + `
const IO_MODULE = (() => {
    function cellVolume(c) { return 1; }   // scoped to the IIFE
    function escapeHtml(s) { return s; }
    return { cellVolume, escapeHtml };
})();
`);
check('names inside an IIFE are not reported', { expectCode: 0 });

/* ---- 9. a regex literal containing braces must not break brace depth --- */
freshCopy();
edit('data_io.js', s => s + `
function ioRegexTrap(t) {
    // If the stripper mistakes this for a division the brace count desyncs
    // and every declaration after it is misclassified.
    const re = /\\d{2,3}[{}]/g;
    const s2 = String(t).replace(re, '');
    return s2;
}
function ioAfterTheTrap() { return 1; }
`);
check('a regex literal with braces does not desync the depth count', { expectCode: 0 });

console.log('');
console.log(`${pass} passed, ${fail} failed`);
fs.rmSync(TMP, { recursive: true, force: true });
process.exit(fail ? 1 : 0);
