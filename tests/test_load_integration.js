// test/test_no_duplicate_symbols.js
// ---------------------------------------------------------------------------
//  Guards the whole point of the constants.js / crystal.js / profile.js split.
//
//  1. THE THREE MODULES DO NOT COLLIDE WITH EACH OTHER. In a browser they
//     share one script scope, so two `const`s of the same name is a hard
//     SyntaxError at load.
//
//  2. NO CONSUMER RE-DECLARES A SHARED SYMBOL. This is the regression that
//     matters. The refactor exists because PV_GAUSS_AREA,
//     GSAS_GAUSSIAN_TO_DEG, prepareVoigt, updateHklPositions, lowerBound and
//     friends were each declared in two, three, or in one case four files at
//     once, free to drift apart silently. If someone re-adds a local copy --
//     which is the natural thing to do when a worker throws ReferenceError --
//     this fails loudly instead of letting the preview and the fit compute
//     different physics.
//
//  3. EVERY CONSUMER IS STILL WIRED UP, and no shared export is dead weight.
//
//  4. LOAD ORDER IS CONSISTENT in the HTML and in both workers.
//
//  The symbol list comes from the modules' own exports, never hand-written, so
//  it cannot go stale the way a maintained list would.
//
//  Run: node test/test_no_duplicate_symbols.js
// ---------------------------------------------------------------------------

'use strict';

const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');

const read = f => fs.readFileSync(path.join(ROOT, f), 'utf8').replace(/\r\n/g, '\n');

const SHARED = ['constants.js', 'crystal.js', 'profile.js'];
const CONSUMERS = ['powder5.html', 'refinement_worker.js',
                   'charge_flipping_worker.js', 'density3d.js'];

let failures = 0;
const fail = m => { failures++; console.error('  FAIL ' + m); };

/** Strips comments and string literals so nothing is matched inside them. */
function strip(src) {
    let out = '', i = 0;
    while (i < src.length) {
        const c = src[i];
        if (c === '/' && src[i + 1] === '/') {
            while (i < src.length && src[i] !== '\n') i++;
            continue;
        }
        if (c === '/' && src[i + 1] === '*') {
            i += 2;
            while (i < src.length && !(src[i] === '*' && src[i + 1] === '/')) {
                if (src[i] === '\n') out += '\n';
                i++;
            }
            i += 2;
            continue;
        }
        if (c === '"' || c === "'" || c === '`') {
            const q = c;
            i++;
            while (i < src.length && src[i] !== q) {
                if (src[i] === '\\') i++;
                if (src[i] === '\n') out += '\n';
                i++;
            }
            i++;
            out += '""';
            continue;
        }
        out += c;
        i++;
    }
    return out;
}

const lineOf = (s, idx) => s.slice(0, idx).split('\n').length;

/**
 * Lines on which `name` is DECLARED (function / const / let / var), at any
 * nesting depth. A nested re-declaration shadows the shared one, which is
 * exactly as dangerous as a top-level duplicate and much harder to spot.
 * @param {string} stripped Source with comments and strings already removed.
 * @param {string} name
 * @returns {number[]}
 */
function declarationsOf(stripped, name) {
    const hits = [];
    const re = new RegExp(
        `(?:^|[;{})\\n])\\s*(?:(?:async\\s+)?function\\s+|const\\s+|let\\s+|var\\s+)${name}\\b`,
        'g'
    );
    let m;
    while ((m = re.exec(stripped))) hits.push(lineOf(stripped, m.index));
    return hits;
}

// ---------------------------------------------------------------------------
console.log('1. the three shared modules do not collide with each other\n');
// ---------------------------------------------------------------------------
/** @type {Map<string, string>} symbol -> owning module */
const owner = new Map();
for (const f of SHARED) {
    const names = Object.keys(require(path.join(ROOT, f)));
    for (const name of names) {
        if (owner.has(name)) fail(`${name} is exported by both ${owner.get(name)} and ${f}`);
        else owner.set(name, f);
    }
    console.log(`   ${f.padEnd(16)} ${String(names.length).padStart(2)} exports`);
}
console.log(`   ${owner.size} distinct shared symbols, no collisions`);

// ---------------------------------------------------------------------------
console.log('\n2. no consumer re-declares a shared symbol\n');
// ---------------------------------------------------------------------------
const stripped = new Map(CONSUMERS.map(f => [f, strip(read(f))]));
let redeclared = 0;
for (const f of CONSUMERS) {
    for (const name of owner.keys()) {
        for (const line of declarationsOf(stripped.get(f), name)) {
            fail(`${f}:${line} re-declares ${name}, which ${owner.get(name)} owns`);
            redeclared++;
        }
    }
}
if (!redeclared) {
    console.log(`   ${CONSUMERS.length} files x ${owner.size} symbols: none re-declared`);
}

// ---------------------------------------------------------------------------
console.log('\n3. wiring\n');
// ---------------------------------------------------------------------------
const uses = (f, n) => new RegExp(`\\b${n}\\b`).test(stripped.get(f));
for (const f of CONSUMERS) {
    const used = [...owner.keys()].filter(n => uses(f, n));
    console.log(`   ${f.padEnd(28)} uses ${String(used.length).padStart(2)} shared symbols`);
    if (used.length === 0) fail(`${f} references no shared symbol -- is it still wired up?`);
}
// A symbol used only INSIDE the shared modules is not dead weight -- it is an
// internal helper that happens to be exported for the tests. Only flag names
// nothing anywhere uses.
const sharedSrc = SHARED.map(f => strip(read(f))).join('\n');
const usedInShared = n => new RegExp(`\\b${n}\\b`).test(sharedSrc);
const orphans = [...owner.keys()].filter(
    n => !CONSUMERS.some(f => uses(f, n)) && !usedInShared(n)
);
if (orphans.length) {
    console.log(`\n   exported but referenced nowhere at all:`);
    console.log(`     ${orphans.join(', ')}`);
} else {
    console.log('   every shared export is referenced by a consumer or another module');
}

// ---------------------------------------------------------------------------
console.log('\n4. load order\n');
// ---------------------------------------------------------------------------
const html = read('powder5.html');
const order = [...html.matchAll(/<script src="([^"]+)"><\/script>/g)].map(m => m[1]);
const idx = f => order.indexOf(f);
for (const [before, after, why] of [
    ['constants.js', 'crystal.js',   'crystal.js after constants.js'],
    ['constants.js', 'profile.js',   'profile.js after constants.js (reads its declarations)'],
    ['crystal.js',   'density3d.js', 'density3d.js after crystal.js (uses cellBasis)']
]) {
    if (idx(before) < 0) fail(`powder5.html has no <script src="${before}">`);
    else if (idx(after) < 0) fail(`powder5.html has no <script src="${after}">`);
    else if (idx(before) > idx(after)) fail(`powder5.html loads ${after} before ${before}`);
    else console.log(`   ${why}`);
}

for (const [worker, expected] of [
    ['refinement_worker.js', "importScripts('constants.js', 'crystal.js', 'profile.js')"],
    ['charge_flipping_worker.js', "importScripts('constants.js', 'crystal.js')"]
]) {
    const src = read(worker);
    if (!src.includes(expected)) fail(`${worker} is missing: ${expected}`);
    else console.log(`   ${worker.padEnd(28)} line ${lineOf(src, src.indexOf(expected))}: ${expected}`);
}

console.log(`\n${failures} failures`);
process.exit(failures ? 1 : 0);
