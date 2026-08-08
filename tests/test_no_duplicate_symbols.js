// test/test_pt_cost_scale.js
// ---------------------------------------------------------------------------
//  The PT energy scale is the one BEHAVIOURAL change in this batch, so it gets
//  its own suite.
//
//  Background. Both PT acceptance tests divide the cost difference by a shared
//  reference scale, making the temperature ladder relative:
//      accept = exp(-(dChi2 / costScale) / T)
//  The scale used to be frozen at the initial cost for the whole run. A
//  successful global search drops chi-square by one to two orders of
//  magnitude, at which point a frozen scale leaves the ladder 10-100x hotter
//  in relative terms than designed: the cold replicas stop being cold.
//
//  What is asserted here:
//    1. The scale is NOT re-quoted for ordinary jitter (detailed balance
//       between refresh boundaries).
//    2. It IS re-quoted for a genuine change of regime.
//    3. With a frozen scale, the cold end of the ladder demonstrably stops
//       discriminating; with rescaling, it does not.
//    4. The ladder is sane at both ends to begin with -- which is the reason
//       NOT to "fix" costScale to 1.0, a change that would freeze every chain.
//
//  Run: node test/test_pt_cost_scale.js
// ---------------------------------------------------------------------------

'use strict';

const fs = require('fs');
const path = require('path');
const C = require('../constants.js');

let failures = 0, checks = 0;
const check = (cond, msg) => {
    checks++;
    if (!cond) { failures++; console.error('  FAIL ' + msg); }
};

// ---------------------------------------------------------------------------
//  The rescaling rule, kept in step with refinement_worker.js. The assertion
//  at the bottom of this file checks the worker still contains it.
// ---------------------------------------------------------------------------
function makeScaler(initialCost) {
    let costScale = Math.max(C.PT_COST_SCALE_FLOOR, Math.abs(initialCost));
    return {
        get value() { return costScale; },
        maybeRescale(reps) {
            let sum = 0, n = 0;
            for (const r of reps) {
                if (r && isFinite(r.cost)) { sum += Math.abs(r.cost); n++; }
            }
            if (n === 0) return false;
            const proposed = Math.max(C.PT_COST_SCALE_FLOOR, sum / n);
            const ratio = proposed / costScale;
            if (ratio > (1 - C.PT_RESCALE_THRESHOLD) &&
                ratio < 1 / (1 - C.PT_RESCALE_THRESHOLD)) return false;
            costScale = proposed;
            return true;
        }
    };
}

const ladder = () => Array.from({ length: C.PT_NUM_REPLICAS }, (_, i) =>
    C.PT_MAX_TEMP * Math.pow(C.PT_MIN_TEMP / C.PT_MAX_TEMP, i / (C.PT_NUM_REPLICAS - 1 || 1)));
const reps = cost => Array.from({ length: C.PT_NUM_REPLICAS }, () => ({ cost }));
const accept = (dChi2, scale, T) => Math.exp(Math.max(-700, -(dChi2 / scale) / T));

// ---------------------------------------------------------------------------
console.log('1. jitter does not move the scale; a change of regime does\n');
// ---------------------------------------------------------------------------
for (const [cost, shouldMove] of [
    [1.00e6, false], [9.0e5, false], [7.0e5, false], [5.1e5, false],
    [4.9e5, true],  [1.0e5, true],  [1.0e4, true],
    [1.9e6, false], [2.1e6, true]
]) {
    const s = makeScaler(1e6);
    const moved = s.maybeRescale(reps(cost));
    const verb = moved ? 'rescaled' : 'held    ';
    console.log(`   chi2 ${cost.toExponential(1)}  ->  ${verb}  (scale ${s.value.toExponential(1)})`);
    check(moved === shouldMove,
        `chi2 ${cost.toExponential(1)}: expected ${shouldMove ? 'rescale' : 'hold'}`);
}

// ---------------------------------------------------------------------------
console.log('\n2. NaN / empty replica sets never corrupt the scale\n');
// ---------------------------------------------------------------------------
{
    const s = makeScaler(1e6);
    check(s.maybeRescale([]) === false, 'empty replica list must not rescale');
    check(s.maybeRescale([{ cost: NaN }, { cost: Infinity }]) === false,
          'all-non-finite costs must not rescale');
    check(s.value === 1e6, 'scale unchanged after degenerate input');

    // A mix: only the finite ones count.
    s.maybeRescale([{ cost: NaN }, { cost: 1e4 }, { cost: 1e4 }]);
    check(Math.abs(s.value - 1e4) < 1e-9, `mixed input should give 1e4, got ${s.value}`);
    console.log('   empty, all-NaN and mixed inputs handled');

    const z = makeScaler(0);
    check(z.value === C.PT_COST_SCALE_FLOOR, 'zero initial cost must fall back to the floor');
    console.log('   zero initial cost falls back to PT_COST_SCALE_FLOOR');
}

// ---------------------------------------------------------------------------
console.log('\n3. the cold end of the ladder keeps discriminating\n');
// ---------------------------------------------------------------------------
{
    const T = ladder();
    const coldT = T[T.length - 1];
    const initial = 1e6;

    // A search that has succeeded: chi-square down two decades.
    const nowChi2 = 1e4;
    // A move that worsens the fit by 0.1% -- a cold chain must reject it.
    const dChi2 = 0.001 * nowChi2;

    const frozen = accept(dChi2, initial, coldT);
    const s = makeScaler(initial);
    s.maybeRescale(reps(nowChi2));
    const rescaled = accept(dChi2, s.value, coldT);

    console.log(`   after a 100x descent, coldest replica accepts a +0.1% move with p =`);
    console.log(`     frozen scale    ${frozen.toExponential(2)}   <- not cold any more`);
    console.log(`     rescaled        ${rescaled.toExponential(2)}`);

    check(frozen > 0.01, 'the frozen-scale pathology should be visible (p should be large)');
    check(rescaled < 1e-6, `rescaled cold replica should reject: p = ${rescaled}`);
    check(rescaled < frozen * 1e-6, 'rescaling should change the cold end by many orders');
}

// ---------------------------------------------------------------------------
console.log('\n4. the ladder is sane at both ends (why costScale != 1)\n');
// ---------------------------------------------------------------------------
{
    const T = ladder();
    const chi2 = 1e6;
    const scale = chi2;

    // Hot end: a move that DOUBLES chi-square should still be accepted often,
    // or the chain cannot explore.
    const hot = accept(chi2, scale, T[0]);
    // Cold end: a 1% worsening should be rejected outright.
    const cold = accept(0.01 * chi2, scale, T[T.length - 1]);

    console.log(`   hot  replica (T=${T[0].toExponential(1)}): doubling chi2 accepted with p = ${hot.toFixed(3)}`);
    console.log(`   cold replica (T=${T[T.length - 1].toExponential(1)}): +1% chi2 accepted with p = ${cold.toExponential(2)}`);
    check(hot > 0.2 && hot < 0.8, `hot replica should explore freely, p = ${hot}`);
    check(cold < 1e-6, `cold replica should be near-greedy, p = ${cold}`);

    // The review proposed costScale = 1.0 with the raw chi-square difference.
    // Show what that does: every chain, hot included, freezes solid.
    const naiveHot = accept(chi2, 1.0, T[0]);
    console.log(`\n   with costScale = 1 (as a code review proposed):`);
    console.log(`     hot replica accepts doubling chi2 with p = ${naiveHot.toExponential(2)}`);
    check(naiveHot < 1e-100,
        'the costScale=1 variant should be demonstrably frozen, not merely worse');
    console.log('     -> the hot chain cannot move at all; PT degenerates to greedy descent.');
}

// ---------------------------------------------------------------------------
console.log('\n5. rescaling is rare, and only at refresh boundaries\n');
// ---------------------------------------------------------------------------
{
    const s = makeScaler(1e6);
    let n = 0;
    // A two-decade descent sampled at 40 refresh boundaries.
    for (let k = 0; k < 40; k++) {
        if (s.maybeRescale(reps(1e6 * Math.pow(0.1, k / 20)))) n++;
    }
    console.log(`   ${n} rescalings over 40 refresh boundaries spanning a 100x descent`);
    check(n >= 3 && n <= 10, `expected a handful of rescalings, got ${n}`);

    // And the call site really is the refresh boundary, not the inner loop.
    const worker = fs.readFileSync(
        path.join(__dirname, '..', 'refinement_worker.js'), 'utf8').replace(/\r\n/g, '\n');
    check(worker.includes('maybeRescaleCost(replicas);'),
        'refinement_worker.js should call maybeRescaleCost');
    const callIdx = worker.indexOf('maybeRescaleCost(replicas);');
    const refreshIdx = worker.indexOf('iter % LEBAIL_PT_REFRESH_INTERVAL === 0');
    check(refreshIdx >= 0 && callIdx > refreshIdx && callIdx - refreshIdx < 1200,
        'maybeRescaleCost must sit inside the Le Bail refresh block');
    // It must NOT appear anywhere near the per-move Metropolis test.
    const moveIdx = worker.indexOf('-delta_cost / (costScale * replica.temp)');
    check(moveIdx < 0 || Math.abs(moveIdx - callIdx) > 500,
        'maybeRescaleCost must not be called from the per-move acceptance path');
    console.log('   worker calls it from the refresh block only');
}

console.log(`\n${checks} checks, ${failures} failures`);
process.exit(failures ? 1 : 0);
