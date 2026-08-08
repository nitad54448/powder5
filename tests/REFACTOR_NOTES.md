# Refactor notes — shared modules, PT energy scale, GPU buffer lifetime

Work done against the code review. Everything here is verified by
`node test/run_all.js` (five suites, ~420k assertions, all passing).

---

## New files

| File | Contents |
|---|---|
| `constants.js` | Every physical and numerical constant, each with its justification |
| `crystal.js` | Metric tensors, d-spacings, cell geometry, monoclinic convention |
| `profile.js` | Peak shift, profile widths, `prepareVoigt` / `evalVoigt`, areas |
| `test/*.js` | Five suites plus `run_all.js` |

**Deployment:** the three shared modules must sit next to `powder5.html`,
alongside `sg_engine.js` and the space-group JSON. `importScripts` resolves
relative to the worker's URL, so the workers need them in the same directory.
No build step; they are plain classic scripts.

---

## What was actually done

### 1. Duplication removed (review items 16, 17, 18, 19, 20)

The review found four duplicated things. There were **five**, and two had
already silently diverged.

| Symbol | Was defined in | Now |
|---|---|---|
| `PV_*`, `GSAS_*`, `STEPHENS_*`, window multipliers | powder5.html **and** refinement_worker.js | `constants.js` |
| `prepareVoigt`, `evalVoigt`, width calculators (~200 lines) | powder5.html **and** refinement_worker.js | `profile.js` |
| `updateHklPositions` (byte-identical copies) | powder5.html **and** refinement_worker.js | `crystal.js` |
| `metricTensor`, `fracDistance`, `cellVolume`, `cellBasis`, `fracToCart` | charge_flipping_worker.js **and** density3d.js | `crystal.js` |
| `lowerBound` | powder5.html **and** refinement_worker.js | `constants.js` |
| `buildInvDsq` (the reciprocal metric, a **fifth** copy) | powder5.html | `crystal.js` as `buildInvDsqEvaluator` |

**The two that had already diverged**, which is the concrete evidence the
refactor was worth doing:

- the worker's `evalVoigt` carried a reciprocal-precompute optimisation
  (`delta * inv_H_G` rather than `delta / H_G`) that the main thread never
  received;
- the main thread carried `getPseudoVoigtArea`, which the worker never had.

Neither changed results *this time*. Nothing in the program could have told you
if it had — the preview would have drawn one peak shape, the fit would have
minimised another, and the report would have been perfectly self-consistent
about the wrong answer.

The fifth copy is worth singling out. `buildInvDsq` in `powder5.html` carried a
comment promising it used *"exactly the same metric as `updateHklPositions()`
so that generation and the downstream 2-theta filter can never disagree."*
That promise was being kept by hand, across two files. It is now structural:
`buildInvDsqEvaluator` is built from the same `buildInvDsq` that
`updateHklPositions` uses, and the two are asserted to agree to 4.4e-16 across
every crystal system and all three monoclinic settings.

**Magic numbers (#20).** Every literal now carries the reason it has the value
it has, not just a label. `PROFILE_WINDOW_MAX_DEG = 10` gets both its cost
argument and its physical one. `FD_BASE_FRACTION = 1e-4` gets the
sqrt(eps_eff) analysis and an explanation of why it deliberately sits *above*
the theoretical optimum. `LEBAIL_REVIVAL_FRACTION = 1e-5` gets a two-sided
bound showing the value is not critical — which is more useful than a
justification for the exact number, because it tells you what you may safely
change it to.

**Types (#19).** JSDoc typedefs on the shared numerical surface: `ProfileParams`,
`Reflection`, `Widths`, `VoigtPrep`, `Cell`, `InvDsqForm`, `InvDsqEvaluator`,
plus `@param`/`@returns` on every exported function. Enough for editor
tooling and `checkJs` to catch the number-vs-Float64Array confusion the review
was worried about, without a build step.

### 2. GPU buffer lifetime (review item 5) — real bug, wrong diagnosis

The review said the leak happened when `cfAcquireGPU()` failed. It cannot:
`cfAcquireGPU` runs before any buffer exists.

The real gap: buffers were created at lines 1237–1311 but `try` did not open
until 1399. A throw in that window leaked everything already allocated, with no
reference left to destroy it. Now every allocation goes through a `track()`
helper that registers it at creation, and `try` opens before the first one.
Verified: no `device.createBuffer` call outside `track` remains.

### 3. PT energy scale — the real issue under review item 9

The review's diagnosis was wrong and its fix (`costScale = 1.0`) would have
frozen every chain in the ladder solid. But there was a genuine bug underneath:
`costScale` was frozen at the *initial* cost for the whole run.

A successful global search drops χ² by one to two orders of magnitude, and a
frozen scale then leaves the ladder 10–100× hotter in relative terms than
designed. Measured — after a 100× descent, the **coldest** replica accepts a
+0.1% χ² move with:

```
frozen scale    p = 3.68e-01     <- has stopped being cold
rescaled        p = 3.72e-44
```

The scale cannot be re-quoted continuously — a denominator that moves under the
chain's feet on every step breaks detailed balance. So it is re-quoted only at
the periodic Le Bail refresh boundaries, where the objective is *already* being
redefined and every replica cost re-quoted anyway, and only when χ² has moved
by more than a factor of two. Between two refreshes the acceptance criterion is
exactly constant, which is what the Metropolis argument needs. Fires ~5 times
over a two-decade descent, not once per refresh.

### 4. `Math.hypot` (review item 14) — done, and it is bigger than claimed

The review estimated 3–5×. Measured 16× over 4e6 iterations, with identical
checksums. Four hot sites in `charge_flipping_worker.js`; the one in the origin
search runs over the whole N³ grid (2.1e6 calls at 128³).

### 5. FFT twiddle tables (review item 7) — **implemented, measured, reverted**

This is the one that did not go as planned, and I think the negative result is
the useful part.

I built the precomputed-twiddle version the review asked for. Then I measured
it against a reference DFT and against the original:

| | rel. error vs DFT | 64³ forward transform |
|---|---|---|
| recurrence (original) | 2.9e-14 @ n=128, 2.7e-13 @ n=1024 | 23.4 ms |
| table | 2.9e-14 @ n=128, 2.7e-13 @ n=1024 | 26.8 ms |
| table, per-direction (no sign multiply) | — | 27.3 ms |

**Accuracy is identical to two significant figures at every size.** The error
is dominated by the O(log n) accumulation through the butterflies, not by the
twiddle recurrence, so removing the recurrence removes nothing measurable. At
3e-14 relative it is twelve orders of magnitude below the counting statistics
on the intensities being transformed — "visible phase noise at 128³" does not
survive contact with a measurement.

**And the table is 12–15% slower.** The recurrence keeps its two doubles in
registers across the inner loop; the table costs two loads per butterfly. This
is the opposite of the usual C intuition, which is exactly why it needed
measuring rather than assuming.

So I reverted it. What survives is a comment above `fft1d` carrying the table
of measurements and the instruction to run `node test/test_fft.js` before
changing anything, and a test suite that keeps the alternative implementation
alive purely so those numbers stay reproducible. The next person to read this
code and think "recursive twiddles, that's a bug" will find the answer instead
of redoing the work.

---

## Not done, and why

Items 1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 15 were false positives — see the
earlier analysis. Three of them (#3 sigma scaling, #8 Surface Nets, #9
`costScale = 1.0`) would have introduced real bugs. Two guards now sit in the
code against re-breaking them:

- `crystal.js` warns that the general monoclinic branch is correct and has been
  "fixed" into being wrong before;
- `constants.js` documents `LEBAIL_RESET_EACH_ITER` as **must stay true**, with
  the measured Rwp numbers.

---

## Verification

```
$ node test/run_all.js

test_refactor_equivalence.js    420,414 checks, 0 failures
test_no_duplicate_symbols.js    0 failures
test_fft.js                     0 failures
test_pt_cost_scale.js           24 checks, 0 failures
test_load_integration.js        40 checks, 0 failures

All 5 suites passed.
```

**`test_refactor_equivalence.js`** is the load-bearing one. It extracts the
pre-refactor code into a reference module and compares bit-for-bit (`Object.is`,
so ±0 and NaN are distinguished) across all three profile types, all eight
crystal systems, all three monoclinic settings, 61 sample points per profile
window including both exact edges, and degenerate inputs (zero cell edge,
179.99° angles, zero wavelength). A refactor of numerical code is only safe if
it is provably a no-op; this is that proof.

**`test_no_duplicate_symbols.js`** derives its symbol list from the modules'
own exports, so it cannot go stale. It found three duplicates I had missed,
including the fifth metric copy.

**`test_load_integration.js`** matters more than it looks. The other suites use
`require()`, which is *not* how these files load in production. This one
concatenates them the way a browser script scope and `importScripts` do, in one
`vm` realm, so a duplicate `const` is a SyntaxError here exactly as it would be
on first paint. It also verifies the cross-module link is live — change
`MIN_PROFILE_FWHM_DEG` in constants.js and the width floor in profile.js
follows, which a `require()`-based test cannot check because there the two
modules would hold separate bindings.

---

## Suggested next step

`test_refactor_equivalence.js` compares against a snapshot of the old code, so
it proves the refactor changed nothing — but it says nothing about whether the
physics was right to begin with. The natural follow-up is the integration test
the review suggested: a synthetic pattern from a known structure, refined from
a perturbed starting cell, asserting the refinement recovers the parameters.
That would be the first test in this codebase that could catch a wrong answer
rather than a changed one.
