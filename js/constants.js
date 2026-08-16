// constants.js
// version 157, 16 aug 2026
// ---------------------------------------------------------------------------
//  SINGLE SOURCE OF TRUTH for the physical and numerical constants shared by
//  the main thread (powder5.html), refinement_worker.js and
//  charge_flipping_worker.js.
//
//  Why this file exists
//  --------------------
//  PV_GAUSS_AREA, PV_LORENTZ_AREA, GSAS_GAUSSIAN_TO_DEG and STEPHENS_PREFACTOR
//  used to be declared independently in powder5.html and in
//  refinement_worker.js. Two copies of a physical constant is not a style
//  problem: if one is edited and the other is not, the preview and the fit are
//  computing DIFFERENT profile functions, and nothing in the program can
//  detect it. The refinement converges, the report is self-consistent, and the
//  numbers are wrong.
//
//  Loading
//  -------
//    browser :  script src="constants.js"
//                 (written without angle brackets on purpose: a literal
//                  closing script tag inside this file would terminate the
//                  page early if anyone ever inlined it into powder5.html)   (before everything else)
//    worker  :  importScripts('constants.js')
//    node    :  require('./constants.js')              (for the unit tests)
//
//  In the first two cases every declaration below lands in the realm's shared
//  script scope, so `PV_G0` etc. are simply visible. The module.exports footer
//  is only for node.
//
//  EVERY numeric literal in this file carries its justification. A constant
//  without a reason is a constant nobody can safely change.
// ---------------------------------------------------------------------------

'use strict';

// ===========================================================================
//  1. Profile truncation window
// ===========================================================================
//  A Gaussian is dead by ~4 FWHM; a Lorentzian is not. At 8 FWHM a pure
//  Lorentzian still loses ~4% of its area to truncation and another ~4% to the
//  pedestal subtraction in prepareVoigt. That in itself is survivable, because
//  the Le Bail extraction measures the area of the profile as actually drawn --
//  but the loss depends on eta, and eta moves as the widths refine, so a FIXED
//  window leaks truncation error straight into the intensity/width
//  correlation. The window is therefore scaled with the Lorentzian content,
//  and the exact analytic area of the UNtruncated shape is carried on the prep
//  object (prep.analyticArea) so integrated intensities can be reported free
//  of the truncation bias.
// ---------------------------------------------------------------------------

/**
 * Window half-width for a PURE GAUSSIAN, in multiples of the FWHM.
 * WHY 8: a Gaussian at 4 FWHM from centre is at exp(-ln2*64) ~ 1e-19 of peak,
 * i.e. below double-precision relevance. 8 is 4 with a factor-of-two margin so
 * that the eta-interpolated window never dips below the Gaussian requirement.
 * @type {number}
 */
const CALCULATION_WINDOW_MULTIPLIER_G = 8.0;

/**
 * Window half-width for a PURE LORENTZIAN, in multiples of the FWHM.
 * WHY 25: the truncated fraction of a Lorentzian beyond n FWHM is
 * (2/pi)*atan(1/(2n)) ~ 1/(pi*n). At n = 8 that is 4.0%; at n = 25 it is 1.3%,
 * and the residual is removed exactly by prep.analyticArea. Going further buys
 * little and costs linearly in evaluation time -- the profile loop is the
 * hottest path in the program.
 * @type {number}
 */
const CALCULATION_WINDOW_MULTIPLIER_L = 25.0;

/**
 * Hard cap on the profile half-window, in DEGREES 2-theta.
 * WHY 10: two purposes. (a) Cost. The window is a multiple of the FWHM, so a
 * refinement that transiently drives a width to a nonsensical value would
 * otherwise make one reflection touch the entire pattern, turning an
 * O(n_peaks * points_per_peak) evaluation into O(n_peaks * n_points) -- a
 * three-orders-of-magnitude stall exactly when the fit is already in trouble.
 * (b) Physics. A genuine Bragg reflection with a half-width beyond 10 deg 2th
 * is not a peak any more; capping there means the cost function stops
 * rewarding "explain the background with one enormous reflection".
 * @type {number}
 */
const PROFILE_WINDOW_MAX_DEG = 10.0;

/** Legacy alias, retained for older call sites. @type {number} */
const CALCULATION_WINDOW_MULTIPLIER = CALCULATION_WINDOW_MULTIPLIER_G;

// ===========================================================================
//  2. Minimum resolvable profile FWHM
// ===========================================================================
//  Set from the real data step (setMinProfileFwhmFromAxis) as soon as the
//  pattern is known.
//
//  A peak narrower than a couple of data steps is not measurable, and letting
//  the refiner go there is exactly how a fit with a bad cell "improves" Rwp:
//  a spike in the wrong place costs less than a broad peak in the wrong place,
//  so the profile collapses instead of the cell correcting.
// ---------------------------------------------------------------------------

/**
 * Minimum resolvable profile FWHM, degrees 2-theta. Mutable: overwritten by
 * setMinProfileFwhmFromAxis once the real axis is known. The 1e-3 seed is a
 * placeholder for the case where a profile is evaluated before any data is
 * loaded (tooltips, defaults), not a physically meaningful value.
 * @type {number}
 */
let MIN_PROFILE_FWHM_DEG = 1e-3;

/**
 * Sets MIN_PROFILE_FWHM_DEG from the measured step of a 2-theta axis.
 *
 * WHY 2 x step: Nyquist. A peak whose FWHM equals one step is sampled once and
 * carries no shape information at all; two steps is the coarsest FWHM at which
 * the profile has any width the data can constrain.
 * WHY the 1e-4 clamp: guards against a degenerate or single-point axis
 * producing a zero or absurdly small floor, which would re-open the collapse
 * the floor exists to prevent.
 *
 * @param {ArrayLike<number>} tthAxis Ascending 2-theta axis, degrees.
 * @returns {void}
 */
function setMinProfileFwhmFromAxis(tthAxis) {
    if (!tthAxis || tthAxis.length < 3) return;
    const step = (tthAxis[tthAxis.length - 1] - tthAxis[0]) / (tthAxis.length - 1);
    if (isFinite(step) && step > 0) MIN_PROFILE_FWHM_DEG = Math.max(1e-4, 2 * step);
}

/** @returns {number} The current minimum resolvable FWHM, degrees 2-theta. */
function getMinProfileFwhmDeg() { return MIN_PROFILE_FWHM_DEG; }

// ===========================================================================
//  3. Smooth positive branch
// ===========================================================================

/**
 * Smooth positive branch: >= floor, strictly increasing, derivative never 0.
 *     softPositive(x, f) = (x + sqrt(x^2 + 4 f^2)) / 2
 * f at x = 0, ~x for x >> f, ~f^2/|x| for x << -f, d/dx > 0 everywhere.
 *
 * FIX: replaces the `if (sigma_sq_cdeg2 > 0) { ... }` guards in the width
 * calculators. Those guards were a hard dead zone -- once a width combination
 * went non-positive the width stuck at its default value, the one-sided
 * finite-difference column came out identically zero, and the parameter was
 * frozen there for the rest of the fit. It could not recover even after the
 * cell and zero-point had been corrected.
 *
 * @param {number} x     Value to floor.
 * @param {number} floor Strictly positive floor.
 * @returns {number}
 */
function softPositive(x, floor) {
    if (!isFinite(x)) return floor;
    return 0.5 * (x + Math.sqrt(x * x + 4 * floor * floor));
}

// ===========================================================================
//  4. Pseudo-Voigt normalisation
// ===========================================================================

/**
 * Peak height of the UNIT-AREA Gaussian, in units of 1/FWHM. = 0.9394372787
 * Exact: 2*sqrt(ln2/pi). Not a fitted number.
 * @type {number}
 */
const PV_G0 = 2 * Math.sqrt(Math.LN2 / Math.PI);

/**
 * Peak height of the UNIT-AREA Lorentzian, in units of 1/FWHM. = 0.6366197724
 * Exact: 2/pi. Not a fitted number.
 * @type {number}
 */
const PV_L0 = 2 / Math.PI;

/** Area of the UNIT-HEIGHT Gaussian, in units of FWHM. = 1.0644670194 @type {number} */
const PV_GAUSS_AREA = 1 / PV_G0;

/** Area of the UNIT-HEIGHT Lorentzian, in units of FWHM. = 1.5707963268 @type {number} */
const PV_LORENTZ_AREA = 1 / PV_L0;

/**
 * Mixing coefficients for a pseudo-Voigt evaluated on UNIT-HEIGHT components.
 *
 * FIX (eta convention). The pseudo-Voigt is defined as eta*L + (1-eta)*G with
 * L and G normalised to unit AREA -- that is what makes eta a mixing fraction,
 * and it is the convention the Thompson-Cox-Hastings eta polynomial was fitted
 * for. This code carries a peak HEIGHT per reflection, so the components are
 * evaluated unit-height and the area normalisation is folded into these
 * coefficients, which are then renormalised so the mixture is still 1 at the
 * peak centre.
 *
 * Mixing unit-height components directly (the previous behaviour) is a
 * DIFFERENT function: L(0)/G(0) = 0.637/0.939 = 0.678 at equal FWHM, so a TCH
 * eta of 0.578 (Gamma_L/Gamma = 0.5) was rendered as a shape whose true area
 * mixing fraction is 0.48. For simple_pvoigt, where eta is refined, that was
 * only a mislabelling; for tch_aniso, where eta is COMPUTED from the
 * polynomial, the modelled Lorentzian content was systematically too high and
 * X / Y had to absorb the difference.
 *
 * @param {number} eta Lorentzian AREA fraction, clamped to [0, 1].
 * @param {number} hg  Gaussian FWHM, degrees 2-theta.
 * @param {number} hl  Lorentzian FWHM, degrees 2-theta.
 * @returns {{cL: number, cG: number}} Weights for unit-HEIGHT components.
 */
function pvMixCoeffs(eta, hg, hl) {
    const e = Math.max(0, Math.min(1, eta));
    const wL = e * PV_L0 / Math.max(1e-12, hl);
    const wG = (1 - e) * PV_G0 / Math.max(1e-12, hg);
    const s = wL + wG;
    if (!(s > 0) || !isFinite(s)) return { cL: e, cG: 1 - e };
    return { cL: wL / s, cG: wG / s };
}

// ===========================================================================
//  5. GSAS unit conversions and Stephens anisotropic strain
// ===========================================================================
/*  ============================================================
 *  GSAS profile-function conventions (CW types 3 & 4)
 *  ============================================================
 *  Units of the profile parameters follow GSAS exactly:
 *
 *    Gaussian (Caglioti):   sigma^2[cdeg^2]/(8*ln2)
 *                          = GU*tan^2(th) + GV*tan(th) + GW + GP/cos^2(th)
 *      -> H_G[deg] = sqrt(8*ln2 * sigma^2) / 100
 *
 *    Lorentzian:            gamma_L[cdeg]
 *                          = LX/cos(th) + LY*tan(th)          (TCH)
 *                          = LX/cos(th)                       (Simple/Split)
 *      -> gamma_L[deg] = gamma_L[cdeg] / 100
 *
 *    Conversions (Kaduk & Reid, Powder Diffr. 26, 88, 2011):
 *      GSAS_GU = 1803.4 * FullProf_U          (1803.4 = 100^2/(8*ln2))
 *      GSAS_LX = 100   * FullProf_Y_size
 *      GSAS_LY = 100   * FullProf_X_strain
 *
 *  In this UI:
 *    - Simple pVoigt  GU,GV,GW,GP,LX           = GSAS GU,GV,GW,GP,LX
 *    - TCH            U,V,W                    = GSAS GU,GV,GW
 *                     X = GSAS LY (strain),  Y = GSAS LX (size)
 *      (matches Thompson-Cox-Hastings 1987 notation, which GSAS adopted.)
 *    - Split pVoigt   GU_L/R, GV_L/R, GW_L/R, LX_L/R = GSAS units
 *
 *  Stephens (1999) anisotropic strain (TOPAS / GSAS-II form):
 *      M_HKL = SUM S_HKL * h^H k^K l^L            (degree-4 polynomial)
 *      pp    = d^2 * sqrt(max(M_HKL, 0)) / 1000
 *      gamma_L_aniso[deg] = (1.8/pi) * pp * eta     * tan(th)
 *      gamma_G_aniso[deg] = (1.8/pi) * pp * (1-eta) * tan(th)
 *
 *  All width helpers in profile.js return widths in DEGREES 2-theta.
 */

/**
 * sqrt(8*ln2)/100: converts sqrt(sigma^2[cdeg^2]/(8*ln2)) -> FWHM in degrees.
 * The sqrt(8 ln 2) is the exact sigma -> FWHM factor for a Gaussian; the /100
 * is the centidegree -> degree conversion GSAS uses for its width parameters.
 * @type {number}
 */
const GSAS_GAUSSIAN_TO_DEG = Math.sqrt(8 * Math.LN2) / 100.0;

/** 1/100: converts gamma_L[cdeg] -> gamma_L[deg]. @type {number} */
const GSAS_LORENTZIAN_TO_DEG = 1.0 / 100.0;

/**
 * Stephens prefactor, 1.8/pi = 0.5730.
 * WHY this value: it is the canonical constant used by GSAS-II and TOPAS in
 * their implementation of Stephens (1999), J. Appl. Cryst. 32, 281. It folds
 * together the sigma -> FWHM conversion and the 1/1000 scaling built into the
 * definition of pp. It is a CONVENTION, not a derived quantity -- changing it
 * makes the refined S_HKL values incomparable with values from those programs,
 * which is the whole point of matching them.
 * @type {number}
 */
const STEPHENS_PREFACTOR = 1.8 / Math.PI;

/**
 * Default Gaussian/Lorentzian split of the Stephens broadening when the user
 * supplies no explicit eta_aniso. 1.0 = purely LORENTZIAN.
 * WHY pure Lorentzian: microstrain broadening from a distribution of lattice
 * parameters produces a peak shape whose tails fall off far more slowly than a
 * Gaussian, and GSAS-II defaults to the Lorentzian-dominated form for the same
 * reason. If your sample's strain broadening is closer to Gaussian, set
 * params.eta_aniso explicitly rather than editing this default.
 * @type {number}
 */
const STEPHENS_DEFAULT_ETA = 1.0;

// ===========================================================================
//  6. Le Bail decomposition
// ===========================================================================
//  Le Bail revival floor.
//
//  BUG FIX (absorbing zero). calculatePattern skips reflections with
//  intensity <= 0, and the decomposition weight is
//      f_j(i) = I_j * phi_j(i) / SUM_k I_k * phi_k(i),
//  so a reflection sitting at exactly zero gets f = 0 at every point, is
//  handed no intensity, and stays at zero forever. Since the extracted area is
//  clipped at zero and the net observed profile is (correctly) allowed to go
//  negative, a reflection whose window happens to sit over background -- which
//  is precisely what happens while the cell or the zero-point is still wrong --
//  had a ~50% chance of dying on the FIRST pass, including during the seeding
//  cycles, and could never come back once the cell was corrected.
//
//  The heights used for modelling are therefore floored at a tiny fraction of
//  the mean, which keeps f_j > 0 so the partition can revive, while
//  contributing nothing measurable to the pattern. The MEASURED quantities --
//  I_integrated and I_sigma -- are reported unclipped, so a genuinely absent
//  reflection is still detectable as I/sigma ~ 0 (or negative) rather than
//  being hidden behind a floor.
// ---------------------------------------------------------------------------

/**
 * Revival floor as a fraction of the mean peak height.
 * WHY 1e-5: it has to clear two bars from opposite directions.
 *   Lower bound  -- it must survive f_j = I_j*phi_j / SUM in double precision
 *                   against neighbours of order 1, so anything below ~1e-12 is
 *                   indistinguishable from the zero it exists to avoid.
 *   Upper bound  -- it must be invisible in the modelled pattern. Counting
 *                   statistics on a 10^4-count peak give ~1% noise; 1e-5 of
 *                   the mean is three orders of magnitude under that, so a
 *                   truly absent reflection contributes nothing measurable.
 * Anything in 1e-8 .. 1e-4 satisfies both; 1e-5 sits in the middle on a log
 * scale. The value is NOT critical -- if a fit is sensitive to it, the problem
 * is elsewhere.
 * @type {number}
 */
const LEBAIL_REVIVAL_FRACTION = 1e-5;

/**
 * Decomposition passes per re-extraction.
 * WHY 6: the Le Bail partition is a fixed-point iteration. Starting from flat
 * intensities it settles to a relative change below ~1e-4 in 4-6 passes on
 * ordinary patterns; 6 buys the margin for heavily overlapped ones. The cost
 * is small next to building n_params Jacobian columns, which is what the
 * surrounding LM iteration does anyway.
 * @type {number}
 */
const LEBAIL_PASSES_PER_ITER = 6;

/**
 * Hard cap on decomposition passes when the extraction is run to convergence.
 *
 * WHY A CONVERGENCE TEST AT ALL: LEBAIL_PASSES_PER_ITER is a fixed count, and
 * a fixed count is either wasteful or wrong depending on the pattern. On a
 * well-separated pattern the partition settles in two passes and the other
 * four are pure cost; on a heavily overlapped one -- the case that actually
 * needs the method -- six is not enough, and the LM then fits its profile
 * parameters against intensities that are still moving. Iterating to a
 * relative change below LEBAIL_PASS_TOL fixes both ends.
 *
 * WHY 24: six times the old fixed count. Nothing that has not settled by then
 * is going to; it is a runaway guard, not a target.
 * @type {number}
 */
const LEBAIL_MAX_PASSES = 24;

/**
 * Relative change in the extracted intensities below which the decomposition
 * is considered settled.
 *
 * WHY 1e-4: the fixed six-pass schedule was chosen because it reached "a
 * relative change below ~1e-4 on ordinary patterns" (see above), so this is
 * the same target, now measured instead of assumed.
 * @type {number}
 */
const LEBAIL_PASS_TOL = 1e-4;

/**
 * Restart each re-extraction from flat intensities?
 *
 * MUST stay true. A single pass does not settle, so the intensities carry
 * HISTORY: the cost at a given set of parameters depends on how the fit
 * arrived there. That makes the Le Bail objective path-dependent, and a
 * path-dependent objective cannot be minimised meaningfully -- the
 * decomposition keeps adapting to cover mispositioned peaks, so the cost can
 * fall while the cell drifts away. Measured on a synthetic cubic pattern, the
 * SETTLED objective is emphatic (Rwp 9.4% at the true cell against 91.3% at
 * the aliased cell+zero-point minimum), but with one pass per iteration the
 * refiner slid into the aliased minimum anyway and reported 81%.
 *
 * Restarting from flat makes the extraction a deterministic function of the
 * parameters, which is what the outer monotonicity guard needs in order to
 * mean anything.
 * @type {boolean}
 */
const LEBAIL_RESET_EACH_ITER = true;

/**
 * Flat seed height used when LEBAIL_RESET_EACH_ITER restarts the partition.
 * WHY 100 and why the exact value does not matter: the first decomposition
 * pass computes f_j = I_j*phi_j / SUM_k I_k*phi_k, which is invariant under a
 * common scaling of every I_j. Starting all reflections EQUAL is the whole
 * point (it is what makes the extraction path-independent); the magnitude
 * cancels. 100 is simply a readable number in the middle of the range the
 * extracted intensities occupy.
 * @type {number}
 */
const LEBAIL_FLAT_SEED = 100.0;

/**
 * How often PT re-extracts the Le Bail intensities, in iterations.
 *
 * BUG FIX. Holding the intensities fixed for a whole LM iteration is correct
 * (see the long note in refineParametersLM), but that argument does NOT carry
 * over to a Monte Carlo run of thousands of moves. Intensities frozen at the
 * INITIAL parameters are the decomposition of the observed profile against the
 * starting cell, so the objective's minimum sits exactly where the search
 * started -- which makes a global search over the cell unable to find anything
 * else, the one job PT exists to do. Re-extract periodically, per replica, at
 * that replica's own parameters.
 *
 * WHY 40: a full re-extraction is LEBAIL_PASSES_PER_ITER pattern evaluations
 * per replica, against one evaluation per proposed move. At 40 the refresh
 * costs ~15% of the run, and the objective is stale for at most 40 moves --
 * short compared with the ~10^2-10^3 moves a replica needs to travel anywhere
 * interesting in cell space.
 * @type {number}
 */
const LEBAIL_PT_REFRESH_INTERVAL = 40;

// ===========================================================================
//  7. Levenberg-Marquardt tuning
// ===========================================================================
//  Finite differences. The escalation ladder is BOUNDED: an older version
//  multiplied the probe by 100 up to five times (1e8 x the base step). A width
//  parameter sitting in a numerical dead zone could be driven to a nonsensical
//  value, the pattern then changed wildly, `responded` flipped true, and the
//  resulting secant -- across eight orders of magnitude -- was accepted as a
//  derivative.
// ---------------------------------------------------------------------------

/**
 * First finite-difference probe, as a fraction of the parameter's scale.
 * WHY 1e-4: for a one-sided difference the total error is
 * O(h * f'') + O(eps_eff / h), minimised at h ~ sqrt(eps_eff). The cost here
 * is a SUM over ~10^3-10^4 points, so its effective relative noise is nearer
 * 1e-12 than machine epsilon, putting the optimum around 1e-6..1e-5. 1e-4 sits
 * deliberately above that: the extra truncation error is tolerable, and the
 * larger probe is far more likely to escape the flat spots the profile
 * functions genuinely have (a width sitting on its soft floor, a Stephens term
 * at zero).
 * @type {number}
 */
const FD_BASE_FRACTION = 1e-4;

/**
 * Hard cap on the escalated probe, as a fraction of the parameter's own
 * magnitude.
 * WHY 2%: large enough to escape a numerical dead zone, small enough that the
 * secant is still a derivative. Above a few percent the profile response to a
 * cell edge or a width is visibly non-linear, and the "derivative" that comes
 * back is a chord across curvature -- which is what produced the eight-orders-
 * of-magnitude garbage the ladder cap exists to prevent.
 * @type {number}
 */
const FD_MAX_PROBE_FRACTION = 0.02;

/**
 * Maximum probe escalations before the column is declared empty.
 * WHY 4: the ladder multiplies by 100 each time, so from FD_BASE_FRACTION 1e-4
 * four steps already reach FD_MAX_PROBE_FRACTION. More attempts cannot go
 * anywhere the cap allows.
 * @type {number}
 */
const FD_MAX_ATTEMPTS = 4;

/**
 * Relative tolerance below which a Jacobian column counts as empty, i.e. the
 * parameter is not determined by the data at this point. Compared against the
 * LARGEST diagonal of JtJ, so it is scale-free.
 * WHY 1e-12: measured JtJ diagonals here span ten orders of magnitude purely
 * because the parameters have different units, so the tolerance has to sit
 * below that spread while staying above double-precision accumulation noise on
 * a 10^4-term sum (~1e-13 relative). 1e-12 is the one decade that satisfies
 * both.
 * @type {number}
 */
const JTJ_RANK_TOL = 1e-12;

/**
 * Relative cost drop below which an LM iteration counts as converged.
 * WHY 1e-7 and not smaller: in Le Bail mode the intensity re-extraction at the
 * top of the next iteration moves the cost by far more than 1e-9, so the old
 * 1e-9 threshold could never fire and LM always burned maxIter. 1e-7 is a
 * relative change no longer distinguishable from the re-extraction's own
 * jitter.
 * @type {number}
 */
const LM_COST_TOL = 1e-7;

/**
 * Consecutive iterations under LM_COST_TOL required to declare convergence.
 * WHY 2: one flat step is common mid-fit (a rejected trial, a lambda
 * adjustment); two in a row is not.
 * @type {number}
 */
const LM_CONVERGED_REPEATS = 2;

/**
 * Relative slack before an outer (extract-then-refine) iteration is judged to
 * have made the TRUE, re-extracted cost worse.
 * WHY 1e-6: an order of magnitude looser than LM_COST_TOL, so a step that
 * merely converged is never mistaken for a step that regressed.
 * @type {number}
 */
const LM_OUTER_TOL = 1e-6;

/**
 * Outer rejections tolerated before giving up.
 * WHY 6: each rejection halves the trust region, so six take the step cap down
 * by 2^-6 ~ 1.6%. Below that the step is smaller than the finite-difference
 * probe and further iterations cannot learn anything new.
 * @type {number}
 */
const LM_MAX_OUTER_REJECTS = 6;

/**
 * Initial trust-region scaling on the per-parameter step caps.
 * WHY 0.05: the caps in getParameterMapping are PHYSICAL limits (half an
 * angstrom on a cell edge, ~2 deg on a cell angle). Starting at 5% of those
 * means the first LM step is ~0.025 A on a cell edge -- big enough to make
 * real progress, small enough that a bad initial Jacobian cannot throw the
 * cell somewhere unrecoverable.
 * @type {number}
 */
const LM_TRUST_INITIAL = 0.05;

/**
 * Trust-region growth factor after an iteration that improved the
 * re-extracted cost.
 * WHY 1.6: geometric growth reaches the full physical cap in ~6 good
 * iterations. A larger factor (2 or more) overshoots into the rejection
 * branch too readily; a smaller one wastes iterations creeping.
 * @type {number}
 */
const LM_TRUST_GROWTH = 1.6;

// ===========================================================================
//  8. Parallel tempering
// ===========================================================================

/**
 * Number of PT replicas.
 * WHY 8: the swap acceptance between neighbours falls off as the temperature
 * ratio grows. Spanning PT_MAX_TEMP..PT_MIN_TEMP (five decades) with 8
 * replicas gives a ratio of 10^(5/7) ~ 5.2 per rung, which keeps neighbouring
 * chains overlapping enough to exchange while staying affordable -- the cost
 * is linear in replica count and every replica evaluates the full pattern.
 * @type {number}
 */
const PT_NUM_REPLICAS = 8;

/**
 * Hottest replica temperature, in units of the RELATIVE cost scale (see
 * ptCostScale). At T = 1 a move that doubles chi-squared is accepted with
 * probability e^-1 ~ 0.37, which is the definition of a chain that explores
 * freely.
 * @type {number}
 */
const PT_MAX_TEMP = 1.0;

/**
 * Coldest replica temperature, same units. At T = 1e-5 a move that worsens
 * chi-squared by 0.01% is already rejected with probability 1 - e^-10, i.e.
 * the chain is effectively a greedy minimiser -- which is what the cold end of
 * a PT ladder is for.
 * @type {number}
 */
const PT_MIN_TEMP = 1e-5;

/**
 * Iterations between replica-exchange attempts.
 * WHY 10: swaps must be rare enough that a replica decorrelates between them
 * (otherwise the ladder just shuffles identical states) and frequent enough to
 * carry a good configuration down the ladder within a run.
 * @type {number}
 */
const PT_SWAP_INTERVAL = 10;

/**
 * Floor on the PT cost scale, guarding a division by zero for a fit that
 * starts at a numerically perfect cost.
 * @type {number}
 */
const PT_COST_SCALE_FLOOR = 1e-9;

/**
 * Fractional cost change below which the PT energy scale is NOT rescaled.
 *
 * WHY this exists. Both PT acceptance tests divide the cost difference by a
 * reference energy scale, which makes the temperature ladder RELATIVE:
 * exp(-(dChi2/scale)/T). The scale used to be frozen at the initial cost. That
 * is fine for the first stretch of a run, but a successful global search drops
 * chi-squared by one to two orders of magnitude, and a frozen scale then means
 * the ladder is silently 10-100x hotter in relative terms than it was
 * designed to be: the cold replicas stop being cold, and PT loses the
 * exploit end of its explore/exploit balance exactly when it is closing in on
 * a solution.
 *
 * Rescaling on every refresh would break detailed balance outright (the
 * acceptance criterion would change under the chain's feet on every step), so
 * the scale is instead re-quoted only at the periodic Le Bail refresh
 * boundaries -- points at which the objective function is ALREADY being
 * redefined and every replica cost re-quoted -- and only when the cost has
 * moved by more than this fraction.
 *
 * WHY 0.5: a factor-of-two change in chi-squared is a genuine change of regime
 * for the search, well outside Monte Carlo jitter, and rescaling at that
 * granularity means a two-decade descent triggers ~7 rescalings across the
 * whole run rather than one per refresh.
 * @type {number}
 */
const PT_RESCALE_THRESHOLD = 0.5;

// ===========================================================================
//  9. Small shared numeric utilities
// ===========================================================================

/**
 * First index i with arr[i] >= target, on an ASCENDING array. O(log n).
 * PERF FIX: replaces `let s=0; while (s<n && arr[s]<target) s++;`, which made
 * every pattern evaluation O(n_peaks * n_points) instead of
 * O(n_peaks * points_per_peak).
 *
 * @param {ArrayLike<number>} arr Ascending array.
 * @param {number} target
 * @returns {number} Index in [0, arr.length].
 */
function lowerBound(arr, target) {
    let lo = 0, hi = arr.length;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (arr[mid] < target) lo = mid + 1; else hi = mid;
    }
    return lo;
}

// ---------------------------------------------------------------------------
//  Node export footer. MIN_PROFILE_FWHM_DEG is exposed as a live getter so a
//  node consumer sees mutations from setMinProfileFwhmFromAxis; a plain value
//  would snapshot the 1e-3 placeholder.
// ---------------------------------------------------------------------------
if (typeof module !== 'undefined' && module.exports) {
    const _exports = {
        CALCULATION_WINDOW_MULTIPLIER_G, CALCULATION_WINDOW_MULTIPLIER_L,
        PROFILE_WINDOW_MAX_DEG, CALCULATION_WINDOW_MULTIPLIER,
        setMinProfileFwhmFromAxis, getMinProfileFwhmDeg, softPositive,
        PV_G0, PV_L0, PV_GAUSS_AREA, PV_LORENTZ_AREA, pvMixCoeffs,
        GSAS_GAUSSIAN_TO_DEG, GSAS_LORENTZIAN_TO_DEG,
        STEPHENS_PREFACTOR, STEPHENS_DEFAULT_ETA,
        LEBAIL_REVIVAL_FRACTION, LEBAIL_PASSES_PER_ITER,
        LEBAIL_MAX_PASSES, LEBAIL_PASS_TOL,
        LEBAIL_RESET_EACH_ITER, LEBAIL_FLAT_SEED, LEBAIL_PT_REFRESH_INTERVAL,
        FD_BASE_FRACTION, FD_MAX_PROBE_FRACTION, FD_MAX_ATTEMPTS,
        JTJ_RANK_TOL, LM_COST_TOL, LM_CONVERGED_REPEATS,
        LM_OUTER_TOL, LM_MAX_OUTER_REJECTS, LM_TRUST_INITIAL, LM_TRUST_GROWTH,
        PT_NUM_REPLICAS, PT_MAX_TEMP, PT_MIN_TEMP, PT_SWAP_INTERVAL,
        PT_COST_SCALE_FLOOR, PT_RESCALE_THRESHOLD,
        lowerBound
    };
    Object.defineProperty(_exports, 'MIN_PROFILE_FWHM_DEG', {
        enumerable: true,
        get: () => MIN_PROFILE_FWHM_DEG
    });
    module.exports = _exports;
}

const CF_DEFAULTS = {
    gridSize: 32,
    thresholdSigma: 1.10,
    maxIterations: 1000,
    randomStarts: 3,
    peakMerge: 0.9,
    symLambda: 0.5,
    overlapTolTth: 0.05
};
