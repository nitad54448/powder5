const LP_DEG2RAD = Math.PI / 180;
const POLARIZATION_DEFAULTS = { mode: 'lab', monoTth: 0, polFrac: 0.95 };

function normalizePolarization(pol) {
    const p = (pol && typeof pol === 'object') ? pol : {};
    let mode = String(p.mode == null ? POLARIZATION_DEFAULTS.mode : p.mode).toLowerCase();
    if (mode !== 'none' && mode !== 'synchrotron') mode = 'lab';

    let monoTth = Number(p.monoTth);
    if (!isFinite(monoTth) || monoTth < 0 || monoTth >= 180) monoTth = 0;

    let polFrac = Number(p.polFrac);
    if (!isFinite(polFrac)) polFrac = POLARIZATION_DEFAULTS.polFrac;
    polFrac = Math.min(1, Math.max(0, polFrac));

    const K = sharkoPolarizationK({ mode, monoTth, polFrac });
    return { mode, monoTth, polFrac, K };
}

// ---------------------------------------------------------------------------
//  EVERY Lp IN THIS PROGRAM COMES FROM crystal.js.
//
//  These are named wrappers over sharko* and nothing else -- no arithmetic
//  lives in this file. The formula existed in one place already; what did not
//  was the DEGENERATE-ANGLE POLICY. lorentzPolarization() here returned NaN at
//  a 2-theta of zero or 180, while a function of the same name in each of the
//  two workers delegated to sharkoLp() and returned 1. One name, two answers,
//  decided by which file the caller happened to be in.
//
//  Both policies were correct for their own callers, so neither was removed:
//  they are now sharkoLp (lenient, for reducing a pattern) and sharkoLpStrict
//  (NaN, for a report or an export that must not invent a correction), and the
//  page takes the strict one, which is what it has always used.
// ---------------------------------------------------------------------------

/** Lorentz factor alone. NaN where the geometry is degenerate. */
function lorentzFactor(tthDeg) {
    return sharkoLorentzFactor(tthDeg);
}

/** Polarisation factor alone, from a descriptor. */
function polarizationFactor(tthDeg, pol) {
    return sharkoPolarizationFactor(tthDeg, normalizePolarization(pol).K);
}

/**
 * Lp for the page: reporting, exports and structure-factor conversion.
 *
 * NaN at a degenerate angle, deliberately -- see sharkoLpStrict. Callers here
 * test for it and write a CIF "?" or drop the row; substituting 1 would put a
 * finite, wrong |F|^2 into a file someone else will read as a measurement.
 */
function lorentzPolarization(tthDeg, pol) {
    return sharkoLpStrict(tthDeg, normalizePolarization(pol));
}

function getPolarizationSettings() {
    const modeEl = document.getElementById('polarization-mode');
    const monoEl = document.getElementById('mono-2theta');
    const fracEl = document.getElementById('pol-fraction');
    return normalizePolarization({
        mode:    modeEl ? modeEl.value : POLARIZATION_DEFAULTS.mode,
        monoTth: monoEl ? parseFloat(monoEl.value) : POLARIZATION_DEFAULTS.monoTth,
        polFrac: fracEl ? parseFloat(fracEl.value) : POLARIZATION_DEFAULTS.polFrac
    });
}

function polarizationFromParams(params) {
    if (params && params.polarization) return normalizePolarization(params.polarization);
    return getPolarizationSettings();
}

function polarizationDescription(pol) {
    const P = normalizePolarization(pol);
    if (P.mode === 'none') return 'None (Lorentz factor only, P = 1)';
    if (P.mode === 'synchrotron')
        return `Synchrotron (perpendicular polarised fraction f = ${P.polFrac.toFixed(3)})`;
    return P.monoTth > 0
        ? `Lab X-ray tube with monochromator, 2th_M = ${P.monoTth.toFixed(2)} deg`
        : 'Lab X-ray tube, no monochromator (unpolarised beam)';
}

function polarizationFormulaLines(pol) {
    const P = normalizePolarization(pol);
    return [
        'Lp(2th) = L(2th) * P(2th)',
        '  L(2th) = 1 / ( sin^2(th) * cos(th) )',
        `  P(2th) = ( 1 + K * cos^2(2th) ) / ( 1 + K ),   K = ${P.K.toFixed(6)}`
    ];
}

/**
 * Lp for the solution paths, where a bad angle must not drop a reflection.
 *
 * The lenient policy, reached through the strict one so that "not finite" is
 * decided in exactly one place. Equivalent to sharkoLp(); kept as a name
 * because charge flipping, the Wyckoff launcher and the space-group test all
 * call it and the name says which policy they are choosing.
 */
function cfLorentzPolarization(tthDeg, pol) {
    const lp = lorentzPolarization(tthDeg, pol);
    return (Number.isFinite(lp) && lp > 0) ? lp : 1;
}

/**
 * Complementary error function, Numerical Recipes' erfcc.
 *
 * THE FIFTH COEFFICIENT IS -0.18628806, NOT -1.18628806.
 *
 * It was the latter here, which is not an inaccuracy but a different function:
 * erfc(0) must be exactly 1 and came out as 0.368, an error of 63%. The whole
 * approximation is wrong across the range that matters -- 34% at x = 0.5, 18%
 * at x = 1 -- and it is worst in the middle, where a weak reflection's
 * posterior actually lives.
 *
 * It never fired, because powder5.html carried a byte-identical copy of this
 * function apart from that digit, and its inline script is evaluated after
 * this file. Load order was the only thing standing between the French-Wilson
 * correction and a broken normal CDF: reorder the script tags, or delete the
 * inline copy as the redundant duplicate it appears to be, and every weak
 * reflection silently gets the wrong posterior intensity. The duplicate has
 * now been removed from powder5.html and this is the single definition, so the
 * digit below is no longer covered by anything.
 */
function cfErfc(x) {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const ans = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
        t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
        t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? ans : 2 - ans;
}

function cfNormPdf(t) { return Math.exp(-0.5 * t * t) / Math.sqrt(2 * Math.PI); }
function cfNormCdf(t) { return 0.5 * cfErfc(-t / Math.SQRT2); }

function cfMillsRatio(t) {
    if (t > -6) {
        const c = cfNormCdf(t);
        return c > 1e-300 ? cfNormPdf(t) / c : -t;
    }
    const t2 = t * t;
    return -t * (1 + 1 / t2 - 2 / (t2 * t2));
}

function cfPosteriorAcentric(E, sigma, Sigma) {
    const mu = E - sigma * sigma / Sigma;
    const t = mu / sigma;
    return sigma * (t + cfMillsRatio(t));
}

function cfPosteriorCentric(E, sigma, Sigma) {
    const hi = Math.max(E + 8 * sigma, 8 * sigma);
    if (!(hi > 0)) return 0;
    const scale = Math.max(1e-12, Math.min(sigma, Sigma));
    const n = Math.min(20000, Math.max(400, Math.ceil(hi / (scale / 8))));
    const h = hi / n;
    let maxLog = -Infinity;
    const logw = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        const J = (i + 0.5) * h;
        const lw = -0.5 * Math.log(J) - J / (2 * Sigma)
                 - (E - J) * (E - J) / (2 * sigma * sigma);
        logw[i] = lw;
        if (lw > maxLog) maxLog = lw;
    }
    let num = 0, den = 0;
    for (let i = 0; i < n; i++) {
        const w = Math.exp(logw[i] - maxLog);
        num += w * (i + 0.5) * h;
        den += w;
    }
    return den > 0 ? num / den : Math.max(0, E);
}

function cfIsCentric(h, k, l, symops) {
    if (!Array.isArray(symops) || symops.length === 0) return false;
    for (const op of symops) {
        const r = op.r;
        if (!r) continue;
        if (h * r[0] + k * r[3] + l * r[6] === -h &&
            h * r[1] + k * r[4] + l * r[7] === -k &&
            h * r[2] + k * r[5] + l * r[8] === -l) return true;
    }
    return false;
}

function frenchWilsonCorrect(list, symops, pol) {
    const usable = list.filter(r => Number.isFinite(r.sigma) && r.sigma > 0);
    if (usable.length < 8) {
        return { applied: false, reason: 'no usable standard uncertainties', corrected: 0 };
    }

    const polModel = normalizePolarization(pol);
    for (const r of list) {
        const lp = Number.isFinite(r.lp) && r.lp > 0 ? r.lp : cfLorentzPolarization(r.tth, polModel);
        const m = (Number.isFinite(r.multiplicity) && r.multiplicity > 0) ? r.multiplicity : 1;
        r._w = m * lp;
        r._E = r.area / r._w;
        r._sE = Number.isFinite(r.sigma) ? r.sigma / r._w : NaN;
        r._q = Number.isFinite(r.d) && r.d > 0 ? 1 / (r.d * r.d) : 0;
    }

    const order = list.slice().sort((a, b) => a._q - b._q);
    const perShell = Math.max(20, Math.ceil(order.length / 12));
    let corrected = 0, negatives = 0;

    for (let start = 0; start < order.length; start += perShell) {
        const shell = order.slice(start, Math.min(start + perShell, order.length));
        if (shell.length === 0) continue;
        let sum = 0;
        for (const r of shell) sum += r._E;
        const Sigma = sum / shell.length;
        if (!(Sigma > 0)) continue;

        for (const r of shell) {
            if (!(r._sE > 0)) continue;
            if (r._E < 0) negatives++;
            if (r._E > 3 * r._sE) continue;
            const post = cfIsCentric(r.h, r.k, r.l, symops)
                ? cfPosteriorCentric(r._E, r._sE, Sigma)
                : cfPosteriorAcentric(r._E, r._sE, Sigma);
            if (Number.isFinite(post) && post >= 0) {
                r.area = post * r._w;
                corrected++;
            }
        }
    }

    for (const r of list) { delete r._w; delete r._E; delete r._sE; delete r._q; }
    return { applied: true, corrected, negatives, weak: corrected };
}