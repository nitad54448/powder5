// test/test_fft.js
// ---------------------------------------------------------------------------
//  Validates the table-driven radix-2 FFT in charge_flipping_worker.js:
//    1. against a reference DFT, forward and inverse;
//    2. round-trip identity on a 3D grid;
//    3. quantified against the OLD recurrence-based twiddles, so the accuracy
//       claim in the code comment is a measurement and not an assertion.
//
//  Run: node test/test_fft.js
// ---------------------------------------------------------------------------

'use strict';

const fs = require('fs');
const path = require('path');

// Pull fft1d / fftTwiddles / fft3d out of the worker without executing the
// rest of it (the worker's top level touches fetch/importScripts).
const src = fs.readFileSync(
    path.join(__dirname, '..', 'charge_flipping_worker.js'), 'utf8'
).replace(/\r\n/g, '\n');
const a = src.indexOf('function fft1d(');
if (a < 0) throw new Error('fft1d not found -- has the worker been restructured?');
const fft3dStart = src.indexOf('function fft3d', a);
const fft3dEnd = src.indexOf('\n}\n', src.indexOf('if (inverse) {', fft3dStart)) + 3;
const chunk = src.slice(a, fft3dEnd);
// Evaluate in an isolated scope and hand back the two entry points.
const { fft1d, fft3d } = new Function(chunk + '\nreturn { fft1d, fft3d };')();

// The precomputed-table alternative, kept ONLY so the comparison documented in
// the comment above fft1d in charge_flipping_worker.js stays reproducible. It
// is deliberately NOT what ships: it measured 12-15% slower with identical
// accuracy. If you change the production FFT, re-run this and update that
// comment with the new numbers.
const _twCache = new Map();
function twiddles(n, inverse) {
    const key = n * 2 + (inverse ? 1 : 0);
    let t = _twCache.get(key);
    if (t) return t;
    const cosT = new Float64Array(Math.max(1, n));
    const sinT = new Float64Array(Math.max(1, n));
    const sgn = inverse ? 1 : -1;
    for (let len = 2; len <= n; len <<= 1) {
        const half = len >> 1, off = half - 1, ang = sgn * 2 * Math.PI / len;
        for (let k = 0; k < half; k++) {
            cosT[off + k] = Math.cos(ang * k);
            sinT[off + k] = Math.sin(ang * k);
        }
    }
    t = { cosT, sinT };
    _twCache.set(key, t);
    return t;
}
function fft1d_table(re, im, base, stride, n, inverse) {
    for (let i = 1, j = 0; i < n; i++) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            const p = base + i * stride, q = base + j * stride;
            let t = re[p]; re[p] = re[q]; re[q] = t;
            t = im[p]; im[p] = im[q]; im[q] = t;
        }
    }
    const { cosT, sinT } = twiddles(n, inverse);
    for (let len = 2; len <= n; len <<= 1) {
        const half = len >> 1, off = half - 1;
        for (let i = 0; i < n; i += len) {
            for (let k = 0; k < half; k++) {
                const cr = cosT[off + k], ci = sinT[off + k];
                const iu = base + (i + k) * stride;
                const iv = base + (i + k + half) * stride;
                const uRe = re[iu], uIm = im[iu];
                const vRe = re[iv] * cr - im[iv] * ci;
                const vIm = re[iv] * ci + im[iv] * cr;
                re[iu] = uRe + vRe; im[iu] = uIm + vIm;
                re[iv] = uRe - vRe; im[iv] = uIm - vIm;
            }
        }
    }
}

function makeSignal(n, seed) {
    let t = seed >>> 0;
    const rnd = () => {
        t = (t + 0x6D2B79F5) | 0;
        let r = Math.imul(t ^ (t >>> 15), 1 | t);
        r = (r + Math.imul(r ^ (r >>> 7), 61 | r)) ^ r;
        return ((r ^ (r >>> 14)) >>> 0) / 4294967296 - 0.5;
    };
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let i = 0; i < n; i++) { re[i] = rnd(); im[i] = rnd(); }
    return { re, im };
}

/** Reference DFT in double precision. */
function dft(re, im, inverse) {
    const n = re.length;
    const R = new Float64Array(n), I = new Float64Array(n);
    const sgn = inverse ? 1 : -1;
    for (let k = 0; k < n; k++) {
        let sr = 0, si = 0;
        for (let t = 0; t < n; t++) {
            const ang = sgn * 2 * Math.PI * k * t / n;
            const c = Math.cos(ang), s = Math.sin(ang);
            sr += re[t] * c - im[t] * s;
            si += re[t] * s + im[t] * c;
        }
        R[k] = sr; I[k] = si;
    }
    return { R, I };
}

function relErr(re, im, R, I) {
    let err = 0, norm = 0;
    for (let k = 0; k < re.length; k++) {
        err = Math.max(err, Math.hypot(re[k] - R[k], im[k] - I[k]));
        norm = Math.max(norm, Math.hypot(R[k], I[k]));
    }
    return err / norm;
}

let failures = 0;
const TOL = 1e-12;

console.log('1. forward and inverse against a reference DFT\n');
console.log('        n   production   table alt   ratio');
for (const n of [64, 128, 256, 512, 1024, 2048]) {
    for (const inverse of [false, true]) {
        const sig = makeSignal(n, 12345 + n);
        const ref = dft(sig.re, sig.im, inverse);

        const t1 = sig.re.slice(), t2 = sig.im.slice();
        fft1d(t1, t2, 0, 1, n, inverse);
        const eNew = relErr(t1, t2, ref.R, ref.I);

        const r1 = sig.re.slice(), r2 = sig.im.slice();
        fft1d_table(r1, r2, 0, 1, n, inverse);
        const eOld = relErr(r1, r2, ref.R, ref.I);   // 'Old' = the table alternative

        if (!inverse) {
            console.log(`  ${String(n).padStart(7)}   ${eNew.toExponential(2)}    ` +
                        `${eOld.toExponential(2)}     ${(eOld / eNew).toFixed(1)}x`);
        }
        if (!(eNew < TOL)) {
            failures++;
            console.error(`  FAIL n=${n} inverse=${inverse}: rel err ${eNew}`);
        }
        // The whole point of the measurement: the two must agree to within a
        // small factor. If the table ever becomes MUCH better, the comment
        // above fft1d in the worker is out of date and should be revisited.
        if (!(eNew <= eOld * 3)) {
            failures++;
            console.error(`  FAIL n=${n}: production FFT is materially less accurate ` +
                          `than the table alternative (${eNew} vs ${eOld}) -- ` +
                          `revisit the note above fft1d.`);
        }
    }
}

console.log('\n2. strided access (the 3D transform uses stride N and N*N)\n');
for (const n of [64, 256]) {
    for (const stride of [1, 3, 7]) {
        const flatRe = new Float64Array(n * stride + 5);
        const flatIm = new Float64Array(n * stride + 5);
        const sig = makeSignal(n, 999 + stride);
        const base = 2;
        for (let i = 0; i < n; i++) {
            flatRe[base + i * stride] = sig.re[i];
            flatIm[base + i * stride] = sig.im[i];
        }
        fft1d(flatRe, flatIm, base, stride, n, false);
        const got = new Float64Array(n), gotI = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            got[i] = flatRe[base + i * stride];
            gotI[i] = flatIm[base + i * stride];
        }
        const ref = dft(sig.re, sig.im, false);
        const e = relErr(got, gotI, ref.R, ref.I);
        console.log(`  n=${String(n).padStart(4)} stride=${stride}   rel err ${e.toExponential(2)}`);
        if (!(e < TOL)) { failures++; console.error(`  FAIL stride ${stride}`); }
        // Elements outside the transform must be untouched.
        for (let i = 0; i < flatRe.length; i++) {
            const inSpan = (i >= base) && ((i - base) % stride === 0) && ((i - base) / stride < n);
            if (!inSpan && (flatRe[i] !== 0 || flatIm[i] !== 0)) {
                failures++;
                console.error(`  FAIL stride ${stride}: wrote outside its span at ${i}`);
                break;
            }
        }
    }
}

console.log('\n3. 3D round trip on an N^3 grid\n');
for (const N of [16, 32, 64]) {
    const N3 = N * N * N;
    const sig = makeSignal(N3, 4242 + N);
    const re = sig.re.slice(), im = sig.im.slice();
    fft3d(re, im, N, false);
    fft3d(re, im, N, true);
    let worst = 0, scale = 0;
    for (let i = 0; i < N3; i++) {
        worst = Math.max(worst, Math.hypot(re[i] - sig.re[i], im[i] - sig.im[i]));
        scale = Math.max(scale, Math.hypot(sig.re[i], sig.im[i]));
    }
    const e = worst / scale;
    console.log(`  ${N}^3 forward+inverse   rel err ${e.toExponential(2)}`);
    if (!(e < TOL)) { failures++; console.error(`  FAIL round trip N=${N}`); }
}

console.log('\n4. Parseval (energy conservation), 3D\n');
for (const N of [16, 32]) {
    const N3 = N * N * N;
    const sig = makeSignal(N3, 77 + N);
    const re = sig.re.slice(), im = sig.im.slice();
    let eIn = 0;
    for (let i = 0; i < N3; i++) eIn += sig.re[i] ** 2 + sig.im[i] ** 2;
    fft3d(re, im, N, false);
    let eOut = 0;
    for (let i = 0; i < N3; i++) eOut += re[i] ** 2 + im[i] ** 2;
    const ratio = eOut / (eIn * N3);
    console.log(`  ${N}^3   sum|F|^2 / (N^3 sum|f|^2) = ${ratio.toFixed(15)}`);
    if (Math.abs(ratio - 1) > 1e-12) { failures++; console.error(`  FAIL Parseval N=${N}`); }
}

console.log('\n5. throughput, 64^3 forward transform\n');
{
    const N = 64, N3 = N * N * N;
    const sig = makeSignal(N3, 5);
    const warm = { re: sig.re.slice(), im: sig.im.slice() };
    fft3d(warm.re, warm.im, N, false);

    const time = (fn) => {
        const re = sig.re.slice(), im = sig.im.slice();
        const t0 = process.hrtime.bigint();
        fn(re, im);
        return Number(process.hrtime.bigint() - t0) / 1e6;
    };
    const tNew = time((re, im) => fft3d(re, im, N, false));

    // fft3d wired to the old inner routine, for a like-for-like comparison.
    const fft3dTable = (re, im, N) => {
        const N2 = N * N;
        for (let l = 0; l < N; l++) for (let k = 0; k < N; k++) fft1d_table(re, im, k * N + l * N2, 1, N, false);
        for (let l = 0; l < N; l++) for (let h = 0; h < N; h++) fft1d_table(re, im, h + l * N2, N, N, false);
        for (let k = 0; k < N; k++) for (let h = 0; h < N; h++) fft1d_table(re, im, h + k * N, N2, N, false);
    };
    const tOld = time((re, im) => fft3dTable(re, im, N));
    console.log(`  production (recurrence)  ${tNew.toFixed(1)} ms`);
    console.log(`  table alternative        ${tOld.toFixed(1)} ms   ` +
                `(${(tOld / tNew).toFixed(2)}x production)`);
    console.log(tOld > tNew
        ? '  -> the table is slower, as documented above fft1d. Leave it alone.'
        : '  -> the table now wins on this engine; revisit the note above fft1d.');
}

console.log(`\n${failures} failures`);
process.exit(failures ? 1 : 0);
