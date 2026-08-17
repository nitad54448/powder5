// Stubs for the globals data_io.js reaches for.
global.polarizationFromParams = () => ({ K: 1 });
global.lorentzPolarization = (tth, pol) => {
    const th = tth * Math.PI / 360, t2 = tth * Math.PI / 180;
    return (1 / (Math.sin(th) ** 2 * Math.cos(th))) *
           ((1 + pol.K * Math.cos(t2) ** 2) / (1 + pol.K));
};
global.calculatePeakShift = () => 0;
global.integratedPeakArea = (hkl) => hkl.__area;

const io = require('/home/claude/data_io.js');

const params = {
    a: 8.48180, b: 5.39980, c: 6.96140, alpha: 90, beta: 90, gamma: 90,
    lambda: 1.540560, lambda2: 1.544390, ratio: 0.497, zeroShift: -0.00450,
    profileType: 'tch_aniso', system: 'orthorhombic',
    U: 0.0123, V: -0.0045, W: 0.0067, X: 0.0891, Y: 0.0, SL: 0.0021, HL: 0.0013,
    S400: 1.5e-5, S040: 2.5e-5, S004: 3.5e-5, S220: 0, S202: 0, S022: 0,
    polarization: { model: 'lab', fraction: 0.5 }
};
const hklList = [
    { h_orig:0, k_orig:1, l_orig:1, d:4.26667, tth:20.7739, multiplicity:4,
      intensity:11260.628, __area:2072.4, hkl_list:['(0,1,1)'] },
    { h_orig:0, k_orig:0, l_orig:3, d:2.29529, tth:39.2,    multiplicity:2,
      intensity:-90.0,     __area:-15.3,  hkl_list:['(0,0,3)'] },
];
const bundle = {
    tth: [10, 10.02, 10.04, 10.06], obs: [100, 110, 120, 130],
    calc: [101, 111, 119, 131], bkg: [10, 10, 10, 10],
    params, spaceGroup: { symbol: 'Pnma', number: 62 }, hklList,
    stats: { rwp: 5.43, r_p: 4.01, chi2: 0.052 },
    sourceName: 'roundtrip', scaleFactor: 1,
    polarisationUi: { mode: 'lab', monoTth: 0, fraction: 0.5 }
};

const text = io.writePdCIF(bundle);
const back = io.parsePdCifFile(text);

let bad = 0;
const chk = (name, want, got, tol = 1e-6) => {
    const ok = (typeof want === 'number')
        ? (Number.isFinite(got) && Math.abs(got - want) <= tol)
        : String(got) === String(want);
    if (!ok) bad++;
    console.log(`  ${ok ? 'PASS' : 'FAIL'}  ${name.padEnd(26)} wrote ${String(want).padStart(11)}   read ${String(got).padStart(11)}`);
};

console.log('cell');
for (const k of ['a','b','c','alpha','beta','gamma']) chk(k, params[k], back.cell && back.cell[k], 1e-5);
console.log('space group');
chk('symbol', 'Pnma', back.spaceGroup && back.spaceGroup.symbol);
chk('IT number', 62, back.spaceGroup && back.spaceGroup.number);
console.log('radiation');
chk('lambda1', params.lambda, back.wavelength, 1e-6);
chk('lambda2', params.lambda2, back.wavelength2, 1e-6);
chk('ratio I2/I1', params.ratio, back.ratio21, 1e-5);
console.log('zero point');
chk('zeroShift', params.zeroShift, back.zeroShift, 1e-5);
console.log('profile');
chk('profileType', 'tch_aniso', back.profileType);
for (const k of ['U','V','W','X','SL','HL','S400','S040','S004'])
    chk('param ' + k, params[k], back.profileParams[k], 1e-6);
console.log('polarisation');
chk('mode', 'lab', back.polarisation && back.polarisation.mode);
chk('mono 2theta', 0, back.polarisation && back.polarisation.monoTth, 1e-6);
chk('fraction', 0.5, back.polarisation && back.polarisation.fraction, 1e-6);
console.log('profile data');
chk('n points', 4, back.intensity.length);
chk('first 2theta', 10, back.tth[0], 1e-6);

console.log('\nreflection loop as written:');
text.split('\n').filter(l => /^  [-\d]/.test(l) && l.split(/\s+/).length === 8)
    .forEach(l => console.log('   ', l.trim()));

// The doublet must vanish entirely at ratio 0.
const mono = io.parsePdCifFile(io.writePdCIF({ ...bundle,
    params: { ...params, ratio: 0 } }));
console.log('\nratio = 0 export:');
chk('lambda2 absent', 'null', String(mono.wavelength2));
chk('ratio absent', 'null', String(mono.ratio21));

console.log(bad ? `\n${bad} FAILURE(S)` : '\nround trip complete: every field survived');
process.exit(bad ? 1 : 0);
