// The real pipeline: parsers.js -> detectAndParseFile -> normalizeParsedData,
// exercised on the actual file the user exported.
global.polarizationFromParams = () => ({ K: 1 });
global.lorentzPolarization = () => 1;
global.calculatePeakShift = () => 0;
global.integratedPeakArea = () => 1;
const io = require('./data_io.js');
global.parsePdCifFile = io.parsePdCifFile;

// parsers.js declares its helpers with const inside a block, so evaluate it and
// reach the two entry points the loader uses.
const src = require('fs').readFileSync('./parsers.js', 'utf8');
const factory = new Function('parsePdCifFile', src + `
    return { detectAndParseFile, normalizeParsedData };`);
const { detectAndParseFile, normalizeParsedData } = factory(io.parsePdCifFile);

const cif = require('fs').readFileSync(require('path').join(__dirname, 'pdcif_pipeline_sample.cif'), 'utf8');
const rawParsed = detectAndParseFile('PBSO4.cif', cif);
const parsedData = normalizeParsedData(rawParsed, 'PBSO4.cif');

let bad = 0;
const chk = (name, got, want) => {
    const ok = JSON.stringify(got) === JSON.stringify(want);
    if (!ok) bad++;
    console.log(`  ${ok ? 'PASS' : 'FAIL'}  ${name.padEnd(14)} ${JSON.stringify(got)}`);
};
console.log('what reaches applyPdCifMetadata(parsedData):');
chk('cell', parsedData.cell, { a:8.48181, b:5.39977, c:6.96142, alpha:90, beta:90, gamma:90 });
chk('spaceGroup', parsedData.spaceGroup, { symbol:'Pnma', number:62 });
chk('wavelength', parsedData.wavelength, 1.54056);
chk('wavelength2', parsedData.wavelength2, 1.54439);
chk('ratio21', parsedData.ratio21, 0.497);
chk('zeroShift', parsedData.zeroShift, -0.02783);
chk('profileType', parsedData.profileType, 'split_pvoigt');
chk('profileParams', Object.keys(parsedData.profileParams).length, 11);
chk('polarisation', parsedData.polarisation, { mode:'lab', monoTth:0, fraction:0.95 });
chk('points', parsedData.tth.length, 3);
console.log('\nprofile parameters recovered:');
for (const [k, v] of Object.entries(parsedData.profileParams)) {
    const id = 'param-' + k.toLowerCase();
    console.log(`   ${k.padEnd(12)} ${String(v).padStart(12)}   -> #${id}`);
}
console.log(bad ? `\n${bad} FAILURE(S)` : '\nnothing lost between the reader and the loader');
process.exit(bad ? 1 : 0);
