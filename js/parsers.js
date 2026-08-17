// ---------------------------------------------------------------------------
//  ZIP CONTAINERS
//
//  .brml (Bruker) and .rasx (Rigaku) are not XML files. They are ZIP archives
//  with XML and text INSIDE them, and the loader read every file as text -- so
//  a .brml arrived as the raw deflate stream, its extension matched the BRML
//  rule, DOMParser was handed binary and threw, and the file fell through to
//  the two-column fallback which found nothing. The registry entry for .rasx
//  additionally routed to the BRML parser, which would not have worked even on
//  an unzipped one: RASX keeps its intensities in a plain text profile, not in
//  the XML at all.
//
//  Requires lib/jszip.min.js. Without it the archive formats are refused with
//  a message that says so, rather than failing as a parse error.
// ---------------------------------------------------------------------------

/** ZIP local file header magic, "PK\x03\x04". */
const isZipBuffer = (buf) => {
    const b = new Uint8Array(buf, 0, Math.min(4, buf.byteLength));
    return b.length >= 4 && b[0] === 0x50 && b[1] === 0x4B && b[2] === 0x03 && b[3] === 0x04;
};

/**
 * Bytes to text, without mangling the legacy formats.
 *
 * Most instrument text files are ASCII, but the older ones carry degree signs
 * and micro symbols in a single-byte codepage. Strict UTF-8 would reject those
 * bytes and lenient UTF-8 turns them into replacement characters, so a failed
 * UTF-8 decode falls back to windows-1252 rather than keeping the damage.
 */
const decodeTextBuffer = (buf) => {
    let s = new TextDecoder('utf-8', { fatal: false }).decode(buf);
    if (s.indexOf('\uFFFD') >= 0) {
        try { s = new TextDecoder('windows-1252').decode(buf); } catch (e) { /* keep utf-8 */ }
    }
    return s;
};

/**
 * Pull a wavelength out of whatever XML an archive happens to carry.
 *
 * Every vendor spells it differently and the surrounding schema is large, so
 * this looks for the number rather than the structure. Bounded to 0.3-3.0 A:
 * that covers every laboratory anode and neutron instrument, and rejects the
 * unrelated numbers -- goniometer radii, slit widths, timestamps -- that a
 * looser pattern would pick up.
 */
const sniffWavelengthFromXml = (xml) => {
    const pats = [
        /kAlpha1[^0-9\-]{0,40}([0-9]*\.[0-9]+)/i,
        /<Wavelength[^>]*>\s*([0-9]*\.[0-9]+)/i,
        /Wavelength[^0-9\-]{0,40}([0-9]*\.[0-9]+)/i
    ];
    for (const p of pats) {
        const m = xml.match(p);
        if (m) {
            let v = parseFloat(m[1]);
            // Rigaku states it in nm in some files.
            if (v > 0.03 && v < 0.3) v *= 10;
            if (v >= 0.3 && v <= 3.0) return v;
        }
    }
    return null;
};

/**
 * Parse an unzipped archive, given its entries as {name, text}.
 *
 * Tried in order of how specifically the entry names identify the vendor, so a
 * container that happens to hold several readable things is read as the format
 * it actually is rather than as whichever entry sorted first.
 *
 * @param {{name: string, text: string}[]} entries
 * @param {string} fileName
 */
const parseZipEntries = (entries, fileName) => {
    const xml = entries.filter(e => /\.xml$/i.test(e.name));
    const problems = [];

    //   Bruker BRML: the scan lives in Experiment*/RawData*.xml
    const raw = entries.filter(e => /rawdata.*\.xml$/i.test(e.name));
    for (const e of raw) {
        try {
            const r = parseBrukerBrmlFile(e.text);
            if (r && r.tth && r.tth.length) {
                if (!r.wavelength) {
                    for (const x of xml) {
                        const wl = sniffWavelengthFromXml(x.text);
                        if (wl) { r.wavelength = wl; break; }
                    }
                }
                return r;
            }
        } catch (err) { problems.push(`${e.name}: ${err.message}`); }
    }

    //   Rigaku RASX: intensities in Data*/Profile*.txt, conditions in XML
    const prof = entries.filter(e => /profile.*\.txt$/i.test(e.name));
    for (const e of prof) {
        try {
            const r = parseDataFile(e.text, `${fileName}:${e.name}`);
            if (r && r.tth && r.tth.length) {
                for (const x of xml) {
                    const wl = sniffWavelengthFromXml(x.text);
                    if (wl) { r.wavelength = wl; break; }
                }
                return r;
            }
        } catch (err) { problems.push(`${e.name}: ${err.message}`); }
    }

    //   An XRDML or anything else the ordinary detector recognises.
    for (const e of entries) {
        // Skip the metadata files every container carries; running the
        // two-column fallback over them yields a few stray numbers that look
        // like a one-point pattern.
        if (/(^|\/)(\[Content_Types\]|_rels|docProps)/i.test(e.name)) continue;
        try {
            const r = detectAndParseFile(e.name, e.text);
            if (r && r.tth && r.tth.length > 1) return r;
        } catch (err) { problems.push(`${e.name}: ${err.message}`); }
    }

    throw new Error(
        `${fileName} is a ZIP archive, but none of its ${entries.length} entries ` +
        `held a readable powder scan.` +
        (problems.length ? ` Tried: ${problems.slice(0, 3).join('; ')}` : ''));
};

/**
 * The loader's entry point: bytes in, parsed pattern out.
 *
 * Async because unzipping is. Everything that is not an archive takes the same
 * synchronous path it always did, decoded from the buffer instead of being
 * handed to FileReader.readAsText.
 *
 * @param {string} fileName
 * @param {ArrayBuffer} buffer
 * @returns {Promise<object>}
 */
const detectAndParseBuffer = async (fileName, buffer) => {
    if (!isZipBuffer(buffer)) {
        return detectAndParseFile(fileName, decodeTextBuffer(buffer));
    }
    if (typeof JSZip === 'undefined') {
        throw new Error(
            `${fileName} is a ZIP archive (.brml and .rasx both are), and the unzip ` +
            `library is not loaded. Add lib/jszip.min.js, or unzip the file and open ` +
            `the RawData XML or Profile text inside it.`);
    }
    const zip = await JSZip.loadAsync(buffer);
    const entries = [];
    const names = Object.keys(zip.files).filter(n => !zip.files[n].dir);
    for (const n of names) {
        // Only the text members are of interest; an embedded image or a
        // thumbnail decoded as text is megabytes of noise for the parsers to
        // walk through.
        if (!/\.(xml|txt|csv|dat|asc|raw|xrdml|ras|rels)$/i.test(n) && !/rawdata/i.test(n)) continue;
        entries.push({ name: n, text: await zip.files[n].async('string') });
    }
    if (!entries.length) {
        throw new Error(`${fileName} is a ZIP archive with no text entries to read.`);
    }
    return parseZipEntries(entries, fileName);
};

/**
         * Smart file detector. It checks for known headers and extensions
         * and falls back to a generic 2-column parser.
         */
        const detectAndParseFile = (fileName, fileContent) => {
            const name = fileName.toLowerCase();
            const lines = fileContent.trim().split(/\r?\n/);
            const firstLine = lines.length > 0 ? lines[0].trim() : '';
            const upperContent = fileContent.substring(0, 500).toUpperCase(); // Check first 500 chars

            //   Parser Registry  
            const PARSER_REGISTRY = [
                { // pdCIF -- must be tested BEFORE the generic fallback, and
                  // before the XML tests, since a CIF is neither XML nor two
                  // columns and the fallback would read its loop headers as data.
                    test: (name, content) => name.endsWith('.cif') ||
                                             /\bdata_/.test(content) &&
                                             /_pd_(meas|proc)_/.test(content),
                    // RETURNED WHOLE, NOT WHITELISTED.
                    //
                    // This used to rebuild the result as
                    //     { tth, intensity, wavelength, zeroShift }
                    // which silently discarded everything else the reader had
                    // extracted -- the cell, the space group, the Ka2
                    // wavelength and its ratio, the profile function and its
                    // parameters, the polarisation model. A pdCIF written by
                    // this program therefore round-tripped its pattern, its
                    // wavelength and its zero point, and threw away the entire
                    // refinement that had produced them, with nothing anywhere
                    // reporting a loss.
                    //
                    // A whitelist here is the wrong shape of code: it has to be
                    // edited every time the reader learns a new tag, and the
                    // failure when someone forgets is invisible. parsePdCifFile
                    // already returns exactly the fields it means to return, so
                    // this hands them on untouched. normalizeParsedData() below
                    // spreads the object rather than rebuilding it, so they
                    // survive the rest of the way.
                    parser: (content) => parsePdCifFile(content)
                },
                { // XRDML
                    test: (name, content) => name.endsWith('.xrdml') || (content.includes('<?xml') && content.includes('<xrdMeasurement')),
                    parser: parseXrdmlFile
                },
                { // BRML
                    test: (name, content) => name.endsWith('.brml') || (content.includes('<?xml') && content.includes('<RawDataFile')),
                    parser: parseBrukerBrmlFile
                },
                { // Rigaku RASX, already unzipped by hand.
                  //
                  // This used to route to parseBrukerBrmlFile, which looks for
                  // <RawDataFile> and Bruker's TwoTheta axis attributes -- tags
                  // a RASX does not contain. RASX keeps its intensities in a
                  // plain two-column Profile text member, so an unzipped one
                  // reaching this point is either that text or an XML that only
                  // carries conditions. Both are handled by falling through.
                    test: (name, content) => name.endsWith('.rasx') && !content.includes('<?xml'),
                    parser: parseDataFile
                },
                { // UXD
                    test: (name, content, firstLine) => name.endsWith('.uxd') || firstLine.startsWith('_FILEVERSION'),
                    parser: parseUxdFile
                },
                { // Rigaku RAS
                    test: (name, content, firstLine, upper) => name.endsWith('.ras') || upper.includes('*RAS_HEADER_START'),
                    parser: parseRigakuRasFile
                },
                { // Philips UDF/RD/SD
                    test: (name) => name.endsWith('.udf') || name.endsWith('.rd') || name.endsWith('.sd'),
                    parser: parsePhilipsUdfFile
                },
                { // GSAS ESD/XRA
                    test: (name, content, firstLine, upper, allLines) => allLines.some(line => line.trim().toUpperCase().startsWith('BANK')),
                    parser: (content, allLines) => {
                        const bankLine = allLines.find(line => line.trim().toUpperCase().startsWith('BANK'));
                        if (bankLine.toUpperCase().includes('STD')) {
                            return parseGsasXraFile(content);
                        }
                        return parseGsasEsdFile(content);
                    }
                },
                { // Jade MDI (treat as 2-column)
                     test: (name, content, firstLine, upper) => name.endsWith('.mdi') && (upper.includes('2-THETA, INTENSITY') || upper.startsWith('(SAMPLE')),
                     parser: parseDataFile
                }
            ];
            
            //   Iterate Registry  
            for (const rule of PARSER_REGISTRY) {
                try {
                    if (rule.test(name, fileContent, firstLine, upperContent, lines)) {
                        // Pass 'content' to parser, but 'lines' to the special GSAS one
                        if (rule.parser.length > 1) {
                             return rule.parser(fileContent, lines); // For GSAS parser
                        }
                        return rule.parser(fileContent);
                    }
                } catch (e) {
                    console.warn(`Parser ${rule.parser.name} failed, trying next...`, e.message);
                }
            }

            //   Fallback for all other 2-column-like formats  
            // This will attempt to parse: .xy, .csv, .txt, .dat, .asc, etc.
            return parseDataFile(fileContent, fileName);
        };
    
        /**
         * Final gate every parser's output passes through before it becomes
         * fullExperimentalData.
         *
         * Two things were previously only warned about in the console and then
         * used anyway:
         *   - Non-finite points. A single NaN 2-theta poisons the axis range and
         *     every downstream binary search.
         *   - Descending 2-theta. lowerBound() is a binary search that requires
         *     an ASCENDING array, and the loader takes tth[0] / tth[n-1] as the
         *     slider min / max. A descending file therefore produced a negative
         *     slider step and silently wrong pattern slicing rather than an error.
         */
        const normalizeParsedData = (parsed, fileName = "") => {
            if (!parsed || !parsed.tth || !parsed.intensity) {
                throw new Error(`Parser returned no data for ${fileName}.`);
            }
            const n = Math.min(parsed.tth.length, parsed.intensity.length);
            if (parsed.tth.length !== parsed.intensity.length) {
                console.warn(`Data File (${fileName}) Warning: 2-theta and intensity counts differ ` +
                             `(${parsed.tth.length} vs ${parsed.intensity.length}); using the first ${n}.`);
            }

            const pts = [];
            let dropped = 0;
            for (let i = 0; i < n; i++) {
                const x = parsed.tth[i], y = parsed.intensity[i];
                if (isFinite(x) && isFinite(y)) pts.push({ x, y }); else dropped++;
            }
            if (dropped > 0) {
                console.warn(`Data File (${fileName}) Warning: discarded ${dropped} non-finite point(s).`);
            }
            if (pts.length === 0) {
                throw new Error(`No usable numeric data points were found in ${fileName}.`);
            }

            let sorted = true;
            for (let i = 1; i < pts.length; i++) {
                if (pts[i].x < pts[i - 1].x) { sorted = false; break; }
            }
            if (!sorted) {
                pts.sort((a, b) => a.x - b.x);
                console.warn(`Data File (${fileName}) Notice: 2-theta was not ascending; the data has been sorted.`);
            }

            // Collapse exact duplicate 2-theta values: they give a zero-width
            // interval that the spline interpolator cannot handle.
            const outX = [], outY = [];
            let duplicates = 0;
            for (let i = 0; i < pts.length; i++) {
                if (i > 0 && pts[i].x === pts[i - 1].x) { duplicates++; continue; }
                outX.push(pts[i].x); outY.push(pts[i].y);
            }
            if (duplicates > 0) {
                console.warn(`Data File (${fileName}) Notice: removed ${duplicates} duplicate 2-theta point(s).`);
            }
            if (outX.length < 2) {
                throw new Error(`Only ${outX.length} usable data point(s) in ${fileName}; at least two are required.`);
            }

            return { ...parsed, tth: outX, intensity: outY };
        };

        /**
         * Generic 2-column parser. This is the fallback for most text files.
         * Includes validation logic for 2-theta (X) and step size (dX).
         */
        const parseDataFile = (text, fileName = "") => {
            const lines = text.trim().split(/\r?\n/);
            const tth = [], intensity = [];
            let last_x = -Infinity;
            let suspicious_steps = 0;
            let positive_x_values = 0;
            let negative_steps = 0;
            let headerLines = 0;
            let dataStarted = false;

            lines.forEach(line => {
                const trimmed = line.trim();

                // Skip commented or empty lines
                if (trimmed === '' || trimmed.startsWith('#') || trimmed.startsWith('//') ||
                    trimmed.startsWith('!') || trimmed.startsWith(';')) {
                    if (!dataStarted) headerLines++;
                    return;
                }

                // FIX: header detection used to reject any line containing a
                // letter BEFORE attempting to parse it. The "E" of scientific
                // notation is a letter, so a file written as "10.02  1.5E+03"
                // had every single line classified as a header, leaving nothing
                // to plot and failing with "Could not parse any 2-column data".
                // Try to read the line as a number pair first; only fall back to
                // treating it as a header when that fails. A genuine header such
                // as "2Theta  Intensity" still fails to parse and is skipped.
                const parts = trimmed.split(/[\s,;]+/).filter(Boolean);
                const x = parts.length >= 2 ? parseFloat(parts[0]) : NaN;
                const y = parts.length >= 2 ? parseFloat(parts[1]) : NaN;

                if (!isFinite(x) || !isFinite(y)) {
                    if (!dataStarted) headerLines++; // Still in the header
                    return;
                }
                
                dataStarted = true; // First valid numeric pair found

                //   vérif
                if (x > 0) positive_x_values++;

                if (last_x !== -Infinity) {
                    const dX = x - last_x;
                    if (dX < 0) {
                        negative_steps++; // Data is descending
                    } else if (dX > 0 && (dX < 0.0001 || dX > 0.2)) { 
                        suspicious_steps++; // Step size is weird
                    }
                }
                last_x = x;
                

                tth.push(x);
                intensity.push(y);
            });

            //   Final checks (log warnings to console)  
            if (tth.length > 10) { 
                if (positive_x_values / tth.length < 0.5) {
                    console.warn(`Data File (${fileName}) Warning: Most 2-theta (X) values are zero or negative. This is unusual for XRD data.`);
                }
                if (negative_steps / tth.length > 0.8) {
                     console.warn(`Data File (${fileName}) Warning: Data appears to be sorted in descending 2-theta order.`);
                }
                if (suspicious_steps / tth.length > 0.2) {
                    console.warn(`Data File (${fileName}) Warning: Many data points have a step size outside the typical range (0.0001° - 0.2°). Check file format.`);
                }
            } else if (tth.length === 0) {
                 throw new Error(`Could not parse any 2-column data from ${fileName}. File may be binary or have an unknown header.`);
            }

            return { tth, intensity };
        };

        const parseXrdmlFile = (xmlString) => { const parser = new DOMParser(); const xmlDoc = parser.parseFromString(xmlString, "application/xml"); if (xmlDoc.querySelector("parsererror")) { throw new Error("Error parsing XRDML file."); } let wavelength = null; const kAlpha1Node = xmlDoc.querySelector("kAlpha1"); if (kAlpha1Node?.textContent) wavelength = parseFloat(kAlpha1Node.textContent); const intensityNode = xmlDoc.querySelector("intensities") || xmlDoc.querySelector("counts"); if (!intensityNode) throw new Error("Could not find <intensities> or <counts> in XRDML file."); const intensity = intensityNode.textContent.trim().split(/\s+/).map(Number); const positionsNode = xmlDoc.querySelector('positions[axis="2Theta"]'); if (!positionsNode) throw new Error("Could not find <positions> in XRDML file."); const startPosNode = positionsNode.querySelector("startPosition"); const endPosNode = positionsNode.querySelector("endPosition"); if (!startPosNode || !endPosNode) throw new Error("Could not find start/end positions in XRDML."); const startPos = parseFloat(startPosNode.textContent); const endPos = parseFloat(endPosNode.textContent); if (!isFinite(startPos) || !isFinite(endPos)) throw new Error("XRDML start/end positions are not numeric."); if (intensity.length < 2) throw new Error("XRDML file contains fewer than two data points."); /* FIX: (length - 1) was an unguarded divisor -> Infinity for a 1-point scan. */ const step = (endPos - startPos) / (intensity.length - 1); const tth = Array.from({ length: intensity.length }, (_, i) => startPos + i * step); return { tth, intensity, wavelength }; };
        const parseBrukerBrmlFile = (xmlString) => { const parser = new DOMParser(); const xmlDoc = parser.parseFromString(xmlString, "application/xml"); if (xmlDoc.querySelector("parsererror")) { throw new Error("Error parsing BRML file."); } let wavelength = null; const wlNode = xmlDoc.querySelector('usedWavelength'); if (wlNode) { const kAlpha1 = wlNode.getAttribute('kAlpha1'); if (kAlpha1) wavelength = parseFloat(kAlpha1); } const intensityNode = xmlDoc.querySelector("dataPoints > counts"); if (!intensityNode) throw new Error("No <counts> data found in BRML file."); const intensity = intensityNode.textContent.trim().split(/\s+/).map(Number); const startPosNode = xmlDoc.querySelector('startPosition[axis="TwoTheta"]'); const stepSizeNode = xmlDoc.querySelector('increment[axis="TwoTheta"]'); if (!startPosNode || !stepSizeNode) throw new Error("Could not find scan parameters in BRML file."); const startPos = parseFloat(startPosNode.textContent); const stepSize = parseFloat(stepSizeNode.textContent); const tth = Array.from({ length: intensity.length }, (_, i) => startPos + i * stepSize); return { tth, intensity, wavelength }; };
        const parseRigakuRasFile = (text) => { const lines = text.trim().split(/\r?\n/); const tth = [], intensity = []; let inDataSection = false; let wavelength = null; for (const line of lines) { const upperLine = line.toUpperCase(); if (upperLine.startsWith('*WAVE_LENGTH') || upperLine.startsWith('*MEAS_COND_XG_WAVE_LENGTH')) { const parts = line.trim().split(/\s+/); if (parts.length > 1) { const wl = parseFloat(parts[1]); if (!isNaN(wl)) wavelength = wl; } } if (upperLine.startsWith('*RAS_INT_START')) { inDataSection = true; continue; } if (upperLine.startsWith('*RAS_INT_END')) break; if (inDataSection) { const parts = line.trim().split(/[\s,]+/); if (parts.length >= 2) { const x = parseFloat(parts[0]); const y = parseFloat(parts[1]); if (!isNaN(x) && !isNaN(y)) { tth.push(x); intensity.push(y); } } } } if (tth.length === 0) throw new Error("No data found in RAS file data section."); return { tth, intensity, wavelength }; };
        const parseGsasEsdFile = (text) => { const lines = text.trim().split(/\r?\n/); let wavelength = null; let startTth, stepSize; let dataStartIndex = -1; lines.forEach((line, index) => { const upperLine = line.toUpperCase(); if (upperLine.includes('WAVELENGTH')) { const match = line.match(/wavelength\s+([0-9.]+)/i); if (match && match[1]) wavelength = parseFloat(match[1]); } if (upperLine.startsWith('BANK')) { const parts = line.trim().split(/\s+/); /* FIX: was >= 6 while reading parts[6]; a 6-token BANK line produced stepSize = NaN, which passed the `undefined` guard below and made every 2-theta NaN. */ if (parts.length >= 7 && parts[4].toUpperCase() === 'CONST') { startTth = parseFloat(parts[5]) / 100.0; stepSize = parseFloat(parts[6]) / 100.0; dataStartIndex = index + 1; } } }); if (!isFinite(startTth) || !isFinite(stepSize) || stepSize === 0) throw new Error("GSAS Parse Error: Could not find a valid 'BANK' line with CONST scan parameters."); if (dataStartIndex !== -1 && lines[dataStartIndex]?.toUpperCase().includes('STD')) dataStartIndex++; if (dataStartIndex === -1 || dataStartIndex >= lines.length) throw new Error("GSAS Parse Error: Found scan parameters but no subsequent data lines."); const intensity = []; for (let i = dataStartIndex; i < lines.length; i++) { const parts = lines[i].trim().split(/\s+/); for (let j = 1; j < parts.length; j += 2) { const val = parseFloat(parts[j]); if (!isNaN(val)) intensity.push(val); } } if (intensity.length === 0) throw new Error("GSAS Parse Error: No intensity data could be parsed."); const tth = Array.from({ length: intensity.length }, (_, i) => startTth + i * stepSize); return { tth, intensity, wavelength }; };
        
        const parseGsasXraFile = (text) => {
    const lines = text.trim().split(/\r?\n/);
    let wavelength = null;
    let startTth, stepSize;
    let dataStartIndex = -1;
    lines.forEach((line, index) => {
        const upperLine = line.toUpperCase();
        if (upperLine.includes('WAVELENGTH')) {
            const match = line.match(/wavelength\s+([0-9.]+)/i);
            if (match && match[1]) wavelength = parseFloat(match[1]);
        }
        if (upperLine.startsWith('BANK')) {
            const parts = line.trim().split(/\s+/);
            if (parts.length >= 7 && parts[4].toUpperCase() === 'CONST') {
                startTth = parseFloat(parts[5]) / 100.0;
                stepSize = parseFloat(parts[6]) / 100.0;
                dataStartIndex = index + 1;
            }
        }
    });

    // FIX: NaN !== undefined, so a malformed BANK line slipped through here too.
    if (!isFinite(startTth) || !isFinite(stepSize) || stepSize === 0) throw new Error("GSAS XRA Parse Error: Could not find a valid 'BANK' line with CONST scan parameters.");
    if (dataStartIndex === -1 || dataStartIndex >= lines.length) throw new Error("GSAS XRA Parse Error: Found scan parameters but no subsequent data lines.");

    const intensity = [];
    for (let i = dataStartIndex; i < lines.length; i++) {
        if (lines[i].trim() === '') continue;
        const parts = lines[i].trim().split(/\s+/);
        for (let j = 0; j < parts.length; j++) {
            const val = parseFloat(parts[j]);
            if (!isNaN(val)) intensity.push(val);
        }
    }
    if (intensity.length === 0) throw new Error("GSAS XRA Parse Error: No intensity data could be parsed.");
    const tth = Array.from({ length: intensity.length }, (_, i) => startTth + i * stepSize);
    return { tth, intensity, wavelength };
};
        
        const parseUxdFile = (text) => { const lines = text.trim().split(/\r?\n/); const intensity = []; let startTth, stepSize, wavelength; let inDataSection = false; for (const line of lines) { const trimmedLine = line.trim(); if (inDataSection) { const parts = trimmedLine.split(/\s+/); parts.forEach(part => { const val = parseFloat(part); if (!isNaN(val)) intensity.push(val); }); } else { if (trimmedLine.toUpperCase().startsWith('_START=')) startTth = parseFloat(trimmedLine.substring(7)); else if (trimmedLine.toUpperCase().startsWith('_STEPSIZE=')) stepSize = parseFloat(trimmedLine.substring(10)); else if (trimmedLine.toUpperCase().startsWith('_WL1=')) wavelength = parseFloat(trimmedLine.substring(5)); else if (trimmedLine.toUpperCase() === '_COUNTS') inDataSection = true; } } if (!isFinite(startTth) || !isFinite(stepSize) || stepSize === 0) throw new Error("Could not find valid _START and _STEPSIZE in UXD file."); if (intensity.length === 0) throw new Error("No intensity data found after _COUNTS in UXD file."); const tth = Array.from({ length: intensity.length }, (_, i) => startTth + i * stepSize); return { tth, intensity, wavelength }; };
        const parsePhilipsUdfFile = (text) => { const lines = text.trim().split(/\r?\n/); const tth = [], intensity = []; let inDataSection = false; let wavelength = null; for (const line of lines) { const trimmedLine = line.trim(); if (trimmedLine.toUpperCase().startsWith('LAMBDA')) { const parts = trimmedLine.split('='); if (parts.length > 1) wavelength = parseFloat(parts[1]); } if (trimmedLine.toUpperCase() === '[DATA]') { inDataSection = true; continue; } if (trimmedLine.startsWith('[') && trimmedLine.toUpperCase() !== '[DATA]') inDataSection = false; if (inDataSection) { const parts = trimmedLine.split(/,/).map(p => p.trim()); if(parts.length >= 2) { const x = parseFloat(parts[0]); const y = parseFloat(parts[1]); if (!isNaN(x) && !isNaN(y)) { tth.push(x); intensity.push(y); } } } } if (tth.length === 0) throw new Error("No [Data] section found in UDF file."); return { tth, intensity, wavelength }; };