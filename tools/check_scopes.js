#!/usr/bin/env node
/* ===========================================================================
 *  check_scopes.js  --  find names that exist in two places at once.
 *
 *  WHAT THIS IS FOR
 *
 *  Powder 5 loads fifteen classic <script> files into ONE shared global scope,
 *  and then powder5.html declares another ~1000 names inside a single
 *  DOMContentLoaded closure that spans eight thousand lines. Nothing in that
 *  arrangement objects when the same name is defined twice. The later
 *  definition silently wins, or the closure silently shadows the global, and
 *  the losing copy sits there looking maintained.
 *
 *  That is not hypothetical. In this codebase it has happened at least four
 *  times, and three of the four had already drifted into DIFFERENT BEHAVIOUR
 *  by the time anyone looked:
 *
 *    cfErfc                  polarization.js vs powder5.html -- one digit of
 *                            one coefficient differed, giving erfc(0) = 0.368
 *                            where it must be 1. Load order was the only thing
 *                            keeping the correct copy in play.
 *    runFit                  refinement_controller.js vs powder5.html -- 47
 *                            lines apart. A fix applied to the controller copy
 *                            did nothing at all, because the buttons call the
 *                            other one.
 *    generateMasterHklList   same pair, same cause.
 *    drawTheoreticalPreview  charting.js vs powder5.html -- differing guards.
 *
 *  A duplicate that is identical today is not safe, it is merely not yet
 *  wrong. This tool reports both, and says which.
 *
 *  USAGE
 *      node tools/check_scopes.js [rootDir]      # default: parent of tools/
 *      node tools/check_scopes.js --verbose      # list every name checked
 *
 *  Exit code 0 when clean, 1 when anything is reported, 2 on a usage error.
 *  Intended for a pre-commit hook or CI step.
 *
 *  WHAT IT DELIBERATELY DOES NOT DO
 *
 *  It is not a linter and it does not parse JavaScript properly -- there is no
 *  dependency to install, which is the point. It strips comments, strings,
 *  template literals and regex literals, then counts braces. That is enough to
 *  tell a top-level declaration from a local one in code formatted like this,
 *  and any construct it cannot read makes it report MORE, never less.
 * ======================================================================== */

'use strict';
const fs = require('fs');
const path = require('path');

const args = process.argv.slice(2);
const VERBOSE = args.includes('--verbose');
const ROOT = path.resolve(args.find(a => !a.startsWith('--')) ||
                          path.join(__dirname, '..'));

/* ---------------------------------------------------------------------------
 *  1. STRIPPING
 *
 *  Replace every comment, string, template literal and regex literal with
 *  same-length filler, so that offsets and line numbers survive but no brace,
 *  quote or keyword inside them can be mistaken for code.
 *
 *  The regex-literal case is the one that matters: a pattern like /\d{2}/
 *  carries braces that would otherwise unbalance the depth count for the rest
 *  of the file. A '/' begins a regex only where a VALUE cannot appear, which
 *  is decided by the previous significant character.
 * ------------------------------------------------------------------------ */
function strip(src) {
    const out = Array.from(src);
    const blank = (from, to, keepNewlines = true) => {
        for (let i = from; i < to && i < out.length; i++) {
            if (!(keepNewlines && out[i] === '\n')) out[i] = ' ';
        }
    };
    let i = 0, prev = '';           // prev = last significant char seen
    while (i < src.length) {
        const c = src[i], d = src[i + 1];

        if (c === '/' && d === '/') {
            let j = src.indexOf('\n', i); if (j < 0) j = src.length;
            blank(i, j); i = j; continue;
        }
        if (c === '/' && d === '*') {
            let j = src.indexOf('*/', i + 2); j = j < 0 ? src.length : j + 2;
            blank(i, j); i = j; continue;
        }
        if (c === '"' || c === "'" || c === '`') {
            let j = i + 1;
            while (j < src.length) {
                if (src[j] === '\\') { j += 2; continue; }
                if (src[j] === c) { j++; break; }
                j++;
            }
            // A template literal's ${...} is brace-balanced, so blanking the
            // whole span leaves the depth count intact.
            blank(i + 1, j - 1, c === '`');
            i = j; prev = 'x'; continue;
        }
        if (c === '/') {
            // Regex, or division? A regex may only start where a value may.
            const canBeRegex = prev === '' || '(,=:[!&|?{};+-*%~^<>'.includes(prev) ||
                               /[\s]/.test(prev);
            if (canBeRegex) {
                let j = i + 1, inClass = false, ok = false;
                while (j < src.length) {
                    const ch = src[j];
                    if (ch === '\\') { j += 2; continue; }
                    if (ch === '\n') break;                 // unterminated: not a regex
                    if (ch === '[') inClass = true;
                    else if (ch === ']') inClass = false;
                    else if (ch === '/' && !inClass) { ok = true; j++; break; }
                    j++;
                }
                if (ok) {
                    while (j < src.length && /[a-z]/.test(src[j])) j++;   // flags
                    blank(i, j, false);
                    i = j; prev = 'x'; continue;
                }
            }
        }
        if (!/\s/.test(c)) prev = c;
        i++;
    }
    return out.join('');
}

/* ---------------------------------------------------------------------------
 *  2. DECLARATIONS, WITH THEIR BRACE DEPTH
 *
 *  Depth 0 is the top level of the file (or of an inline <script>). A classic
 *  script's depth-0 declarations are the ones that land in the shared global
 *  scope and can therefore collide with another file's.
 * ------------------------------------------------------------------------ */
const DECL = /(?:^|[^.\w$])(?:(function)\s+([A-Za-z_$][\w$]*)|(const|let|var)\s+([A-Za-z_$][\w$]*)\s*[=;])/g;

function declarations(src, lineOffset = 0) {
    const clean = strip(src);
    // Precompute brace depth at every offset.
    const depth = new Int32Array(clean.length + 1);
    let d = 0;
    for (let i = 0; i < clean.length; i++) {
        depth[i] = d;
        if (clean[i] === '{') d++;
        else if (clean[i] === '}') d--;
    }
    depth[clean.length] = d;

    const lineAt = (off) => {
        let n = lineOffset + 1;
        for (let i = 0; i < off; i++) if (clean[i] === '\n') n++;
        return n;
    };
    // One pass to build a newline index instead of rescanning per match.
    const nl = [];
    for (let i = 0; i < clean.length; i++) if (clean[i] === '\n') nl.push(i);
    const fastLine = (off) => {
        let lo = 0, hi = nl.length;
        while (lo < hi) { const m = (lo + hi) >> 1; if (nl[m] < off) lo = m + 1; else hi = m; }
        return lineOffset + lo + 1;
    };

    const out = [];
    let m;
    DECL.lastIndex = 0;
    while ((m = DECL.exec(clean)) !== null) {
        const kind = m[1] || m[3];
        const name = m[2] || m[4];
        const at = m.index + m[0].indexOf(kind);
        out.push({ name, kind, depth: depth[at], line: fastLine(at), offset: at });
        DECL.lastIndex = m.index + m[0].length - 1;   // allow adjacent matches
    }
    return { decls: out, clean };
}

/* Body of a declaration, normalised for comparison: comments and blank lines
 * gone, whitespace collapsed. Two copies that differ only in formatting or in
 * their comments count as identical. */
function bodyOf(clean, offset) {
    const open = clean.indexOf('{', offset);
    if (open < 0) return null;
    let d = 0, i = open;
    for (; i < clean.length; i++) {
        if (clean[i] === '{') d++;
        else if (clean[i] === '}') { d--; if (d === 0) { i++; break; } }
    }
    return clean.slice(offset, i).replace(/\s+/g, ' ').trim();
}

/* ---------------------------------------------------------------------------
 *  3. THE FILE MAP
 *
 *  Load order matters: among classic scripts sharing the global scope, the
 *  LAST `function` or `var` declaration of a name wins. (Two top-level `const`
 *  or `let` of the same name in two classic scripts is a hard SyntaxError, so
 *  it cannot survive in working code -- but it is worth saying so if seen.)
 * ------------------------------------------------------------------------ */
function readHtml(file) {
    const src = fs.readFileSync(file, 'utf8');
    const scripts = [];
    const re = /<script\b([^>]*)>([\s\S]*?)<\/script>/gi;
    let m;
    while ((m = re.exec(src)) !== null) {
        const attrs = m[1];
        const srcAttr = /\bsrc\s*=\s*["']([^"']+)["']/i.exec(attrs);
        const line = src.slice(0, m.index).split('\n').length;
        if (srcAttr) scripts.push({ type: 'src', path: srcAttr[1], line });
        else scripts.push({ type: 'inline', body: m[2],
                            line: line + m[0].slice(0, m[0].indexOf('>')).split('\n').length - 1 });
    }
    return { src, scripts };
}

function workerImports(file) {
    if (!fs.existsSync(file)) return [];
    const clean = strip(fs.readFileSync(file, 'utf8'));
    const raw = fs.readFileSync(file, 'utf8');
    const names = [];
    const re = /importScripts\s*\(([^)]*)\)/g;
    let m;
    while ((m = re.exec(raw)) !== null) {
        for (const s of m[1].split(',')) {
            const t = s.trim().replace(/^['"`]|['"`]$/g, '');
            if (t && /\.js$/.test(t)) names.push(t);
        }
    }
    return names;
}

/* ------------------------------------------------------------------------ */
const findings = [];
const note = (severity, title, lines) => findings.push({ severity, title, lines });

const htmlPath = path.join(ROOT, 'powder5.html');
if (!fs.existsSync(htmlPath)) {
    console.error(`check_scopes: no powder5.html in ${ROOT}`);
    process.exit(2);
}
const html = readHtml(htmlPath);

/* --- the page's classic scripts, in load order ---------------------------
 *
 *  The src attribute is resolved against ROOT first, then against ROOT with
 *  the leading directory dropped. That second attempt is what lets the tool
 *  run against a FLATTENED copy of the tree -- a scratch directory where
 *  everything sits side by side -- as well as against the deployed layout
 *  with its js/ subdirectory. Both are real working arrangements here and a
 *  checker that only understood one of them would be run in neither.
 * ------------------------------------------------------------------------ */
function resolveScript(rel) {
    const candidates = [
        path.join(ROOT, rel),
        path.join(ROOT, path.basename(rel)),
        path.join(ROOT, 'js', path.basename(rel))
    ];
    return candidates.find(c => fs.existsSync(c)) || null;
}

const pageFiles = [];
for (const s of html.scripts) {
    if (s.type !== 'src') continue;
    if (/^https?:/i.test(s.path)) continue;
    // Third-party bundles are minified into a single line and are not ours to
    // police; skipping them silently keeps the report about this codebase.
    if (/(^|\/)lib\//.test(s.path)) continue;
    const found = resolveScript(s.path);
    if (found) pageFiles.push({ label: s.path, path: found });
    else note('warn', 'Script tag points at a file that is not present',
              [`  ${s.path}  (referenced from powder5.html line ${s.line})`,
               `  Looked in ${ROOT}, ${path.join(ROOT,'js')} and the root by basename.`]);
}

/* --- gather page-scope declarations -------------------------------------- */
const pageGlobals = new Map();     // name -> [{file, kind, line, body}]
for (const f of pageFiles) {
    const src = fs.readFileSync(f.path, 'utf8');
    const { decls, clean } = declarations(src);
    for (const dcl of decls) {
        if (dcl.depth !== 0) continue;
        if (!pageGlobals.has(dcl.name)) pageGlobals.set(dcl.name, []);
        pageGlobals.get(dcl.name).push({
            file: f.label, kind: dcl.kind, line: dcl.line,
            body: dcl.kind === 'function' ? bodyOf(clean, dcl.offset) : null
        });
    }
}

/* inline scripts: depth 0 is global; deeper is the closure */
const inlineGlobals = [];
const closureDecls = [];
for (const s of html.scripts) {
    if (s.type !== 'inline') continue;
    const { decls, clean } = declarations(s.body, s.line);
    for (const dcl of decls) {
        const rec = { file: 'powder5.html', kind: dcl.kind, line: dcl.line,
                      depth: dcl.depth,
                      body: dcl.kind === 'function' ? bodyOf(clean, dcl.offset) : null };
        if (dcl.depth === 0) inlineGlobals.push([dcl.name, rec]);
        else if (dcl.depth === 1) closureDecls.push([dcl.name, rec]);
    }
}
for (const [name, rec] of inlineGlobals) {
    if (!pageGlobals.has(name)) pageGlobals.set(name, []);
    pageGlobals.get(name).push(rec);
}

/* --- FINDING A: two definitions in the shared page scope ----------------- */
const dupLines = [];
for (const [name, defs] of [...pageGlobals].sort()) {
    if (defs.length < 2) continue;
    const kinds = new Set(defs.map(d => d.kind));
    const bodies = defs.map(d => d.body).filter(Boolean);
    const drifted = bodies.length >= 2 && new Set(bodies).size > 1;
    const winner = defs[defs.length - 1];
    dupLines.push(`  ${name}`);
    for (const d of defs) {
        const mark = (d === winner) ? '  <-- this one wins (loaded last)' : '';
        dupLines.push(`      ${d.kind.padEnd(8)} ${d.file}:${d.line}${mark}`);
    }
    if (drifted) {
        dupLines.push('      STATUS: the copies DIFFER. One of them is being ignored,');
        dupLines.push('              and whichever it is has already stopped matching.');
    } else if (bodies.length >= 2) {
        dupLines.push('      STATUS: identical today. That is not safe, only not yet wrong --');
        dupLines.push('              an edit to the losing copy will do nothing, silently.');
    }
    if (kinds.has('const') || kinds.has('let')) {
        dupLines.push('      NOTE: a repeated top-level const/let across two classic scripts');
        dupLines.push('            is a SyntaxError at load, not a silent override.');
    }
}
if (dupLines.length) note('error', 'Same name defined twice in the shared page scope', dupLines);

/* --- FINDING B: closure declarations shadowing a page global ------------- */
const shadowLines = [];
for (const [name, rec] of closureDecls.sort((a, b) => a[0].localeCompare(b[0]))) {
    const outer = pageGlobals.get(name);
    if (!outer || !outer.length) continue;
    if (outer.every(o => o.file === 'powder5.html')) continue;   // same file, not a shadow
    shadowLines.push(`  ${name}`);
    shadowLines.push(`      ${rec.kind.padEnd(8)} powder5.html:${rec.line}   (inside the closure -- THIS is what runs)`);
    for (const o of outer) {
        shadowLines.push(`      ${o.kind.padEnd(8)} ${o.file}:${o.line}   (shadowed; unreachable from the closure)`);
    }
    const bodies = [rec.body, ...outer.map(o => o.body)].filter(Boolean);
    if (bodies.length >= 2 && new Set(bodies).size > 1) {
        shadowLines.push('      STATUS: the copies DIFFER -- a fix applied to the outer one does nothing.');
    }
}
if (shadowLines.length) {
    note('error',
         'A closure declaration shadows a global from a module file',
         shadowLines);
}

/* --- FINDING C: duplicates inside each worker's own scope ---------------- */
const workerDirs = [ROOT, path.join(ROOT, 'js')].filter(d => fs.existsSync(d));
const workers = [];
for (const dir of workerDirs) {
    for (const f of fs.readdirSync(dir)) {
        if (/_worker\.js$/.test(f)) workers.push(path.join(dir, f));
    }
}
for (const w of workers) {
    const label = path.basename(w);
    const set = new Map();
    // A worker's importScripts paths are relative to the WORKER's directory.
    const wdir = path.dirname(w);
    const files = [...workerImports(w).map(n => {
        const direct = path.join(wdir, n);
        return fs.existsSync(direct) ? direct : (resolveScript(n) || direct);
    }), w];
    for (const f of files) {
        if (!fs.existsSync(f)) {
            note('warn', 'A worker imports a file that is not present',
                 [`  ${label} -> importScripts('${path.basename(f)}')`]);
            continue;
        }
        const { decls, clean } = declarations(fs.readFileSync(f, 'utf8'));
        for (const dcl of decls) {
            if (dcl.depth !== 0) continue;
            if (!set.has(dcl.name)) set.set(dcl.name, []);
            set.get(dcl.name).push({ file: path.basename(f), kind: dcl.kind, line: dcl.line,
                                     body: dcl.kind === 'function' ? bodyOf(clean, dcl.offset) : null });
        }
    }
    const wl = [];
    for (const [name, defs] of [...set].sort()) {
        if (defs.length < 2) continue;
        wl.push(`  ${name}`);
        for (const d of defs) wl.push(`      ${d.kind.padEnd(8)} ${d.file}:${d.line}`);
        const bodies = defs.map(d => d.body).filter(Boolean);
        if (bodies.length >= 2 && new Set(bodies).size > 1) {
            wl.push('      STATUS: the copies DIFFER.');
        }
    }
    if (wl.length) note('error', `Same name defined twice inside ${label}'s scope`, wl);
}

/* --- report -------------------------------------------------------------- */
const errors = findings.filter(f => f.severity === 'error');
const warns  = findings.filter(f => f.severity === 'warn');

console.log('check_scopes  ' + ROOT);
console.log('  page scripts : ' + pageFiles.length +
            '   inline blocks : ' + html.scripts.filter(s => s.type === 'inline').length +
            '   workers : ' + workers.length);
console.log('  page-scope names : ' + pageGlobals.size +
            '   closure names : ' + closureDecls.length);
console.log('');

for (const f of [...errors, ...warns]) {
    console.log((f.severity === 'error' ? 'FAIL  ' : 'warn  ') + f.title);
    for (const l of f.lines) console.log(l);
    console.log('');
}

if (VERBOSE) {
    console.log('--- every page-scope name ---');
    for (const [n, d] of [...pageGlobals].sort()) {
        console.log('  ' + n.padEnd(38) + d.map(x => x.file + ':' + x.line).join(', '));
    }
    console.log('');
}

if (!errors.length) {
    console.log('OK  no duplicate or shadowed names.');
    if (warns.length) console.log('    (' + warns.length + ' warning(s) above)');
}
process.exit(errors.length ? 1 : 0);
