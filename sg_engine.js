// sg_engine.js
// ------------------------------------------------------------------
// Single space-group engine, driven by the cctbx JSON database
// (sg_database.js). Loaded by both the main thread (via <script>)
// and the refinement worker (via importScripts).
//
// Public API exposed as SG_ENGINE on `window` / `self`:
//
//   SG_ENGINE.resolve(input)                -> normalised SG record or null
//   SG_ENGINE.list()                        -> array of all standard symbols (autocomplete)
//   SG_ENGINE.isReflectionAllowed(h,k,l,sg) -> boolean
//   SG_ENGINE.getMultiplicity(h,k,l,laue)   -> {multiplicity, canonical_hkl_obj}
//
// A resolved SG record has the shape:
//   {
//     number, symbol, setting_description, hall,
//     system,           // 'cubic' | 'tetragonal' | 'orthorhombic' | ...
//     laue_class,       // 'm-3m' | 'm-3' | '6/mmm' | ... (standard 11 Laue classes)
//     point_group,
//     centering,        // 'P' | 'I' | 'F' | 'C' | 'A' | 'B' | 'R'
//     centrosymmetric,
//     reflection_conditions  // copy from the cctbx setting
//   }
// ------------------------------------------------------------------

(function (root) {
    'use strict';

    if (!root.SG_DATABASE) {
        throw new Error('SG_ENGINE: SG_DATABASE must be loaded before sg_engine.js');
    }
    const DB = root.SG_DATABASE.space_groups;

    // ------------------------------------------------------------------
    // 1. Point-group -> Laue-class map (the 11 centrosymmetric Laue groups).
    //    Source: International Tables Vol. A, Table 3.2.1.4.
    // ------------------------------------------------------------------
    const POINT_GROUP_TO_LAUE = {
        '1':    '-1',   '-1':   '-1',
        '2':    '2/m',  'm':    '2/m',  '2/m':  '2/m',
        '222':  'mmm',  'mm2':  'mmm',  'mmm':  'mmm',
        '4':    '4/m',  '-4':   '4/m',  '4/m':  '4/m',
        '422':  '4/mmm','4mm':  '4/mmm','-42m': '4/mmm','4/mmm':'4/mmm',
        '3':    '-3',   '-3':   '-3',
        '32':   '-3m',  '3m':   '-3m',  '-3m':  '-3m',
        '6':    '6/m',  '-6':   '6/m',  '6/m':  '6/m',
        '622':  '6/mmm','6mm':  '6/mmm','-62m': '6/mmm','6/mmm':'6/mmm',
        '23':   'm-3',  'm-3':  'm-3',
        '432':  'm-3m', '-43m': 'm-3m', 'm-3m': 'm-3m'
    };

    // The 'system' in the main thread's old code used 'rhombohedral' as a
    // distinct bucket; the cctbx JSON keeps them all under 'trigonal'. The
    // downstream HKL generator branches on system name, so preserve the
    // distinction for R-centred groups.
    function deriveSystem(record) {
        const cs = record.crystal_system;
        if (cs !== 'trigonal') return cs;
        const sym = record.standard_symbol || '';
        return sym.trim().charAt(0).toUpperCase() === 'R' ? 'rhombohedral' : 'trigonal';
    }

    // ------------------------------------------------------------------
    // 2. Symbol normalisation.
    //    The user may type "P 21 21 21", "P212121", "p 2_1 2_1 2_1",
    //    "Fd-3m", "Fd3m", "P63/mmc", "R -3 c :H", etc. Strip everything
    //    that doesn't carry meaning and uppercase the first char of each
    //    "word" for matching purposes.
    // ------------------------------------------------------------------
    function normSymbol(s) {
        if (typeof s != 'string') return '';
        return s
            .replace(/_/g, '')             // 2_1 -> 21
            .replace(/\s+/g, '')           // drop spaces
            .replace(/:[a-zA-Z0-9]+$/, '') // drop :H / :R / :1 / :2 origin tags
            .toUpperCase();
    }

    // "Short" Hermann-Mauguin symbol: drop positions that are just "1".
    //   "P 1 21/c 1"   -> "P 21/c"   -> normalised "P21/C"
    //   "P 1 2 1"      -> "P 2"      -> "P2"
    //   "Fd-3m"        -> "Fd-3m"    (no spaces; returned unchanged)
    // Input may or may not have spaces; we work from whatever we've got.
    function shortSymbol(raw) {
        if (typeof raw !== 'string') return '';
        const s = raw.replace(/_/g, '').replace(/:[a-zA-Z0-9]+$/, '');
        // If the symbol has interior spaces, treat them as position delimiters
        // and drop any position that is exactly "1".
        if (/\s/.test(s.trim())) {
            const parts = s.trim().split(/\s+/);
            const kept = [parts[0]];  // centring letter
            for (let i = 1; i < parts.length; i++) {
                if (parts[i] !== '1') kept.push(parts[i]);
            }
            return kept.join('').toUpperCase();
        }
        // No spaces: nothing to shorten safely.
        return s.toUpperCase();
    }

    // The cctbx standard_symbol has embedded spaces (e.g. "P 21 21 21");
    // keep a pre-computed index for fast lookups.
    const BY_NUMBER = {};          // number -> primary setting record
    const BY_SYMBOL = new Map();   // normalised symbol -> setting record
    const BY_HALL   = new Map();   // normalised hall -> setting record
    const ALL_STANDARD_SYMBOLS = [];

    function pickPrimarySetting(sgEntry) {
        // Prefer a setting marked 'standard'. For monoclinic, none are; in
        // that case prefer the 'b' (unique axis b) setting, which is what
        // the HKL generator assumes. Otherwise fall back to the first.
        const settings = sgEntry.settings || [];
        let s = settings.find(x => x.description === 'standard');
        if (s) return s;
        s = settings.find(x => x.description === 'b');
        if (s) return s;
        return settings[0];
    }

    function buildRecord(sgEntry, setting) {
        const laue = POINT_GROUP_TO_LAUE[sgEntry.point_group] || null;
        if (!laue) {
            // Should not happen for the 230 standard groups, but be defensive.
            console.warn('SG_ENGINE: unknown point group', sgEntry.point_group,
                         'for SG', sgEntry.number);
        }
        return {
            number: sgEntry.number,
            symbol: setting.symbol,
            standard_symbol: sgEntry.standard_symbol,
            setting_description: setting.description,
            hall: setting.hall,
            system: deriveSystem(sgEntry),
            laue_class: laue,
            point_group: sgEntry.point_group,
            centering: (setting.symbol || '').charAt(0).toUpperCase(),
            centrosymmetric: !!sgEntry.centrosymmetric,
            reflection_conditions: setting.reflection_conditions || {}
        };
    }

    (function buildIndices() {
        Object.keys(DB).forEach(key => {
            const sgEntry = DB[key];
            const primary = pickPrimarySetting(sgEntry);
            if (!primary) return;
            const primaryRec = buildRecord(sgEntry, primary);
            BY_NUMBER[sgEntry.number] = primaryRec;
            ALL_STANDARD_SYMBOLS.push(sgEntry.standard_symbol);

            // Index every setting of every group by its symbol, Hall name,
            // AND the short-form symbol ("P 1 21/c 1" -> "P21/C").
            (sgEntry.settings || []).forEach(setting => {
                const rec = buildRecord(sgEntry, setting);
                const sKey = normSymbol(setting.symbol);
                if (sKey && !BY_SYMBOL.has(sKey)) BY_SYMBOL.set(sKey, rec);
                const shortKey = normSymbol(shortSymbol(setting.symbol));
                if (shortKey && !BY_SYMBOL.has(shortKey)) BY_SYMBOL.set(shortKey, rec);
                const hKey = normSymbol(setting.hall);
                if (hKey && !BY_HALL.has(hKey)) BY_HALL.set(hKey, rec);
            });
            // Also map the canonical "standard_symbol" (and its short form)
            // to the primary record. Short form works because the
            // standard_symbol is stored with spaces, e.g. "P 1 21/c 1".
            const stdKey = normSymbol(sgEntry.standard_symbol);
            if (stdKey && !BY_SYMBOL.has(stdKey)) BY_SYMBOL.set(stdKey, primaryRec);
            const stdShort = normSymbol(shortSymbol(sgEntry.standard_symbol));
            if (stdShort && !BY_SYMBOL.has(stdShort)) BY_SYMBOL.set(stdShort, primaryRec);
        });
    })();

    function resolve(input) {
        if (input == null) return null;
        // Number, or numeric string
        if (typeof input === 'number' && Number.isFinite(input)) {
            return BY_NUMBER[input] || null;
        }
        if (typeof input === 'string') {
            const trimmed = input.trim();
            if (/^\d+$/.test(trimmed)) {
                const n = parseInt(trimmed, 10);
                return BY_NUMBER[n] || null;
            }
            // First try the user's input verbatim (normalised).
            const key = normSymbol(trimmed);
            const direct = BY_SYMBOL.get(key) || BY_HALL.get(key);
            if (direct) return direct;
            // Then try the short form (drops stand-alone "1" positions).
            const shortKey = normSymbol(shortSymbol(trimmed));
            return BY_SYMBOL.get(shortKey) || null;
        }
        return null;
    }

    // ------------------------------------------------------------------
    // 3. Reflection-condition evaluator.
    //    The cctbx JSON uses a tiny grammar:
    //      "<expr>=<d>n"   (e.g. "h+k=2n", "-h+k+l=3n", "2h+l=4n")
    //    plus the family key ("hkl", "h00", "0k0", "00l", "h0l", "0kl",
    //    "hk0", "hhl"). A reflection is ALLOWED iff, for every family
    //    it belongs to, ALL conditions listed for that family hold.
    //
    //    Conditions in the JSON are cumulative, not alternative: e.g.
    //    F-centring yields three rules on 'hkl' ('h+k=2n', 'h+l=2n',
    //    'k+l=2n'), all three of which must be satisfied.
    // ------------------------------------------------------------------
    function inFamily(h, k, l, family) {
        switch (family) {
            case 'hkl': return true;
            case 'h00': return k === 0 && l === 0;
            case '0k0': return h === 0 && l === 0;
            case '00l': return h === 0 && k === 0;
            case 'h0l': return k === 0;
            case '0kl': return h === 0;
            case 'hk0': return l === 0;
            case 'hhl': return h === k;
            default:    return false;
        }
    }

    // Parses once, returns a closure over (h,k,l). Cached per SG record.
    const CONDITION_RE = /^\s*(-?\d*[a-z](?:\s*[+\-]\s*\d*[a-z])*)\s*=\s*(\d+)\s*n\s*$/i;

    function parseExpression(expr) {
        // Return a function (h,k,l) -> integer
        // Accept tokens like: h, k, l, 2h, -h, +h, 2l
        const clean = expr.replace(/\s+/g, '');
        // Split into signed terms: split on '+' or '-' keeping the sign
        const tokens = clean.match(/[+\-]?\d*[hkl]/g);
        if (!tokens) return null;
        const parts = tokens.map(tok => {
            const m = tok.match(/^([+\-]?)(\d*)([hkl])$/);
            if (!m) return null;
            const sign = m[1] === '-' ? -1 : 1;
            const coeff = m[2] === '' ? 1 : parseInt(m[2], 10);
            return { sign, coeff, var: m[3] };
        });
        if (parts.some(p => p === null)) return null;
        return function (h, k, l) {
            let total = 0;
            for (const p of parts) {
                const v = p.var === 'h' ? h : p.var === 'k' ? k : l;
                total += p.sign * p.coeff * v;
            }
            return total;
        };
    }

    function parseCondition(str) {
        const m = str.match(CONDITION_RE);
        if (!m) {
            console.warn('SG_ENGINE: unparsable reflection condition:', str);
            return null;
        }
        const fn = parseExpression(m[1]);
        const divisor = parseInt(m[2], 10);
        if (!fn || !divisor) return null;
        return { fn, divisor, raw: str };
    }

    // Compile reflection_conditions into parsed form on first use; cache
    // on the record itself to avoid reparsing on every HKL.
    function compileConditions(sg) {
        if (sg.__compiled_conditions) return sg.__compiled_conditions;
        const out = {};
        const src = sg.reflection_conditions || {};
        for (const fam in src) {
            const parsed = [];
            for (const rule of src[fam]) {
                const p = parseCondition(rule);
                if (p) parsed.push(p);
            }
            out[fam] = parsed;
        }
        Object.defineProperty(sg, '__compiled_conditions', {
            value: out, enumerable: false, writable: false
        });
        return out;
    }

    function isReflectionAllowed(h, k, l, sg) {
        if (!sg) return true;
        const conds = compileConditions(sg);
        // For every family the hkl is a member of, ALL conditions must hold.
        for (const fam in conds) {
            if (!inFamily(h, k, l, fam)) continue;
            const rules = conds[fam];
            for (const r of rules) {
                const val = r.fn(h, k, l);
                if (val % r.divisor !== 0) return false;
            }
        }
        return true;
    }

    // ------------------------------------------------------------------
    // 4. Multiplicity for the 11 Laue classes. Identical logic to what
    //    was in powder5.html, kept here so both main thread and worker
    //    share a single definition. Operates on unsigned indices and the
    //    standard conventions for each Laue group.
    // ------------------------------------------------------------------
    function getMultiplicity(h, k, l, laue_class) {
        if (h === 0 && k === 0 && l === 0) {
            return { multiplicity: 1, canonical_hkl_obj: [0, 0, 0] };
        }
        const abs_h = Math.abs(h);
        const abs_k = Math.abs(k);
        const abs_l = Math.abs(l);
        const sorted = [abs_h, abs_k, abs_l].sort((a, b) => b - a);
        const h_p = sorted[0], k_p = sorted[1], l_p = sorted[2];
        let m = 0;

        switch (laue_class) {
            case 'm-3m':
                // Special forms first (they must take precedence over the
                // general h>k>=l branches, which would otherwise swallow
                // them because 0 >= 0 is true).
                if (h_p > 0 && k_p === 0 && l_p === 0) m = 6;            // {h00}
                else if (h_p === k_p && l_p === 0 && h_p > 0) m = 12;    // {hh0}
                else if (h_p === k_p && k_p === l_p && l_p > 0) m = 8;   // {hhh}
                else if (h_p > k_p && k_p > 0 && l_p === 0) m = 24;      // {hk0}
                else if (h_p > k_p && k_p > l_p && l_p > 0) m = 48;      // {hkl}, h>k>l>0
                else if (h_p === k_p && k_p > l_p && l_p > 0) m = 24;    // {hhl}, h=k>l>0
                else if (h_p > k_p && k_p === l_p && l_p > 0) m = 24;    // {hkk}, h>k=l>0
                else m = 1;
                break;
            case 'm-3':
                if (h_p > 0 && k_p === 0 && l_p === 0) m = 6;            // {h00}
                else if (h_p === k_p && l_p === 0 && h_p > 0) m = 12;    // {hh0}
                else if (h_p === k_p && k_p === l_p && l_p > 0) m = 8;   // {hhh}
                else if (h_p > k_p && k_p > 0 && l_p === 0) m = 12;      // {hk0} - note 12, not 24
                else if (h_p > k_p && k_p > l_p && l_p > 0) m = 24;      // {hkl}
                else if ((h_p === k_p && k_p > l_p && l_p > 0) ||
                         (h_p > k_p && k_p === l_p && l_p > 0)) m = 12;
                else m = 1;
                break;
            case '6/mmm':
                if (l_p > 0) {
                    if (abs_h === 0 && abs_k === 0) m = 2;
                    else if (abs_h > 0 && abs_k === 0) m = 12;
                    else if (abs_h === abs_k && abs_k > 0) m = 12;
                    else if (abs_h > abs_k && abs_k >= 0) m = 24;
                    else m = 24;
                } else {
                    if (abs_h === 0 && abs_k === 0) m = 1;
                    else if (abs_h > 0 && abs_k === 0) m = 6;
                    else if (abs_h === abs_k && abs_k > 0) m = 6;
                    else if (abs_h > abs_k && abs_k >= 0) m = 12;
                    else m = 12;
                }
                break;
            case '6/m':
                if (l_p > 0) m = (abs_h > 0 || abs_k > 0) ? 12 : 2;
                else m = (abs_h > 0 || abs_k > 0) ? 6 : 1;
                break;
            case '-3m':
                if (l_p !== 0) {
                    if (abs_h === 0 && abs_k === 0) m = 2;
                    else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) m = 12;
                    else m = 24;
                } else {
                    if (abs_h === 0 && abs_k === 0) m = 1;
                    else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) m = 6;
                    else m = 12;
                }
                break;
            case '-3':
                if (abs_h === 0 && abs_k === 0 && l_p === 0) m = 1;
                else if (abs_h === 0 && abs_k === 0) m = 2;
                else m = 6;
                break;
            case '4/mmm':
                if (l_p > 0) {
                    if (abs_h === 0 && abs_k === 0) m = 2;
                    else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) m = 8;
                    else m = 16;
                } else {
                    if (abs_h === 0 && abs_k === 0) m = 1;
                    else if (abs_h === 0 || abs_k === 0 || abs_h === abs_k) m = 4;
                    else m = 8;
                }
                break;
            case '4/m':
                if (l_p > 0) m = (abs_h > 0 || abs_k > 0) ? 8 : 2;
                else m = (abs_h > 0 || abs_k > 0) ? 4 : 1;
                break;
            case 'mmm':
                if (abs_h > 0 && abs_k > 0 && l_p > 0) m = 8;
                else if ((abs_h > 0 && abs_k > 0 && l_p === 0) ||
                         (abs_h > 0 && abs_k === 0 && l_p > 0) ||
                         (abs_h === 0 && abs_k > 0 && l_p > 0)) m = 4;
                else if (abs_h > 0 || abs_k > 0 || l_p > 0) m = 2;
                else m = 1;
                break;
            case '2/m':
                // Unique axis b: the k axis is special.
                if (abs_k > 0) m = 4;
                else if (abs_k === 0 && (abs_h !== 0 || l_p !== 0)) m = 2;
                else m = 1;
                break;
            case '-1':
                if (abs_h === 0 && abs_k === 0 && l_p === 0) m = 1;
                else m = 2;
                break;
            default:
                console.warn('SG_ENGINE: unknown Laue class', laue_class, '-> m=1');
                m = 1;
        }
        return { multiplicity: m, canonical_hkl_obj: [h, k, l] };
    }

    // ------------------------------------------------------------------
    // 5. Helpers exposed for the UI: list of all 230 standard symbols
    //    for <datalist> autocomplete, and an "is supported" probe.
    // ------------------------------------------------------------------
    function list() {
        return ALL_STANDARD_SYMBOLS.slice();
    }

    // Return every setting of a given space group. Merges pairs that are
    // physically equivalent for powder diffraction (same HM symbol, same
    // description, identical reflection conditions — the typical "origin
    // choice 1 vs 2" redundancy found in SG 59, 68, 70, 85-88, 125, 129,
    // 130, 133, 141, 142, 201, 203, 222, 224, 227, 228, ...). Genuinely
    // different settings (axis permutations, hex vs. rhombohedral axes on
    // R-centred groups) are kept separate and get a clarifying suffix
    // added to the `display_label`.
    function listSettings(input) {
        const head = resolve(input);
        if (!head) return [];
        const sgEntry = DB[String(head.number)];
        if (!sgEntry) return [];
        const rawSettings = sgEntry.settings || [];

        // Group by (symbol, description) -- these are the entries the user
        // currently sees as duplicates.
        const buckets = new Map();
        rawSettings.forEach(s => {
            const key = s.symbol + '|' + (s.description || '');
            if (!buckets.has(key)) buckets.set(key, []);
            buckets.get(key).push(s);
        });

        const out = [];
        const condSig = (s) => JSON.stringify(s.reflection_conditions || {}, Object.keys(s.reflection_conditions || {}).sort());
        const isRhomb = (s) => (s.hall || '').indexOf('3*') !== -1;

        for (const [, group] of buckets) {
            // Partition by whether the Hall symbol marks rhombohedral axes.
            const hex = group.filter(s => !isRhomb(s));
            const rho = group.filter(s =>  isRhomb(s));

            // Within each axis choice, dedupe by reflection-condition
            // signature: identical conditions => same powder pattern =>
            // only emit one record (arbitrarily the first). Keep all
            // Hall symbols so advanced callers can still pick one.
            const emitFromPartition = (partition, axisLabel) => {
                if (partition.length === 0) return;
                const seen = new Map();
                partition.forEach(s => {
                    const sig = condSig(s);
                    if (!seen.has(sig)) seen.set(sig, []);
                    seen.get(sig).push(s);
                });
                seen.forEach(dupes => {
                    const rec = buildRecord(sgEntry, dupes[0]);
                    rec.display_label = rec.symbol
                        + (rec.setting_description && rec.setting_description !== 'standard'
                            ? `  (${rec.setting_description})` : '')
                        + (axisLabel ? `  [${axisLabel}]` : '');
                    rec.alternate_halls = dupes.map(d => d.hall);
                    out.push(rec);
                });
            };
            // Only tag an axis choice when both exist (otherwise it'd be noise).
            if (hex.length && rho.length) {
                emitFromPartition(hex, 'hex. axes');
                emitFromPartition(rho, 'rhomb. axes');
            } else {
                emitFromPartition(hex.length ? hex : rho, null);
            }
        }
        return out;
    }

    root.SG_ENGINE = {
        resolve:              resolve,
        isReflectionAllowed:  isReflectionAllowed,
        getMultiplicity:      getMultiplicity,
        list:                 list,
        listSettings:         listSettings,
        // Exposed for tests only.
        _parseCondition:      parseCondition,
        _inFamily:            inFamily,
        _POINT_GROUP_TO_LAUE: POINT_GROUP_TO_LAUE
    };
})(typeof self !== 'undefined' ? self : globalThis);
