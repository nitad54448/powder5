# cctbx_generate_sg_harker_v7.py
# 1 Aug 2026 - adds the data needed to build a STRUCTURE, not just index peaks.
#
# NOTE ON LINEAGE: this file is v5 plus the additions listed below. It does NOT
# contain whatever changes went into your own v6 -- if that version altered the
# absence analysis, the Harker geometry or the output schema, merge those in
# here before running. The additions below are confined to three new helper
# functions and the block marked "NEW in v7" inside generate_all_space_groups();
# nothing else in v5 was touched, so a merge should be mechanical.
#
# Changes from v5:
#
# 1. NEW get_symmetry_operations(). The old get_rotations() kept only the
#    rotation part of each operator AND deduplicated on it, which loses
#    everything a structure solution needs:
#       - the translation parts (screw axes, glide planes) were thrown away,
#         so Pnma and Pmmm produced an identical list;
#       - deduplication on the rotation collapsed the centring operators, so
#         an I- or F-centred group lost half or three quarters of its
#         operators outright.
#    The new function emits all sg.order_z() operators with their translations,
#    as EXACT rationals (numerator/denominator) so no 1/3 ever arrives in the
#    browser as 0.333333.
#    "rotations" is still emitted, unchanged, so the existing JS keeps working.
#
# 2. NEW centring translations (n_ltr, ltr) and the conventional centring
#    letter, so a JS consumer can expand a primitive operator list.
#
# 3. NEW Wyckoff table per setting: letter, multiplicity, site symmetry and
#    the special-position operator. This is what turns "a peak at 0.12, 0.25,
#    0.31" into "Wyckoff 4c, site symmetry .m., x,1/4,z".
#
# 4. NEW laue_class and centering at the top level of each space group. The JS
#    (powder5.html) already reads currentSG.laue_class and currentSG.centering,
#    but v5 never wrote them -- sg_engine.js was deriving them.
#
# 5. NEW z_prime / order_z and is_chiral / is_centric flags.
#
# Everything from v5 is preserved byte for byte unless listed above.
#
# --- v5 header retained below ---
# 17 Jul 2026 - Rule-string ordering fix for the JS (brutus) parser
# 1. Diagonal glide rules now always emit the DOUBLED term first
#    ("2*l+h=4n" instead of "h+2*l=4n").
# 2. Removed the redundant axial rule alongside a diagonal 4n rule.
# - Fixed "Space group symbol not recognized" via the "Hall:" prefix.
# - Strict logic for d-glides (4n) and hhl zones.

import json
import os
from collections import OrderedDict, defaultdict
from cctbx import sgtbx
from fractions import Fraction
import sys


def get_rotations(sg_info):
    """Extracts all real-space rotational matrices for the space group.

    UNCHANGED from v5, and deliberately so: powder5.html may still read it.
    Do NOT use this for building a structure -- see get_symmetry_operations().
    """
    sg = sg_info.group()
    rots = []
    for i in range(sg.order_z()):
        r = tuple(sg(i).r().as_double())
        if r not in rots:
            rots.append(r)
    return rots


def _rational_triplet(num, den):
    """(numerators, denominator) -> list of 'p/q' strings in lowest terms."""
    out = []
    for n in num:
        f = Fraction(int(n), int(den))
        out.append(f"{f.numerator}/{f.denominator}" if f.denominator != 1 else str(f.numerator))
    return out


def get_symmetry_operations(sg_info):
    """Full operator list: rotation, translation, and the xyz triplet.

    Returned per operator:
        xyz    : 'x,y,z' style string, exactly as cctbx formats it
        r      : 9 integers, row-major, x' = r . x + t
        t_num  : 3 integers, numerators of the translation
        t_den  : single integer denominator (cctbx uses 12, or 24 for some
                 d-glides); t = t_num / t_den exactly
        t      : the same translation as floats, for convenience

    Every operator of the group is listed, including the centring ones, so a
    consumer can apply the list directly with no further expansion.
    """
    sg = sg_info.group()
    ops = []
    for i in range(sg.order_z()):
        op = sg(i)
        r = op.r()
        t = op.t()
        r_num = [int(v) for v in r.num()]
        r_den = int(r.den())
        if r_den != 1:
            # cctbx always uses r_den == 1 for space-group rotations; guard
            # anyway so a silent non-integer matrix can never slip through.
            r_num = [v // r_den for v in r_num]
        t_num = [int(v) for v in t.num()]
        t_den = int(t.den())
        ops.append(OrderedDict([
            ("xyz", str(op)),
            ("r", r_num),
            ("t_num", t_num),
            ("t_den", t_den),
            ("t", [v / t_den for v in t_num]),
            ("t_frac", _rational_triplet(t_num, t_den)),
        ]))
    return ops


def get_centring(sg_info):
    """Centring translations as exact rationals, plus the conventional letter.

    Derived from the OPERATOR LIST rather than from sg.ltr(): the centring
    translations are exactly the translation parts of the operators whose
    rotation is the identity. The indexed accessor sg.ltr(i) is not available
    in every cctbx build ("_.ltr() takes 1 positional argument but 2 were
    given"), and there is no need for it -- this derivation is exact and works
    everywhere.
    """
    sg = sg_info.group()
    ltr = []
    seen = set()
    identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)

    for i in range(sg.order_z()):
        op = sg(i)
        r = op.r()
        r_num = tuple(int(v) for v in r.num())
        r_den = int(r.den())
        if r_den != 1:
            r_num = tuple(v // r_den for v in r_num)
        if r_num != identity:
            continue
        t = op.t()
        num = [int(x) for x in t.num()]
        den = int(t.den())
        key = tuple(Fraction(x, den) % 1 for x in num)
        if key in seen:
            continue
        seen.add(key)
        ltr.append(OrderedDict([
            ("t_num", num),
            ("t_den", den),
            ("t", [x / den for x in num]),
            ("t_frac", _rational_triplet(num, den)),
        ]))

    letter = "?"
    for getter in (
        lambda: sg.conventional_centring_type_symbol(),
        lambda: sg_info.type().lookup_symbol().strip().split(":")[0].strip()[0],
    ):
        try:
            value = getter()
            if value:
                letter = str(value).strip()[0]
                break
        except Exception:
            continue

    return letter, ltr


# cctbx does not expose the Wyckoff table at one fixed place. Depending on the
# build it is a method on space_group_info, or a submodule of sgtbx that
# "from cctbx import sgtbx" does NOT import for you (which is what produced
# "module 'cctbx.sgtbx' has no attribute 'wyckoff'" -- nothing to pip install,
# just an import that was never done). Try each in turn and remember which one
# worked so the probing happens once.
_WYCKOFF_ACCESSOR = "unprobed"


def _get_wyckoff_table_object(sg_info):
    global _WYCKOFF_ACCESSOR

    def _via_info(info):
        return info.wyckoff_table()

    def _via_submodule(info):
        from cctbx.sgtbx import wyckoff as _w
        return _w.table(info)

    def _via_import(info):
        import importlib
        _w = importlib.import_module("cctbx.sgtbx.wyckoff")
        return _w.table(info)

    def _via_attr(info):
        return sgtbx.wyckoff.table(info)

    def _via_flat(info):
        return sgtbx.wyckoff_table(info)

    routes = [
        ("space_group_info.wyckoff_table()", _via_info),
        ("from cctbx.sgtbx import wyckoff", _via_submodule),
        ("importlib cctbx.sgtbx.wyckoff", _via_import),
        ("sgtbx.wyckoff.table()", _via_attr),
        ("sgtbx.wyckoff_table()", _via_flat),
    ]

    # Once a route is known to work, stop probing.
    if _WYCKOFF_ACCESSOR not in ("unprobed", "none"):
        for name, fn in routes:
            if name == _WYCKOFF_ACCESSOR:
                try:
                    return fn(sg_info)
                except Exception:
                    break
    if _WYCKOFF_ACCESSOR == "none":
        return None

    errors = []
    for name, fn in routes:
        try:
            table = fn(sg_info)
            if table is not None:
                if _WYCKOFF_ACCESSOR == "unprobed":
                    print(f"\n[i] Wyckoff tables found via {name}.")
                    sys.stdout.flush()
                _WYCKOFF_ACCESSOR = name
                return table
        except Exception as e:
            errors.append(f"{name}: {type(e).__name__}: {e}")

    _WYCKOFF_ACCESSOR = "none"
    print("\n[!] No Wyckoff accessor worked in this cctbx build. This is OPTIONAL:")
    print("    symops, centring and everything else are still written, and the")
    print("    structure builder derives site multiplicity from the operators")
    print("    itself. Only the ITA letters (4c, 8d, ...) will be missing.")
    for e in errors:
        print(f"      tried {e}")
    sys.stdout.flush()
    return None


def get_wyckoff_table(sg_info):
    """Wyckoff positions: letter, multiplicity, site symmetry, special op.

    The special_op xyz string is the constraint on a site at that position --
    'x,1/4,z' for a mirror site, '0,0,0' for a fixed one. Projecting a peak
    through it is how you decide whether a maximum is really on a special
    position or just close to one.
    """
    positions = []
    table = _get_wyckoff_table_object(sg_info)
    if table is None:
        return positions

    def _safe(fn, default):
        try:
            return fn()
        except Exception:
            return default

    for i in range(table.size()):
        try:
            pos = table.position(i)
        except Exception:
            continue
        # point_group_type() lives directly on the position in some builds and
        # behind site_symmetry_group() in others; try both before giving up.
        site_sym = _safe(lambda: str(pos.point_group_type()), None)
        if site_sym is None:
            site_sym = _safe(lambda: str(pos.site_symmetry_group().point_group_type()), "?")
        positions.append(OrderedDict([
            ("letter", _safe(lambda: pos.letter(), "?")),
            ("multiplicity", _safe(lambda: int(pos.multiplicity()), 0)),
            ("site_symmetry", site_sym),
            ("special_op", _safe(lambda: str(pos.special_op()), "")),
        ]))
    return positions


# ---------------------------------------------------------------------------
#  Everything below this line is unchanged from v5 except where marked.
# ---------------------------------------------------------------------------

def get_harker_geometry(sg_info):
    sg = sg_info.group()
    sections = []
    coords = ['u', 'v', 'w']
    vars_str = ['x', 'y', 'z']

    # Safely iterate through all symmetry operations
    for i in range(1, sg.order_z()):
        op = sg(i)
        r = op.r().as_double()
        t = op.t().as_double()

        # Construct R - I matrix
        r_minus_i = [
            [r[0]-1, r[1],   r[2]  ],
            [r[3],   r[4]-1, r[5]  ],
            [r[6],   r[7],   r[8]-1]
        ]

        # Identify zero rows (which define the Harker planes/lines)
        zero_rows = [idx for idx in range(3) if all(abs(val) < 1e-5 for val in r_minus_i[idx])]
        if not zero_rows:
            continue

        # Identify the submatrix to invert
        nz_rows = [idx for idx in range(3) if idx not in zero_rows]
        nz_cols = [idx for idx in range(3) if any(abs(r_minus_i[r_idx][idx]) > 1e-5 for r_idx in range(3))]

        solver = {"x": "?", "y": "?", "z": "?"}

        # Algebraically invert the submatrix to solve for x, y, z
        if len(nz_rows) == len(nz_cols) and len(nz_rows) > 0:
            n = len(nz_rows)
            A = [[r_minus_i[row][col] for col in nz_cols] for row in nz_rows]
            invA = None

            if n == 1:
                if abs(A[0][0]) > 1e-5:
                    invA = [[1.0 / A[0][0]]]
            elif n == 2:
                det = A[0][0]*A[1][1] - A[0][1]*A[1][0]
                if abs(det) > 1e-5:
                    invA = [[A[1][1]/det, -A[0][1]/det], [-A[1][0]/det, A[0][0]/det]]
            elif n == 3:
                det = (A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
                     - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
                     + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]))
                if abs(det) > 1e-5:
                    invA = [
                        [(A[1][1]*A[2][2] - A[1][2]*A[2][1])/det, (A[0][2]*A[2][1] - A[0][1]*A[2][2])/det, (A[0][1]*A[1][2] - A[0][2]*A[1][1])/det],
                        [(A[1][2]*A[2][0] - A[1][0]*A[2][2])/det, (A[0][0]*A[2][2] - A[0][2]*A[2][0])/det, (A[0][2]*A[1][0] - A[0][0]*A[1][2])/det],
                        [(A[1][0]*A[2][1] - A[1][1]*A[2][0])/det, (A[0][1]*A[2][0] - A[0][0]*A[2][1])/det, (A[0][0]*A[1][1] - A[0][1]*A[1][0])/det]
                    ]

            # Format the output strings for the JS evaluator
            if invA:
                for c_idx, col_name in enumerate(nz_cols):
                    terms = []
                    for r_idx, row_name in enumerate(nz_rows):
                        coeff = invA[c_idx][r_idx]
                        if abs(coeff) > 1e-5:
                            terms.append(f"{coeff:g}*({coords[row_name]}-{t[row_name]:g})")
                    if terms:
                        solver[vars_str[col_name]] = " + ".join(terms).replace("+ -", "- ")

        sec_type = "plane" if len(zero_rows) == 1 else "line"

        for z_idx in zero_rows:
            val = t[z_idx] % 1.0
            if not any(s['coordinate'] == coords[z_idx] and abs(s['value'] - val) < 1e-5 for s in sections):
                sections.append({
                    "type": sec_type,
                    "coordinate": coords[z_idx],
                    "value": val,
                    "solver": solver
                })

    return sections


def evaluate_rule(h, k, l, rule_str):
    """Evaluates if a reflection satisfies a textual rule safely."""
    try:
        if "=" not in rule_str:
            return False

        lhs, rhs = rule_str.split("=")

        # Map the right-hand side to its modulus integer
        mod_map = {"2n": 2, "3n": 3, "4n": 4, "6n": 6}
        if rhs not in mod_map:
            return False

        mod_val = mod_map[rhs]

        # Evaluate the mathematical expression on the left-hand side natively
        context = {'h': h, 'k': k, 'l': l}
        lhs_val = eval(lhs, {}, context)

        return (lhs_val % mod_val) == 0
    except Exception:
        return False


def get_condition_string(present_values, axis_name):
    """Deduces axial conditions (2n, 4n, 6n)."""
    if not present_values:
        return None

    vals = sorted(list(set(present_values)))
    if not vals:
        return None

    # Strongest conditions first
    if all(v % 6 == 0 for v in vals): return f"{axis_name}=6n"
    if all(v % 4 == 0 for v in vals): return f"{axis_name}=4n"
    if all(v % 3 == 0 for v in vals): return f"{axis_name}=3n"
    if all(v % 2 == 0 for v in vals): return f"{axis_name}=2n"

    return None


def check_zonal_from_points(present_tuples, idx1, idx2):
    """
    Deduces zonal conditions, returning a list of all coexisting independent rules.
    """
    if not present_tuples:
        return []

    names = {0: 'h', 1: 'k', 2: 'l'}
    n1, n2 = names[idx1], names[idx2]

    v1 = [p[idx1] for p in present_tuples]
    v2 = [p[idx2] for p in present_tuples]
    sums = [p[idx1] + p[idx2] for p in present_tuples]
    diffs = [p[idx1] - p[idx2] for p in present_tuples]

    rules = []

    # --- 1. PRIORITY: Diamond/Complex Glides (4n) ---
    diag_forced_axis = None
    if all((2*p[idx1] + p[idx2]) % 4 == 0 for p in present_tuples):
        rules.append(f"2*{n1}+{n2}=4n")
        diag_forced_axis = n2
    elif all((2*p[idx1] - p[idx2]) % 4 == 0 for p in present_tuples):
        rules.append(f"2*{n1}-{n2}=4n")
        diag_forced_axis = n2
    elif all((p[idx1] + 2*p[idx2]) % 4 == 0 for p in present_tuples):
        rules.append(f"2*{n2}+{n1}=4n")   # doubled term first
        diag_forced_axis = n1

    # --- 2. PRIORITY: Standard Glide Sums (4n) ---
    has_4n_sum = False
    if all(s % 4 == 0 for s in sums):
        rules.append(f"{n1}+{n2}=4n")
        has_4n_sum = True
    elif all(d % 4 == 0 for d in diffs):
        rules.append(f"{n1}-{n2}=4n")

    # --- 3. PRIORITY: Axial Conditions (2n) ---
    has_n1_2n = all(v % 2 == 0 for v in v1)
    has_n2_2n = all(v % 2 == 0 for v in v2)

    if has_n1_2n and diag_forced_axis != n1:
        rules.append(f"{n1}=2n")
    if has_n2_2n and diag_forced_axis != n2:
        rules.append(f"{n2}=2n")

    # --- 4. PRIORITY: Weak Centering Sums (3n, 2n) ---
    if all(s % 3 == 0 for s in sums):
        rules.append(f"{n1}+{n2}=3n")
    elif all(d % 3 == 0 for d in diffs):
        rules.append(f"{n1}-{n2}=3n")

    if not (has_4n_sum or (has_n1_2n and has_n2_2n)):
        if all(s % 2 == 0 for s in sums):
            rules.append(f"{n1}+{n2}=2n")
        elif all(d % 2 == 0 for d in diffs):
            rules.append(f"{n1}-{n2}=2n")

    return rules


def is_redundant(specific_rule, general_rules, axes_indices):
    """Checks if a specific rule is already covered by general HKL rules."""
    if not general_rules:
        return False

    valid_points_under_general = []
    is_diagonal = axes_indices.get('diagonal', False)
    active_axes = [key for key, v in axes_indices.items() if v is None and key != 'diagonal']

    for i in range(1, 25):
        if len(active_axes) == 1 and not is_diagonal:  # Axial check
            h, k, l = 0, 0, 0
            if active_axes[0] == 'h': h = i
            elif active_axes[0] == 'k': k = i
            elif active_axes[0] == 'l': l = i

            if all(evaluate_rule(h, k, l, r) for r in general_rules):
                valid_points_under_general.append((h, k, l))

        elif len(active_axes) >= 2 or is_diagonal:  # Zonal check
            for j in range(1, 25):
                h, k, l = 0, 0, 0
                if is_diagonal:  # hhl
                    h, k, l = i, i, j
                elif axes_indices['l'] == 0:  # hk0
                    h, k = i, j
                elif axes_indices['k'] == 0:  # h0l
                    h, l = i, j
                elif axes_indices['h'] == 0:  # 0kl
                    k, l = i, j

                if all(evaluate_rule(h, k, l, r) for r in general_rules):
                    valid_points_under_general.append((h, k, l))

    if not valid_points_under_general:
        return False

    if len(active_axes) == 1 and not is_diagonal:
        idx_map = {'h': 0, 'k': 1, 'l': 2}
        axis_char = active_axes[0]
        vals = [p[idx_map[axis_char]] for p in valid_points_under_general]
        return get_condition_string(vals, axis_char) == specific_rule

    else:
        idx_map = {'h': 0, 'k': 1, 'l': 2}
        idx1 = 0 if is_diagonal else idx_map[active_axes[0]]
        idx2 = 2 if is_diagonal else idx_map[active_axes[1]]
        return specific_rule in check_zonal_from_points(valid_points_under_general, idx1, idx2)


def analyze_systematic_absences(space_group_info):
    sg = space_group_info.group()

    # --- 1. HKL Conditions (Centering) ---
    present_hkl = []
    for h in range(1, 8):
        for k in range(1, 8):
            for l in range(1, 8):
                if not sg.is_sys_absent((h, k, l)):
                    present_hkl.append((h, k, l))

    hkl_rules = []
    if present_hkl:
        if all((h+k) % 2 == 0 for h, k, l in present_hkl): hkl_rules.append("h+k=2n")
        if all((h+l) % 2 == 0 for h, k, l in present_hkl): hkl_rules.append("h+l=2n")
        if all((k+l) % 2 == 0 for h, k, l in present_hkl): hkl_rules.append("k+l=2n")
        if all((h+k+l) % 2 == 0 for h, k, l in present_hkl): hkl_rules.append("h+k+l=2n")

        # Rhombohedral check (Obverse vs Reverse)
        if all((-h+k+l) % 3 == 0 for h, k, l in present_hkl): hkl_rules.append("-h+k+l=3n")
        if all((h-k+l) % 3 == 0 for h, k, l in present_hkl): hkl_rules.append("h-k+l=3n")

        # F-centering cleanup (only if strictly F and not I)
        if "h+k=2n" in hkl_rules and "h+l=2n" in hkl_rules and "k+l=2n" in hkl_rules:
            if "h+k+l=2n" not in hkl_rules:
                hkl_rules = ["h+k=2n", "h+l=2n", "k+l=2n"]

    conditions = OrderedDict()
    if hkl_rules:
        conditions['hkl'] = hkl_rules

    # --- 2. Zonal Conditions ---
    zones = [
        ('0kl', {'h': 0, 'k': None, 'l': None}, 1, 2),
        ('h0l', {'h': None, 'k': 0, 'l': None}, 0, 2),
        ('hk0', {'h': None, 'k': None, 'l': 0}, 0, 1),
        ('hhl', {'h': None, 'k': None, 'l': None, 'diagonal': True}, 0, 2)
    ]

    for zone_name, axes_map, idx1, idx2 in zones:
        points = []
        for i in range(1, 20):
            for j in range(1, 20):
                hkl = [0, 0, 0]
                if zone_name == 'hhl':
                    hkl[0] = i
                    hkl[1] = i
                    hkl[2] = j
                else:
                    hkl[idx1] = i
                    hkl[idx2] = j

                if not sg.is_sys_absent(tuple(hkl)):
                    points.append(tuple(hkl))

        rules = check_zonal_from_points(points, idx1, idx2)
        if rules:
            valid_rules = [r for r in rules if not is_redundant(r, hkl_rules, axes_map)]
            if valid_rules:
                conditions[zone_name] = valid_rules

    # --- 3. Axial Conditions ---
    axes_defs = [
        ('h00', {'h': None, 'k': 0, 'l': 0}, 0, 'h'),
        ('0k0', {'h': 0, 'k': None, 'l': 0}, 1, 'k'),
        ('00l', {'h': 0, 'k': 0, 'l': None}, 2, 'l')
    ]

    for axis_name, axes_map, idx, label in axes_defs:
        points = []
        for i in range(1, 30):
            hkl = [0, 0, 0]
            hkl[idx] = i
            if not sg.is_sys_absent(tuple(hkl)):
                points.append(hkl[idx])

        rule = get_condition_string(points, label)
        if rule:
            is_covered = False
            if is_redundant(rule, hkl_rules, axes_map):
                is_covered = True

            if not is_covered:
                parents = []
                if axis_name == 'h00': parents = ['hk0', 'h0l']
                elif axis_name == '0k0': parents = ['hk0', '0kl']
                elif axis_name == '00l': parents = ['h0l', '0kl', 'hhl']

                for p in parents:
                    if p in conditions:
                        if is_redundant(rule, conditions[p], axes_map):
                            is_covered = True
                            break

            if not is_covered:
                conditions[axis_name] = [rule]

    return conditions


_WARNED = set()


def _warn_once(what, err):
    """Report an unavailable cctbx accessor once, not 700 times."""
    key = f"{what}:{type(err).__name__}:{err}"
    if key in _WARNED:
        return
    _WARNED.add(key)
    print(f"\n[!] {what} unavailable in this cctbx build ({err}). "
          f"That field will be empty; everything else is still written.")
    sys.stdout.flush()


def generate_all_space_groups():
    all_data = defaultdict(lambda: {
        "number": 0, "standard_symbol": "", "crystal_system": "",
        "point_group": "", "laue_class": "", "centrosymmetric": False,
        "chiral": False, "settings": []
    })

    iterator = sgtbx.space_group_symbol_iterator()
    processed_settings = defaultdict(set)

    print("Processing settings...")
    count = 0

    while True:
        try:
            symbols = iterator.next()
        except StopIteration:
            break

        if symbols.number() == 0:
            break

        try:
            sg_num = symbols.number()
            hm_symbol = symbols.hermann_mauguin()
            hall = symbols.hall().strip()

            uid = hall
            if uid in processed_settings[sg_num]:
                continue
            processed_settings[sg_num].add(uid)

            sg_info = sgtbx.space_group_info(symbol=f"Hall: {hall}")
            sg = sg_info.group()

            if all_data[str(sg_num)]["number"] == 0:
                all_data[str(sg_num)]["number"] = sg_num
                all_data[str(sg_num)]["standard_symbol"] = sg_info.type().lookup_symbol()
                all_data[str(sg_num)]["crystal_system"] = str(sg.crystal_system()).lower()
                all_data[str(sg_num)]["point_group"] = str(sg.point_group_type())
                all_data[str(sg_num)]["centrosymmetric"] = sg.is_centric()
                # NEW in v7 -- powder5.html already reads currentSG.laue_class
                try:
                    all_data[str(sg_num)]["laue_class"] = str(sg.laue_group_type())
                except Exception as e:
                    _warn_once("laue_group_type", e)
                    all_data[str(sg_num)]["laue_class"] = ""
                try:
                    all_data[str(sg_num)]["chiral"] = bool(sg.is_chiral())
                except Exception:
                    all_data[str(sg_num)]["chiral"] = False

            conditions = analyze_systematic_absences(sg_info)
            harker_data = get_harker_geometry(sg_info)
            rotations = get_rotations(sg_info)

            # --- NEW in v7 ---
            # Each of these is isolated. The enclosing try/except skips the
            # WHOLE setting on any exception, so an API difference in one
            # accessor used to silently drop every space group in the file
            # (all 230 of them, leaving an empty JSON). A failure here now
            # costs only the field it belongs to.
            try:
                symops = get_symmetry_operations(sg_info)
            except Exception as e:
                _warn_once("symops", e)
                symops = []
            try:
                centring_letter, centring_ltr = get_centring(sg_info)
            except Exception as e:
                _warn_once("centring", e)
                centring_letter, centring_ltr = "?", []
            try:
                wyckoff = get_wyckoff_table(sg_info)
            except Exception as e:
                _warn_once("wyckoff", e)
                wyckoff = []
            try:
                order_z = int(sg.order_z())
            except Exception:
                order_z = len(symops)
            try:
                order_p = int(sg.order_p())
            except Exception:
                order_p = 0

            desc = symbols.qualifier() if symbols.qualifier() else "standard"
            clean_sym = hm_symbol.replace(" ", "")

            all_data[str(sg_num)]["settings"].append({
                "symbol": clean_sym,
                "description": desc,
                "hall": hall,
                "reflection_conditions": conditions,
                "harker_sections": harker_data,
                "rotations": rotations,
                # --- NEW in v7 ---
                "symops": symops,
                "order_z": order_z,
                "order_p": order_p,
                "centering": centring_letter,
                "centring_translations": centring_ltr,
                "wyckoff": wyckoff,
            })

            count += 1
            if count % 50 == 0:
                print(f"Processed {count} settings...", end="\r")
                sys.stdout.flush()

        except Exception as e:
            print(f"\n[!] Error processing SG {sg_num}: {e}")
            sys.stdout.flush()
            continue

    return all_data


# The filename the app fetches. It appears in exactly two other places:
#   powder5.html          (the DOMContentLoaded fetch)
#   refinement_worker.js  (the initPromise fetch)
# All three must match. Both fetches use a RELATIVE path, so the file has to
# sit next to powder5.html -- not next to this script.
OUTPUT_NAME = 'cctbx_space_groups_all_settings_v7.json'


def resolve_output_path():
    """Work out where to write, and say so in absolute terms.

    Running this from a scripts/ subfolder used to drop the JSON there, while
    the dev server serves the project root -- giving a 404 for a file that had
    just been written successfully. The app's own location decides now.

    Order: an explicit path on the command line, then whichever nearby folder
    actually contains powder5.html, then the current directory.
    """
    if len(sys.argv) > 1:
        target = os.path.abspath(sys.argv[1])
        if os.path.isdir(target):
            target = os.path.join(target, OUTPUT_NAME)
        print(f"Output path given on the command line: {target}")
        return target

    here = os.path.abspath(os.getcwd())
    script_dir = os.path.abspath(os.path.dirname(__file__))
    for folder in (here, script_dir, os.path.dirname(script_dir), os.path.dirname(here)):
        if folder and os.path.isfile(os.path.join(folder, 'powder5.html')):
            target = os.path.join(folder, OUTPUT_NAME)
            if os.path.abspath(folder) != here:
                print(f"\n[i] powder5.html found in {folder}")
                print("    Writing the JSON there, not into the current folder, because")
                print("    the app fetches it with a relative path from its own location.")
            return target

    print("\n[!] powder5.html was not found nearby, so the JSON goes in the current")
    print(f"    folder ({here}). If the browser reports a 404 for {OUTPUT_NAME},")
    print("    move it next to powder5.html, or re-run as:")
    print(f"        python {os.path.basename(__file__)} <path-to-folder-with-powder5.html>")
    return os.path.join(here, OUTPUT_NAME)


def main():
    print("=" * 60)
    print("Space Group Generator v7 - symops, centring and Wyckoff added")
    print("=" * 60)

    data = generate_all_space_groups()

    total_settings = sum(len(x['settings']) for x in data.values())
    print(f"\nFinal processing complete.")
    print(f"Total space group numbers: {len(data)}")
    print(f"Total settings: {total_settings}")

    sorted_data = OrderedDict()
    for k in sorted(data.keys(), key=int):
        sorted_data[k] = data[k]
        sorted_data[k]["settings"].sort(key=lambda x: x["symbol"])

    output_wrapper = OrderedDict()
    output_wrapper["schema_version"] = 7
    output_wrapper["space_groups"] = sorted_data

    outfile = resolve_output_path()
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(output_wrapper, f, indent=2, ensure_ascii=False)

    print(f"\nSaved to {outfile}")

    # --- self-check: the operator list must be a closed group -------------
    print("\nVerifying operator lists form closed groups...")
    bad = 0
    checked = 0
    missing_symops = 0
    missing_wyckoff = 0
    for num, entry in sorted_data.items():
        for setting in entry["settings"]:
            ops = setting["symops"]
            if not setting.get("wyckoff"):
                missing_wyckoff += 1
            if not ops:
                missing_symops += 1
                continue
            if len(ops) != setting["order_z"]:
                print(f"  [!] SG {num} {setting['symbol']}: {len(ops)} ops but order_z={setting['order_z']}")
                bad += 1
                continue
            keys = set()
            for op in ops:
                keys.add((tuple(op["r"]), tuple(Fraction(n, op["t_den"]) % 1 for n in op["t_num"])))
            if len(keys) != len(ops):
                print(f"  [!] SG {num} {setting['symbol']}: duplicate operators")
                bad += 1
                continue
            # closure: every product of two operators must already be present
            closed = True
            for a in ops:
                for b in ops:
                    ra, ta, da = a["r"], a["t_num"], a["t_den"]
                    rb, tb, db = b["r"], b["t_num"], b["t_den"]
                    rc = [sum(ra[3 * i + k] * rb[3 * k + j] for k in range(3)) for i in range(3) for j in range(3)]
                    tc = [Fraction(ta[i], da) + sum(Fraction(ra[3 * i + k] * tb[k], db) for k in range(3))
                          for i in range(3)]
                    key = (tuple(rc), tuple(t % 1 for t in tc))
                    if key not in keys:
                        closed = False
                        break
                if not closed:
                    break
            if not closed:
                print(f"  [!] SG {num} {setting['symbol']}: operator list is not closed")
                bad += 1
            checked += 1
    print(f"Checked {checked} settings, {bad} problem(s).")
    if missing_symops:
        print(f"[!] {missing_symops} setting(s) have NO symops -- structure building "
              f"will not work for those. Scroll up for the accessor that failed.")
    if missing_wyckoff:
        print(f"[i] {missing_wyckoff} setting(s) have no Wyckoff table (optional field).")
    if not missing_symops:
        print("All settings carry a full, closed operator list.")
    print("=" * 60)


if __name__ == "__main__":
    main()
