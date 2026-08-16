// density3d.js
// version 140, 2 august 2026
//
// Isosurface view of the charge-flipping electron density, drawn with
// three.js. Companion to the 2D section viewer: same map, same absolute
// scale, three dimensions instead of one slice.
//
// ---------------------------------------------------------------------------
//  Why Surface Nets and not Marching Cubes
// ---------------------------------------------------------------------------
//  Marching Cubes needs a 256-entry triangle table of about 2500 hand-entered
//  integers. A single wrong number produces holes that only show up on
//  particular topologies, and there is no way to check the table short of
//  reproducing it from a trusted source.
//
//  Surface Nets (Gibson 1998) needs no such table. Its two lookup tables are
//  GENERATED at load time from the geometry of the cube, in eleven lines, and
//  the whole thing is verifiable: the test in the header comment of
//  surfaceNets() checks that the extracted mesh is closed, which is the
//  property a broken table destroys first.
//
//  It also produces a better surface for this particular job. Surface Nets
//  places one vertex per cell, at the centroid of the edge crossings, which
//  gives a smoother mesh with far fewer degenerate slivers than Marching
//  Cubes on noisy data -- and a charge-flipping map at 1.5 sigma is noisy.
//  The cost is that sharp features get rounded, which electron density does
//  not have.
//
// ---------------------------------------------------------------------------
//  Requires THREE (r147-r160, UMD build) as a global. Degrades to a message
//  in the panel if it is absent, so the rest of the app is unaffected.
// ---------------------------------------------------------------------------

(function (global) {
'use strict';

// ---------------------------------------------------------------------------
//  Lattice geometry comes from crystal.js, which the document loads
//  immediately before this file. cellBasis and fracToCart used to be defined
//  here as well as in powder5.html and (in metric-tensor form) in
//  charge_flipping_worker.js.
//
//  In a browser crystal.js's declarations are already global and the scope
//  chain finds them; under node they are module-scoped, so pull them in.
// ---------------------------------------------------------------------------
if (typeof cellBasis === 'undefined'
        && typeof require === 'function' && typeof module !== 'undefined') {
    Object.assign(global, require('./crystal.js'));
}

// ===========================================================================
//  Lookup tables, generated rather than transcribed.
//
//  cube_edges holds the 12 edges of the unit cube as pairs of corner indices,
//  where a corner index is a 3-bit number whose bits are its x, y, z
//  coordinates. Two corners share an edge exactly when their indices differ
//  in one bit, which is what the XOR below enumerates.
//
//  edge_table[m] is then the 12-bit mask of edges crossed by the surface when
//  the corner sign pattern is m: an edge is crossed when its two corners fall
//  on opposite sides.
// ===========================================================================
const cube_edges = new Int32Array(24);
const edge_table = new Int32Array(256);

(function initTables() {
    let k = 0;
    for (let i = 0; i < 8; ++i) {
        for (let j = 1; j <= 4; j <<= 1) {
            const p = i ^ j;
            if (i <= p) { cube_edges[k++] = i; cube_edges[k++] = p; }
        }
    }
    for (let i = 0; i < 256; ++i) {
        let em = 0;
        for (let j = 0; j < 24; j += 2) {
            const a = !!(i & (1 << cube_edges[j]));
            const b = !!(i & (1 << cube_edges[j + 1]));
            em |= (a !== b) ? (1 << (j >> 1)) : 0;
        }
        edge_table[i] = em;
    }
})();

// ===========================================================================
//  surfaceNets(data, dims)
//
//  data  : Float32Array of dims[0]*dims[1]*dims[2] samples, x fastest.
//          NEGATIVE is inside the surface.
//  dims  : [nx, ny, nz] sample counts (cells = dims-1 along each axis).
//
//  Returns { positions: Float32Array (3 per vertex, in sample coordinates),
//            indices:   Uint32Array  (3 per triangle) }.
//
//  Verified against an analytic sphere: every mesh edge is shared by exactly
//  two triangles (closed manifold) and every vertex lies within half a sample
//  spacing of the true radius.
// ===========================================================================
function surfaceNets(data, dims) {
    const verts = [];
    const faces = [];
    const dx = dims[0], dy = dims[1], dz = dims[2];

    const R = new Int32Array([1, dx + 1, (dx + 1) * (dy + 1)]);
    const grid = new Float64Array(8);
    const buffer = new Int32Array(R[2] * 2);
    const x = new Int32Array(3);
    const v = [0, 0, 0];

    let buf_no = 1;
    let n = 0;

    for (x[2] = 0; x[2] < dz - 1; ++x[2], n += dx, buf_no ^= 1, R[2] = -R[2]) {
        // Pointer into the two-slab ring buffer that remembers which vertex
        // each cell produced, so faces can be stitched to the neighbours
        // already visited.
        let m = 1 + (dx + 1) * (1 + buf_no * (dy + 1));

        for (x[1] = 0; x[1] < dy - 1; ++x[1], ++n, m += 2)
        for (x[0] = 0; x[0] < dx - 1; ++x[0], ++n, ++m) {

            // The eight corner values of this cell, and the sign mask.
            let mask = 0, g = 0, idx = n;
            for (let k = 0; k < 2; ++k, idx += dx * (dy - 2))
            for (let j = 0; j < 2; ++j, idx += dx - 2)
            for (let i = 0; i < 2; ++i, ++g, ++idx) {
                const p = data[idx];
                grid[g] = p;
                mask |= (p < 0) ? (1 << g) : 0;
            }

            // Entirely inside or entirely outside: no surface here.
            if (mask === 0 || mask === 0xff) continue;

            const edge_mask = edge_table[mask];
            v[0] = 0; v[1] = 0; v[2] = 0;
            let e_count = 0;

            for (let i = 0; i < 12; ++i) {
                if (!(edge_mask & (1 << i))) continue;
                ++e_count;
                const e0 = cube_edges[i << 1];
                const e1 = cube_edges[(i << 1) + 1];
                const g0 = grid[e0], g1 = grid[e1];
                let t = g0 - g1;
                if (Math.abs(t) > 1e-10) t = g0 / t; else { --e_count; continue; }
                for (let j = 0, k = 1; j < 3; ++j, k <<= 1) {
                    const a = e0 & k, b = e1 & k;
                    if (a !== b) v[j] += a ? 1.0 - t : t;
                    else        v[j] += a ? 1.0 : 0.0;
                }
            }
            if (e_count === 0) continue;

            const s = 1.0 / e_count;
            buffer[m] = verts.length / 3;
            verts.push(x[0] + s * v[0], x[1] + s * v[1], x[2] + s * v[2]);

            // Each of the three axis-aligned edges through the cell's origin
            // corner that the surface crosses closes a quad with the three
            // cells already visited around that edge.
            for (let i = 0; i < 3; ++i) {
                if (!(edge_mask & (1 << i))) continue;
                const iu = (i + 1) % 3, iv = (i + 2) % 3;
                if (x[iu] === 0 || x[iv] === 0) continue;
                const du = R[iu], dv = R[iv];
                const a = buffer[m], b = buffer[m - du],
                      c = buffer[m - du - dv], d = buffer[m - dv];
                if (mask & 1) {
                    faces.push(a, b, c, a, c, d);
                } else {
                    faces.push(a, d, c, a, c, b);
                }
            }
        }
    }

    return {
        positions: Float32Array.from(verts),
        indices: Uint32Array.from(faces),
        vertexCount: verts.length / 3,
        triangleCount: faces.length / 3
    };
}

// ===========================================================================
//  Crystallography
// ===========================================================================

// MOVED. cellBasis and fracToCart now live in crystal.js (see the loader at
// the top of this IIFE), alongside metricTensor, cellVolume and the d-spacing
// helpers. They are used here exactly as before.

// Whole-map statistics, so the isolevel can be quoted in sigma exactly as the
// 2D section caption does.
function mapStats(map) {
    let lo = Infinity, hi = -Infinity, sum = 0, sum2 = 0;
    for (let i = 0; i < map.length; i++) {
        const v = map[i];
        if (v < lo) lo = v;
        if (v > hi) hi = v;
        sum += v; sum2 += v * v;
    }
    const n = map.length || 1;
    const mean = sum / n;
    return { lo, hi, mean, sigma: Math.sqrt(Math.max(0, sum2 / n - mean * mean)) };
}

// ===========================================================================
//  Build the isosurface for one unit cell.
//
//  The grid is periodic, so the sample block is built with an extra plane on
//  each axis (index N is a copy of index 0, positioned at N). Marching the
//  wrapped grid directly instead would stitch the cell's far face back to its
//  near face and drag triangles right across the box.
//
//  `stride` decimates the grid before extraction. The UI no longer exposes it
//  -- the full grid is fast enough that offering a coarse mode only invited
//  people to look at a worse surface than they had to -- but the parameter is
//  kept because a 128^3 map is two million cells and a future caller on a
//  weaker machine may want it.
// ===========================================================================
function buildIsosurface(map, N, level, stride) {
    const st = Math.max(1, Math.floor(stride || 1));
    const M = Math.max(2, Math.floor(N / st));         // cells per axis
    const S = M + 1;                                   // samples per axis
    const N2 = N * N;

    const data = new Float32Array(S * S * S);
    let p = 0;
    for (let k = 0; k < S; k++) {
        const zk = ((k * st) % N) * N2;
        for (let j = 0; j < S; j++) {
            const yj = ((j * st) % N) * N;
            for (let i = 0; i < S; i++) {
                // Negative inside, so the sign mask marks the dense region.
                data[p++] = level - map[((i * st) % N) + yj + zk];
            }
        }
    }

    const mesh = surfaceNets(data, [S, S, S]);
    return { mesh, samplesPerAxis: S, scale: st / N };  // sample index -> fractional
}

// ===========================================================================
//  Viewer
// ===========================================================================
// One colour per distinct site so a heavy atom and its light neighbours can be
// told apart. Sixteen entries rather than eight: a ten-site solution used to
// wrap around and give sites 9 and 10 the same colours as sites 1 and 2, which
// is precisely the ambiguity the legend exists to remove.
const SITE_PALETTE = [
    0x5ac8fa, 0xffd60a, 0x30d158, 0xff9f0a,
    0xbf5af2, 0xff6482, 0x64d2ff, 0xa2845e,
    0xff453a, 0x32ade6, 0xffd426, 0x66d4cf,
    0xda8fff, 0xac8e68, 0x8e8cd8, 0x4cd964
];

// ---------------------------------------------------------------------------
//  Per-site colour overrides, rank -> 0xRRGGBB.
//
//  The palette above is the DEFAULT. A user who has just solved a structure
//  knows things the program does not -- which lobe is the heavy atom, which
//  two sites are the same element -- and the only way to express that in a
//  picture is to colour the spheres. Overrides live here rather than in
//  powder5.html so that siteColor(), buildSiteGroup() and the legend swatch
//  cannot disagree about what colour a site is, which was the whole reason
//  siteColor() was exposed in the first place.
//
//  Keyed by rank, so they survive a re-render of the same structure. They are
//  NOT meaningful across a new solution -- rank 3 of one solution has nothing
//  to do with rank 3 of the next -- so the caller clears them when a fresh
//  structure arrives.
// ---------------------------------------------------------------------------
/** @type {Map<number, number>} */
const SITE_COLOR_OVERRIDES = new Map();

/** Default palette colour for a rank, as 0xRRGGBB. @param {number} rank */
function defaultSiteColorInt(rank) {
    const i = (Number.isFinite(rank) ? rank : 1) - 1;
    const n = SITE_PALETTE.length;
    return SITE_PALETTE[((i % n) + n) % n];
}

/** Effective colour for a rank, override first. @param {number} rank */
function siteColorInt(rank) {
    const o = SITE_COLOR_OVERRIDES.get(rank);
    return (o === undefined) ? defaultSiteColorInt(rank) : o;
}

/**
 * Parses '#rrggbb' (or 'rrggbb') to an integer, or null if unparseable.
 * @param {string} hex
 * @returns {number|null}
 */
function parseHexColor(hex) {
    if (typeof hex !== 'string') return null;
    const m = hex.trim().match(/^#?([0-9a-fA-F]{6})$/);
    return m ? parseInt(m[1], 16) : null;
}

// Spheres for one expanded site list. Shared geometry and material per rank,
// so a fifty-position cell costs two GPU resources per distinct site rather
// than two per sphere.
function buildSiteGroup(t, basis, sites) {
    const THREE = t.THREE;
    if (!Array.isArray(sites) || !sites.length) return null;

    const group = new THREE.Group();
    const scaleRef = Math.cbrt(Math.abs(
        basis.av[0] * (basis.bv[1] * basis.cv[2] - basis.bv[2] * basis.cv[1]) +
        basis.bv[0] * (basis.cv[1] * basis.av[2] - basis.cv[2] * basis.av[1]) +
        basis.cv[0] * (basis.av[1] * basis.bv[2] - basis.av[2] * basis.bv[1])));
    const baseRad = Math.max(0.04, scaleRef * 0.030);
    const geoms = new Map(), mats = new Map();

    for (const st of sites) {
        const rank = Number.isFinite(st.rank) ? st.rank : 1;
        const rel = Number.isFinite(st.relative) ? Math.max(0.12, st.relative) : 1;
        if (!geoms.has(rank)) {
            // Radius follows the cube root of peak height, roughly how ionic
            // radius tracks Z.
            geoms.set(rank, new THREE.SphereGeometry(baseRad * (0.55 + 0.75 * Math.cbrt(rel)), 16, 12));
            mats.set(rank, new THREE.MeshPhongMaterial({
                color: siteColorInt(rank), shininess: 60
            }));
        }
        const c = fracToCart(basis, st.x, st.y, st.z);
        const m = new THREE.Mesh(geoms.get(rank), mats.get(rank));
        m.position.set(c[0], c[1], c[2]);
        group.add(m);
    }
    group.userData.sharedGeometries = [...geoms.values()];
    group.userData.sharedMaterials = [...mats.values()];
    // Kept so setSiteColor can repaint one rank in place. Recolouring by
    // rebuilding the group would drop and reallocate every sphere for what is
    // a single uniform change -- visible as a hitch while dragging a colour
    // picker, which fires continuously.
    group.userData.materialsByRank = mats;
    return group;
}

const Density3D = {
    _three: null,
    _state: null,

    isAvailable() {
        return typeof global.THREE !== 'undefined' && !!global.THREE.WebGLRenderer;
    },

    // Sets up the renderer once. Safe to call repeatedly.
    init(canvas) {
        if (!this.isAvailable() || !canvas) return false;
        if (this._three && this._three.canvas === canvas) return true;

        const THREE = global.THREE;
        let renderer;
        try {
            renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
        } catch (e) {
            console.warn('Density3D: WebGL unavailable:', e);
            return false;
        }
        renderer.setPixelRatio(Math.min(2, global.devicePixelRatio || 1));

        const scene = new THREE.Scene();
        const camera = new THREE.PerspectiveCamera(38, 1, 0.05, 5000);
        const root = new THREE.Group();
        scene.add(root);

        // Two lights and a weak ambient: enough to read the shape of a lobe
        // without the far side going pure black. Intensities are set by
        // _applyTheme, because what reads well on black washes out on white.
        const key = new THREE.DirectionalLight(0xffffff, 2.0);
        key.position.set(1, 1.2, 1.4);
        const fill = new THREE.DirectionalLight(0x88aaff, 0.7);
        fill.position.set(-1, -0.8, -0.6);
        const ambient = new THREE.AmbientLight(0xffffff, 0.45);
        scene.add(key, fill, ambient);

        this._three = {
            THREE, canvas, renderer, scene, camera, root, key, fill, ambient,
            surface: null, cell: null, sites: null, basis: null,
            // Spherical camera pose. Written directly by the pointer handlers
            // rather than through a controls addon, which keeps the dependency
            // to the single UMD build file.
            theta: -1.05, phi: 1.15, radius: 3, target: new THREE.Vector3(),
            needsRender: true
        };

        this._attachControls();

        // The results panel can be display:none when the surface is first
        // built, in which case clientWidth is 0 and the renderer would size
        // itself to a fallback and then never correct. Watching the element
        // covers both that and ordinary window resizing.
        if (typeof global.ResizeObserver === 'function') {
            const ro = new global.ResizeObserver(() => { if (this._three) this._three.needsRender = true; });
            ro.observe(canvas);
            this._three.resizeObserver = ro;
        }

        this._loop();
        return true;
    },

    _attachControls() {
        const t = this._three;
        const el = t.canvas;
        let dragging = false, lastX = 0, lastY = 0;

        el.addEventListener('pointerdown', (e) => {
            dragging = true; lastX = e.clientX; lastY = e.clientY;
            el.setPointerCapture(e.pointerId);
        });
        el.addEventListener('pointerup', (e) => {
            dragging = false;
            try { el.releasePointerCapture(e.pointerId); } catch (_) { /* already gone */ }
        });
        el.addEventListener('pointercancel', () => { dragging = false; });
        el.addEventListener('pointermove', (e) => {
            if (!dragging) return;
            t.theta -= (e.clientX - lastX) * 0.008;
            t.phi   -= (e.clientY - lastY) * 0.008;
            // Stop short of the poles, where the up vector flips and the view
            // snaps through 180 degrees.
            t.phi = Math.max(0.05, Math.min(Math.PI - 0.05, t.phi));
            lastX = e.clientX; lastY = e.clientY;
            t.needsRender = true;
        });
        el.addEventListener('wheel', (e) => {
            e.preventDefault();
            t.radius *= Math.exp(e.deltaY * 0.0012);
            t.radius = Math.max(0.4, Math.min(50, t.radius));
            t.needsRender = true;
        }, { passive: false });
    },

    _loop() {
        const t = this._three;
        if (!t) return;
        const step = () => {
            if (!this._three) return;
            if (t.needsRender) { this._draw(); t.needsRender = false; }
            global.requestAnimationFrame(step);
        };
        global.requestAnimationFrame(step);
    },

    _draw() {
        const t = this._three;
        const canvas = t.canvas;
        const w = canvas.clientWidth || 400, h = canvas.clientHeight || 400;
        if (canvas.width !== w || canvas.height !== h) {
            t.renderer.setSize(w, h, false);
            t.camera.aspect = w / Math.max(1, h);
            t.camera.updateProjectionMatrix();
        }
        const r = t.radius * (t.frameRadius || 1);

        // C IS VERTICAL.
        //
        // Three.js takes +Y as up by default, and fracToCart puts a along x,
        // b along y and c along z -- so the default camera stood the cell on
        // its b axis. Every crystallographic drawing convention puts c up, and
        // a structure shown on the wrong axis is genuinely hard to read
        // against a published figure.
        //
        // Fixed at the CAMERA rather than by rotating the scene: the geometry
        // stays in the crystallographic Cartesian frame, so coordinates,
        // distances and the a/b/c labels all keep meaning exactly what they
        // say. Rotating the model would have made the axis letters lie.
        //
        // phi is now the polar angle from +c and theta the azimuth in the a-b
        // plane, which is also what the drag handler already assumes.
        t.camera.up.set(0, 0, 1);
        t.camera.position.set(
            t.target.x + r * Math.sin(t.phi) * Math.cos(t.theta),
            t.target.y + r * Math.sin(t.phi) * Math.sin(t.theta),
            t.target.z + r * Math.cos(t.phi)
        );
        t.camera.lookAt(t.target);
        t.renderer.render(t.scene, t.camera);
    },

    invalidate() { if (this._three) this._three.needsRender = true; },

    // -----------------------------------------------------------------------
    //  THEME
    //
    //  The canvas is created with alpha:true and no scene background, so the
    //  page shows through and the viewer inherits the light or dark surface it
    //  sits on. That is the right choice, but it means every colour in the
    //  scene was tuned against black and none of it was re-tuned against
    //  white: a 0.45 white ambient plus a 2.0 key washes a light orange
    //  isosurface out to near-white, and the grey cell edges disappear
    //  altogether. The fix is not a darker background -- it is to light the
    //  scene for the surface it is drawn on.
    //
    //  Called on init and again whenever the theme is switched, so it never
    //  needs the viewer to be rebuilt.
    // -----------------------------------------------------------------------
    _themeIsDark() {
        try {
            return (global.document.documentElement.getAttribute('data-theme') || 'dark') !== 'light';
        } catch (e) { return true; }
    },

    /**
     * Re-light and re-colour the scene for the current theme.
     * @param {boolean} [dark]  Defaults to reading data-theme off <html>.
     */
    applyTheme(dark) {
        const t = this._three;
        if (!t) return;
        const isDark = (dark === undefined) ? this._themeIsDark() : !!dark;
        t.isDark = isDark;

        // Less ambient and a weaker key on white: on a light page the eye has
        // no dark reference, so the same intensities read as glare.
        t.ambient.intensity = isDark ? 0.45 : 0.28;
        t.key.intensity     = isDark ? 2.00 : 1.25;
        t.fill.intensity    = isDark ? 0.70 : 0.45;
        // A blue fill on a white page tints the shadows lilac; neutral there.
        t.fill.color.setHex(isDark ? 0x88aaff : 0xbfcbdd);

        if (t.surface && t.surface.material) {
            // Deeper and more saturated in light mode. The dark-mode amber has
            // roughly the same lightness as the page it would sit on.
            t.surface.material.color.setHex(isDark ? 0xffb648 : 0xd2761a);
            t.surface.material.specular.setHex(isDark ? 0x332211 : 0x6b5236);
            t.surface.material.opacity = isDark ? 0.92 : 0.95;
            t.surface.material.needsUpdate = true;
        }
        if (t.cell) {
            // The nine non-origin edges are the grey ones; the three coloured
            // axis edges stay as they are, since red/green/blue read on both.
            t.cell.traverse(o => {
                if (o.isLineSegments && o.material) {
                    o.material.opacity = isDark ? 0.9 : 1.0;
                    o.material.needsUpdate = true;
                }
            });
        }
        t.needsRender = true;
    },


    // -----------------------------------------------------------------------
    //  Rebuild everything from a map. Returns a small report for the caption.
    // -----------------------------------------------------------------------
    // The colour a site is drawn in, as a CSS hex string. Exposed so a legend
    // can label the spheres without keeping a second copy of the palette that
    // could drift out of step with this one.
    siteColor(rank) {
        return '#' + siteColorInt(rank).toString(16).padStart(6, '0');
    },

    /** The palette default for a rank, ignoring any override. */
    defaultSiteColor(rank) {
        return '#' + defaultSiteColorInt(rank).toString(16).padStart(6, '0');
    },

    /** True if this rank has been recoloured away from the palette default. */
    hasSiteColorOverride(rank) {
        return SITE_COLOR_OVERRIDES.has(rank);
    },

    /**
     * Recolours one site. Repaints the existing material in place, so this is
     * cheap enough to call from a colour picker's continuous `input` event.
     *
     * @param {number} rank  Site rank.
     * @param {string|null} hex '#rrggbb', or null to restore the palette default.
     * @returns {boolean} True if the colour was applied.
     */
    setSiteColor(rank, hex) {
        if (!Number.isFinite(rank)) return false;
        if (hex === null || hex === undefined || hex === '') {
            SITE_COLOR_OVERRIDES.delete(rank);
        } else {
            const v = parseHexColor(hex);
            if (v === null) return false;
            SITE_COLOR_OVERRIDES.set(rank, v);
        }
        return this._repaintRank(rank);
    },

    /** Drops every override and repaints. Call when a new structure arrives. */
    resetSiteColors() {
        const ranks = [...SITE_COLOR_OVERRIDES.keys()];
        SITE_COLOR_OVERRIDES.clear();
        for (const r of ranks) this._repaintRank(r);
        return ranks.length;
    },

    /** Pushes the current colour for `rank` into the live material, if any. */
    _repaintRank(rank) {
        const t = this._three;
        const mats = t && t.sites && t.sites.userData && t.sites.userData.materialsByRank;
        if (!mats) return true;              // nothing drawn yet; colour applies on next build
        const m = mats.get(rank);
        if (m && m.color) {
            m.color.setHex(siteColorInt(rank));
            t.needsRender = true;
        }
        return true;
    },

    // Redraw only the atom spheres. Toggling one site off used to go through
    // show(), which re-ran marching cubes over the whole grid for a change
    // that touches nothing but a handful of spheres.
    updateSites(sites) {
        const t = this._three;
        if (!t || !t.basis) return false;
        this._disposeObject(t.sites);
        t.sites = buildSiteGroup(t, t.basis, sites);
        if (t.sites) t.root.add(t.sites);
        t.needsRender = true;
        return true;
    },

    show(opts) {
        if (!this._three) return null;
        const t = this._three, THREE = t.THREE;
        const { map, gridSize, cell, sigmaLevel, stride, showCell, showSites, sites } = opts;
        // showSurface:false draws the cell and the atoms with NO density. The
        // isosurface is by far the most expensive thing here, so it is skipped
        // outright rather than built and hidden.
        const showSurface = opts.showSurface !== false;
        if (!map || !gridSize || !cell) return null;

        const stats = mapStats(map);
        const level = (sigmaLevel || 3) * (stats.sigma || 1);

        const t0 = (global.performance && performance.now) ? performance.now() : Date.now();
        const iso = showSurface ? buildIsosurface(map, gridSize, level, stride) : null;
        const t1 = (global.performance && performance.now) ? performance.now() : Date.now();

        const basis = cellBasis(cell);

        // Sample index -> fractional -> Cartesian, in place.
        const pos = iso ? iso.mesh.positions : null;
        if (iso) {
            const s = iso.scale;
            for (let i = 0; i < pos.length; i += 3) {
                const c = fracToCart(basis, pos[i] * s, pos[i + 1] * s, pos[i + 2] * s);
                pos[i] = c[0]; pos[i + 1] = c[1]; pos[i + 2] = c[2];
            }
        }

        // --- surface ---
        this._disposeObject(t.surface);
        t.surface = null;
        if (iso && iso.mesh.triangleCount > 0) {
            const geo = new THREE.BufferGeometry();
            geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
            geo.setIndex(new THREE.BufferAttribute(iso.mesh.indices, 1));
            geo.computeVertexNormals();
            const mat = new THREE.MeshPhongMaterial({
                color: 0xffb648, specular: 0x332211, shininess: 22,
                side: THREE.DoubleSide, flatShading: false,
                transparent: true, opacity: 0.92
            });
            t.surface = new THREE.Mesh(geo, mat);
            t.root.add(t.surface);
        }

        // --- unit cell wireframe, with RGB axes ---
        //
        // The nine edges that do not touch the origin stay grey; the three that
        // do are coloured a = red, b = green, c = blue and carry a letter, so
        // the orientation is readable at a glance from any viewpoint. Without
        // that there is no way to tell which axis you are looking down, and a
        // structure viewed along the wrong axis is easy to misread.
        this._disposeObject(t.cell);
        t.cell = null;
        if (showCell !== false) {
            const corners = [];
            for (let k = 0; k < 2; k++) for (let j = 0; j < 2; j++) for (let i = 0; i < 2; i++) {
                corners.push(fracToCart(basis, i, j, k));
            }
            // Vertex colours, so these are fixed when the geometry is built
            // rather than adjusted later by applyTheme. The axis triplet reads
            // on either background; the neutral edges do not -- 0.44 grey is
            // nearly invisible on a white page, so it darkens in light mode.
            const lightMode = !this._themeIsDark();
            const AXIS_COLORS = lightMode
                ? [[0.80, 0.18, 0.16], [0.16, 0.60, 0.20], [0.16, 0.38, 0.80]]
                : [[0.95, 0.30, 0.28], [0.42, 0.85, 0.38], [0.36, 0.62, 0.98]];
            const GREY = lightMode ? [0.34, 0.37, 0.42] : [0.44, 0.47, 0.51];

            const edges = [], colors = [];
            for (let A = 0; A < 8; A++) for (let B = A + 1; B < 8; B++) {
                // Corners are adjacent when their bit patterns differ in one place.
                const diff = A ^ B;
                if (diff !== 1 && diff !== 2 && diff !== 4) continue;
                edges.push(...corners[A], ...corners[B]);
                // Corner 0 is the origin; diff identifies which axis the edge
                // runs along (1 = a, 2 = b, 4 = c).
                const c = (A === 0) ? AXIS_COLORS[diff === 1 ? 0 : diff === 2 ? 1 : 2] : GREY;
                colors.push(...c, ...c);
            }

            const cellGroup = new THREE.Group();
            const g = new THREE.BufferGeometry();
            g.setAttribute('position', new THREE.Float32BufferAttribute(edges, 3));
            g.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
            const lines = new THREE.LineSegments(g,
                new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.9 }));
            cellGroup.add(lines);

            // Axis letters, drawn into a small canvas and hung on a sprite so
            // they always face the camera.
            const labelDefs = lightMode
                ? [['a', [1, 0, 0], '#c22a24'],
                   ['b', [0, 1, 0], '#2a992f'],
                   ['c', [0, 0, 1], '#2a5fc9']]
                : [['a', [1, 0, 0], '#f24d47'],
                   ['b', [0, 1, 0], '#6bd961'],
                   ['c', [0, 0, 1], '#5c9efa']];
            const labelScale = Math.hypot(...fracToCart(basis, 1, 1, 1)) * 0.09;
            for (const [text, f, col] of labelDefs) {
                const cv = document.createElement('canvas');
                cv.width = cv.height = 64;
                const cx = cv.getContext('2d');
                cx.font = 'bold 46px system-ui, sans-serif';
                cx.fillStyle = col;
                cx.textAlign = 'center';
                cx.textBaseline = 'middle';
                cx.fillText(text, 32, 34);
                const tex = new THREE.CanvasTexture(cv);
                const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false }));
                const p = fracToCart(basis, f[0] * 1.12, f[1] * 1.12, f[2] * 1.12);
                sp.position.set(p[0], p[1], p[2]);
                sp.scale.set(labelScale, labelScale, labelScale);
                cellGroup.add(sp);
            }

            t.cell = cellGroup;
            t.root.add(cellGroup);
        }

        // The surface and the cell have just been rebuilt, so re-apply the
        // theme to the new materials. Without this a scene built while light
        // mode was already active kept the dark-mode lighting until the user
        // happened to toggle the theme.
        this.applyTheme();

        // --- atom positions ---
        //
        // The caller passes positions ALREADY EXPANDED by the symmetry
        // operators. The site list itself is one representative per orbit, so
        // drawing it raw put a single sphere next to one of the four lobes of
        // a multiplicity-4 site and left the other three bare, which looks
        // exactly like a bug in the map.
        this._disposeObject(t.sites);
        t.basis = basis;
        t.sites = showSites ? buildSiteGroup(t, basis, sites) : null;
        if (t.sites) t.root.add(t.sites);

        // Centre on the cell body-diagonal midpoint and frame it.
        const mid = fracToCart(basis, 0.5, 0.5, 0.5);
        t.target.set(mid[0], mid[1], mid[2]);
        const diag = Math.hypot(...fracToCart(basis, 1, 1, 1));
        t.frameRadius = Math.max(1e-3, diag * 0.75);
        t.needsRender = true;

        return {
            level,
            sigma: stats.sigma,
            maxRho: stats.hi,
            surfaceDrawn: !!iso,
            triangles: iso ? iso.mesh.triangleCount : 0,
            vertices: iso ? iso.mesh.vertexCount : 0,
            cells: iso ? (iso.samplesPerAxis - 1) ** 3 : 0,
            ms: Math.round(t1 - t0)
        };
    },

    setVisibility(what, on) {
        const t = this._three;
        if (!t) return;
        if (what === 'cell' && t.cell) t.cell.visible = !!on;
        if (what === 'sites' && t.sites) t.sites.visible = !!on;
        t.needsRender = true;
    },

    resetView() {
        const t = this._three;
        if (!t) return;
        t.theta = -1.05; t.phi = 1.15; t.radius = 3;
        t.needsRender = true;
    },

    _disposeObject(obj) {
        const t = this._three;
        if (!obj || !t) return;
        t.root.remove(obj);
        // Geometries and materials shared across many meshes are held on the
        // group's userData and disposed once; disposing them per-mesh would
        // free the same resource repeatedly.
        const shared = new Set();
        const ud = obj.userData || {};
        for (const list of [ud.sharedGeometries, ud.sharedMaterials]) {
            if (Array.isArray(list)) list.forEach(r => { shared.add(r); r.dispose?.(); });
        }
        if (obj.traverse) {
            obj.traverse(o => {
                if (o.geometry && !shared.has(o.geometry)) o.geometry.dispose?.();
                if (o.material && !shared.has(o.material)) {
                    if (o.material.map) o.material.map.dispose?.();
                    o.material.dispose?.();
                }
            });
        }
    },

    dispose() {
        const t = this._three;
        if (!t) return;
        if (t.resizeObserver) t.resizeObserver.disconnect();
        [t.surface, t.cell, t.sites].forEach(o => this._disposeObject(o));
        t.renderer.dispose();
        this._three = null;
    },

    // Exposed for the unit test.
    _internals: { surfaceNets, cellBasis, fracToCart, mapStats, buildIsosurface }
};

global.Density3D = Density3D;
if (typeof module !== 'undefined' && module.exports) module.exports = Density3D;

})(typeof globalThis !== 'undefined' ? globalThis : this);
