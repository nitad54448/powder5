// swarm_density.wgsl
// Real-space electron density fitness for Charge Flipping

struct Params {
    o0: vec4<f32>, o1: vec4<f32>, o2: vec4<f32>,
    r0: vec4<f32>, r1: vec4<f32>, r2: vec4<f32>,
    nTot: f32, maxSites: f32, numParticles: f32, nGroupsActive: f32, 
    nElem: f32, nBondRules: f32, rMinOff: f32, ruleOff: f32,
    fTabOff: f32, nRefl: f32, penClash: f32, penBond: f32,
    penCoord: f32, penScale: f32, centro: f32, pad0: f32,
    symOff: f32, tablesOff: f32, projOff: f32, pad1: f32,
};

struct ParticleState {
    fit: f32,
    cc: f32,
    pen: f32,
    stepSize: f32,
    temp: f32,
    seed: u32,
    accepted: u32,
    pad: u32,

    // Best-ever state for this particle's chain, kept so the CPU can batch
    // several generations of dispatches without a chain's discovery being
    // lost when a later Metropolis step accepts a worse move and walks away
    // from it. bestFit is a raw comparison value taken at whatever penScale
    // was active when it was written and is NOT valid to compare directly
    // against a later f_new computed at a different penScale - see the
    // rescoring in the accept block below, which mirrors recordBest() on
    // the JS side. -1e30 in bestFit means "nothing recorded yet".
    bestFit: f32,
    bestCC: f32,
    bestPen: f32,
    bestPad: f32,
};

@group(0) @binding(0) var<storage, read_write> particles: array<f32>;
@group(0) @binding(1) var<storage, read> particleAssign: array<u32>;
@group(0) @binding(2) var<storage, read> genPack: array<u32>;
@group(0) @binding(3) var<storage, read> staticFloats: array<f32>;
@group(0) @binding(4) var<storage, read> reflPack: array<u32>;
@group(0) @binding(5) var<storage, read> densityMap: array<f32>;
@group(0) @binding(6) var<storage, read_write> mcmcState: array<ParticleState>;
@group(0) @binding(7) var<uniform> params: Params;

const WG: u32 = 64u;
const MAX_GEN_ATOMS: u32 = 384u; //__MAX_GEN_ATOMS__
const MIN_IMAGE_SHELL: i32 = 0;  //__MIN_IMAGE_SHELL__
const MAX_BOND_RULES: u32 = 8u;
const RULE_STRIDE: u32 = 6u;
const MAX_COORD_SLOTS: u32 = 16u;
const NO_PARTNER_CHARGE: f32 = 5.0;
const COORD_SOFT_A: f32 = 0.25;

// HARD CAPACITIES OF THIS KERNEL. Enforced on the host too - wyckoff_assign.js
// throws in buildAssignmentTables() and buildRestraintTables() before a single
// dispatch - and clamped again here, because WGSL does not trap an
// out-of-range index: it clamps it silently. An overflow therefore shows up
// not as a crash but as two atoms scoring each other's restraints, which is
// indistinguishable from a bad structure.
//
//   MAX_ELEM  - distinct elements in the composition. genPack carries 4 bits
//               of type, so the packing tolerates 16; the binding constraint
//               is local_rMin, which is MAX_ELEM^2 floats of workgroup memory.
//   MAX_SITES - independent Wyckoff sites in one assignment, i.e. the width of
//               prop_sites. The proposal step gives one lane to each site, so
//               this can never usefully exceed the workgroup size of 64.
//
//               32 rather than 16 since a 66-atom molecular structure in a
//               group of order 4 needs 17. The cost is 256 bytes of workgroup
//               memory, 1.6% of the guaranteed 16 kB, and it does not reduce
//               MAX_GEN_ATOMS: swMaxGenAtoms() budgets a flat 20 bytes an atom
//               against 85% of the limit, which already leaves ~3.6 kB spare
//               at the atom cap. Raising it further is cheap in memory and
//               expensive in enumeration: candidate assignments grow steeply
//               with the site budget, and the swarm gains three free
//               coordinates to fit per extra site.
const MAX_ELEM: u32 = 8u;
const MAX_SITES: u32 = 32u;

const TWO_PI: f32 = 6.28318530718;
const NO_NEIGHBOUR: f32 = 1.0e9;

// Murmur-style finalizer, used only to decorrelate the per-lane RNG streams.
//
// The seeding was `seed + lid * 1013u`: an additive offset into one LCG, so
// adjacent lanes ran streams differing by a fixed constant and stayed
// correlated through the permutation. Worse, the seed only advances by the
// number of draws LANE 0 makes, so after a few generations lane 1's stream has
// walked into the region lane 0 already used. Hashing the (seed, lid) pair
// gives each lane an independent-looking stream for three multiplies, once per
// dispatch.
fn hash_u32(x: u32) -> u32 {
    var h = x;
    h ^= h >> 16u;  h *= 0x7feb352du;
    h ^= h >> 15u;  h *= 0x846ca68bu;
    h ^= h >> 16u;
    return h;
}

var<private> rng_state: u32;

fn rand_float() -> f32 {
    rng_state = rng_state * 747796405u + 2891336453u;
    let word = ((rng_state >> ((rng_state >> 28u) + 4u)) ^ rng_state) * 277803737u;
    let res = (word >> 22u) ^ word;
    return f32(res) / 4294967296.0;
}

fn rand_gaussian() -> f32 {
    let u1 = max(rand_float(), 1e-7);
    let u2 = rand_float();
    let r = sqrt(-2.0 * log(u1));
    let theta = TWO_PI * u2;
    return r * cos(theta);
}

fn project_onto_wyckoff(pos: vec3<f32>, p_mat: array<f32, 9>, t_vec: vec3<f32>) -> vec3<f32> {
    var out: vec3<f32>;
    out.x = p_mat[0] * pos.x + p_mat[1] * pos.y + p_mat[2] * pos.z + t_vec.x;
    out.y = p_mat[3] * pos.x + p_mat[4] * pos.y + p_mat[5] * pos.z + t_vec.y;
    out.z = p_mat[6] * pos.x + p_mat[7] * pos.y + p_mat[8] * pos.z + t_vec.z;
    return out - floor(out);
}

var<workgroup> prop_sites: array<vec3<f32>, 32>;   // MAX_SITES
var<workgroup> gx: array<f32, MAX_GEN_ATOMS>;
var<workgroup> gy: array<f32, MAX_GEN_ATOMS>;
var<workgroup> gz: array<f32, MAX_GEN_ATOMS>;
var<workgroup> gT: array<u32, MAX_GEN_ATOMS>;

var<workgroup> rFit: array<f32, WG>;
var<workgroup> rPen: array<f32, WG>;
// Total scattering weight of the atoms this workgroup summed, so the fitness
// can be normalised to a weighted mean rather than left as a bare sum.
var<workgroup> rWgt: array<f32, WG>;
var<workgroup> local_rMin: array<f32, 64>;         // MAX_ELEM * MAX_ELEM

// Cartesian lengths of the 27 shell translations, indexed (a+1)*9+(b+1)*3+(c+1).
//
// Constant for the whole cell, so they are computed once per dispatch rather
// than re-derived for every atom pair. They exist to let the counting shell be
// SKIPPED cheaply: an image displaced by t is at least ||t|| - d away, where d
// is the nearest-image distance already in hand, so one subtract and compare
// rejects a translation without the matrix multiply and square root that
// finding its true distance would cost.
var<workgroup> shellLen: array<f32, 27>;


// COUNTING NEIGHBOURS NEEDS EVERY IMAGE, NOT THE NEAREST ONE.
//
// cartDist below returns a single scalar: the distance to the CLOSEST image of
// atom j. That is the right answer for a clash floor and for "how far is the
// nearest partner", and it is the wrong answer for "how many partners are
// there". Two lattice translates of the same j can both sit inside a
// coordination shell whenever the repeat is short enough -- a cubic perovskite
// B cation has its six oxygens as three symmetry images and their three
// translates, so a "B O 6" rule counted three and charged the structure for the
// other three forever. The host reference in contacts.js has always enumerated
// the neighbouring cells and counted each image separately; this makes the
// kernel agree with it.
//
// COUNT_SHELL is applied ONLY to the rule-matching pairs, so the 27x cost falls
// on the handful of A-B pairs a constraint actually names rather than on the
// whole contact loop. Injected like MIN_IMAGE_SHELL; setting it to 0 restores
// the previous nearest-image-only counting.
const COUNT_SHELL: i32 = 1;  //__COUNT_SHELL__

// The fractional separation, in the reduced basis, wrapped to the nearest
// image. Split out of cartDist so the shell can be walked without recomputing
// the basis transform for every translation.
fn fracDelta(p1: vec3<f32>, p2: vec3<f32>) -> vec3<f32> {
    let d = p1 - p2;
    let q = vec3<f32>(
        params.r0.x*d.x + params.r0.y*d.y + params.r0.z*d.z,
        params.r1.x*d.x + params.r1.y*d.y + params.r1.z*d.z,
        params.r2.x*d.x + params.r2.y*d.y + params.r2.z*d.z);
    return q - round(q);
}

/** Cartesian length of a fractional vector in the reduced basis. */
fn cartLen(t: vec3<f32>) -> f32 {
    let cx = params.o0.x*t.x + params.o0.y*t.y + params.o0.z*t.z;
    let cy = params.o1.x*t.x + params.o1.y*t.y + params.o1.z*t.z;
    let cz = params.o2.x*t.x + params.o2.y*t.y + params.o2.z*t.z;
    return sqrt(cx*cx + cy*cy + cz*cz);
}

fn cartDist(p1: vec3<f32>, p2: vec3<f32>) -> f32 {
    let d = p1 - p2;
    var q = vec3<f32>(
        params.r0.x*d.x + params.r0.y*d.y + params.r0.z*d.z,
        params.r1.x*d.x + params.r1.y*d.y + params.r1.z*d.z,
        params.r2.x*d.x + params.r2.y*d.y + params.r2.z*d.z);
    q = q - round(q);

    var best: f32 = 1.0e30;
    for (var i = -MIN_IMAGE_SHELL; i <= MIN_IMAGE_SHELL; i = i + 1) {
        for (var j = -MIN_IMAGE_SHELL; j <= MIN_IMAGE_SHELL; j = j + 1) {
            for (var k = -MIN_IMAGE_SHELL; k <= MIN_IMAGE_SHELL; k = k + 1) {
                let t = q + vec3<f32>(f32(i), f32(j), f32(k));
                let cx = params.o0.x*t.x + params.o0.y*t.y + params.o0.z*t.z;
                let cy = params.o1.x*t.x + params.o1.y*t.y + params.o1.z*t.z;
                let cz = params.o2.x*t.x + params.o2.y*t.y + params.o2.z*t.z;
                best = min(best, cx*cx + cy*cy + cz*cz);
            }
        }
    }
    return sqrt(best);
}

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wgId: vec3<u32>,
        @builtin(local_invocation_index) lid: u32) {

    let pIdx = wgId.x;
    if (pIdx >= u32(params.numParticles)) { return; }

    let nTot   = min(u32(params.nTot), MAX_GEN_ATOMS);
    let nRules = min(u32(params.nBondRules), MAX_BOND_RULES);
    let ruleOff = u32(params.ruleOff);

    // TWO DIFFERENT NUMBERS, AND THEY MUST NOT BE CONFLATED.
    //
    // siteStride is the stride the HOST allocated `particles` and `siteProj`
    // with, so every index derived from it has to use the raw uniform:
    // clamping it would slide every particle's slice of the buffer.
    //
    // maxSites is what THIS KERNEL can hold in prop_sites. The host refuses an
    // assignment with more than MAX_SITES sites, so on any run that starts the
    // two are equal; the min exists only so that a host/kernel mismatch
    // truncates instead of writing past the end of a workgroup array.
    let siteStride = u32(params.maxSites);
    let maxSites   = min(siteStride, MAX_SITES);

        let A     = particleAssign[pIdx];
        let gBase = A * nTot;
        let pBase = pIdx * siteStride * 3u;

        // `particles` is JS-allocated at double length: the first half is
        // the current, live coordinates (what pBase already indexes into);
        // the second half is a best-ever cache, one slot per particle, that
        // only this kernel ever writes and only a Replica Exchange sync
        // reads back.
        let currentStride  = siteStride * 3u;
        let bestBaseOffset = u32(params.numParticles) * currentStride;
        let bestBase       = bestBaseOffset + pIdx * currentStride;

        let sd = mcmcState[pIdx].stepSize;
        rng_state = hash_u32(mcmcState[pIdx].seed ^ (lid * 0x9E3779B9u));
        
        if (lid < maxSites) {
            let s = lid;
            let c = pBase + s * 3u;
            var pos = vec3<f32>(particles[c], particles[c + 1u], particles[c + 2u]);
            
            if (params.pad0 > 0.5) {
                pos.x = pos.x + rand_gaussian() * sd;
                pos.y = pos.y + rand_gaussian() * sd;
                pos.z = pos.z + rand_gaussian() * sd;
            }
            
            let p_off = u32(params.projOff);
            // siteStride, NOT maxSites. siteProj was allocated by the host at
            // (nAssign x T.maxSites x 12) floats, so the row for assignment A
            // starts at A * siteStride; using the kernel-clamped value would
            // read another assignment's projector the moment the two differ.
            let b = (A * siteStride + s) * 12u;
            var p_mat: array<f32, 9>;
            for (var i = 0u; i < 9u; i = i + 1u) { p_mat[i] = staticFloats[p_off + b + i]; }
            let t_vec = vec3<f32>(staticFloats[p_off + b + 9u], staticFloats[p_off + b + 10u], staticFloats[p_off + b + 11u]);
            
            prop_sites[s] = project_onto_wyckoff(pos, p_mat, t_vec);
        }
        
        if (lid == 0u) {
            mcmcState[pIdx].seed = rng_state;
        }
        workgroupBarrier();

        for (var g = lid; g < nTot; g = g + WG) {
            let pk = genPack[gBase + g];
            let s  = pk & 0xFFu;
            let o  = (pk >> 8u) & 0xFFFu;
            let ty = (pk >> 20u) & 0xFu;

            let p = prop_sites[s];

            // Symmetry operators are packed 12 floats each: r[9] then t[3].
            // The operator base offset `ob` MUST appear in all twelve reads.
            // It previously appeared in exactly one (the r[0] term), so every
            // symmetry image was built from operator 0's matrix and
            // translation - i.e. the identity - and all N images of a site
            // collapsed onto very nearly the same point. That silently turned
            // the whole search into a fit of one unexpanded asymmetric unit,
            // and it rewarded piling atoms on top of each other, because
            // coincident images each collected the same density.
            let s_off = u32(params.symOff);
            let ob = o * 12u;
            let rb = s_off + ob;
            var n = vec3<f32>(
                p.x*staticFloats[rb+0u] + p.y*staticFloats[rb+1u] + p.z*staticFloats[rb+2u] + staticFloats[rb+9u],
                p.x*staticFloats[rb+3u] + p.y*staticFloats[rb+4u] + p.z*staticFloats[rb+5u] + staticFloats[rb+10u],
                p.x*staticFloats[rb+6u] + p.y*staticFloats[rb+7u] + p.z*staticFloats[rb+8u] + staticFloats[rb+11u]
            );
            n = n - floor(n);





        gx[g] = n.x; gy[g] = n.y; gz[g] = n.z;
        // Clamped ONCE, here, rather than at each of the places gT is read.
        // genPack carries 4 bits of type, so a malformed table could present
        // 15; local_rMin holds MAX_ELEM^2 entries and the scattering table has
        // a row of nElem, so either reader would run off the end.
        gT[g] = min(ty, MAX_ELEM - 1u);
    }
    workgroupBarrier();

    // Z-WEIGHTED, NORMALISED DENSITY OVERLAP.
    //
    // This was a bare sum of the interpolated density over every generated
    // atom, which had three faults that between them produced structures with
    // no bonds in them at all:
    //
    //   - The atom TYPE was loaded into gT and then never consulted, so a
    //     light atom counted exactly as much as a heavy one. In PbSO4 that
    //     let 16 oxygens outvote 4 leads 4:1 in an electron density map where
    //     the lead peak is an order of magnitude taller than an oxygen one.
    //
    //   - Nothing normalised the total, so the score grew with the number of
    //     atoms in the cell and ran 0..nTot. The penalty weights below - 0.05
    //     per clash, 0.10 per Angstrom - were chosen against a correlation
    //     coefficient in [0, 1], so against a sum over two dozen atoms they
    //     were roughly twenty times too weak to overrule anything. A 0.43 A
    //     Pb-O contact cost about 1.2 against a fitness of order 20.
    //
    //   - Coincident atoms each collected the FULL peak height, so stacking
    //     was positively rewarded.
    //
    // Weighting by Z and dividing by the total weight fixes all three: the
    // result is a weighted mean density per unit scattering power, bounded in
    // [0, 1] exactly like the correlation the reflection kernel returns, so
    // the penalties are back on the scale they were tuned for and the two
    // kernels are directly comparable. zByType sits at offset 0 of the
    // restraint tables (see buildRestraintTables), so no new binding is
    // needed to reach it.
    var sFit: f32 = 0.0;
    var sWeight: f32 = 0.0;
    let z_off = u32(params.tablesOff);
    let N_grid = u32(params.nGroupsActive); 
    let N_f = f32(N_grid);

    for (var gi = lid; gi < nTot; gi = gi + WG) {
        var fx = gx[gi] * N_f;
        var fy = gy[gi] * N_f;
        var fz = gz[gi] * N_f;
        
        fx = fx - floor(fx / N_f) * N_f;
        fy = fy - floor(fy / N_f) * N_f;
        fz = fz - floor(fz / N_f) * N_f;
        
        let ix = u32(fx) % N_grid;
        let iy = u32(fy) % N_grid;
        let iz = u32(fz) % N_grid;
        
        let tx = fx - floor(fx);
        let ty = fy - floor(fy);
        let tz = fz - floor(fz);
        
        var acc: f32 = 0.0;
        for (var dz = 0u; dz < 2u; dz = dz + 1u) {
            let wz = select(1.0 - tz, tz, dz == 1u);
            let jz = (iz + dz) % N_grid;
            for (var dy = 0u; dy < 2u; dy = dy + 1u) {
                let wy = select(1.0 - ty, ty, dy == 1u);
                let jy = (iy + dy) % N_grid;
                for (var dx = 0u; dx < 2u; dx = dx + 1u) {
                    let wx = select(1.0 - tx, tx, dx == 1u);
                    let jx = (ix + dx) % N_grid;
                    let idx = jx + jy * N_grid + jz * N_grid * N_grid;
                    acc = acc + wx * wy * wz * densityMap[idx];
                }
            }
        }
        // Weight by the atom's atomic number. A heavier atom is expected to
        // sit on a taller peak, and this is what makes the objective say so.
        let zw = max(staticFloats[z_off + gT[gi]], 1.0);
        sFit = sFit + zw * acc;
sWeight = sWeight + zw;
    }

    let t_off = u32(params.tablesOff);
    // nElem is the ROW STRIDE of the r_min table as the host packed it, so it
    // is not clamped; the fill length is, because local_rMin is MAX_ELEM^2.
    let nElem   = u32(params.nElem);
    let nElemSq = min(nElem * nElem, MAX_ELEM * MAX_ELEM);
    for (var idx = lid; idx < nElemSq; idx = idx + WG) {
        local_rMin[idx] = staticFloats[t_off + u32(params.rMinOff) + idx];
    }
    if (lid < 27u) {
        let sa = i32(lid / 9u) - 1;
        let sb = i32((lid / 3u) % 3u) - 1;
        let sc = i32(lid % 3u) - 1;
        shellLen[lid] = cartLen(vec3<f32>(f32(sa), f32(sb), f32(sc)));
    }
    workgroupBarrier();

    var ruleA: array<u32, MAX_BOND_RULES>;
    var ruleB: array<u32, MAX_BOND_RULES>;
    var ruleMin: array<f32, MAX_BOND_RULES>;
    var ruleMax: array<f32, MAX_BOND_RULES>;
    var ruleN: array<f32, MAX_BOND_RULES>;
    var ruleMode: array<u32, MAX_BOND_RULES>;
    for (var k = 0u; k < nRules; k = k + 1u) {
        let rb = ruleOff + k * RULE_STRIDE;
        ruleA[k] = u32(staticFloats[t_off + rb]);
        ruleB[k] = u32(staticFloats[t_off + rb + 1u]);
        ruleMin[k] = staticFloats[t_off + rb + 2u];
        ruleMax[k] = staticFloats[t_off + rb + 3u];
        ruleN[k] = staticFloats[t_off + rb + 4u];
        ruleMode[k] = u32(staticFloats[t_off + rb + 5u]);
    }

    // slotOff[k] is where rule k's nearest-neighbour list starts inside `slot`,
    // and slotN[k] is how many entries it actually got. The two are separate
    // because the coordination numbers can sum past MAX_COORD_SLOTS: the host
    // refuses that up front, but if the two ever drift apart this truncates the
    // last rule rather than letting it write over the first one's distances --
    // an overflow WGSL clamps in silence and that looks exactly like a
    // structure which cannot satisfy its own restraints.
    var slotOff: array<u32, MAX_BOND_RULES>;
    var slotN:   array<u32, MAX_BOND_RULES>;
    {
        var acc: u32 = 0u;
        for (var k = 0u; k < nRules; k = k + 1u) {
            slotOff[k] = acc;
            var nk: u32 = 0u;
            if (ruleMode[k] != 0u) {
                nk = min(u32(ruleN[k]), MAX_COORD_SLOTS - acc);
            }
            slotN[k] = nk;
            acc = acc + nk;
        }
    }

    var pen: f32 = 0.0;
    for (var i = lid; i < nTot; i = i + WG) {
        let ti = gT[i];
        let pi = vec3<f32>(gx[i], gy[i], gz[i]);

        var myRules: u32 = 0u;
        for (var k = 0u; k < nRules; k = k + 1u) {
            if (ruleA[k] == ti) { myRules = myRules | (1u << k); }
        }
        var nearest: array<f32, MAX_BOND_RULES>;
        var inWindow: array<f32, MAX_BOND_RULES>;
        for (var k = 0u; k < MAX_BOND_RULES; k = k + 1u) {
            nearest[k] = NO_NEIGHBOUR;
            inWindow[k] = 0.0;
        }
        var slot: array<f32, MAX_COORD_SLOTS>;
        for (var m = 0u; m < MAX_COORD_SLOTS; m = m + 1u) { slot[m] = NO_NEIGHBOUR; }

        for (var j = 0u; j < nTot; j = j + 1u) {
            let tj = gT[j];
            let pj = vec3<f32>(gx[j], gy[j], gz[j]);
            let d = cartDist(pi, pj);

            // The clash floor is a statement about the CLOSEST approach, so it
            // stays on the nearest image. j == i is excluded because an atom
            // cannot collide with itself; its lattice translates are handled
            // by the counting pass below, which can legitimately treat them as
            // neighbours in a short-repeat structure.
            if (j != i) {
                let dmin = local_rMin[ti * nElem + tj];
                if (j > i && d < dmin) {
                    let overlap = (dmin - d) / max(dmin, 1e-3);
                    pen = pen + params.penClash * (1.0 + 9.0 * overlap);
                }
            }

            if (myRules == 0u) { continue; }

            // Computed once per PAIR, not once per rule. Two rules naming the
            // same A-B pair were each redoing the basis transform below.
            let qf = fracDelta(pi, pj);
            // The length of the t = 0 image, which is what the triangle
            // inequality in the shell below needs: |qf + t| >= ||t|| - |qf|.
            //
            // NOT `d`. d is cartDist, the minimum over the MIN_IMAGE_SHELL
            // search, so d <= |qf| and using it INFLATES the lower bound -- the
            // skip then rejects translations that are genuinely inside the
            // window. Identical to |qf| whenever MIN_IMAGE_SHELL is 0, which is
            // why this never showed up on a well-shaped cell; the thin cells
            // that set it to 1 are exactly the ones that need the shell.
            let dq = cartLen(qf);

            for (var k = 0u; k < nRules; k = k + 1u) {
                if ((myRules & (1u << k)) == 0u) { continue; }
                if (ruleB[k] != tj) { continue; }

                let nk = slotN[k];
                let counted = (nk > 0u);

                // The nearest partner is exactly what cartDist already returns,
                // so it is taken from d rather than re-derived in the shell.
                // That keeps the mode-0 test exact even when the shell below
                // prunes away the distant images.
                if (j != i) { nearest[k] = min(nearest[k], d); }

                // EVERY IMAGE OF j THAT CAN MATTER, COUNTED SEPARATELY.
                //
                // Distances outside the window still enter the slot list: the
                // deficit is graded on how far each of the N nearest partners
                // misses by, and dropping the far ones would collapse "partner
                // at 2.5 A, wanted 1.65" and "no partner at all" onto one flat
                // penalty, taking the gradient the search descends with it.
                //
                // "Can matter" is the cutoff: an image must either reach the
                // window or beat the worst slot currently held. An image
                // displaced by t is at least ||t|| - d away, so a translation
                // whose length exceeds cutoff + d cannot qualify and is
                // rejected on a subtract, without the matrix multiply and
                // square root its true distance would cost. In a cell wide
                // enough for the host to have compiled COUNT_SHELL to 0 this
                // loop is a single pass; in the thin cells that need the shell
                // it typically keeps two or three translations out of 27.
                var cutoff = ruleMax[k];
                if (counted) { cutoff = max(cutoff, slot[slotOff[k] + nk - 1u]); }
                for (var sa = -COUNT_SHELL; sa <= COUNT_SHELL; sa = sa + 1) {
                for (var sb = -COUNT_SHELL; sb <= COUNT_SHELL; sb = sb + 1) {
                for (var sc = -COUNT_SHELL; sc <= COUNT_SHELL; sc = sc + 1) {
                    let si = u32((sa + 1) * 9 + (sb + 1) * 3 + (sc + 1));
                    if (shellLen[si] - dq > cutoff) { continue; }
                    let dd = cartLen(qf + vec3<f32>(f32(sa), f32(sb), f32(sc)));
                    if (dd < 1.0e-4) { continue; }   // the atom itself
                    if (dd >= ruleMin[k] && dd <= ruleMax[k]) {
                        inWindow[k] = inWindow[k] + 1.0;
                    }
                    if (counted) {
                        let lo = slotOff[k];
                        let hi = lo + nk;
                        if (dd < slot[hi - 1u]) {
                            var m = hi - 1u;
                            loop {
                                if (m > lo && slot[m - 1u] > dd) {
                                    slot[m] = slot[m - 1u];
                                    m = m - 1u;
                                } else { break; }
                            }
                            slot[m] = dd;
                            cutoff = max(ruleMax[k], slot[hi - 1u]);
                        }
                    }
                }}}
            }
        }

        for (var k = 0u; k < nRules; k = k + 1u) {
            if ((myRules & (1u << k)) == 0u) { continue; }

            if (ruleMode[k] == 0u) {
                let nk = nearest[k];
                if (nk < NO_NEIGHBOUR && nk > ruleMax[k]) {
                    pen = pen + params.penBond * (nk - ruleMax[k]);
                }
            } else {
                let nk = slotN[k];
                var deficit: f32 = 0.0;
                for (var m = 0u; m < nk; m = m + 1u) {
                    let dm = slot[slotOff[k] + m];
                    if (dm >= NO_NEIGHBOUR) {
                        deficit = deficit + NO_PARTNER_CHARGE;
                    } else {
                        let miss = max(0.0, dm - ruleMax[k]) + max(0.0, ruleMin[k] - dm);
                        deficit = deficit + miss * (1.0 + miss / COORD_SOFT_A);
                    }
                }
                pen = pen + params.penBond * deficit;

                if (ruleMode[k] == 1u) {
                    pen = pen + params.penCoord * max(0.0, inWindow[k] - ruleN[k]);
                }
            }
        }
    }

    pen = pen * params.penScale;

    rFit[lid] = sFit; 
    rPen[lid] = pen;
    rWgt[lid] = sWeight;
    workgroupBarrier();

    for (var stride = WG / 2u; stride > 0u; stride = stride >> 1u) {
        if (lid < stride) {
            rFit[lid] = rFit[lid] + rFit[lid + stride];
            rPen[lid] = rPen[lid] + rPen[lid + stride];
            rWgt[lid] = rWgt[lid] + rWgt[lid + stride];
        }
        workgroupBarrier();
    }

if (lid == 0u) {
        // Normalise to a weighted MEAN density, in [0, 1] because the host
        // scales the map so its maximum is 1. This is the quantity the
        // penalties are weighted against, and it is comparable between
        // assignments with different numbers of atoms - a bare sum was not,
        // and quietly favoured whichever assignment put more atoms in the
        // cell.
        let cc = rFit[0] / max(rWgt[0], 1e-9);
        let pen = rPen[0];
        let f_new = cc - pen;
        
        var st = mcmcState[pIdx];

        // RE-SCORE THE INCUMBENT AT THE CURRENT PENALTY WEIGHT.
        //
        // st.fit was measured at whatever penScale was live when this chain
        // last moved, and the host ramps penScale EVERY generation
        // (penRampStart 0.2 -> penRampEnd 4.0). A straight comparison
        // therefore tests a proposal charged at today's weight against an
        // incumbent charged at an older, gentler one, and the gap grows
        // monotonically with the ramp: by the second half of a run essentially
        // nothing is accepted and the ladder stops sampling.
        //
        // st.cc and st.pen are stored raw (weight 1) for exactly this reason --
        // the same convention exchangeSweep(), recordBest() and the
        // bestRescored comparison below all already use.
        //
        // The guard covers a chain's FIRST evaluation, where the host seeds
        // st.fit with -inf and st.cc with NaN. Without it the NaN propagates
        // and the opening proposal of every restart is rejected.
        var f_old = st.fit;
        if (f_old > -1.0e29) { f_old = st.cc - st.pen * params.penScale; }
        
        if (params.pad0 > 0.5) {
            var accept = false;
            if (f_new >= f_old) {
                accept = true;
            } else {
                rng_state = st.seed;
                let prob = exp((f_new - f_old) / max(st.temp, 1e-9));
                if (rand_float() < prob) {
                    accept = true;
                }
                st.seed = rng_state;
            }
            
            if (accept) {
                st.fit = f_new;
                st.cc = cc;
                st.pen = pen / max(params.penScale, 1e-9);
                st.accepted = st.accepted + 1u;
                
                for (var s = 0u; s < maxSites; s = s + 1u) {
                    let c = pBase + s * 3u;
                    particles[c] = prop_sites[s].x;
                    particles[c + 1u] = prop_sites[s].y;
                    particles[c + 2u] = prop_sites[s].z;
                }

                // Historical best for this chain, at weight 1 (raw), so it
                // can be re-scored against whatever penScale is current at
                // read-back time - the same trick recordBest() uses on the
                // JS side to compare archived entries fairly against a
                // moving penalty ramp. Comparing f_new directly to a stale
                // st.bestFit here would be wrong: an easy structure found
                // early under a gentle ramp would look unbeatable forever.
                let hasBest = st.bestFit > -1.0e29;
                let bestRescored = select(-1.0e30, st.bestCC - st.bestPen * params.penScale, hasBest);
                if (f_new > bestRescored) {
                    st.bestFit = f_new;
                    st.bestCC  = cc;
                    st.bestPen = pen / max(params.penScale, 1e-9);

                    for (var s = 0u; s < maxSites; s = s + 1u) {
                        let c = bestBase + s * 3u;
                        particles[c]      = prop_sites[s].x;
                        particles[c + 1u] = prop_sites[s].y;
                        particles[c + 2u] = prop_sites[s].z;
                    }
                }
            }
            
            let STEP_UP = 1.057596;
            let STEP_DOWN = 0.9762857;
            if (accept) {
                st.stepSize = clamp(st.stepSize * STEP_UP, 2e-4, 0.5);
            } else {
                st.stepSize = clamp(st.stepSize * STEP_DOWN, 2e-4, 0.5);
            }
            
            mcmcState[pIdx] = st;
        } else {
            st.fit = f_new;
            st.cc = cc;
            st.pen = pen / max(params.penScale, 1e-9);
            mcmcState[pIdx] = st;
        }
    }
}