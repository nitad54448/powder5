// swarm_reflection.wgsl
// Reciprocal-space |Fcalc|^2 fitness for the Wyckoff search.
//
// ---------------------------------------------------------------------------
//  WHAT THIS IS, AND WHY IT IS THE ONLY SWARM KERNEL
// ---------------------------------------------------------------------------
//  There used to be a sibling, swarm_density.wgsl, which scored a candidate by
//  how much CHARGE FLIPPING DENSITY sat under its atoms. It is DELETED, and
//  the objective with it -- see the reasoning at the top of wyRunWithDevice()
//  in wyckoff_worker.js. In short: point-sampling a map couples this method to
//  the charge-flipping result it exists to be independent of, and the sample
//  is biased in the wrong direction, since a heavy atom's own
//  series-termination ripples dig a trough underneath it and rank it BELOW a
//  light one.
//
//  What that means for anyone editing this file: the shared machinery -- the
//  RNG, the Wyckoff projection, the symmetry expansion, the Metropolis step,
//  the bond-rule penalties, the reduction -- no longer has a second copy to be
//  kept in step with. A fix here is the whole fix. If you resurrect a second
//  objective, resurrect the diffing discipline with it.
//
//  WHY THE DIRECT-SPACE OBJECTIVE IS WORTH HAVING. The density objective can
//  only be as good as the map, and on a heavy-atom structure the map is the
//  weak link: in PbSO4 the lead ripples run at the same amplitude as the
//  oxygen peaks, so the light atoms are not merely hard to see, they are
//  underneath the noise floor. This objective never forms a map. It compares
//  |Fcalc|^2 against |Fobs|^2 directly, so it has no phase problem, no origin
//  or hand ambiguity, and no ripples -- and because the model carries the
//  composition, it cannot invent an atom where there is no scattering to
//  explain it.
//
//  It is also better conditioned. PbSO4 in Pnma has eleven free positional
//  parameters against several hundred observations. Charge flipping solves for
//  thousands of grid amplitudes and phases from the same data.
//
// ---------------------------------------------------------------------------
//  BINDING 5 IS THE ONLY INTERFACE DIFFERENCE
// ---------------------------------------------------------------------------
//  The density kernel binds the map there. This one binds `groupData`, which
//  the host packs as two regions in one buffer:
//
//      [0 .. fTabOff)              per-group metadata, GROUP_STRIDE floats each
//      [fTabOff .. end)            scattering factors, nRefl * nElem floats
//
//  The scattering factors are ALREADY f(s) * exp(-B s^2 / 4), tabulated per
//  (reflection, element) by swScatteringTable on the host. That is why this
//  kernel needs no d-spacings, no B, and no form-factor evaluation: the entire
//  s-dependence has been folded in before the buffer was written.
// ---------------------------------------------------------------------------

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
// [group metadata | scattering factor table]; see the header.
@group(0) @binding(5) var<storage, read> groupData: array<f32>;
@group(0) @binding(6) var<storage, read_write> mcmcState: array<ParticleState>;
@group(0) @binding(7) var<uniform> params: Params;

const WG: u32 = 64u;
const MAX_GEN_ATOMS: u32 = 384u; //__MAX_GEN_ATOMS__
const MIN_IMAGE_SHELL: i32 = 0;  //__MIN_IMAGE_SHELL__
// Floats per group in the metadata region.
//   3 -> start, count, Iobs                  (every group weighted equally)
//   4 -> start, count, Iobs, weight          (weight is typically 1/sigma^2)
// Injected by the host exactly like MAX_GEN_ATOMS, so the weighted and
// unweighted packers can share this kernel. It is a compile-time constant, so
// the `weight` branch below folds away entirely when it is 3.
const GROUP_STRIDE: u32 = 3u; //__GROUP_STRIDE__
const MAX_BOND_RULES: u32 = 8u;
const RULE_STRIDE: u32 = 6u;
// `slot` is a FUNCTION-SCOPE array, not workgroup storage, so this constant
// costs nothing against the 16 kB budget that MAX_GEN_ATOMS is competing for.
// It is per-lane scratch: 64 f32 is 256 bytes a thread, and because
// slot[slotOff[k] + m] is dynamically indexed the array was already spilling
// out of registers at 16 -- raising it resizes an existing spill rather than
// creating one.
//
// 64 = MAX_BOND_RULES * 8, so the sum-of-coordination-numbers check on the
// host stops binding for any chemistry that fits in eight rules at all. The
// per-atom initialisation below runs to slotTotal, the slots the rules
// actually claimed, so a run needing 20 does not pay for the other 44.
const MAX_COORD_SLOTS: u32 = 64u;
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

var<workgroup> rPen: array<f32, WG>;
// Accumulators for a WEIGHTED PEARSON CORRELATION between Icalc and Iobs over
// the active groups. Six partial sums rather than one, because the correlation
// needs both means and both variances, and each thread only sees the groups it
// happened to be handed.
//
// Pearson rather than a plain cosine similarity, and rather than an R factor:
//   - It is scale-free, so there is NO SCALE FACTOR to fit. That removes a
//     parameter, and more importantly removes a degeneracy the swarm would
//     otherwise have to explore.
//   - Subtracting the means makes it discriminating. Icalc and Iobs are both
//     non-negative and share a large positive offset; a cosine similarity is
//     dominated by that offset and reads high for almost any structure.
//   - It lands in [0, 1] after clamping, which is the range the bond and
//     clash penalties below were tuned against. An R factor would have needed
//     every penalty weight retuned.
var<workgroup> rSw:  array<f32, WG>;   // sum w
var<workgroup> rSc:  array<f32, WG>;   // sum w * Icalc
var<workgroup> rSo:  array<f32, WG>;   // sum w * Iobs
var<workgroup> rScc: array<f32, WG>;   // sum w * Icalc^2
var<workgroup> rSoo: array<f32, WG>;   // sum w * Iobs^2
var<workgroup> rSco: array<f32, WG>;   // sum w * Icalc * Iobs
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
    // Hoisted above the generation loop so gT can be clamped against the
    // COMPOSITION rather than only against MAX_ELEM. Both tables gT indexes
    // are nElem wide, not MAX_ELEM wide: local_rMin is filled to nElem^2 and
    // the scattering row for a reflection is nElem long. A type of 7 in a
    // three-element composition was therefore inside the old clamp and still
    // read past the filled region of one and into the next reflection's row of
    // the other. Unreachable while genType comes from the host's demand index,
    // which is what this clamp is for.
    let nElem      = u32(params.nElem);
    let nElemClamp = min(nElem, MAX_ELEM);
    let maxType    = max(nElemClamp, 1u) - 1u;
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
        // RAW params.nTot, not the clamped `nTot`, for exactly the reason
        // siteStride is raw above: this is the stride the HOST allocated
        // genPack with (nAssign x nTot entries), so clamping it would slide
        // every assignment's row rather than shortening it. The loop below
        // still runs to the clamped `nTot`, so a host/kernel mismatch
        // truncates the cell -- which the reduction and the restraints both
        // survive -- instead of reading another assignment's atoms, which
        // they do not.
        let gBase = A * u32(params.nTot);
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
        gT[g] = min(ty, maxType);
    }
    workgroupBarrier();

    // ----------------------------------------------------------------------
    //  |Fcalc|^2 AGAINST |Fobs|^2, GROUP BY GROUP
    //
    //  Work is split over GROUPS, not over atoms: one thread owns a whole
    //  group and walks its reflections serially. Every thread therefore reads
    //  the full generated cell out of gx/gy/gz/gT, which is a workgroup-memory
    //  broadcast rather than a gather, and no thread has to combine a partial
    //  structure factor with another thread's.
    //
    //  A "group" is a set of reflections the powder pattern cannot separate.
    //  The host has already summed their observed intensities into one number,
    //  so the comparison is made where the measurement actually exists --
    //  summing the calculated intensities the same way. This is the same
    //  reason the report prints overlap clusters rather than trusting the
    //  individual Pawley intensities of a heavily overlapped block.
    //
    //  Iobs is packed as sum(mult * Fo^2) over the group, so Icalc is
    //  accumulated as sum(mult * |F|^2) to match it term for term.
    // ----------------------------------------------------------------------
    // nElem is declared at the top of main() so the generation loop's type
    // clamp can use it; it is the row stride of both the r_min table and the
    // scattering table, and is deliberately NOT clamped where it is used as a
    // stride - see the fill length at nElemSq below.
    let fTabOff = u32(params.fTabOff);
    let nGroups = u32(params.nGroupsActive);
    let nReflTot = u32(params.nRefl);
    let isCentro = params.centro > 0.5;

    var sw: f32 = 0.0;
    var sc: f32 = 0.0;
    var so: f32 = 0.0;
    var scc: f32 = 0.0;
    var soo: f32 = 0.0;
    var sco: f32 = 0.0;

    for (var g = lid; g < nGroups; g = g + WG) {
        let gb    = g * GROUP_STRIDE;
        let start = u32(groupData[gb]);
        let count = u32(groupData[gb + 1u]);
        let iObs  = groupData[gb + 2u];

        // Compile-time: with GROUP_STRIDE == 3 this whole branch disappears.
        var wgt: f32 = 1.0;
        if (GROUP_STRIDE >= 4u) { wgt = groupData[gb + 3u]; }
        if (!(wgt > 0.0)) { continue; }

        var iCalc: f32 = 0.0;
        for (var m = 0u; m < count; m = m + 1u) {
            let r = start + m;
            if (r >= nReflTot) { break; }

            // hkl is bit-packed ten bits per index, biased by 512 so negative
            // indices survive. Unpack with a mask, then remove the bias.
            let pk = reflPack[r * 2u];
            let h = f32(i32(pk & 0x3FFu) - 512);
            let k = f32(i32((pk >> 10u) & 0x3FFu) - 512);
            let l = f32(i32((pk >> 20u) & 0x3FFu) - 512);
            let mult = f32(reflPack[r * 2u + 1u]);

            let fBase = fTabOff + r * nElem;

            var fr: f32 = 0.0;
            var fi: f32 = 0.0;
            if (isCentro) {
                // The cell contents are closed under inversion through the
                // origin, so every atom at r has a partner at -r of the same
                // type and the sine terms cancel exactly. F is real, and the
                // hottest loop in the program loses half its transcendentals.
                //
                // The host TESTS this rather than assuming it from the space
                // group: a centrosymmetric group whose inversion centre is not
                // at the origin carries (-I | t) with t nonzero, for which F
                // is complex.
                for (var a = 0u; a < nTot; a = a + 1u) {
                    // Reduce to one turn BEFORE scaling by 2*pi. At h,k,l up to
                    // +/-511 the unreduced argument reaches ~3200 rad, outside
                    // the range WGSL guarantees for cos/sin and outside what
                    // several backends' hardware instructions reduce well.
                    let q  = h * gx[a] + k * gy[a] + l * gz[a];
                    let ph = TWO_PI * (q - floor(q));
                    fr = fr + groupData[fBase + gT[a]] * cos(ph);
                }
            } else {
                for (var a = 0u; a < nTot; a = a + 1u) {
                    // See the centrosymmetric branch: reduced to one turn
                    // before scaling, so cos/sin never see a large argument.
                    let q  = h * gx[a] + k * gy[a] + l * gz[a];
                    let ph = TWO_PI * (q - floor(q));
                    let f  = groupData[fBase + gT[a]];
                    fr = fr + f * cos(ph);
                    fi = fi + f * sin(ph);
                }
            }
            iCalc = iCalc + mult * (fr * fr + fi * fi);
        }

        sw  = sw  + wgt;
        sc  = sc  + wgt * iCalc;
        so  = so  + wgt * iObs;
        scc = scc + wgt * iCalc * iCalc;
        soo = soo + wgt * iObs  * iObs;
        sco = sco + wgt * iCalc * iObs;
    }

    let t_off = u32(params.tablesOff);
    // nElem (declared with the reflection block above) is the ROW STRIDE of the
    // r_min table as the host packed it, so it is not clamped; the fill length
    // is, because local_rMin is MAX_ELEM^2.
    let nElemSq = min(nElem * nElem, MAX_ELEM * MAX_ELEM);
    for (var idx = lid; idx < nElemSq; idx = idx + WG) {
        local_rMin[idx] = staticFloats[t_off + u32(params.rMinOff) + idx];
    }
    if (lid < 27u) {
        // ta/tb/tc, not sa/sb/sc: `sc` is the running sum(w * Icalc) in this
        // function's scope. WGSL lets the inner declaration shadow it, so the
        // old spelling compiled and worked -- and would have stopped working,
        // silently, the moment either declaration moved.
        let ta = i32(lid / 9u) - 1;
        let tb = i32((lid / 3u) % 3u) - 1;
        let tc = i32(lid % 3u) - 1;
        shellLen[lid] = cartLen(vec3<f32>(f32(ta), f32(tb), f32(tc)));
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
    // How many slots the rules between them actually claimed. Hoisted out of
    // the block below because the per-atom loop needs it: `slot` is re-declared
    // on every iteration of that loop and has to be re-initialised each time,
    // and initialising all MAX_COORD_SLOTS entries would make the cost of the
    // constant proportional to the constant rather than to the demand. Entries
    // at or above slotTotal are never read - no slotOff[k] + m reaches them -
    // so stopping there is not an approximation.
    var slotTotal: u32 = 0u;
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
        slotTotal = acc;
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
        for (var m = 0u; m < slotTotal; m = m + 1u) { slot[m] = NO_NEIGHBOUR; }

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

    rPen[lid] = pen;
    rSw[lid]  = sw;
    rSc[lid]  = sc;
    rSo[lid]  = so;
    rScc[lid] = scc;
    rSoo[lid] = soo;
    rSco[lid] = sco;
    workgroupBarrier();

    for (var stride = WG / 2u; stride > 0u; stride = stride >> 1u) {
        if (lid < stride) {
            rPen[lid] = rPen[lid] + rPen[lid + stride];
            rSw[lid]  = rSw[lid]  + rSw[lid + stride];
            rSc[lid]  = rSc[lid]  + rSc[lid + stride];
            rSo[lid]  = rSo[lid]  + rSo[lid + stride];
            rScc[lid] = rScc[lid] + rScc[lid + stride];
            rSoo[lid] = rSoo[lid] + rSoo[lid + stride];
            rSco[lid] = rSco[lid] + rSco[lid + stride];
        }
        workgroupBarrier();
    }

if (lid == 0u) {
        // Weighted Pearson correlation between Icalc and Iobs, clamped into
        // [0, 1] so it occupies the same range the density kernel returns and
        // the penalty weights below keep their tuning.
        //
        // A zero variance on either side means the comparison carries no
        // information: every group calculated the same intensity (all atoms
        // stacked on one point, or no atoms at all), or the observations are
        // flat. Returning 0 there is correct and is what the host expects --
        // random starting positions legitimately score 0.
        let wSum = max(rSw[0], 1e-9);
        let mc = rSc[0] / wSum;
        let mo = rSo[0] / wSum;
        let varC = rScc[0] / wSum - mc * mc;
        let varO = rSoo[0] / wSum - mo * mo;
        let cov  = rSco[0] / wSum - mc * mo;
        var cc: f32 = 0.0;
        let denom = varC * varO;
        if (denom > 1e-20) {
            cc = clamp(cov / sqrt(denom), 0.0, 1.0);
        }
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