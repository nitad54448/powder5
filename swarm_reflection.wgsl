// swarm_reflection.wgsl
// Reciprocal-space |Fcalc|^2 fitness for the Wyckoff search.
//
// ---------------------------------------------------------------------------
//  WHAT THIS IS, AND WHY IT IS A SEPARATE FILE
// ---------------------------------------------------------------------------
//  swarm_density.wgsl scores a candidate structure by how much CHARGE FLIPPING
//  DENSITY sits under its atoms. This file scores the same candidate by how
//  well its calculated intensities match the PAWLEY INTENSITIES. Everything
//  else -- the RNG, the Wyckoff projection, the symmetry expansion, the
//  Metropolis step, the bond-rule penalties, the reduction -- is identical and
//  deliberately kept identical, so the two objectives stay comparable and a
//  fix to the shared machinery can be applied to both by diffing them.
//
//  Two objectives, two files, rather than one file with a mode switch: the
//  scoring sections share no arithmetic at all, and a `if (mode)` in the
//  hottest loop in the program would cost both paths for the benefit of
//  neither.
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
const MAX_COORD_SLOTS: u32 = 16u;
const NO_PARTNER_CHARGE: f32 = 5.0;
const COORD_SOFT_A: f32 = 0.25;
const MAX_ELEM: u32 = 16u;
const TWO_PI: f32 = 6.28318530718;
const NO_NEIGHBOUR: f32 = 1.0e9;

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

var<workgroup> prop_sites: array<vec3<f32>, 64>;
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
var<workgroup> local_rMin: array<f32, 256>;

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

    let nTot     = min(u32(params.nTot), MAX_GEN_ATOMS);
    let maxSites = u32(params.maxSites);
    let nRules   = min(u32(params.nBondRules), MAX_BOND_RULES);
    let rMinOff  = u32(params.rMinOff);
    let ruleOff  = u32(params.ruleOff);

        let A     = particleAssign[pIdx];
        let gBase = A * nTot;
        let pBase = pIdx * maxSites * 3u;

        // `particles` is JS-allocated at double length: the first half is
        // the current, live coordinates (what pBase already indexes into);
        // the second half is a best-ever cache, one slot per particle, that
        // only this kernel ever writes and only a Replica Exchange sync
        // reads back.
        let currentStride  = maxSites * 3u;
        let bestBaseOffset = u32(params.numParticles) * currentStride;
        let bestBase        = bestBaseOffset + pIdx * currentStride;

       let sd = mcmcState[pIdx].stepSize;
        rng_state = mcmcState[pIdx].seed + lid * 1013u;
        
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
            let b = (A * maxSites + s) * 12u;
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
        gT[g] = ty;
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
    let nElem   = u32(params.nElem);
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
                    let ph = TWO_PI * (h * gx[a] + k * gy[a] + l * gz[a]);
                    fr = fr + groupData[fBase + gT[a]] * cos(ph);
                }
            } else {
                for (var a = 0u; a < nTot; a = a + 1u) {
                    let ph = TWO_PI * (h * gx[a] + k * gy[a] + l * gz[a]);
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
    let nElemSq = u32(params.nElem) * u32(params.nElem);
    for (var idx = lid; idx < nElemSq; idx = idx + WG) {
        local_rMin[idx] = staticFloats[t_off + u32(params.rMinOff) + idx];
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

    var slotOff: array<u32, MAX_BOND_RULES>;
    {
        var acc: u32 = 0u;
        for (var k = 0u; k < nRules; k = k + 1u) {
            slotOff[k] = acc;
            if (ruleMode[k] != 0u) { acc = acc + u32(ruleN[k]); }
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
            if (j == i) { continue; }
            let tj = gT[j];
            let d = cartDist(pi, vec3<f32>(gx[j], gy[j], gz[j]));

            let dmin = local_rMin[ti * u32(params.nElem) + tj];
            if (j > i && d < dmin) {
                let overlap = (dmin - d) / max(dmin, 1e-3);
                pen = pen + params.penClash * (1.0 + 9.0 * overlap);
            }

            if (myRules != 0u) {
                for (var k = 0u; k < nRules; k = k + 1u) {
                    if ((myRules & (1u << k)) == 0u) { continue; }
                    if (ruleB[k] != tj) { continue; }
                    nearest[k] = min(nearest[k], d);
                    if (d >= ruleMin[k] && d <= ruleMax[k]) {
                        inWindow[k] = inWindow[k] + 1.0;
                    }
                    let nk = u32(ruleN[k]);
                    if (ruleMode[k] != 0u && nk > 0u) {
                        let lo = slotOff[k];
                        let hi = lo + nk;              
                        if (d < slot[hi - 1u]) {
                            var m = hi - 1u;
                            loop {
                                if (m > lo && slot[m - 1u] > d) {
                                    slot[m] = slot[m - 1u];
                                    m = m - 1u;
                                } else { break; }
                            }
                            slot[m] = d;
                        }
                    }
                }
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
                let nk = u32(ruleN[k]);
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
        let f_old = st.fit;
        
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