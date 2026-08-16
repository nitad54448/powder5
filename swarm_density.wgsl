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

var<workgroup> rFit: array<f32, WG>;
var<workgroup> rPen: array<f32, WG>;
// Total scattering weight of the atoms this workgroup summed, so the fitness
// can be normalised to a weighted mean rather than left as a bare sum.
var<workgroup> rWgt: array<f32, WG>;
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