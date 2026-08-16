// charge_flipping.wgsl
// version 138, 2 august 2026
//
// WebGPU kernels for dual-space charge flipping on powder data, with the
// space group applied DIRECTLY inside the iteration loop.
//
// ---------------------------------------------------------------------------
//  What changed relative to v137, and why
// ---------------------------------------------------------------------------
//  * BINDING COLLISION FIXED. v137 declared obsIdx/obsAmp at @binding(4)/(5)
//    and then declared groups/groupIdx at the SAME two binding points. WGSL
//    permits one resource per binding point per module, so createShaderModule
//    rejected the whole file, runChargeFlippingGPU always threw, and the host
//    silently fell back to a CPU path that was itself broken. Nothing ever
//    ran. The dead obsIdx/obsAmp declarations are gone.
//
//  * TWO BIND GROUPS. Group 0 carries the per-pass state (uniform + the two
//    complex grids + the reduction scratch); group 1 carries the static
//    symmetry tables. The FFT pipelines are built with a layout containing
//    group 0 only, so no single pipeline ever exceeds the default limit of
//    8 storage buffers per stage.
//
//  * PING-PONG INSTEAD OF COPY. v137 ran a copy_b_to_a pass after every one
//    of the 3*log2(N) Stockham stages, doubling the memory traffic of every
//    transform. Here the host alternates two bind groups whose src/dst are
//    swapped, so the copy happens at most once per transform (and only when
//    3*log2(N) is odd).
//
//  * SYMMETRY. Three new kernels put the space group into the loop:
//        symmetrize   - imposes F(hR) = F(h)*exp(-2*pi*i*h.t) over each orbit
//        constrain    - forces EQUAL |F| across a symmetry orbit
//        zero_absent  - forces F = 0 on systematically absent reflections
//    plus `repartition`, which shares the observed intensity between the
//    DISTINCT orbits that overlap in 2-theta (the genuinely powder-specific
//    ambiguity). v137 had this exactly backwards: it repartitioned intensity
//    freely among members of one symmetry orbit -- which are required by
//    symmetry to be equal -- and did nothing at all about real overlaps.
//
//  * THE BEST ITERATE IS CHOSEN ON THE GPU (in its own bind group, group 2;
//    see the declaration for why it cannot live in group 0). v138 copied partials[0..3] back
//    and mapped it EVERY cycle, purely so the host could decide whether the
//    R factor had improved and, if so, issue a buffer copy. That is a full
//    pipeline flush per iteration: at 64^3 the GPU work in one cycle is small
//    enough that the round-trip latency dominated the whole run. keep_best
//    and commit_best do the comparison and the snapshot on the device, so the
//    host never has to see R in order to make progress. It still reads the R
//    history back for the convergence plot, but in batches, and being a few
//    iterations behind costs nothing there.
//
//  * Orbit LOST ITS target_i FIELD. It was written by the host and read by
//    nobody: repartition overwrites orbitTarget[] from the calculated
//    intensities on the first cycle and constrain reads only that. Two
//    sources of truth for the same quantity, one of them dead.
//
//  * NO PER-ITERATION UNIFORM WRITES. Every field of Params is static for a
//    trial. The flip threshold is computed on the GPU by finalize_moments and
//    left in partials[0], so the host no longer has to read sigma back,
//    compute delta and write it into a uniform between passes.
//
// ---------------------------------------------------------------------------
//  Data layout
// ---------------------------------------------------------------------------
//  Complex grids are interleaved f32: 2*i = real, 2*i+1 = imaginary, with
//  i = h + k*N + l*N*N and Miller indices wrapped into 0..N-1.
//
//  partials[] is one buffer with a fixed header:
//        [0] delta      - flip threshold, written by finalize_moments
//        [1] R numerator
//        [2] R denominator
//        [3] best R so far, initialised to +inf by the host
//        [4 + 2*w], [5 + 2*w]  - per-workgroup reduction results
//
//  orbitIdx[] packs a grid index in bits 0..30 and a "this member stores the
//  Friedel mate, so the stored value is the complex conjugate" flag in bit 31.
//  The largest supported grid is 128^3 = 2097152 < 2^31, so the packing is
//  safe with a wide margin.
// ---------------------------------------------------------------------------

struct Params {
    n            : u32,   // grid edge
    axis         : u32,   // 0 = x, 1 = y, 2 = z  (FFT only)
    stage_len    : u32,   // L: half-size of the current sub-transform
    inverse      : u32,   // 1 for the inverse transform

    delta_sigma  : f32,   // flip threshold in units of the map's own sigma
    n_orbits     : u32,
    n_clusters   : u32,
    absent_start : u32,   // offset into orbitIdx of the absent-reflection list

    absent_count : u32,
    sym_lambda   : f32,   // 0 = pure P1, 1 = strict space-group symmetry
    seed         : u32,
    wg_full      : u32,   // workgroups used by reduce_moments

    wg_cluster   : u32,   // workgroups used by reduce_rfactor
    _pad0        : u32,
    _pad1        : u32,
    _pad2        : u32,
};

@group(0) @binding(0) var<uniform> P : Params;
@group(0) @binding(1) var<storage, read_write> src      : array<f32>;
@group(0) @binding(2) var<storage, read_write> dst      : array<f32>;
@group(0) @binding(3) var<storage, read_write> partials : array<f32>;
// (the best-iterate grid is declared with the other group-2 resource below)

// A symmetry orbit: `count` consecutive entries of orbitIdx starting at
// `start`. Its share of the observed intensity lives in orbitTarget[g], which
// repartition writes and constrain reads; the orbit itself carries no copy.
struct Orbit {
    start    : u32,
    count    : u32,
};

// A set of orbits that overlap in 2-theta and were measured as one peak.
// obs_i is the observed total; the split between the member orbits is free.
struct Cluster {
    orbit_start : u32,
    n_orbits    : u32,
    obs_i       : f32,
};

@group(1) @binding(0) var<storage, read>       orbits      : array<Orbit>;
@group(1) @binding(1) var<storage, read>       orbitIdx    : array<u32>;
@group(1) @binding(2) var<storage, read>       phase       : array<f32>;  // 2 per member
@group(1) @binding(3) var<storage, read>       clusters    : array<Cluster>;
@group(1) @binding(4) var<storage, read_write> orbitTarget : array<f32>;

// The best iterate so far. Written only by keep_best.
//
// IT HAS ITS OWN BIND GROUP, and that is not cosmetic. The per-stage storage
// buffer limit is validated against the PIPELINE LAYOUT, not against what an
// entry point actually reads. Group 0 carries three storage buffers and group
// 1 carries five, so the symmetry pipelines -- which bind both -- already sit
// exactly on the guaranteed limit of eight. Putting `best` in group 0 made
// them nine and createPipelineLayout threw, even though symmetrize never
// touches it. Only keep_best binds group 2, so only keep_best pays for it:
// 3 + 0 + 1 = 4.
//
// Eight is the floor every WebGPU adapter must support. Plenty of hardware
// reports 16 and we could ask for it in requiredLimits, but that would refuse
// the device outright on anything older, which is a bad trade for one buffer.
@group(2) @binding(0) var<storage, read_write> best : array<f32>;

const PI : f32 = 3.141592653589793;

const P_DELTA : u32 = 0u;
const P_RNUM  : u32 = 1u;
const P_RDEN  : u32 = 2u;
const P_BESTR : u32 = 3u;
const P_BASE  : u32 = 4u;

const IDX_MASK : u32 = 0x7fffffffu;
const CONJ_BIT : u32 = 0x80000000u;

fn voxel_offset(x : u32, y : u32, z : u32, n : u32) -> u32 {
    return 2u * ((z * n + y) * n + x);
}

// Cheap integer hash (Murmur-style finalizer) used only to re-randomise the
// phase of an orbit whose calculated intensity has collapsed to zero. Doing
// this on the GPU avoids shipping a random buffer across, which would have
// cost a ninth storage binding.
fn hash_u32(x : u32) -> u32 {
    var h = x;
    h ^= h >> 16u;  h *= 0x7feb352du;
    h ^= h >> 15u;  h *= 0x846ca68bu;
    h ^= h >> 16u;
    return h;
}
fn rand01(x : u32) -> f32 {
    return f32(hash_u32(x) & 0x00ffffffu) / 16777216.0;
}

// ===========================================================================
//  FFT: batched 1D Stockham radix-2, three passes (x, then y, then z).
//  Reads `src`, writes `dst`. The host swaps the two through alternating
//  bind groups, so successive stages ping-pong with no copy in between.
// ===========================================================================
@compute @workgroup_size(64)
fn fft_stage(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let n = P.n;
    let half = n >> 1u;
    let total = half * n * n;
    let gid = gid3.x;
    if (gid >= total) { return; }

    let L = P.stage_len;
    let j    = gid % half;
    let line = gid / half;

    let k   = j % L;
    let blk = j / L;
    let in0 = blk * L + k;
    let in1 = in0 + half;
    let out0 = blk * (2u * L) + k;
    let out1 = out0 + L;

    var ang = -PI * f32(k) / f32(L);
    if (P.inverse == 1u) { ang = -ang; }
    let wr = cos(ang);
    let wi = sin(ang);

    let axis = P.axis;
    var ax : u32 = 0u; var ay : u32 = 0u; var az : u32 = 0u;
    if (axis == 0u)      { ay = line % n; az = line / n; }
    else if (axis == 1u) { ax = line % n; az = line / n; }
    else                 { ax = line % n; ay = line / n; }

    var iO0 : u32; var iO1 : u32; var oO0 : u32; var oO1 : u32;
    if (axis == 0u) {
        iO0 = voxel_offset(in0, ay, az, n);  iO1 = voxel_offset(in1, ay, az, n);
        oO0 = voxel_offset(out0, ay, az, n); oO1 = voxel_offset(out1, ay, az, n);
    } else if (axis == 1u) {
        iO0 = voxel_offset(ax, in0, az, n);  iO1 = voxel_offset(ax, in1, az, n);
        oO0 = voxel_offset(ax, out0, az, n); oO1 = voxel_offset(ax, out1, az, n);
    } else {
        iO0 = voxel_offset(ax, ay, in0, n);  iO1 = voxel_offset(ax, ay, in1, n);
        oO0 = voxel_offset(ax, ay, out0, n); oO1 = voxel_offset(ax, ay, out1, n);
    }

    let ar = src[iO0];  let ai = src[iO0 + 1u];
    let br = src[iO1];  let bi = src[iO1 + 1u];
    let tr = wr * br - wi * bi;
    let ti = wr * bi + wi * br;
    dst[oO0] = ar + tr;   dst[oO0 + 1u] = ai + ti;
    dst[oO1] = ar - tr;   dst[oO1 + 1u] = ai - ti;
}

// Straight copy src -> dst. Used once per transform, only when the stage
// count is odd and the result has landed in the wrong buffer.
//
// FOUR FLOATS PER INVOCATION, deliberately. One float per invocation needs
// 2*N^3/64 workgroups, which at N = 128 is 65536 -- one over the WebGPU limit
// of 65535 workgroups per dimension. That single dispatch was enough to make
// the whole command buffer invalid, and because WebGPU reports validation
// errors asynchronously the run did not throw: it just produced a grid that
// had never been written and an R factor of infinity. Every other kernel here
// dispatches at most N^3/64 = 32768 groups, so this was the only one over the
// line, but the host now checks all of them against the device limit too.
const COPY_PER_THREAD : u32 = 4u;

@compute @workgroup_size(64)
fn copy_buf(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let total = 2u * P.n * P.n * P.n;
    let base = gid3.x * COPY_PER_THREAD;
    if (base >= total) { return; }
    for (var k = 0u; k < COPY_PER_THREAD; k = k + 1u) {
        let i = base + k;
        if (i < total) { dst[i] = src[i]; }
    }
}

@compute @workgroup_size(64)
fn scale_inverse(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let total = P.n * P.n * P.n;
    let i = gid3.x;
    if (i >= total) { return; }
    let s = 1.0 / f32(total);
    src[2u * i]      = src[2u * i]      * s;
    src[2u * i + 1u] = src[2u * i + 1u] * s;
}

// ===========================================================================
//  Flip threshold: mean and mean-of-squares of the real part, reduced in two
//  stages so the host never has to see sigma.
// ===========================================================================
var<workgroup> wg_a : array<f32, 64>;
var<workgroup> wg_b : array<f32, 64>;

fn wg_reduce(lid : u32) {
    var s : u32 = 32u;
    loop {
        if (s == 0u) { break; }
        if (lid < s) {
            wg_a[lid] = wg_a[lid] + wg_a[lid + s];
            wg_b[lid] = wg_b[lid] + wg_b[lid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
}

@compute @workgroup_size(64)
fn reduce_moments(@builtin(global_invocation_id) gid3 : vec3<u32>,
                  @builtin(local_invocation_id)  lid3 : vec3<u32>,
                  @builtin(workgroup_id)         wid3 : vec3<u32>) {
    let total = P.n * P.n * P.n;
    let i = gid3.x;
    let lid = lid3.x;

    var v : f32 = 0.0;
    if (i < total) { v = src[2u * i]; }
    wg_a[lid] = v;
    wg_b[lid] = v * v;
    workgroupBarrier();
    wg_reduce(lid);
    if (lid == 0u) {
        partials[P_BASE + 2u * wid3.x]      = wg_a[0];
        partials[P_BASE + 2u * wid3.x + 1u] = wg_b[0];
    }
}

@compute @workgroup_size(64)
fn finalize_moments(@builtin(local_invocation_id) lid3 : vec3<u32>) {
    let lid = lid3.x;
    let total = f32(P.n * P.n * P.n);

    // The loop bound depends only on `base`, which is the same in every lane,
    // so the barrier inside wg_reduce is provably in uniform control flow.
    var sa : f32 = 0.0;
    var sb : f32 = 0.0;
    var base : u32 = 0u;
    loop {
        if (base >= P.wg_full) { break; }
        let w = base + lid;
        if (w < P.wg_full) {
            sa += partials[P_BASE + 2u * w];
            sb += partials[P_BASE + 2u * w + 1u];
        }
        base += 64u;
    }
    wg_a[lid] = sa;
    wg_b[lid] = sb;
    workgroupBarrier();
    wg_reduce(lid);

    if (lid == 0u) {
        let mean = wg_a[0] / total;
        let varr = max(0.0, wg_b[0] / total - mean * mean);
        partials[P_DELTA] = P.delta_sigma * sqrt(varr);
    }
}

@compute @workgroup_size(64)
fn flip(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let total = P.n * P.n * P.n;
    let i = gid3.x;
    if (i >= total) { return; }
    let delta = partials[P_DELTA];
    var re = src[2u * i];
    if (re < delta) { re = -re; }
    src[2u * i]      = re;
    src[2u * i + 1u] = 0.0;    // the density is real by construction
}

// ===========================================================================
//  Space-group symmetrisation.
//
//  For a real-space operator x' = xR + t the structure factors of a correctly
//  positioned structure obey
//
//        F(hR) = F(h) * exp(-2*pi*i * h.t)
//
//  Each orbit stores, for every member, the phase factor
//
//        Ebar = exp(+2*pi*i * h0.t)
//
//  where h0 is the orbit's representative. Reducing a member to h0 is a
//  multiply by Ebar (after conjugating if the member holds the Friedel mate),
//  and restoring is the reverse. Averaging the reduced values and pushing
//  them back is a projection onto the space group; sym_lambda damps it, so
//  0 reproduces the classical P1 charge flipping of Oszlanyi & Suto and 1
//  imposes the group exactly.
//
//  Doing this every cycle also fixes the origin: the symmetrised solution
//  comes out in the standard setting instead of shifted by an unknown vector.
// ===========================================================================
@compute @workgroup_size(64)
fn symmetrize(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let g = gid3.x;
    if (g >= P.n_orbits) { return; }
    if (P.sym_lambda <= 0.0) { return; }

    let orb = orbits[g];
    if (orb.count < 2u) { return; }

    var sr : f32 = 0.0;
    var si : f32 = 0.0;
    for (var i = 0u; i < orb.count; i = i + 1u) {
        let p    = orb.start + i;
        let word = orbitIdx[p];
        let idx  = word & IDX_MASK;
        var a = src[2u * idx];
        var b = src[2u * idx + 1u];
        if ((word & CONJ_BIT) != 0u) { b = -b; }
        let cr = phase[2u * p];
        let ci = phase[2u * p + 1u];
        sr += a * cr - b * ci;
        si += a * ci + b * cr;
    }
    let inv = 1.0 / f32(orb.count);
    sr = sr * inv;
    si = si * inv;

    let lam = clamp(P.sym_lambda, 0.0, 1.0);
    for (var i = 0u; i < orb.count; i = i + 1u) {
        let p    = orb.start + i;
        let word = orbitIdx[p];
        let idx  = word & IDX_MASK;
        let cr = phase[2u * p];
        let ci = phase[2u * p + 1u];
        // multiply the orbit mean by conj(Ebar)
        var a = sr * cr + si * ci;
        var b = si * cr - sr * ci;
        if ((word & CONJ_BIT) != 0u) { b = -b; }
        src[2u * idx]      = mix(src[2u * idx],      a, lam);
        src[2u * idx + 1u] = mix(src[2u * idx + 1u], b, lam);
    }
}

// ===========================================================================
//  Powder intensity sharing.
//
//  Only the TOTAL intensity of a 2-theta cluster is observed. Within the
//  cluster the split between distinct orbits is unknown, so it is taken from
//  the current calculated intensities -- the Le Bail step. Symmetry
//  equivalents are NOT involved here: they live inside one orbit and are
//  handled by `constrain`.
// ===========================================================================
@compute @workgroup_size(64)
fn repartition(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let c = gid3.x;
    if (c >= P.n_clusters) { return; }
    let cl = clusters[c];
    if (cl.n_orbits == 0u) { return; }

    var total : f32 = 0.0;
    for (var o = 0u; o < cl.n_orbits; o = o + 1u) {
        let orb = orbits[cl.orbit_start + o];
        var s : f32 = 0.0;
        for (var i = 0u; i < orb.count; i = i + 1u) {
            let idx = orbitIdx[orb.start + i] & IDX_MASK;
            let a = src[2u * idx];
            let b = src[2u * idx + 1u];
            s += a * a + b * b;
        }
        orbitTarget[cl.orbit_start + o] = s;   // scratch: calculated intensity
        total += s;
    }

    if (total > 1e-20) {
        let k = cl.obs_i / total;
        for (var o = 0u; o < cl.n_orbits; o = o + 1u) {
            orbitTarget[cl.orbit_start + o] = orbitTarget[cl.orbit_start + o] * k;
        }
    } else {
        // Nothing calculated yet: share equally, weighted by orbit size so
        // every reflection in the cluster starts from the same |F|.
        var wsum : f32 = 0.0;
        for (var o = 0u; o < cl.n_orbits; o = o + 1u) {
            wsum += f32(orbits[cl.orbit_start + o].count);
        }
        if (wsum <= 0.0) { return; }
        for (var o = 0u; o < cl.n_orbits; o = o + 1u) {
            orbitTarget[cl.orbit_start + o] =
                cl.obs_i * f32(orbits[cl.orbit_start + o].count) / wsum;
        }
    }
}

// ===========================================================================
//  Amplitude constraint.
//
//  Members of one orbit are symmetry equivalent, so their amplitudes are
//  EQUAL: |F| = sqrt(orbitTarget[g] / count). Only the phases survive from the
//  flipped density. v137 rescaled the orbit by a single factor, which
//  preserved whatever amplitude imbalance the flip had introduced -- i.e. it
//  let the map violate the point group. This does not.
// ===========================================================================
@compute @workgroup_size(64)
fn constrain(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let g = gid3.x;
    if (g >= P.n_orbits) { return; }
    let orb = orbits[g];
    if (orb.count == 0u) { return; }

    let tgt = max(0.0, orbitTarget[g]);
    let amp = sqrt(tgt / f32(orb.count));

    for (var i = 0u; i < orb.count; i = i + 1u) {
        let p   = orb.start + i;
        let idx = orbitIdx[p] & IDX_MASK;
        let a = src[2u * idx];
        let b = src[2u * idx + 1u];
        let m = sqrt(a * a + b * b);
        if (m > 1e-20) {
            let s = amp / m;
            src[2u * idx]      = a * s;
            src[2u * idx + 1u] = b * s;
        } else {
            // Dead reflection: restart it on a random phase rather than
            // parking it on the real axis, which biases the map towards a
            // centrosymmetric solution.
            let ang = 2.0 * PI * rand01(P.seed ^ (p * 2654435761u));
            src[2u * idx]      = amp * cos(ang);
            src[2u * idx + 1u] = amp * sin(ang);
        }
    }
}

// Systematically absent reflections are a hard zero. Leaving them free lets
// the density break the lattice centring, which for an F or I lattice means
// most of the reciprocal grid is unconstrained.
@compute @workgroup_size(64)
fn zero_absent(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let i = gid3.x;
    if (i >= P.absent_count) { return; }
    let idx = orbitIdx[P.absent_start + i] & IDX_MASK;
    src[2u * idx]      = 0.0;
    src[2u * idx + 1u] = 0.0;
}

@compute @workgroup_size(1)
fn zero_dc(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    if (gid3.x != 0u) { return; }
    src[0] = 0.0;
    src[1] = 0.0;
}

// ===========================================================================
//  R factor, on the observable quantity: cluster totals, not individual
//  reflections. A powder pattern cannot distinguish the members of a cluster,
//  so scoring them separately would be scoring noise.
// ===========================================================================
@compute @workgroup_size(64)
fn reduce_rfactor(@builtin(global_invocation_id) gid3 : vec3<u32>,
                  @builtin(local_invocation_id)  lid3 : vec3<u32>,
                  @builtin(workgroup_id)         wid3 : vec3<u32>) {
    let c = gid3.x;
    let lid = lid3.x;
    var d : f32 = 0.0;
    var t : f32 = 0.0;

    if (c < P.n_clusters) {
        let cl = clusters[c];
        var total : f32 = 0.0;
        for (var o = 0u; o < cl.n_orbits; o = o + 1u) {
            let orb = orbits[cl.orbit_start + o];
            for (var i = 0u; i < orb.count; i = i + 1u) {
                let idx = orbitIdx[orb.start + i] & IDX_MASK;
                let a = src[2u * idx];
                let b = src[2u * idx + 1u];
                total += a * a + b * b;
            }
        }
        t = sqrt(max(0.0, cl.obs_i));
        d = abs(sqrt(max(0.0, total)) - t);
    }

    wg_a[lid] = d;
    wg_b[lid] = t;
    workgroupBarrier();
    wg_reduce(lid);
    if (lid == 0u) {
        partials[P_BASE + 2u * wid3.x]      = wg_a[0];
        partials[P_BASE + 2u * wid3.x + 1u] = wg_b[0];
    }
}

// ===========================================================================
//  Best-iterate snapshot, on the device.
//
//  keep_best copies src -> best when the R factor just written by
//  finalize_rfactor beats the record in partials[P_BESTR]; commit_best then
//  updates the record. They are separate dispatches because every keep_best
//  invocation must see the OLD record -- a single kernel that both compared
//  and updated would race with itself across workgroups.
//
//  Cost is one extra pass over the grid per cycle, against 3*log2(N) passes
//  for each of the two transforms. That is a few percent, and it buys the
//  removal of a mapAsync stall per iteration.
// ===========================================================================
@compute @workgroup_size(64)
fn keep_best(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    let den = partials[P_RDEN];
    if (!(den > 0.0)) { return; }
    let r = partials[P_RNUM] / den;
    if (!(r < partials[P_BESTR])) { return; }

    let total = 2u * P.n * P.n * P.n;
    let base = gid3.x * COPY_PER_THREAD;
    if (base >= total) { return; }
    for (var k = 0u; k < COPY_PER_THREAD; k = k + 1u) {
        let i = base + k;
        if (i < total) { best[i] = src[i]; }
    }
}

@compute @workgroup_size(1)
fn commit_best(@builtin(global_invocation_id) gid3 : vec3<u32>) {
    if (gid3.x != 0u) { return; }
    let den = partials[P_RDEN];
    if (!(den > 0.0)) { return; }
    let r = partials[P_RNUM] / den;
    if (r < partials[P_BESTR]) { partials[P_BESTR] = r; }
}

@compute @workgroup_size(64)
fn finalize_rfactor(@builtin(local_invocation_id) lid3 : vec3<u32>) {
    let lid = lid3.x;
    var sa : f32 = 0.0;
    var sb : f32 = 0.0;
    var base : u32 = 0u;
    loop {
        if (base >= P.wg_cluster) { break; }
        let w = base + lid;
        if (w < P.wg_cluster) {
            sa += partials[P_BASE + 2u * w];
            sb += partials[P_BASE + 2u * w + 1u];
        }
        base += 64u;
    }
    wg_a[lid] = sa;
    wg_b[lid] = sb;
    workgroupBarrier();
    wg_reduce(lid);
    if (lid == 0u) {
        partials[P_RNUM] = wg_a[0];
        partials[P_RDEN] = wg_b[0];
    }
}
