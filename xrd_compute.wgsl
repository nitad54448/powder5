struct PeakPrep {
    intensity: f32,
    x0: f32,
    max_d: f32,
    pedestal: f32,
    
    p_type: f32,
    asym_param: f32,
    eta: f32,
    fwhm_total: f32,
    
    H_G_L: f32,
    H_L_L: f32,
    cL_L: f32,
    cG_L: f32,
    
    H_G_R: f32,
    H_L_R: f32,
    cL_R: f32,
    cG_R: f32,
};

@group(0) @binding(0) var<storage, read> tth_axis: array<f32>;
@group(0) @binding(1) var<storage, read> peaks: array<PeakPrep>;
@group(0) @binding(2) var<storage, read_write> y_calc: array<f32>;

const CG: f32 = 2.772588722239781; // 4 * ln(2)

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    if (idx >= arrayLength(&tth_axis)) { return; }
    
    let x = tth_axis[idx];
    var total_intensity: f32 = 0.0;
    
    for (var i: u32 = 0u; i < arrayLength(&peaks); i = i + 1u) {
        let p = peaks[i];
        var delta = x - p.x0;
        
        if (delta > p.max_d || delta < -p.max_d) { continue; }
        
        if (abs(p.asym_param) > 0.0 && abs(delta) >= 1e-9) {
            let t = clamp(p.asym_param * delta, -0.95, 0.95);
            delta = delta / (1.0 - t);
        }
        
        var hg: f32; var hl: f32; var cL: f32; var cG: f32;
        
        if (p.p_type > 1.5 && p.p_type < 2.5) { // split_pvoigt
            if (delta < 0.0) {
                hg = p.H_G_L; hl = p.H_L_L; cL = p.cL_L; cG = p.cG_L;
            } else {
                hg = p.H_G_R; hl = p.H_L_R; cL = p.cL_R; cG = p.cG_R;
            }
        } else { // simple_pvoigt and tch_aniso
            hg = p.H_G_L; hl = p.H_L_L; cL = p.cL_L; cG = p.cG_L;
        }
        
        let dg = (delta / hg) * (delta / hg);
        let dl = (delta / hl) * (delta / hl);
        let g = exp(-CG * dg);
        let l = 1.0 / (1.0 + 4.0 * dl);
        
        let v = cL * l + cG * g - p.pedestal;
        if (v > 0.0) {
            total_intensity += p.intensity * v;
        }
    }
    y_calc[idx] = total_intensity;
}