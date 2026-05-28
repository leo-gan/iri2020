
pub fn dregion(z: f32, it: i32, f: f32, vkp: f32, mut f5sw: f32, mut f6wa: f32, elg: &mut [f32; 7]) {
    let a0 = [1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 3.0];
    let a1 = [0.6, 0.8, 1.1, 1.2, 1.3, 1.4, 1.0];
    let a2 = [0.0, 0.0, 0.08, 0.12, 0.05, 0.2, 0.0];
    let a3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    let a4 = [0.0, 0.0, -0.30, 0.10, 0.20, 0.30, 0.15];
    let a5 = [0.0, -0.10, -0.20, -0.25, -0.30, -0.30, 0.0];
    let a6 = [0.0, 0.1, 0.3, 0.6, 1.0, 1.0, 0.7];
    
    let pi = std::f32::consts::PI;
    let f1z = if z <= 45.0 {
        1.0
    } else if z < 90.0 {
        1.1892 * ((z * pi / 180.0).cos()).sqrt()
    } else {
        0.0
    };
    
    let mut f4s = 1.0;
    if it >= 5 && it <= 9 {
        f4s = 0.0;
        f5sw = 0.0;
        f6wa = 0.0;
    }
    if it == 3 || it == 4 || it == 10 || it == 11 {
        f4s = 0.5;
        f5sw = 0.0;
        f6wa = 0.0;
    }
    
    let f2kp = if vkp > 2.0 { 2.0 } else { vkp };
    let f3f = (f - 60.0) / 300.0 * f1z;
    
    for i in 0..7 {
        elg[i] = a0[i] + a1[i] * f1z + a2[i] * f2kp + a3[i] * f3f + a4[i] * f4s + a5[i] * f5sw + a6[i] * f6wa;
    }
}
