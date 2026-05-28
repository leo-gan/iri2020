use std::f32::consts::PI;

const UMR: f32 = PI / 180.0;

/// Calculates foE/MHz by the Edinburgh method.
pub fn foeedi(cov: f32, xhi: f32, mut xhim: f32, xlati: f32) -> f32 {
    let a = 1.0 + 0.0094 * (cov - 66.0);
    let sl = (xlati * UMR).cos();
    let (sm, c) = if xlati < 32.0 {
        (-1.93 + 1.92 * sl, 23.0 + 116.0 * sl)
    } else {
        (0.11 - 0.49 * sl, 92.0 + 35.0 * sl)
    };

    if xhim >= 90.0 {
        xhim = 89.999;
    }
    let b = (xhim * UMR).cos().powf(sm);
    let sp = if xlati > 12.0 { 1.2 } else { 1.31 };

    // adjusted solar zenith angle during nighttime (XHIC)
    let xhic = xhi - 3.0 * (1.0 + ((xhi - 89.98) / 3.0).exp()).ln();
    let d = (xhic * UMR).cos().powf(sp);
    let mut r4foe = a * b * c * d;

    let mut smin = 0.121 + 0.0015 * (cov - 60.0);
    smin = smin * smin;
    if r4foe < smin {
        r4foe = smin;
    }
    r4foe.powf(0.25)
}

/// Calculates electron density of D maximum.
pub fn xmded(xhi: f32, r: f32, yw: f32) -> f32 {
    if xhi >= 90.0 {
        return yw;
    }
    let y = 6.05e8 + 0.088e8 * r;
    let yy = (xhi * UMR).cos();
    let yyy = -0.1 / yy.powf(2.7);
    let ymd = if yyy < -40.0 {
        0.0
    } else {
        y * yyy.exp()
    };
    if ymd < yw {
        yw
    } else {
        ymd
    }
}

/// Calculates E-F valley parameters.
/// Returns `(hvb, vwu, vwa, vdp)`.
pub fn valgul(xhi: f32) -> (f32, f32, f32, f32) {
    let cs = 0.1 + (xhi * UMR).cos();
    let abc = cs.abs();
    let vdp = 0.45 * cs / (0.1 + abc) + 0.55;
    let arl = (0.1 + abc + cs) / (0.1 + abc - cs);
    let zzz = arl.ln();
    let vwu = 45.0 - 10.0 * zzz;
    let vwa = 45.0 - 5.0 * zzz;
    let hvb = 1000.0 / (7.024 + 0.224 * cs + 0.966 * abc);
    (hvb, vwu, vwa, vdp)
}
