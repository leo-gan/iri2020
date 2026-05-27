use crate::igrf::IgrfModel;

pub const ARGMAX: f32 = 87.3;

// --- Trigonometric Wrappers ---

pub fn acos(x: f32) -> f32 {
    x.clamp(-1.0, 1.0).acos()
}

pub fn asin(x: f32) -> f32 {
    x.clamp(-1.0, 1.0).asin()
}

pub fn atan2d(y: f32, x: f32) -> f32 {
    y.atan2(x).to_degrees()
}

// --- Epstein Transition / Profile Step Functions ---

pub fn eptr(x: f32, sc: f32, hx: f32) -> f32 {
    let d1 = (x - hx) / sc;
    if d1.abs() >= ARGMAX {
        if d1 > 0.0 {
            d1
        } else {
            0.0
        }
    } else {
        (1.0 + d1.exp()).ln()
    }
}

pub fn epst(x: f32, sc: f32, hx: f32) -> f32 {
    let d1 = (x - hx) / sc;
    if d1.abs() >= ARGMAX {
        if d1 > 0.0 {
            1.0
        } else {
            0.0
        }
    } else {
        1.0 / (1.0 + (-d1).exp())
    }
}

pub fn epstep(y2: f32, y1: f32, sc: f32, hx: f32, x: f32) -> f32 {
    y1 + (y2 - y1) * epst(x, sc, hx)
}

pub fn epla(x: f32, sc: f32, hx: f32) -> f32 {
    let d1 = (x - hx) / sc;
    if d1.abs() >= ARGMAX {
        0.0
    } else {
        let d0 = d1.exp();
        let d2 = 1.0 + d0;
        d0 / (d2 * d2)
    }
}

pub fn rlay(x: f32, xm: f32, sc: f32, hx: f32) -> f32 {
    let y1 = eptr(x, sc, hx);
    let y1m = eptr(xm, sc, hx);
    let y2m = epst(xm, sc, hx);
    y1 - y1m - (x - xm) * y2m / sc
}

pub fn d1lay(x: f32, xm: f32, sc: f32, hx: f32) -> f32 {
    (epst(x, sc, hx) - epst(xm, sc, hx)) / sc
}

pub fn d2lay(x: f32, _xm: f32, sc: f32, hx: f32) -> f32 {
    epla(x, sc, hx) / (sc * sc)
}

pub fn booker(h: f32, n: usize, ah: &[f32], av: &[f32], d: &[f32]) -> f32 {
    let mut st = vec![0.0; n];
    st[0] = (av[1] - av[0]) / (ah[1] - ah[0]);
    let mut sum = av[0] + st[0] * (h - ah[0]);
    for idx in 0..(n - 2) {
        let aa = eptr(h, d[idx], ah[idx + 1]);
        let bb = eptr(ah[idx], d[idx], ah[idx + 1]);
        st[idx + 1] = (av[idx + 2] - av[idx + 1]) / (ah[idx + 2] - ah[idx + 1]);
        sum += (st[idx + 1] - st[idx]) * (aa - bb) * d[idx];
    }
    sum
}

pub fn booker1(h: f32, m: usize, f1: f32, ah: &[f32], st: &[f32], d: &[f32]) -> f32 {
    let mut sum = f1 + st[0] * (h - ah[0]);
    for idx in 0..m {
        let aa = eptr(h, d[idx], ah[idx + 1]);
        let bb = eptr(ah[0], d[idx], ah[idx + 1]);
        sum += (st[idx + 1] - st[idx]) * (aa - bb) * d[idx];
    }
    sum
}

// --- Solvers & Root Finders ---

pub fn regfa1<F>(
    x11: f32,
    x22: f32,
    fx11: f32,
    fx22: f32,
    eps: f32,
    fw: f32,
    f: F,
) -> (bool, f32)
where
    F: Fn(f32) -> f32,
{
    let mut ep = eps;
    let mut x1 = x11;
    let mut x2 = x22;
    let mut f1 = fx11 - fw;
    let mut f2 = fx22 - fw;
    let mut k = false;
    let mut ng = 2;
    let mut lfd = 0;
    let mut x: f32;

    if f1 * f2 > 0.0 {
        return (true, 0.0);
    }

    let mut l1 = false;
    x = (x1 * f2 - x2 * f1) / (f2 - f1);

    loop {
        let fx = f(x) - fw;
        lfd += 1;
        if lfd > 20 {
            ep *= 10.0;
            lfd = 0;
        }
        let links = f1 * fx > 0.0;
        k = !k;
        if links {
            x1 = x;
            f1 = fx;
        } else {
            x2 = x;
            f2 = fx;
        }

        if (x2 - x1).abs() <= ep {
            break;
        }

        if k {
            l1 = links;
            let mut dx = (x2 - x1) / (ng as f32);
            if !links {
                dx *= (ng - 1) as f32;
            }
            x = x1 + dx;
        } else {
            if (links && !l1) || (!links && l1) {
                ng *= 2;
            }
            x = (x1 * f2 - x2 * f1) / (f2 - f1);
        }
    }

    (false, x)
}

pub fn lnglsn(n: usize, a: &mut [[f32; 5]; 5], b: &mut [f32; 5]) -> bool {
    let mut aus = false;
    let mut azv = [0.0_f32; 10];

    for k in 0..(n - 1) {
        let imax = k;
        let mut l = k;
        let mut izg = 0;
        let amax = a[k][k].abs();

        loop {
            l += 1;
            if l >= n {
                break;
            }
            let hsp = a[l][k].abs();
            if hsp < 1e-8 {
                izg += 1;
            }
            if hsp <= amax {
                continue;
            }
            break;
        }

        if amax >= 1e-10 {
            if imax != k {
                for col in k..n {
                    azv[col + 1] = a[imax][col];
                    a[imax][col] = a[k][col];
                    a[k][col] = azv[col + 1];
                }
                azv[0] = b[imax];
                b[imax] = b[k];
                b[k] = azv[0];
            }
            if izg == (n - (k + 1)) {
                continue;
            }
            let amax_inv = 1.0 / a[k][k];
            azv[0] = b[k] * amax_inv;
            for m in (k + 1)..n {
                azv[m + 1] = a[k][m] * amax_inv;
            }
            for row in (k + 1)..n {
                let row_amax = a[row][k];
                if row_amax.abs() < 1e-8 {
                    continue;
                }
                a[row][k] = 0.0;
                b[row] -= azv[0] * row_amax;
                for m in (k + 1)..n {
                    a[row][m] -= row_amax * azv[m + 1];
                }
            }
        } else {
            aus = true;
            return aus;
        }
    }

    for k_idx in (0..n).rev() {
        let mut amax = 0.0_f32;
        if k_idx < n - 1 {
            for l in (k_idx + 1)..n {
                amax += a[k_idx][l] * a[n - 1][l];
            }
        }
        if a[k_idx][k_idx].abs() < 1e-6 {
            a[n - 1][k_idx] = 0.0;
        } else {
            a[n - 1][k_idx] = (b[k_idx] - amax) / a[k_idx][k_idx];
        }
    }

    aus
}

// --- Index Interpolation ---

pub fn moda(
    in_mode: i32,
    iyear: i32,
    month: &mut i32,
    iday: &mut i32,
    idoy: &mut i32,
    nrdaymo: &mut i32,
) {
    let mut mm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    if iyear % 4 == 0 {
        mm[1] = 29;
    }

    if in_mode <= 0 {
        let m = (*month as usize).clamp(1, 12);
        let mut mosum = 0;
        if m > 1 {
            for i in 0..(m - 1) {
                mosum += mm[i];
            }
        }
        *idoy = mosum + *iday;
        *nrdaymo = mm[m - 1];
    } else {
        let mut imo = 0;
        let mut mobe = 0;
        let mut moold = 0;
        loop {
            imo += 1;
            if imo > 12 {
                break;
            }
            moold = mobe;
            *nrdaymo = mm[imo - 1];
            mobe += *nrdaymo;
            if mobe >= *idoy {
                break;
            }
        }
        *month = imo as i32;
        *iday = *idoy - moold;
    }
}

pub fn tcon(
    yr: i32,
    mm: i32,
    day: i32,
    idn: i32,
    ionoindx: &[f32; 806],
    indrz: &[f32; 806],
    iymst: i32,
    iymend: i32,
    rz: &mut [f32; 3],
    ig: &mut [f32; 3],
    rsn: &mut f32,
    nmonth: &mut i32,
    mess: bool,
) -> Result<(), String> {
    let iytmp = yr * 100 + mm;
    if iytmp < iymst || iytmp > iymend {
        if mess {
            eprintln!(
                "{} ** OUT OF RANGE **\nThe file IG_RZ.DAT which contains the indices Rz12 and IG12 currently only covers the time period (yymm) : {}-{}",
                iytmp, iymst, iymend
            );
            *nmonth = -1;
            return Err(format!("IG_RZ.DAT out of date index range: {}", iytmp));
        }
    }

    let iyst = iymst / 100;
    let imst = iymst - iyst * 100;
    let num = 2 - imst + (yr - iyst) * 12 + mm;

    let num_idx = (num - 1).clamp(0, 805) as usize;
    rz[0] = indrz[num_idx];
    ig[0] = ionoindx[num_idx];

    let mut midm = 15;
    if mm == 2 {
        midm = 14;
    }

    let mut idd1 = 0;
    let mut nrdaym = 0;
    let mut mm_mut = mm;
    let mut midm_mut = midm;
    moda(0, yr, &mut mm_mut, &mut midm_mut, &mut idd1, &mut nrdaym);

    let mut imm2: i32;
    if day < midm {
        imm2 = mm - 1;
        let mut idd2 = 0;
        if imm2 < 1 {
            imm2 = 12;
            idd2 = -16;
        } else {
            let iyy2 = yr;
            let mut midm2 = 15;
            if imm2 == 2 {
                midm2 = 14;
            }
            let mut imm2_mut = imm2;
            let mut midm2_mut = midm2;
            moda(0, iyy2, &mut imm2_mut, &mut midm2_mut, &mut idd2, &mut nrdaym);
        }

        let num_prev_idx = (num - 2).clamp(0, 805) as usize;
        rz[1] = indrz[num_prev_idx];
        ig[1] = ionoindx[num_prev_idx];
        *rsn = (idn - idd2) as f32 / (idd1 - idd2) as f32;
        rz[2] = rz[1] + (rz[0] - rz[1]) * *rsn;
        ig[2] = ig[1] + (ig[0] - ig[1]) * *rsn;
    } else {
        imm2 = mm + 1;
        let mut idd2 = 0;
        if imm2 > 12 {
            let iyy2 = yr + 1;
            idd2 = 380;
            if iyy2 % 4 == 0 {
                idd2 = 381;
            }
        } else {
            let iyy2 = yr;
            let mut midm2 = 15;
            if imm2 == 2 {
                midm2 = 14;
            }
            let mut imm2_mut = imm2;
            let mut midm2_mut = midm2;
            moda(0, iyy2, &mut imm2_mut, &mut midm2_mut, &mut idd2, &mut nrdaym);
        }

        let num_next_idx = num.clamp(0, 805) as usize;
        rz[1] = indrz[num_next_idx];
        ig[1] = ionoindx[num_next_idx];
        *rsn = (idn - idd1) as f32 / (idd2 - idd1) as f32;
        rz[2] = rz[0] + (rz[1] - rz[0]) * *rsn;
        ig[2] = ig[0] + (ig[1] - ig[0]) * *rsn;
    }

    *nmonth = imm2;
    Ok(())
}

// --- Field Line Tracing & Invariant L-Shell computation (FINDB0/SHELLG/STOER) ---

pub const PI: f32 = 3.141592653589793;
pub const UMR: f32 = PI / 180.0;
pub const ERA: f32 = 6371.2;
pub const EREQU: f32 = 6378.16;
pub const ERPOL: f32 = 6356.775;
pub const AQUAD: f32 = EREQU * EREQU;
pub const BQUAD: f32 = ERPOL * ERPOL;

const U: [[f32; 3]; 3] = [
    [0.3511737, -0.9148385, -0.1993679], // Column 1
    [0.9335804, 0.3583680, 0.0000000],  // Column 2
    [0.0714471, -0.1861260, 0.9799247], // Column 3
];

pub fn feldi(model: &IgrfModel, xi: &[f32; 4]) -> [f32; 197] {
    let mut h = [0.0_f32; 197];
    let ihmax = (model.nmax * model.nmax + 1) as isize;
    let last = ihmax + (model.nmax + model.nmax) as isize;
    let imax = (model.nmax + model.nmax - 1) as isize;

    for idx in ihmax..=last {
        h[idx as usize] = model.g[idx as usize];
    }

    for k in (1..=3).step_by(2) {
        let k = k as isize;
        let mut i = imax;
        let mut ih = ihmax;
        loop {
            let il = ih - i;
            let f = 2.0 / (i - k + 2) as f32;
            let x = xi[1] * f;
            let y = xi[2] * f;
            let z = xi[3] * (f + f);
            i -= 2;

            let il_u = il as usize;
            let ih_u = ih as usize;

            if i < 1 {
                h[il_u] = model.g[il_u] + z * h[ih_u] + 2.0 * (x * h[ih_u + 1] + y * h[ih_u + 2]);
            } else if i == 1 {
                h[il_u + 2] = model.g[il_u + 2] + z * h[ih_u + 2] + x * h[ih_u + 4] - y * (h[ih_u + 3] + h[ih_u]);
                h[il_u + 1] = model.g[il_u + 1] + z * h[ih_u + 1] + y * h[ih_u + 4] + x * (h[ih_u + 3] - h[ih_u]);
                h[il_u] = model.g[il_u] + z * h[ih_u] + 2.0 * (x * h[ih_u + 1] + y * h[ih_u + 2]);
            } else {
                let mut m = 3_isize;
                while m <= i {
                    let m_u = m as usize;
                    h[il_u + m_u + 1] = model.g[il_u + m_u + 1] + z * h[ih_u + m_u + 1] + x * (h[ih_u + m_u + 3] - h[ih_u + m_u - 1]) - y * (h[ih_u + m_u + 2] + h[ih_u + m_u - 2]);
                    h[il_u + m_u] = model.g[il_u + m_u] + z * h[ih_u + m_u] + x * (h[ih_u + m_u + 2] - h[ih_u + m_u - 2]) + y * (h[ih_u + m_u + 3] + h[ih_u + m_u - 1]);
                    m += 2;
                }
                h[il_u + 2] = model.g[il_u + 2] + z * h[ih_u + 2] + x * h[ih_u + 4] - y * (h[ih_u + 3] + h[ih_u]);
                h[il_u + 1] = model.g[il_u + 1] + z * h[ih_u + 1] + y * h[ih_u + 4] + x * (h[ih_u + 3] - h[ih_u]);
                h[il_u] = model.g[il_u] + z * h[ih_u] + 2.0 * (x * h[ih_u + 1] + y * h[ih_u + 2]);
            }

            ih = il;
            if i < k {
                break;
            }
        }
    }

    h
}

pub fn stoer(model: &IgrfModel, p: &mut [f32], bq: &mut f32, r: &mut f32) {
    let zm = p[2];
    let fli = p[0] * p[0] + p[1] * p[1] + 1e-15;
    *r = 0.5 * (fli + (fli * fli + (zm + zm) * (zm + zm)).sqrt());
    let rq = *r * *r;
    let wr = r.sqrt();
    let xm = p[0] * wr;
    let ym = p[1] * wr;

    let mut xi = [0.0_f32; 4];
    xi[1] = xm * U[0][0] + ym * U[1][0] + zm * U[2][0];
    xi[2] = xm * U[0][1] + ym * U[1][1] + zm * U[2][1];
    xi[3] = xm * U[0][2] + zm * U[2][2];

    let h = feldi(model, &xi);
    let q = h[1] / rq;
    let dx = h[3] + h[3] + q * xi[1];
    let dy = h[4] + h[4] + q * xi[2];
    let dz = h[2] + h[2] + q * xi[3];

    let dxm = U[0][0] * dx + U[0][1] * dy + U[0][2] * dz;
    let dym = U[1][0] * dx + U[1][1] * dy;
    let dzm = U[2][0] * dx + U[2][1] * dy + U[2][2] * dz;

    let dr = (xm * dxm + ym * dym + zm * dzm) / *r;
    p[3] = (wr * dxm - 0.5 * p[0] * dr) / (*r * dzm);
    p[4] = (wr * dym - 0.5 * p[1] * dr) / (*r * dzm);

    let dsq = rq * (dxm * dxm + dym * dym + dzm * dzm);
    *bq = dsq * rq * rq;
    p[5] = (dsq / (rq + 3.0 * zm * zm)).sqrt();
    p[6] = p[5] * (rq + zm * zm) / (rq * dzm);
}

pub fn shellg(model: &IgrfModel, glat: f32, glon: f32, alt: f32) -> (f32, i32, f32) {
    let rmin = 0.05_f32;
    let rmax = 1.01_f32;
    let mut step = 0.20_f32;
    let steq = 0.03_f32;
    let mut bequ = 1e10_f32;
    let mut iequ = 0;

    let rlat = glat * UMR;
    let ct = rlat.sin();
    let st = rlat.cos();
    let d_denom = (AQUAD - (AQUAD - BQUAD) * ct * ct).sqrt();
    let mut x = [0.0_f32; 4];
    x[1] = (alt + AQUAD / d_denom) * st / ERA;
    x[3] = (alt + BQUAD / d_denom) * ct / ERA;
    let rlon = glon * UMR;
    x[2] = x[1] * rlon.sin();
    x[1] = x[1] * rlon.cos();

    let rq = 1.0 / (x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
    let r3h = (rq * rq.sqrt()).sqrt();

    let mut p = [[0.0_f32; 8]; 3335];
    p[1][0] = (x[1] * U[0][0] + x[2] * U[0][1] + x[3] * U[0][2]) * r3h;
    p[1][1] = (x[1] * U[1][0] + x[2] * U[1][1]) * r3h;
    p[1][2] = (x[1] * U[2][0] + x[2] * U[2][1] + x[3] * U[2][2]) * rq;

    step = -step.copysign(p[1][2]);

    let mut bq2 = 0.0_f32;
    let mut r2 = 0.0_f32;
    stoer(model, &mut p[1], &mut bq2, &mut r2);
    let b0 = bq2.sqrt();

    p[2][0] = p[1][0] + 0.5 * step * p[1][3];
    p[2][1] = p[1][1] + 0.5 * step * p[1][4];
    p[2][2] = p[1][2] + 0.5 * step;

    let mut bq3 = 0.0_f32;
    let mut r3 = 0.0_f32;
    stoer(model, &mut p[2], &mut bq3, &mut r3);

    p[0][0] = p[1][0] - step * (2.0 * p[1][3] - p[2][3]);
    p[0][1] = p[1][1] - step * (2.0 * p[1][4] - p[2][4]);
    p[0][2] = p[1][2] - step;

    let mut bq1 = 0.0_f32;
    let mut r1 = 0.0_f32;
    stoer(model, &mut p[0], &mut bq1, &mut r1);

    p[2][0] = p[1][0] + step * (20.0 * p[2][3] - 3.0 * p[1][3] + p[0][3]) / 18.0;
    p[2][1] = p[1][1] + step * (20.0 * p[2][4] - 3.0 * p[1][4] + p[0][4]) / 18.0;
    p[2][2] = p[1][2] + step;

    stoer(model, &mut p[2], &mut bq3, &mut r3);

    if bq3 > bq1 {
        step = -step;
        r3 = r1;
        bq3 = bq1;
        for i in 0..7 {
            let zz = p[0][i];
            p[0][i] = p[2][i];
            p[2][i] = zz;
        }
    }

    if bq1 < bequ {
        bequ = bq1;
        iequ = 1;
    }
    if bq2 < bequ {
        bequ = bq2;
        iequ = 2;
    }
    if bq3 < bequ {
        bequ = bq3;
        iequ = 3;
    }

    let step12 = step / 12.0;
    let step2 = step + step;
    let steq = steq.copysign(step);
    let mut fi = 0.0_f32;
    let mut icode = 1;
    let mut oradik = 0.0_f32;
    let mut oterm = 0.0_f32;
    let mut stp = r2 * steq;
    let mut z = p[1][2] + stp;
    stp /= 0.75;
    p[0][7] = step2 * (p[0][0] * p[0][3] + p[0][1] * p[0][4]);
    p[1][7] = step2 * (p[1][0] * p[1][3] + p[1][1] * p[1][4]);

    let mut n_broke = 2;
    let mut exit_reason = 0; // 0: completed, 10: broke to 10, 30: broke to 30

    let mut radik = 0.0_f32;
    let mut c0 = 0.0_f32;
    let mut c1 = 0.0_f32;
    let mut c2 = 0.0_f32;
    let mut c3 = 0.0_f32;

    for n in 2..=3332 {
        n_broke = n;
        p[n][0] = p[n - 1][0] + step12 * (5.0 * p[n][3] + 8.0 * p[n - 1][3] - p[n - 2][3]);
        p[n][1] = p[n - 1][1] + step12 * (5.0 * p[n][4] + 8.0 * p[n - 1][4] - p[n - 2][4]);
        p[n][7] = step2 * (p[n][0] * p[n][3] + p[n][1] * p[n][4]);

        c0 = p[n - 1][0] * p[n - 1][0] + p[n - 1][1] * p[n - 1][1];
        c1 = p[n - 1][7];
        c2 = (p[n][7] - p[n - 2][7]) * 0.25;
        c3 = (p[n][7] + p[n - 2][7] - c1 - c1) / 6.0;

        let d0 = p[n - 1][5];
        let d1 = (p[n][5] - p[n - 2][5]) * 0.5;
        let d2 = (p[n][5] + p[n - 2][5] - d0 - d0) * 0.5;

        let e0 = p[n - 1][6];
        let e1 = (p[n][6] - p[n - 2][6]) * 0.5;
        let e2 = (p[n][6] + p[n - 2][6] - e0 - e0) * 0.5;

        let mut inner_iters = 0;
        loop {
            inner_iters += 1;
            if inner_iters > 10000 {
                exit_reason = 10;
                break;
            }
            let t = (z - p[n - 1][2]) / step;
            if t.is_nan() || t > 1.0 {
                break;
            }
            let hli = 0.5 * (((c3 * t + c2) * t + c1) * t + c0);
            let zq = z * z;
            let r_val = hli + (hli * hli + zq).sqrt();
            if r_val <= rmin {
                exit_reason = 30;
                break;
            }
            let rq_val = r_val * r_val;
            let ff = (1.0 + 3.0 * zq / rq_val).sqrt();
            radik = b0 - ((d2 * t + d1) * t + d0) * r_val * rq_val * ff;
            if r_val > rmax {
                icode = 2;
                radik -= 12.0 * (r_val - rmax) * (r_val - rmax);
            }
            if radik + radik <= oradik {
                exit_reason = 10;
                break;
            }
            let term = radik.sqrt() * ff * ((e2 * t + e1) * t + e0) / (rq_val + zq);
            fi += stp * (oterm + term);
            oradik = radik;
            oterm = term;
            stp = r_val * steq;
            z += stp;
        }

        if exit_reason != 0 {
            break;
        }

        p[n + 1][0] = p[n][0] + step12 * (23.0 * p[n][3] - 16.0 * p[n - 1][3] + 5.0 * p[n - 2][3]);
        p[n + 1][1] = p[n][1] + step12 * (23.0 * p[n][4] - 16.0 * p[n - 1][4] + 5.0 * p[n - 2][4]);
        p[n + 1][2] = p[n][2] + step;

        let mut bq3_loop = 0.0_f32;
        let mut r3_loop = 0.0_f32;
        stoer(model, &mut p[n + 1], &mut bq3_loop, &mut r3_loop);
        if bq3_loop < bequ {
            iequ = n + 1;
            bequ = bq3_loop;
        }
    }

    let fl: f32;

    if exit_reason == 30 {
        icode = 3;
        let t = -p[n_broke - 1][2] / step;
        fl = 1.0 / ((((c3 * t + c2) * t + c1) * t + c0).abs() + 1e-15);
    } else {
        let _iequ_bound = if iequ < 2 { 2 } else { iequ };
        if oradik >= 1e-15 {
            let denom = oradik - radik;
            if denom.abs() > 1e-15 {
                fi += stp / 0.75 * oterm * oradik / denom;
            }
        }
        fi = 0.5 * fi.abs() / b0.sqrt() + 1e-12;

        let dimob0 = model.dimo / b0;
        let arg1 = fi.ln();
        let arg2 = dimob0.ln();
        let xx = 3.0 * arg1 - arg2;
        let gg: f32;

        if xx > 23.0 {
            gg = xx - 3.0460681;
        } else if xx > 11.7 {
            gg = (((((2.8212095e-8 * xx - 3.8049276e-6) * xx + 2.170224e-4) * xx - 6.7310339e-3) * xx + 1.2038224e-1) * xx - 1.8461796e-1) * xx + 2.0007187;
        } else if xx > 3.0 {
            gg = ((((((((6.3271665e-10 * xx - 3.958306e-8) * xx + 9.9766148e-7) * xx - 1.2531932e-5) * xx + 7.9451313e-5) * xx - 3.2077032e-4) * xx + 2.1680398e-3) * xx + 1.2817956e-2) * xx + 4.3510529e-1) * xx + 6.222355e-1;
        } else if xx > -3.0 {
            gg = ((((((((2.6047023e-10 * xx + 2.3028767e-9) * xx - 2.1997983e-8) * xx - 5.3977642e-7) * xx - 3.3408822e-6) * xx + 3.8379917e-5) * xx + 1.1784234e-3) * xx + 1.4492441e-2) * xx + 4.3352788e-1) * xx + 6.228644e-1;
        } else if xx > -22.0 {
            gg = ((((((((-8.1537735e-14 * xx + 8.3232531e-13) * xx + 1.0066362e-9) * xx + 8.1048663e-8) * xx + 3.2916354e-6) * xx + 8.2711096e-5) * xx + 1.3714667e-3) * xx + 1.5017245e-2) * xx + 4.3432642e-1) * xx + 6.2337691e-1;
        } else {
            gg = 3.33338e-1 * xx + 3.0062102e-1;
        }

        fl = (((1.0 + gg.exp()) * dimob0).ln() / 3.0).exp();
    }

    (fl, icode, b0)
}

// --- Haines Scaling Function ---

pub fn tbfit(t1: f64, t2: f64, ibf: i32) -> (f64, f64) {
    let tzero = if ibf <= 1 { (t2 + t1) / 2.0 } else { t1 };
    let mut thint = t2 - t1;
    if ibf <= 2 {
        thint /= 2.0;
    }
    (thint, tzero)
}
