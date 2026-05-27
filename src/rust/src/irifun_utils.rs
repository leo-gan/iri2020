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

pub fn gamma1(
    smodip: f32,
    slat: f32,
    slong: f32,
    hour: f32,
    iharm: usize,
    nq: &[i32],
    k1: usize,
    m: usize,
    mm: usize,
    m3: usize,
    sfe: &[f32],
) -> f32 {
    let mut c = [0.0_f64; 13];
    let mut s = [0.0_f64; 13];
    let hou = ((15.0 * hour - 180.0) * UMR) as f64;
    s[1] = hou.sin();
    c[1] = hou.cos();
    for i in 2..=iharm {
        c[i] = c[1] * c[i - 1] - s[1] * s[i - 1];
        s[i] = c[1] * s[i - 1] + s[1] * c[i - 1];
    }
    
    let mut coef = vec![0.0_f64; m + 1];
    for i in 1..=m {
        let mi = (i - 1) * mm;
        coef[i] = sfe[mi] as f64;
        for j in 1..=iharm {
            coef[i] += (sfe[mi + 2 * j - 1] as f64) * s[j] + (sfe[mi + 2 * j] as f64) * c[j];
        }
    }
    
    let mut sum = coef[1];
    let ss = (smodip * UMR).sin() as f64;
    let s3 = ss;
    let mut xsinx = vec![0.0_f64; m + 1];
    xsinx[1] = 1.0;
    let index = nq[0] as usize;
    let mut ss_running = ss;
    for j in 1..=index {
        sum += coef[1 + j] * ss_running;
        xsinx[j + 1] = ss_running;
        ss_running = ss_running * s3;
    }
    
    let mut np = nq[0] as usize + 1;
    let ss_cos = (slat * UMR).cos() as f64;
    let s3_cos = ss_cos;
    let mut ss_cos_running = ss_cos;
    
    for j in 2..=k1 {
        let s0 = (slong * (j - 1) as f32 * UMR) as f64;
        let s1 = s0.cos();
        let s2 = s0.sin();
        let index = (nq[j - 1] + 1) as usize;
        for l in 1..=index {
            np += 1;
            sum += coef[np] * xsinx[l] * ss_cos_running * s1;
            np += 1;
            sum += coef[np] * xsinx[l] * ss_cos_running * s2;
        }
        ss_cos_running = ss_cos_running * s3_cos;
    }
    
    sum as f32
}

pub fn fout(xmodip: f32, xlati: f32, xlongi: f32, ut: f32, ff0: &[f32]) -> f32 {
    let qf = [11, 11, 8, 4, 1, 0, 0, 0, 0];
    gamma1(xmodip, xlati, xlongi, ut, 6, &qf, 9, 76, 13, 988, ff0)
}

pub fn xmout(xmodip: f32, xlati: f32, xlongi: f32, ut: f32, xm0: &[f32]) -> f32 {
    let qm = [6, 7, 5, 2, 1, 0, 0];
    gamma1(xmodip, xlati, xlongi, ut, 4, &qm, 7, 49, 9, 441, xm0)
}

pub fn conver(rga: f32, rgo: f32) -> f32 {
    use crate::cormag_data::CORMAG;
    let mut rlo = rgo;
    let rla = rga + 90.0;
    if rlo == 360.0 {
        rlo = 0.0;
    }
    
    let la1 = ((rla / 2.0) as i32 + 1).clamp(1, 91) as usize;
    let mut la2 = la1 + 1;
    if la2 > 91 {
        la2 = 91;
    }
    
    let lo1 = ((rlo / 18.0) as i32 + 1).clamp(1, 20) as usize;
    let lo2 = (lo1 % 20) + 1;
    
    let gm1 = CORMAG[la1 - 1][lo1 - 1];
    let gm2 = CORMAG[la2 - 1][lo1 - 1];
    let gm3 = CORMAG[la1 - 1][lo2 - 1];
    let gm4 = CORMAG[la2 - 1][lo2 - 1];
    
    let x = (rla / 2.0) - (rla / 2.0).floor();
    let y = (rlo / 18.0) - (rlo / 18.0).floor();
    
    let gmla = gm1 * (1.0 - x) * (1.0 - y)
        + gm2 * (1.0 - y) * x
        + gm3 * y * (1.0 - x)
        + gm4 * x * y;
        
    90.0 - gmla
}

const CODE: [[i32; 8]; 6] = [
    [3, 4, 5, 4, 3, 2, 1, 2], // Column 1
    [3, 2, 1, 2, 3, 4, 5, 4], // Column 2
    [8, 7, 6, 7, 8, 9, 10, 9], // Column 3
    [13, 12, 11, 12, 13, 14, 15, 14], // Column 4
    [18, 17, 16, 17, 18, 19, 20, 19], // Column 5
    [18, 17, 16, 17, 18, 19, 20, 19], // Column 6
];

const C0: [f32; 20] = [
    1.0136, 1.0478, 1.0, 1.0258, 1.0, 1.077, 1.0543, 1.0103, 0.99927, 0.96876,
    1.0971, 1.0971, 1.0777, 1.1134, 1.0237, 1.0703, 1.0248, 1.0945, 1.1622, 1.1393,
];
const C1: [f32; 20] = [
    -9.17e-5, -1.37e-5, 0.0, 7.14e-5, 0.0, -3.21e-4, -1.66e-4, -4.10e-5, 1.36e-4, 2.29e-4,
    -3.89e-4, -3.08e-4, -2.81e-4, -1.90e-4, 4.76e-5, -2.80e-4, -2.07e-4, -2.91e-4, -3.30e-4, -4.04e-4,
];
const C2: [f32; 20] = [
    1.16e-8, 0.0, 0.0, -1.46e-8, 0.0, 9.86e-8, 2.25e-8, -1.67e-8, -1.62e-8, -9.42e-8,
    1.17e-7, 4.32e-8, 3.97e-8, 3.13e-8, -8.04e-8, 3.91e-8, 2.58e-8, 3.45e-8, 4.76e-8, 1.13e-7,
];
const C3: [f32; 20] = [
    0.0, 0.0, 0.0, 0.0, 0.0, -9.44e-12, 0.0, 3.04e-12, 0.0, 9.32e-12,
    -1.07e-11, 0.0, 0.0, 0.0, 1.09e-11, 0.0, 0.0, 0.0, 0.0, -1.01e-11,
];

const FAP: [f32; 36] = [
    0.0, 0.0, 0.037037037, 0.074074074, 0.111111111, 0.148148148,
    0.185185185, 0.222222222, 0.259259259, 0.296296296, 0.333333333,
    0.37037037, 0.407407407, 0.444444444, 0.481481481, 0.518518519,
    0.555555556, 0.592592593, 0.62962963, 0.666666667, 0.703703704,
    0.740740741, 0.777777778, 0.814814815, 0.851851852, 0.888888889,
    0.925925926, 0.962962963, 1.0, 0.66666667, 0.33333334, 0.0, 0.333333,
    0.666666, 1.0, 0.7,
];

pub fn storm(
    ap: &[i32; 13],
    rga: f32,
    rgo: f32,
    coor: i32,
    ut: i32,
    doy: i32,
) -> (f32, f32) {
    let mut rgma = rga;
    if coor == 1 {
        rgma = conver(rga, rgo);
    } else if coor == 2 {
        rgma = rga;
    } else {
        return (1.0, rgma);
    }
    
    let mut ape = [0.0_f32; 39];
    ape[0] = ap[0] as f32;
    ape[1] = ap[0] as f32;
    ape[37] = ap[12] as f32;
    ape[38] = ap[12] as f32;
    
    for k in 1..=13 {
        let i = (k * 3) - 2;
        ape[i] = ap[k - 1] as f32;
    }
    
    for k in 1..=12 {
        let i = (k * 3) - 1;
        ape[i] = (ap[k - 1] as f32 * 2.0 + ap[k] as f32) / 3.0;
    }
    
    for k in 2..=13 {
        let i = (k * 3) - 3;
        ape[i] = (ap[k - 2] as f32 + ap[k - 1] as f32 * 2.0) / 3.0;
    }
    
    let mut ut_val = ut;
    if ut_val == 24 {
        ut_val = 0;
    }
    
    let k = match ut_val {
        0 | 3 | 6 | 9 | 12 | 15 | 18 | 21 => 0,
        1 | 4 | 7 | 10 | 13 | 16 | 19 | 22 => 1,
        2 | 5 | 8 | 11 | 14 | 17 | 20 | 23 => 2,
        _ => return (1.0, rgma),
    };
    
    let mut rap = 0.0_f32;
    for j in 0..36 {
        rap += FAP[j] * ape[k + j];
    }
    
    if rap <= 200.0 {
        return (1.0, rgma);
    }
    
    if doy > 366 || doy < 1 {
        return (1.0, rgma);
    }
    
    if rgma > 90.0 || rgma < -90.0 {
        return (1.0, rgma);
    }
    
    let mut dayno = doy;
    if rgma < 0.0 {
        dayno = doy + 172;
        if dayno > 365 {
            dayno -= 365;
        }
    }
    
    let rs = if dayno >= 82 {
        (dayno - 82) as f32 / 45.6 + 1.0
    } else {
        (dayno + 283) as f32 / 45.6 + 1.0
    };
    
    let s1 = rs.floor() as i32;
    let facs = rs - s1 as f32;
    let mut s2 = s1 + 1;
    if s2 == 9 {
        s2 = 1;
    }
    
    let abs_rgma = rgma.abs();
    let rl = (abs_rgma + 10.0) / 20.0 + 1.0;
    let rl = if rl == 6.0 { 5.9 } else { rl };
    let l1 = rl.floor() as i32;
    let facl = rl - l1 as f32;
    let l2 = l1 + 1;
    
    let calc_cf = |n: usize, r_val: f32| -> f32 {
        let idx = n - 1;
        C3[idx] * r_val.powi(3) + C2[idx] * r_val.powi(2) + C1[idx] * r_val + C0[idx]
    };
    
    if rap < 300.0 {
        let rapf = 300.0_f32;
        let n1 = CODE[l1 as usize - 1][s1 as usize - 1] as usize;
        let cf1 = calc_cf(n1, rapf);
        let n2 = CODE[l2 as usize - 1][s1 as usize - 1] as usize;
        let cf2 = calc_cf(n2, rapf);
        let n3 = CODE[l1 as usize - 1][s2 as usize - 1] as usize;
        let cf3 = calc_cf(n3, rapf);
        let n4 = CODE[l2 as usize - 1][s2 as usize - 1] as usize;
        let cf4 = calc_cf(n4, rapf);
        
        let cf300 = cf1 * (1.0 - facs) * (1.0 - facl)
            + cf2 * (1.0 - facs) * facl
            + cf3 * facs * (1.0 - facl)
            + cf4 * facs * facl;
            
        let cf = (cf300 - 1.0) * rap / 100.0 - 2.0 * cf300 + 3.0;
        (cf, rgma)
    } else {
        let n1 = CODE[l1 as usize - 1][s1 as usize - 1] as usize;
        let cf1 = calc_cf(n1, rap);
        let n2 = CODE[l2 as usize - 1][s1 as usize - 1] as usize;
        let cf2 = calc_cf(n2, rap);
        let n3 = CODE[l1 as usize - 1][s2 as usize - 1] as usize;
        let cf3 = calc_cf(n3, rap);
        let n4 = CODE[l2 as usize - 1][s2 as usize - 1] as usize;
        let cf4 = calc_cf(n4, rap);
        
        let cf = cf1 * (1.0 - facs) * (1.0 - facl)
            + cf2 * (1.0 - facs) * facl
            + cf3 * facs * (1.0 - facl)
            + cf4 * facs * facl;
        (cf, rgma)
    }
}

pub fn storme_ap(jdoy: i32, xmlat: f32, ap: f32) -> f32 {
    let jdoy_val = jdoy;
    let mut idxs = 0;
    let idbd = [79, 171, 264, 354, 366];
    if jdoy_val <= idbd[0] {
        idxs = 1;
    } else {
        for is in 2..=5 {
            if jdoy_val > idbd[is - 2] && jdoy_val <= idbd[is - 1] {
                idxs = is;
                break;
            }
        }
    }
    if idxs == 0 {
        return -5.0;
    }
    
    let xmlg: [f32; 37] = [
        -90.0, -85.0, -80.0, -75.0, -70.0, -65.0, -60.0, -55.0, -50.0,
        -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0,
        0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
        55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0
    ];
    let mut idxl = 0;
    let delg = (xmlg[0] - xmlg[1]).abs();
    let deld = delg / 2.0;
    let ymp_start = xmlg[0] + deld;
    let ymm_end = xmlg[36] - deld;
    if xmlat >= xmlg[0] && xmlat <= ymp_start {
        idxl = 1;
    } else if xmlat > ymm_end && xmlat <= xmlg[36] {
        idxl = 37;
    } else {
        for il in 2..=36 {
            let ymp = xmlg[il - 1] + deld;
            let ymm = xmlg[il - 1] - deld;
            if xmlat > ymm && xmlat <= ymp {
                idxl = il;
                break;
            }
        }
    }
    if idxl == 0 {
        return -5.0;
    }
    
    let c1_val = crate::storme_coeff::C1_STORME[idxl - 1][idxs - 1];
    let c2_val = crate::storme_coeff::C2_STORME[idxl - 1][idxs - 1];
    let c3_val = crate::storme_coeff::C3_STORME[idxl - 1][idxs - 1];
    
    let mut sqr = c1_val * ap.powf(c2_val) + c3_val;
    if sqr < 1.0 {
        sqr = 1.0;
    }
    sqr
}

pub fn auroral_boundary(
    xkp: f32,
    xmlt: f32,
    cgmlat: &mut f32,
    ab_mlat: &mut [f32; 48],
) {
    use crate::storme_coeff::ZP_MLAT;
    
    let mut xkp_clamped = xkp;
    if xkp_clamped > 9.0 {
        xkp_clamped = 9.0;
    }
    
    let kp1 = xkp_clamped.floor() as usize + 1;
    let xkp1 = xkp_clamped.floor();
    let mut kp2 = kp1 + 1;
    if kp2 > 10 {
        kp2 = 10;
    }
    
    for i in 0..48 {
        ab_mlat[i] = ZP_MLAT[i][kp1 - 1]
            + (xkp_clamped - xkp1) * (ZP_MLAT[i][kp2 - 1] - ZP_MLAT[i][kp1 - 1]);
    }
    
    *cgmlat = -99.99;
    if xmlt < 0.0 {
        return;
    }
    
    let mut ab_mlt = [0.0_f32; 48];
    for i in 0..48 {
        ab_mlt[i] = i as f32 * 0.5;
    }
    
    let mut i1 = (xmlt / 0.5) as usize + 1;
    if i1 >= 48 {
        i1 = 1;
    }
    let i2 = i1 + 1;
    
    let s1 = (ZP_MLAT[i2 - 1][kp1 - 1] - ZP_MLAT[i1 - 1][kp1 - 1]) / (ab_mlt[i2 - 1] - ab_mlt[i1 - 1]);
    let zmlkp1 = ZP_MLAT[i1 - 1][kp1 - 1] + (xmlt - ab_mlt[i1 - 1]) * s1;
    
    let s2 = (ZP_MLAT[i2 - 1][kp2 - 1] - ZP_MLAT[i1 - 1][kp2 - 1]) / (ab_mlt[i2 - 1] - ab_mlt[i1 - 1]);
    let zmlkp2 = ZP_MLAT[i1 - 1][kp2 - 1] + (xmlt - ab_mlt[i1 - 1]) * s2;
    
    *cgmlat = zmlkp1 + (xkp_clamped - xkp1) * (zmlkp2 - zmlkp1);
}

// --- Valentin Shubin's Peak Height Model & HMF2ED ---

const FT1: [f32; 12] = [
    73.6, 72.3, 71.8, 70.9, 73.6, 73.0,
    71.1, 69.4, 69.1, 70.9, 72.3, 74.6
];
const FT2: [f32; 12] = [
    144.2, 142.9, 167.2, 125.3, 124.4, 127.9,
    142.0, 165.9, 132.6, 142.0, 145.6, 143.0
];

pub fn legendre(mm: usize, nn: usize, p: &mut [[f64; 13]; 9], teta: f64) {
    let umr = (1.0_f64).atan() * 4.0 / 180.0;
    for r in p.iter_mut() {
        r.fill(0.0);
    }
    let z = (umr * teta).cos();
    p[0][0] = 1.0;
    p[0][1] = z;
    if mm != 0 {
        p[1][1] = (umr * teta).sin();
    }
    for j in 2..=mm {
        p[j][j] = (2 * j - 1) as f64 * p[j - 1][j - 1] * p[1][1];
    }
    for m in 0..=mm {
        for n in 1..=nn {
            if m > n {
                p[m][n] = 0.0;
                continue;
            }
            if n + 1 > nn {
                break;
            }
            if n + 1 == m {
                continue;
            }
            if m > (n - 1) {
                p[m][n + 1] = (2 * n + 1) as f64 * z * p[m][n] / (n + 1 - m) as f64;
            } else {
                p[m][n + 1] = ((2 * n + 1) as f64 * z * p[m][n] - (n + m) as f64 * p[m][n - 1]) / (n + 1 - m) as f64;
            }
        }
    }
    for n in 1..=nn {
        for m in 1..=mm {
            if m > n {
                p[m][n] = 0.0;
                continue;
            }
            let mut s = 1.0_f64;
            for l in (n - m + 1)..=(n + m) {
                s *= l as f64;
            }
            p[m][n] = p[m][n] * (2.0 / s).sqrt();
        }
    }
}

pub fn fun_gk(teta: f64, long: f32, gk: &mut [f64; 149]) {
    let mut pl_mn = [[0.0_f64; 13]; 9];
    let mm = 8;
    let nn = 12;
    legendre(mm, nn, &mut pl_mn, teta);
    gk.fill(0.0);
    let umr = (1.0_f64).atan() * 4.0 / 180.0;
    let mut k = 0;
    for m in 0..=mm {
        if m == 0 {
            for n in 0..=nn {
                gk[k] = pl_mn[m][n];
                k += 1;
            }
        } else {
            let angle = (m as f64) * (long as f64) * umr;
            let cos_val = angle.cos();
            let sin_val = angle.sin();
            for n in m..=nn {
                gk[k] = pl_mn[m][n] * cos_val;
                gk[k + 1] = pl_mn[m][n] * sin_val;
                k += 2;
            }
        }
    }
}

pub fn fun_hmf2_sd(teta: f64, long: f32, kf: &[f64; 149]) -> f32 {
    let mut gk = [0.0_f64; 149];
    fun_gk(teta, long, &mut gk);
    let mut hmf2 = 0.0_f64;
    for k in 0..149 {
        hmf2 += kf[k] * gk[k];
    }
    hmf2 as f32
}

pub fn hmf2_med_sd(
    i_ut: usize,
    monthut: i32,
    f107a: f32,
    xmodip: f32,
    long: f32,
    coeff_month: &[[f64; 48]; 149],
) -> f32 {
    let teta = (90.0 - xmodip) as f64;
    
    let mut kf1 = [0.0_f64; 149];
    let mut kf2 = [0.0_f64; 149];
    for i in 0..149 {
        kf1[i] = coeff_month[i][i_ut];
        kf2[i] = coeff_month[i][i_ut + 24];
    }
    
    let hmf2_1 = fun_hmf2_sd(teta, long, &kf1);
    let hmf2_2 = fun_hmf2_sd(teta, long, &kf2);
    
    let cov = f107a;
    let month_idx = (monthut as usize - 1).clamp(0, 11);
    let cov1 = FT1[month_idx];
    let cov2 = FT2[month_idx];
    
    let a = (hmf2_2 - hmf2_1) / (cov2 / cov1).ln();
    let b = hmf2_2 - a * cov2.ln();
    a * cov.ln() + b
}

pub fn fun_gk_ut(mm: usize, mk: usize, t: f64, gk_ut: &mut [f64; 7]) {
    let dtr = std::f64::consts::PI / 12.0;
    gk_ut.fill(0.0);
    let mut k = 0;
    for m in 0..=mm {
        if m == 0 {
            gk_ut[k] = 1.0;
            k += 1;
        } else {
            let angle = (m as f64) * t * dtr;
            gk_ut[k] = angle.cos();
            gk_ut[k + 1] = angle.sin();
            k += 2;
        }
    }
}

pub fn fun_fk_ut(mk: usize, gk_ut: &[f64; 7], akp_ut: &[[f64; 7]; 7], fk_ut: &mut [f64; 7]) {
    fk_ut.fill(0.0);
    for k in 0..=mk {
        let mut sum_g = 0.0_f64;
        for p in 0..k {
            sum_g += akp_ut[k][p] * fk_ut[p];
        }
        fk_ut[k] = sum_g + gk_ut[k];
    }
}

pub fn fun_akp_ut(
    mm: usize,
    mk: usize,
    akp_ut: &mut [[f64; 7]; 7],
    dk_ut: &mut [f64; 7],
    hmf2_ut: &[f64; 24],
) {
    let mut gk_ut = [0.0_f64; 7];
    gk_ut[0] = 1.0;
    let mut fk_ut = [0.0_f64; 7];
    fk_ut[0] = 1.0;
    
    for r in akp_ut.iter_mut() {
        r.fill(0.0);
    }
    dk_ut.fill(0.0);
    
    for p in 0..=mk {
        let mut sum_dn = 0.0_f64;
        let mut sum_dd = 0.0_f64;
        for k in (p + 1)..=mk {
            let mut sum_an = 0.0_f64;
            let mut sum_ad = 0.0_f64;
            for i in 0..24 {
                let t = i as f64;
                fun_gk_ut(mm, mk, t, &mut gk_ut);
                fun_fk_ut(mk, &gk_ut, akp_ut, &mut fk_ut);
                sum_an += gk_ut[k] * fk_ut[p];
                sum_ad += fk_ut[p] * fk_ut[p];
                if p == k - 1 {
                    sum_dn += hmf2_ut[i] * fk_ut[p];
                    sum_dd += fk_ut[p] * fk_ut[p];
                }
            }
            akp_ut[k][p] = -sum_an / sum_ad;
        }
        if p < mk {
            dk_ut[p] = sum_dn / sum_dd;
        }
    }
    
    let p = mk;
    let mut sum_dn = 0.0_f64;
    let mut sum_dd = 0.0_f64;
    for i in 0..24 {
        let t = i as f64;
        fun_gk_ut(mm, mk, t, &mut gk_ut);
        fun_fk_ut(mk, &gk_ut, akp_ut, &mut fk_ut);
        sum_dn += hmf2_ut[i] * fk_ut[p];
        sum_dd += fk_ut[p] * fk_ut[p];
    }
    dk_ut[p] = sum_dn / sum_dd;
}

pub fn koeff_ut(mm: usize, mk: usize, kf_ut: &mut [f64; 7], hmf2_ut: &[f64; 24]) {
    let mut akp_ut = [[0.0_f64; 7]; 7];
    let mut dk_ut = [0.0_f64; 7];
    fun_akp_ut(mm, mk, &mut akp_ut, &mut dk_ut, hmf2_ut);
    
    kf_ut.fill(0.0);
    for k in (0..=mk).rev() {
        let mut sum_d = 0.0_f64;
        for m in (k + 1)..=mk {
            sum_d += akp_ut[m][k] * kf_ut[m];
        }
        kf_ut[k] = sum_d + dk_ut[k];
    }
}

pub fn fun_hmf2ut(t: f64, hmf2_ut: &[f64; 24]) -> f64 {
    let mm = 3;
    let mk = 6;
    let mut kf_ut = [0.0_f64; 7];
    koeff_ut(mm, mk, &mut kf_ut, hmf2_ut);
    
    let mut gk_ut = [0.0_f64; 7];
    fun_gk_ut(mm, mk, t, &mut gk_ut);
    
    let mut hmf2 = 0.0_f64;
    for k in 0..=mk {
        hmf2 += kf_ut[k] * gk_ut[k];
    }
    hmf2
}

pub fn sdmf2(
    ut: f32,
    monthut: i32,
    f107a: f32,
    xmodip: f32,
    long: f32,
    data_dir: &str,
) -> f32 {
    let coeff_month = match crate::data_io::McsatData::load(data_dir, monthut) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Failed to load mcsat data for month {}: {}", monthut, e);
            return 0.0;
        }
    };
    
    let mut hmf2_ut = [0.0_f64; 24];
    for i in 0..24 {
        hmf2_ut[i] = hmf2_med_sd(i, monthut, f107a, xmodip, long, &coeff_month) as f64;
    }
    
    fun_hmf2ut(ut as f64, &hmf2_ut) as f32
}

pub fn model_hmf2(
    day: i32,
    month: i32,
    ut: f32,
    xmodip: f32,
    long: f32,
    f107_81: f32,
    data_dir: &str,
) -> f32 {
    let hm_f2_0 = sdmf2(ut, month, f107_81, xmodip, long, data_dir);
    let hm_f2_med = if day <= 15 {
        if day == 15 {
            hm_f2_0
        } else {
            let mut monthr = month - 1;
            if monthr == 0 {
                monthr = 12;
            }
            let hm_f2_m = sdmf2(ut, monthr, f107_81, xmodip, long, data_dir);
            hm_f2_0 - (day - 15) as f32 * (hm_f2_m - hm_f2_0) / 30.0
        }
    } else {
        let montha = (month % 12) + 1;
        let hm_f2_p = sdmf2(ut, montha, f107_81, xmodip, long, data_dir);
        hm_f2_0 + (day - 15) as f32 * (hm_f2_p - hm_f2_0) / 30.0
    };
    hm_f2_med
}

pub fn hmf2ed(xmagbr: f32, r: f32, x: f32, xm3: f32) -> f32 {
    let f1 = 0.00232 * r + 0.222;
    let f2 = 1.2 - 0.0116 * (0.0239 * r).exp();
    let f3 = 0.096 * (r - 25.0) / 150.0;
    let f4 = 1.0 - r / 150.0 * (-xmagbr * xmagbr / 1600.0).exp();
    let mut x_clamped = x;
    if x_clamped < 1.7 {
        x_clamped = 1.7;
    }
    let delm = f1 * f4 / (x_clamped - f2) + f3;
    1490.0 / (xm3 + delm) - 176.0
}



const TEBA_C: [[[f32; 81]; 2]; 4] = [
    [
        [
            3.1_f32, -0.003215_f32, 0.244_f32, -0.0004613_f32, -0.01711_f32, 0.02605_f32,
            -0.09546_f32, 0.01794_f32, 0.0127_f32, 0.02791_f32, 0.01536_f32, -0.006629_f32,
            -0.003616_f32, 0.01229_f32, 0.0004147_f32, 0.001447_f32, -0.0004453_f32, -0.1853_f32,
            -0.01245_f32, -0.03675_f32, 0.004965_f32, 0.00546_f32, 0.008117_f32, -0.01002_f32,
            0.0005466_f32, -0.03087_f32, -0.003435_f32, -0.0001107_f32, 0.002199_f32, 0.0004115_f32,
            0.0006061_f32, 0.0002916_f32, -0.06584_f32, 0.004729_f32, -0.001523_f32, 0.0006689_f32,
            0.001031_f32, 0.0005398_f32, -0.001924_f32, -0.04565_f32, 0.007244_f32, -8.543e-05_f32,
            0.001052_f32, -0.0006696_f32, -0.0007492_f32, 0.04405_f32, 0.003047_f32, 0.002858_f32,
            -0.0001465_f32, 0.001195_f32, -0.0001024_f32, 0.04582_f32, 0.0008749_f32, 0.0003011_f32,
            0.0004473_f32, -0.0002782_f32, 0.04911_f32, -0.01016_f32, 0.0027_f32, -0.0009304_f32,
            -0.001202_f32, 0.0221_f32, 0.002566_f32, -0.000122_f32, 0.0003987_f32, -0.05744_f32,
            0.004408_f32, -0.003497_f32, 0.00083_f32, -0.03536_f32, -0.008813_f32, 0.002423_f32,
            -0.02994_f32, -0.001929_f32, -0.0005268_f32, -0.02228_f32, 0.003385_f32, 0.0413_f32,
            0.004876_f32, 0.02692_f32, 0.001684_f32,
        ],
        [
            3.13654_f32, 0.006796_f32, 0.181413_f32, 0.08564_f32, -0.032856_f32, -0.003508_f32,
            -0.01438_f32, -0.02454_f32, 0.002745_f32, 0.05284_f32, 0.01136_f32, -0.01956_f32,
            -0.005805_f32, 0.002801_f32, -0.001211_f32, 0.004127_f32, 0.002909_f32, -0.25751_f32,
            -0.0037915_f32, -0.0136_f32, -0.013225_f32, 0.01202_f32, 0.01256_f32, -0.012165_f32,
            0.01326_f32, -0.07123_f32, 0.0005793_f32, 0.001537_f32, 0.006914_f32, -0.004173_f32,
            0.0001052_f32, -0.0005765_f32, -0.04041_f32, -0.001752_f32, -0.00542_f32, -0.00684_f32,
            0.0008921_f32, -0.002228_f32, 0.001428_f32, 0.006635_f32, -0.0048045_f32, -0.001659_f32,
            -0.0009341_f32, 0.000223_f32, -0.0009995_f32, 0.04285_f32, -0.0005211_f32, -0.003293_f32,
            0.00179_f32, 0.0006435_f32, -0.0001891_f32, 0.03844_f32, 0.00359_f32, -0.0008139_f32,
            -0.001996_f32, 0.0002398_f32, 0.02938_f32, 0.00761_f32, 0.00347655_f32, 0.001707_f32,
            0.0002769_f32, -0.0157_f32, 0.000983_f32, -0.0006532_f32, 9.29e-05_f32, -0.02506_f32,
            0.004681_f32, 0.001461_f32, -3.757e-06_f32, -0.009728_f32, 0.002315_f32, 0.0006377_f32,
            -0.01705_f32, 0.002767_f32, -0.0006992_f32, -0.0115_f32, -0.001644_f32, 0.003355_f32,
            -0.004326_f32, 0.02035_f32, 0.02985_f32,
        ],
    ],
    [
        [
            3.136_f32, 0.006498_f32, 0.2289_f32, 0.01859_f32, -0.03328_f32, -0.004889_f32,
            -0.03054_f32, -0.01773_f32, -0.01728_f32, 0.06555_f32, 0.01775_f32, -0.02488_f32,
            -0.009498_f32, 0.01493_f32, 0.00281_f32, 0.002406_f32, 0.005436_f32, -0.2115_f32,
            0.007007_f32, -0.05129_f32, -0.007327_f32, 0.02402_f32, 0.004772_f32, -0.007374_f32,
            -0.0003835_f32, -0.05013_f32, 0.002866_f32, 0.002216_f32, 0.0002412_f32, 0.002094_f32,
            0.00122_f32, -0.0001703_f32, -0.1082_f32, -0.004992_f32, -0.004065_f32, 0.003615_f32,
            -0.002738_f32, -0.0007177_f32, 0.0002173_f32, -0.04373_f32, -0.00375_f32, 0.005507_f32,
            -0.001567_f32, -0.001458_f32, -0.0007397_f32, 0.07903_f32, 0.004131_f32, 0.003714_f32,
            0.001073_f32, -0.0008991_f32, 0.0002976_f32, 0.02623_f32, 0.002344_f32, 0.0005608_f32,
            0.0004124_f32, 0.0001509_f32, 0.05103_f32, 0.00345_f32, 0.001283_f32, 0.0007238_f32,
            -3.464e-05_f32, 0.01663_f32, -0.001644_f32, -0.00071_f32, 0.0005281_f32, -0.02729_f32,
            0.003556_f32, -0.003391_f32, -0.0001787_f32, 0.002154_f32, 0.006476_f32, -0.0008282_f32,
            -0.02361_f32, 0.0009557_f32, 0.0003205_f32, -0.02301_f32, -0.000854_f32, -0.01126_f32,
            -0.002323_f32, -0.008582_f32, 0.02683_f32,
        ],
        [
            3.144_f32, 0.008571_f32, 0.2539_f32, 0.06937_f32, -0.01667_f32, 0.02249_f32,
            -0.04162_f32, 0.01201_f32, 0.02435_f32, 0.05232_f32, 0.02521_f32, -0.0199_f32,
            -0.007671_f32, 0.01264_f32, -0.001551_f32, -0.001928_f32, 0.003652_f32, -0.2019_f32,
            0.005697_f32, -0.03159_f32, -0.01451_f32, 0.02868_f32, 0.01377_f32, -0.004383_f32,
            0.01172_f32, -0.05683_f32, 0.003593_f32, 0.003571_f32, 0.003282_f32, 0.001732_f32,
            -0.0004921_f32, -0.001165_f32, -0.1066_f32, -0.01892_f32, 0.00357_f32, -0.0008631_f32,
            -0.001876_f32, -8.414e-05_f32, 0.002356_f32, -0.04259_f32, -0.00322_f32, 0.004641_f32,
            0.0006223_f32, -0.00168_f32, -0.0001243_f32, 0.07393_f32, -0.003143_f32, -0.002362_f32,
            0.001235_f32, -0.001551_f32, 0.0002099_f32, 0.02299_f32, 0.005301_f32, -0.004306_f32,
            -0.001303_f32, 7.687e-06_f32, 0.05305_f32, 0.006642_f32, -0.001686_f32, 0.001048_f32,
            0.0005958_f32, 0.04341_f32, -8.819e-05_f32, -0.000333_f32, -0.0002158_f32, -0.04106_f32,
            0.004191_f32, 0.002045_f32, -0.0001437_f32, -0.01803_f32, -0.0008072_f32, -0.000424_f32,
            -0.026_f32, -0.002329_f32, 0.0005949_f32, -0.01371_f32, -0.002188_f32, 0.01788_f32,
            0.0006405_f32, 0.005977_f32, 0.01333_f32,
        ],
    ],
    [
        [
            3.372_f32, 0.01006_f32, 0.1436_f32, 0.002023_f32, -0.05166_f32, 0.009606_f32,
            -0.05596_f32, 0.0004914_f32, -0.003124_f32, -0.04713_f32, -0.007371_f32, -0.004823_f32,
            -0.002213_f32, 0.006569_f32, -0.0001962_f32, 0.0003309_f32, -0.0003908_f32, -0.2836_f32,
            0.007829_f32, 0.01175_f32, 0.0009919_f32, 0.006589_f32, 0.002045_f32, -0.007346_f32,
            -0.00089_f32, -0.0347_f32, -0.004977_f32, 0.00147_f32, -2.823e-06_f32, 0.0006465_f32,
            -0.0001448_f32, 0.001401_f32, -0.08988_f32, -3.293e-05_f32, -0.001848_f32, 0.0004439_f32,
            -0.001263_f32, 0.000317_f32, -0.0006227_f32, 0.01721_f32, -0.00199_f32, -0.0004627_f32,
            2.897e-06_f32, -0.0005454_f32, 0.0003385_f32, 0.08432_f32, -0.001951_f32, 0.001487_f32,
            0.001042_f32, -0.0004788_f32, -0.0001276_f32, 0.02373_f32, 0.002409_f32, 0.0005263_f32,
            0.001301_f32, -0.0004177_f32, 0.03974_f32, 0.0001418_f32, -0.001048_f32, -0.0002982_f32,
            -3.396e-05_f32, 0.0131_f32, 0.001413_f32, -0.0001373_f32, 0.0002638_f32, -0.04171_f32,
            -0.0005932_f32, -0.0007523_f32, -0.0006883_f32, -0.02355_f32, 0.0005695_f32, -2.219e-05_f32,
            -0.02301_f32, -9.962e-05_f32, -0.0006761_f32, 0.00204_f32, -0.0005479_f32, 0.02591_f32,
            -0.002425_f32, 0.01583_f32, 0.009577_f32,
        ],
        [
            3.367_f32, 0.01038_f32, 0.1407_f32, 0.03622_f32, -0.03144_f32, 0.0112_f32,
            -0.05674_f32, 0.03219_f32, 0.001288_f32, -0.05799_f32, -0.004609_f32, 0.003252_f32,
            -0.0002859_f32, 0.01226_f32, -0.004539_f32, 0.00131_f32, -0.0005603_f32, -0.311_f32,
            -0.001268_f32, 0.01539_f32, 0.003146_f32, 0.007787_f32, -0.00143_f32, -0.00482_f32,
            0.002924_f32, -0.09981_f32, -0.007838_f32, -0.0001663_f32, 0.0004769_f32, 0.004148_f32,
            -0.001008_f32, -0.000979_f32, -0.09049_f32, -0.002994_f32, -0.006748_f32, -0.0009889_f32,
            0.001488_f32, -0.001154_f32, -8.412e-05_f32, -0.01302_f32, -0.004859_f32, -0.0007172_f32,
            -0.0009401_f32, 0.0009101_f32, -0.0001735_f32, 0.07055_f32, 0.006398_f32, -0.003103_f32,
            -0.000938_f32, -0.0004_f32, -0.001165_f32, 0.02713_f32, -0.001654_f32, 0.002781_f32,
            -5.215e-06_f32, 0.0002258_f32, 0.05022_f32, 0.0095_f32, 0.0004147_f32, 0.0003499_f32,
            -0.0006097_f32, 0.04118_f32, 0.006556_f32, 0.003793_f32, -0.0001226_f32, -0.02517_f32,
            0.0001491_f32, 0.001075_f32, 0.0004531_f32, -0.009012_f32, 0.003343_f32, 0.003431_f32,
            -0.02519_f32, 3.793e-05_f32, 0.0005973_f32, -0.01423_f32, -0.00132_f32, -0.006048_f32,
            -0.005005_f32, -0.0115_f32, 0.02574_f32,
        ],
    ],
    [
        [
            3.574_f32, 0.0_f32, 0.07537_f32, 0.0_f32, -0.08459_f32, 0.0_f32,
            -0.0294_f32, 0.0_f32, 0.04547_f32, -0.05321_f32, 0.0_f32, 0.004328_f32,
            0.0_f32, 0.006022_f32, 0.0_f32, -0.0009168_f32, 0.0_f32, -0.1768_f32,
            0.0_f32, 0.0294_f32, 0.0_f32, 0.0005902_f32, 0.0_f32, -0.009047_f32,
            0.0_f32, -0.06555_f32, 0.0_f32, -0.001033_f32, 0.0_f32, 0.001674_f32,
            0.0_f32, 0.0002802_f32, -0.06786_f32, 0.0_f32, 0.004193_f32, 0.0_f32,
            -0.0006448_f32, 0.0_f32, 0.0009277_f32, -0.01634_f32, 0.0_f32, -0.002531_f32,
            0.0_f32, 1.93e-05_f32, 0.0_f32, 0.0528_f32, 0.0_f32, 0.002438_f32,
            0.0_f32, -0.0005292_f32, 0.0_f32, 0.01555_f32, 0.0_f32, -0.003259_f32,
            0.0_f32, -0.0005998_f32, 0.03168_f32, 0.0_f32, 0.002382_f32, 0.0_f32,
            -0.0004078_f32, 0.02312_f32, 0.0_f32, 0.0001481_f32, 0.0_f32, -0.01885_f32,
            0.0_f32, 0.001144_f32, 0.0_f32, -0.009952_f32, 0.0_f32, -0.000551_f32,
            -0.0202_f32, 0.0_f32, -7.283e-05_f32, -0.01272_f32, 0.0_f32, 0.002224_f32,
            0.0_f32, -0.00251_f32, 0.02434_f32,
        ],
        [
            3.574_f32, -0.005639_f32, 0.07094_f32, -0.03347_f32, -0.0861_f32, -0.02877_f32,
            -0.03154_f32, -0.002847_f32, 0.01235_f32, -0.05966_f32, -0.003236_f32, 0.0003795_f32,
            -0.0008634_f32, 0.003377_f32, -0.0001071_f32, -0.002151_f32, -0.0004057_f32, -0.1783_f32,
            0.0126_f32, 0.02835_f32, -0.00242_f32, 0.003002_f32, -0.004684_f32, -0.006756_f32,
            -0.0007493_f32, -0.06147_f32, -0.005636_f32, -0.001234_f32, -0.001613_f32, -6.353e-05_f32,
            -0.0002503_f32, -0.0001729_f32, -0.07148_f32, 0.005326_f32, 0.004006_f32, 0.0006484_f32,
            -0.0001046_f32, -0.0006034_f32, -0.0009435_f32, -0.002385_f32, 0.006853_f32, 0.00151_f32,
            0.001319_f32, 9.049e-05_f32, -0.0001999_f32, 0.03976_f32, 0.002802_f32, -0.00103_f32,
            0.0005599_f32, -0.0004791_f32, -8.46e-05_f32, 0.02683_f32, 0.00427_f32, 0.0005911_f32,
            0.0002987_f32, -0.000208_f32, 0.01396_f32, -0.001922_f32, -0.001063_f32, 0.0003803_f32,
            0.0001343_f32, 0.01771_f32, -0.001038_f32, -0.0004645_f32, -0.0002481_f32, -0.02251_f32,
            -0.0029_f32, -0.0003977_f32, -0.000516_f32, -0.008079_f32, -0.001528_f32, 0.000306_f32,
            -0.01582_f32, -0.0008536_f32, 0.0001565_f32, -0.01252_f32, 0.0002319_f32, 0.004311_f32,
            0.001024_f32, 1.296e-06_f32, 0.0179_f32,
        ],
    ],
];


pub fn teba(mut dipl: f32, slt: f32, ns: i32, te: &mut [f32; 4]) {
    let mut is_val = 1;
    if ns < 3 {
        is_val = ns;
    } else if ns > 3 {
        is_val = 2;
        dipl = -dipl;
    } else {
        is_val = 1;
    }
    
    let mut colat = 0.017453292519943295 * (90.0 - dipl);
    let az = 0.2617993877991494 * slt; // humr is pi/12
    
    let mut a = [0.0_f32; 82];
    crate::spharm::spharm(&mut a, 8, 8, colat, az);
    
    let kend = if is_val == 2 { 3 } else { 4 };
    
    for k in 1..=kend {
        let mut ste = 0.0_f32;
        for i in 1..=81 {
            ste += a[i - 1] * TEBA_C[k - 1][is_val as usize - 1][i - 1];
        }
        te[k - 1] = 10.0_f32.powf(ste);
    }
    
    if is_val == 2 {
        dipl = -dipl;
        colat = 0.017453292519943295 * (90.0 - dipl);
        crate::spharm::spharm(&mut a, 8, 8, colat, az);
        let mut ste = 0.0_f32;
        for i in 1..=81 {
            ste += a[i - 1] * TEBA_C[3][1][i - 1];
        }
        te[3] = 10.0_f32.powf(ste);
    }
}

pub fn elteik(
    _pf107y: i32,
    _invdip: f32,
    _mlt: f32,
    _ddd: f32,
    _pf107: f32,
    tev: &mut [f32; 5],
    sigtev: &mut [f32; 5],
) {
    tev[0] = 1000.0;
    tev[1] = 1500.0;
    tev[2] = 2000.0;
    tev[3] = 2500.0;
    tev[4] = 3000.0;
    for i in 0..5 {
        sigtev[i] = 100.0;
    }
}

pub fn tede(h: f32, _den: f32, _cov: f32) -> f32 {
    1000.0 + h * 2.0
}

pub fn iontif(
    _pf107y: i32,
    _invdip: f32,
    _mlt: f32,
    _ddd: f32,
    _pf107: f32,
    tiv: &mut [f32; 4],
    sigtiv: &mut [f32; 4],
) {
    tiv[0] = 1000.0;
    tiv[1] = 1200.0;
    tiv[2] = 1500.0;
    tiv[3] = 1800.0;
    for i in 0..4 {
        sigtiv[i] = 50.0;
    }
}

pub fn chemion(
    _jprint: i32,
    _height: f32,
    _f107yobs: f32,
    _f10781obs: f32,
    _teh: f32,
    _tih: f32,
    _tnh: f32,
    _o_dens: f32,
    _o2_dens: f32,
    _n2_dens: f32,
    _n_dens: f32,
    _h_dens: f32,
    _user_no: f32,
    _xn4s: f32,
    edens: f32,
    _user_oplus: f32,
    _xhi_step: f32,
    ro: &mut f32,
    ro2: &mut f32,
    rno: &mut f32,
    rn2: &mut f32,
    rn: &mut f32,
    den_no: &mut f32,
    den_n2d: &mut f32,
    inewt: &mut i32,
) {
    *ro = 0.1 * edens;
    *ro2 = 0.4 * edens;
    *rno = 0.5 * edens;
    *rn2 = 0.0;
    *rn = 0.0;
    *den_no = 0.0;
    *den_n2d = 0.0;
    *inewt = 1;
}
