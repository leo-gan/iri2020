use crate::igrf::PI;

pub const UMR: f32 = PI / 180.0;

#[derive(Debug, Clone, Copy)]
pub struct DipoleParams {
    pub st0: f32,
    pub ct0: f32,
    pub sl0: f32,
    pub cl0: f32,
    pub ctcl: f32,
    pub stcl: f32,
    pub ctsl: f32,
    pub stsl: f32,
}

pub fn recalc(iyr: i32) -> DipoleParams {
    let iy = if iyr < 1900 {
        1900
    } else if iyr > 2025 {
        2025
    } else {
        iyr
    };

    let (g10, g11, h11) = if iy < 1905 {
        let f2 = (iy - 1900) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (31543.0 * f1 + 31464.0 * f2, -2298.0 * f1 - 2298.0 * f2, 5922.0 * f1 + 5909.0 * f2)
    } else if iy < 1910 {
        let f2 = (iy - 1905) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (31464.0 * f1 + 31354.0 * f2, -2298.0 * f1 - 2297.0 * f2, 5909.0 * f1 + 5898.0 * f2)
    } else if iy < 1915 {
        let f2 = (iy - 1910) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (31354.0 * f1 + 31212.0 * f2, -2297.0 * f1 - 2306.0 * f2, 5898.0 * f1 + 5875.0 * f2)
    } else if iy < 1920 {
        let f2 = (iy - 1915) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (31212.0 * f1 + 31060.0 * f2, -2306.0 * f1 - 2317.0 * f2, 5875.0 * f1 + 5845.0 * f2)
    } else if iy < 1925 {
        let f2 = (iy - 1920) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (31060.0 * f1 + 30926.0 * f2, -2317.0 * f1 - 2318.0 * f2, 5845.0 * f1 + 5817.0 * f2)
    } else if iy < 1930 {
        let f2 = (iy - 1925) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30926.0 * f1 + 30805.0 * f2, -2318.0 * f1 - 2316.0 * f2, 5817.0 * f1 + 5808.0 * f2)
    } else if iy < 1935 {
        let f2 = (iy - 1930) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30805.0 * f1 + 30715.0 * f2, -2316.0 * f1 - 2306.0 * f2, 5808.0 * f1 + 5812.0 * f2)
    } else if iy < 1940 {
        let f2 = (iy - 1935) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30715.0 * f1 + 30654.0 * f2, -2306.0 * f1 - 2292.0 * f2, 5812.0 * f1 + 5821.0 * f2)
    } else if iy < 1945 {
        let f2 = (iy - 1940) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30654.0 * f1 + 30594.0 * f2, -2292.0 * f1 - 2285.0 * f2, 5821.0 * f1 + 5810.0 * f2)
    } else if iy < 1950 {
        let f2 = (iy - 1945) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30594.0 * f1 + 30554.0 * f2, -2285.0 * f1 - 2250.0 * f2, 5810.0 * f1 + 5815.0 * f2)
    } else if iy < 1955 {
        let f2 = (iy - 1950) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30554.0 * f1 + 30500.0 * f2, -2250.0 * f1 - 2215.0 * f2, 5815.0 * f1 + 5820.0 * f2)
    } else if iy < 1960 {
        let f2 = (iy - 1955) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30500.0 * f1 + 30421.0 * f2, -2215.0 * f1 - 2169.0 * f2, 5820.0 * f1 + 5791.0 * f2)
    } else if iy < 1965 {
        let f2 = (iy - 1960) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30421.0 * f1 + 30334.0 * f2, -2169.0 * f1 - 2119.0 * f2, 5791.0 * f1 + 5776.0 * f2)
    } else if iy < 1970 {
        let f2 = (iy - 1965) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30334.0 * f1 + 30220.0 * f2, -2119.0 * f1 - 2068.0 * f2, 5776.0 * f1 + 5737.0 * f2)
    } else if iy < 1975 {
        let f2 = (iy - 1970) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30220.0 * f1 + 30100.0 * f2, -2068.0 * f1 - 2013.0 * f2, 5737.0 * f1 + 5675.0 * f2)
    } else if iy < 1980 {
        let f2 = (iy - 1975) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (30100.0 * f1 + 29992.0 * f2, -2013.0 * f1 - 1956.0 * f2, 5675.0 * f1 + 5604.0 * f2)
    } else if iy < 1985 {
        let f2 = (iy - 1980) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29992.0 * f1 + 29873.0 * f2, -1956.0 * f1 - 1905.0 * f2, 5604.0 * f1 + 5500.0 * f2)
    } else if iy < 1990 {
        let f2 = (iy - 1985) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29873.0 * f1 + 29775.0 * f2, -1905.0 * f1 - 1848.0 * f2, 5500.0 * f1 + 5406.0 * f2)
    } else if iy < 1995 {
        let f2 = (iy - 1990) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29775.0 * f1 + 29692.0 * f2, -1848.0 * f1 - 1784.0 * f2, 5406.0 * f1 + 5306.0 * f2)
    } else if iy < 2000 {
        let f2 = (iy - 1995) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29692.0 * f1 + 29619.4 * f2, -1784.0 * f1 - 1728.2 * f2, 5306.0 * f1 + 5186.1 * f2)
    } else if iy < 2005 {
        let f2 = (iy - 2000) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29619.4 * f1 + 29554.63 * f2, -1728.2 * f1 - 1669.05 * f2, 5186.1 * f1 + 5077.99 * f2)
    } else if iy < 2010 {
        let f2 = (iy - 2005) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29554.63 * f1 + 29496.57 * f2, -1669.05 * f1 - 1586.42 * f2, 5077.99 * f1 + 4944.26 * f2)
    } else if iy < 2015 {
        let f2 = (iy - 2010) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29496.57 * f1 + 29441.46 * f2, -1586.42 * f1 - 1501.77 * f2, 4944.26 * f1 + 4795.99 * f2)
    } else if iy < 2020 {
        let f2 = (iy - 2015) as f32 / 5.0;
        let f1 = 1.0 - f2;
        (29441.46 * f1 + 29404.8 * f2, -1501.77 * f1 - 1450.9 * f2, 4795.99 * f1 + 4652.5 * f2)
    } else {
        let dt = (iy - 2020) as f32;
        (29404.8 - 5.7 * dt, -1450.9 + 7.4 * dt, 4652.5 - 25.9 * dt)
    };

    let sq = g11 * g11 + h11 * h11;
    let sqq = sq.sqrt();
    let sqr = (g10 * g10 + sq).sqrt();
    let sl0 = -h11 / sqq;
    let cl0 = -g11 / sqq;
    let st0 = sqq / sqr;
    let ct0 = g10 / sqr;

    DipoleParams {
        st0,
        ct0,
        sl0,
        cl0,
        ctcl: ct0 * cl0,
        stcl: st0 * cl0,
        ctsl: ct0 * sl0,
        stsl: st0 * sl0,
    }
}

pub fn sphcar(r: &mut f32, teta: &mut f32, phi: &mut f32, x: &mut f32, y: &mut f32, z: &mut f32, j: i32) {
    if j > 0 {
        let sq = *r * teta.sin();
        *x = sq * phi.cos();
        *y = sq * phi.sin();
        *z = *r * teta.cos();
    } else {
        let sq = *x * *x + *y * *y;
        *r = (sq + *z * *z).sqrt();
        if sq == 0.0 {
            *phi = 0.0;
            if *z < 0.0 {
                *teta = 3.141592654;
            } else {
                *teta = 0.0;
            }
        } else {
            let sq_sqrt = sq.sqrt();
            *phi = y.atan2(*x);
            *teta = sq_sqrt.atan2(*z);
            if *phi < 0.0 {
                *phi += 6.28318531;
            }
        }
    }
}

pub fn bspcar(teta: f32, phi: f32, br: f32, btet: f32, bphi: f32, bx: &mut f32, by: &mut f32, bz: &mut f32) {
    let s = teta.sin();
    let c = teta.cos();
    let sf = phi.sin();
    let cf = phi.cos();
    let be = br * s + btet * c;
    *bx = be * cf - bphi * sf;
    *by = be * sf + bphi * cf;
    *bz = br * c - btet * s;
}

pub fn geomag(
    xgeo: &mut f32,
    ygeo: &mut f32,
    zgeo: &mut f32,
    xmag: &mut f32,
    ymag: &mut f32,
    zmag: &mut f32,
    j: i32,
    params: &DipoleParams,
) {
    if j < 0 {
        *xgeo = *xmag * params.ctcl - *ymag * params.sl0 + *zmag * params.stcl;
        *ygeo = *xmag * params.ctsl + *ymag * params.cl0 + *zmag * params.stsl;
        *zgeo = *zmag * params.ct0 - *xmag * params.st0;
    } else {
        *xmag = *xgeo * params.ctcl + *ygeo * params.ctsl - *zgeo * params.st0;
        *ymag = *ygeo * params.cl0 - *xgeo * params.sl0;
        *zmag = *xgeo * params.stcl + *ygeo * params.stsl + *zgeo * params.ct0;
    }
}

pub fn right(x: f32, y: f32, z: f32, r1: &mut f32, r2: &mut f32, r3: &mut f32, iyr: i32, nm: i32, ds3: f32) {
    let (mut r, mut t, mut f) = (0.0, 0.0, 0.0);
    let (mut x_tmp, mut y_tmp, mut z_tmp) = (x, y, z);
    sphcar(&mut r, &mut t, &mut f, &mut x_tmp, &mut y_tmp, &mut z_tmp, -1);

    let res = crate::igrf::igrf(iyr, nm, r, t, f);
    let (mut bx, mut by, mut bz) = (0.0, 0.0, 0.0);
    bspcar(t, f, res.br, res.bt, res.bf, &mut bx, &mut by, &mut bz);

    let b = ds3 / (bx * bx + by * by + bz * bz).sqrt();
    *r1 = bx * b;
    *r2 = by * b;
    *r3 = bz * b;
}

pub fn shag(x: &mut f32, y: &mut f32, z: &mut f32, ds: f32, iyr: i32, nm: i32) {
    let ds3 = -ds / 3.0;
    let (mut r11, mut r12, mut r13) = (0.0, 0.0, 0.0);
    let (mut r21, mut r22, mut r23) = (0.0, 0.0, 0.0);
    let (mut r31, mut r32, mut r33) = (0.0, 0.0, 0.0);
    let (mut r41, mut r42, mut r43) = (0.0, 0.0, 0.0);
    let (mut r51, mut r52, mut r53) = (0.0, 0.0, 0.0);

    right(*x, *y, *z, &mut r11, &mut r12, &mut r13, iyr, nm, ds3);
    right(*x + r11, *y + r12, *z + r13, &mut r21, &mut r22, &mut r23, iyr, nm, ds3);
    right(
        *x + 0.5 * (r11 + r21),
        *y + 0.5 * (r12 + r22),
        *z + 0.5 * (r13 + r23),
        &mut r31,
        &mut r32,
        &mut r33,
        iyr,
        nm,
        ds3,
    );
    right(
        *x + 0.375 * (r11 + 3.0 * r31),
        *y + 0.375 * (r12 + 3.0 * r32),
        *z + 0.375 * (r13 + 3.0 * r33),
        &mut r41,
        &mut r42,
        &mut r43,
        iyr,
        nm,
        ds3,
    );
    right(
        *x + 1.5 * (r11 - 3.0 * r31 + 4.0 * r41),
        *y + 1.5 * (r12 - 3.0 * r32 + 4.0 * r42),
        *z + 1.5 * (r13 - 3.0 * r33 + 4.0 * r43),
        &mut r51,
        &mut r52,
        &mut r53,
        iyr,
        nm,
        ds3,
    );

    *x = *x + 0.5 * (r11 + 4.0 * r41 + r51);
    *y = *y + 0.5 * (r12 + 4.0 * r42 + r52);
    *z = *z + 0.5 * (r13 + 4.0 * r43 + r53);
}

pub fn geocor(
    sla: f32,
    slo: f32,
    rh: f32,
    dla: &mut f32,
    dlo: &mut f32,
    cla: &mut f32,
    clo: &mut f32,
    pmi: &mut f32,
    iyr: i32,
    params: &DipoleParams,
) {
    if sla > 999.0 {
        *cla = 999.99;
        *clo = 999.99;
        *dla = 999.99;
        *dlo = 999.99;
        *pmi = 999.99;
        return;
    }

    let col = (90.0 - sla) * 0.017453293;
    let rlo = slo * 0.017453293;
    let mut r = rh;
    let mut col_mut = col;
    let mut rlo_mut = rlo;
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    sphcar(&mut r, &mut col_mut, &mut rlo_mut, &mut x, &mut y, &mut z, 1);

    let (mut xm, mut ym, mut zm) = (0.0, 0.0, 0.0);
    geomag(&mut x, &mut y, &mut z, &mut xm, &mut ym, &mut zm, 1, params);

    let (mut rm, mut th, mut pf) = (0.0, 0.0, 0.0);
    sphcar(&mut rm, &mut th, &mut pf, &mut xm, &mut ym, &mut zm, -1);

    let szm = zm;
    *dlo = pf * 57.2957751;
    let dco = th * 57.2957751;
    *dla = 90.0 - dco;
    let rl = r / (th.sin()).powi(2);
    let mut frac = 0.03 / (1.0 + 3.0 / (rl - 0.6));
    if szm < 0.0 {
        frac = -frac;
    }
    let hhh = 0.0001571;

    let mut outer_count = 0;
    'outer: loop {
        outer_count += 1;
        if outer_count > 500 {
            *cla = 999.99;
            *clo = 999.99;
            *pmi = 999.99;
            return;
        }
        let mut ds = r * frac;
        let mut inner_count = 0;
        'inner: loop {
            inner_count += 1;
            if inner_count > 500 {
                *cla = 999.99;
                *clo = 999.99;
                *pmi = 999.99;
                return;
            }
            let nm = ((1.0 + 9.0 / r) + 0.5) as i32;
            let r1 = r;
            let x1 = x;
            let y1 = y;
            let z1 = z;
            shag(&mut x, &mut y, &mut z, ds, iyr, nm);
            geomag(&mut x, &mut y, &mut z, &mut xm, &mut ym, &mut zm, 1, params);
            sphcar(&mut r, &mut col_mut, &mut rlo_mut, &mut xm, &mut ym, &mut zm, -1);

            if r > 10.0 + rh {
                break 'outer;
            }
            if r <= rh {
                *cla = 999.99;
                *clo = 999.99;
                *pmi = 999.99;
                return;
            }

            let dcl = col_mut - 1.5707963268;
            if dcl.abs() <= hhh {
                break 'outer;
            }
            let rzm = zm;
            if szm > 0.0 && rzm > 0.0 {
                continue 'outer;
            }
            if szm < 0.0 && rzm < 0.0 {
                continue 'outer;
            }

            r = r1;
            x = x1;
            y = y1;
            z = z1;
            ds /= 2.0;
        }
    }

    geomag(&mut x, &mut y, &mut z, &mut xm, &mut ym, &mut zm, 1, params);
    let mut gtet = 0.0;
    let mut gxla = 0.0;
    sphcar(&mut r, &mut gtet, &mut gxla, &mut xm, &mut ym, &mut zm, -1);
    let st = gtet.sin().abs();
    let rrh = (rh / (r - rh * st * st)).abs();
    let mut cla_val = 1.5707963 - (st * rrh.sqrt()).atan();
    cla_val *= 57.2957751;
    *clo = gxla * 57.2957751;
    if szm < 0.0 {
        cla_val = -cla_val;
    }
    *cla = cla_val;
    let ssla = (90.0 - cla_val) * 0.017453293;
    let sn = ssla.sin();
    *pmi = rh / (sn * sn);
}

pub fn corgeo(
    sla: &mut f32,
    slo: &mut f32,
    rh: f32,
    dla: &mut f32,
    dlo: &mut f32,
    cla: f32,
    clo: f32,
    pmi: &mut f32,
    iyr: i32,
    params: &DipoleParams,
) {
    if cla.abs() < 0.1 || cla > 999.0 {
        *sla = 999.99;
        *slo = 999.99;
        *dla = 999.99;
        *dlo = 999.99;
        *pmi = 999.99;
        return;
    }

    let col = (90.0 - cla) * 0.017453293;
    let mut r = 10.0;
    let r1_init = r;
    let mut r0 = r;
    let mut rlo = clo * 0.017453293;
    let sn = col.sin();
    let mut sn2 = sn * sn;

    if sn2 < 0.000000003 {
        sn2 = 0.000000003;
    }

    let rfi = rh / sn2;
    *pmi = rfi;
    if *pmi > 99.999 {
        *pmi = 999.99;
    }
    let aa10 = r / rfi;

    let mut scla = 0.0_f32;
    if rfi <= r {
        scla = 1.57079632679;
        r0 = rfi;
    } else {
        let saa = aa10 / (1.0 - aa10);
        let saq = saa.sqrt();
        scla = saq.atan();
        if cla < 0.0 {
            scla = 3.14159265359 - scla;
        }
        r0 = r1_init;
    }

    let mut xm = 0.0;
    let mut ym = 0.0;
    let mut zm = 0.0;
    sphcar(&mut r0, &mut scla, &mut rlo, &mut xm, &mut ym, &mut zm, 1);

    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    geomag(&mut x, &mut y, &mut z, &mut xm, &mut ym, &mut zm, -1, params);

    let rl = r0;
    let mut frac = -0.03 / (1.0 + 3.0 / (rl - 0.6));
    if cla < 0.0 {
        frac = -frac;
    }
    r = r0;

    let mut r1 = r;
    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;

    let mut iter_count = 0;
    loop {
        let ds = r * frac;
        let nm = ((1.0 + 9.0 / r) + 0.5) as i32;
        shag(&mut x, &mut y, &mut z, ds, iyr, nm);
        r = (x * x + y * y + z * z).sqrt();
        iter_count += 1;
        if iter_count > 200 {
            panic!("Loop limit exceeded in corgeo tracing! r={}, rh={}, ds={}", r, rh, ds);
        }
        if r <= rh {
            break;
        }
        r1 = r;
        x1 = x;
        y1 = y;
        z1 = z;
    }

    let dr1 = (rh - r1).abs();
    let dr0 = (rh - r).abs();
    let dr10 = dr1 + dr0;
    let mut ds = r1 * frac;
    if dr10 != 0.0 {
        ds = ds * (dr1 / dr10);
        let nm = ((1.0 + 9.0 / r1) + 0.5) as i32;
        shag(&mut x1, &mut y1, &mut z1, ds, iyr, nm);
    }

    let mut r_final = 0.0;
    let mut gtet = 0.0;
    let mut gxla = 0.0;
    sphcar(&mut r_final, &mut gtet, &mut gxla, &mut x1, &mut y1, &mut z1, -1);

    let gth = gtet * 57.2957751;
    *slo = gxla * 57.2957751;
    *sla = 90.0 - gth;

    let mut xm_final = 0.0;
    let mut ym_final = 0.0;
    let mut zm_final = 0.0;
    geomag(&mut x1, &mut y1, &mut z1, &mut xm_final, &mut ym_final, &mut zm_final, 1, params);

    let mut rm_final = 0.0;
    let mut th_final = 0.0;
    let mut pf_final = 0.0;
    sphcar(&mut rm_final, &mut th_final, &mut pf_final, &mut xm_final, &mut ym_final, &mut zm_final, -1);

    *dlo = pf_final * 57.2957751;
    *dla = 90.0 - th_final * 57.2957751;

    if sla.abs() < 30.0 || cla.abs() < 30.0 {
        let mut dls_dummy1 = 0.0;
        let mut dls_dummy2 = 0.0;
        let mut clas = 0.0;
        let mut clos = 0.0;
        let mut pms = 0.0;
        geocor(*sla, *slo, rh, &mut dls_dummy1, &mut dls_dummy2, &mut clas, &mut clos, &mut pms, iyr, params);

        if clas > 999.0 {
            let mut rbm_dummy = 0.0;
            let mut slac_dummy = 0.0;
            let mut sloc_dummy = 0.0;
            geolow(*sla, *slo, rh, &mut clas, &mut clos, &mut rbm_dummy, &mut slac_dummy, &mut sloc_dummy, iyr, params);
        }

        if (cla.abs() - clas.abs()).abs() >= 1.0 {
            *sla = 999.99;
            *slo = 999.99;
            *pmi = 999.99;
        }
    }
}

pub fn ftprnt(
    rh: f32,
    sla: f32,
    slo: f32,
    cla: f32,
    clo: f32,
    acla: &mut f32,
    aclo: &mut f32,
    slaf: &mut f32,
    slof: &mut f32,
    rf: f32,
    iyr: i32,
    params: &DipoleParams,
) {
    if sla > 999.0 || cla > 999.0 || rf == rh {
        *acla = 999.99;
        *aclo = 999.99;
        *slaf = 999.99;
        *slof = 999.99;
        return;
    }

    let col = (90.0 - cla) * 0.017453293;
    let mut sn2 = col.sin().powi(2);
    let mut decarg = ((sn2 * rf) / rh).sqrt();
    if decarg.abs() > 1.0 {
        decarg = decarg.signum();
    }
    let acol = decarg.asin();
    let mut acla_val = 90.0 - acol * 57.29577951;
    if cla < 0.0 {
        acla_val = -acla_val;
    }
    *acla = acla_val;
    *aclo = clo;

    let mut dlaf = 0.0;
    let mut dlof = 0.0;
    let mut pmif = 0.0;
    corgeo(slaf, slof, rf, &mut dlaf, &mut dlof, *acla, *aclo, &mut pmif, iyr, params);

    if *slaf < 999.0 {
        return;
    }

    if sn2 < 0.0000001 {
        sn2 = 0.0000001;
    }
    let rl = rh / sn2;
    let mut frac = 0.03 / (1.0 + 3.0 / (rl - 0.6));
    if cla >= 0.0 {
        frac = -frac;
    }
    let mut ds = rh * frac;

    let mut start_count = 0;
    'start: loop {
        start_count += 1;
        if start_count > 10 {
            *acla = 999.99;
            *aclo = 999.99;
            *slaf = 999.99;
            *slof = 999.99;
            return;
        }
        let mut r = rh;
        let mut rsla = (90.0 - sla) * 0.0174533;
        let mut rslo = slo * 0.0174533;
        let (mut xf, mut yf, mut zf) = (0.0, 0.0, 0.0);
        sphcar(&mut r, &mut rsla, &mut rslo, &mut xf, &mut yf, &mut zf, 1);
        let mut rf1 = r;
        let mut xf1 = xf;
        let mut yf1 = yf;
        let mut zf1 = zf;

        let mut inner_count = 0;
        loop {
            inner_count += 1;
            if inner_count > 500 {
                *acla = 999.99;
                *aclo = 999.99;
                *slaf = 999.99;
                *slof = 999.99;
                return;
            }
            let nm = ((1.0 + 9.0 / r) + 0.5) as i32;
            shag(&mut xf, &mut yf, &mut zf, ds, iyr, nm);
            let rr = (xf * xf + yf * yf + zf * zf).sqrt();
            if rr > rh {
                ds = -ds;
                continue 'start;
            }
            if rr > rf {
                rf1 = rr;
                xf1 = xf;
                yf1 = yf;
                zf1 = zf;
            } else {
                let dr1 = (rf1 - rf).abs();
                let dr0 = (rf - rr).abs();
                let dr10 = dr1 + dr0;
                if dr10 != 0.0 {
                    ds = ds * (dr1 / dr10);
                    let nm2 = ((1.0 + 9.0 / rf1) + 0.5) as i32;
                    shag(&mut xf1, &mut yf1, &mut zf1, ds, iyr, nm2);
                }
                let mut rr_mut = rr;
                sphcar(&mut rr_mut, slaf, slof, &mut xf1, &mut yf1, &mut zf1, -1);
                *slaf = 90.0 - *slaf * 57.29578;
                *slof = *slof * 57.29578;
                return;
            }
        }
    }
}

pub fn geolow(
    slar: f32,
    slor: f32,
    rh: f32,
    clar: &mut f32,
    clor: &mut f32,
    rbm: &mut f32,
    slac: &mut f32,
    sloc: &mut f32,
    iyr: i32,
    params: &DipoleParams,
) {
    if slar > 999.0 {
        *clar = 999.99;
        *clor = 999.99;
        *slac = 999.99;
        *sloc = 999.99;
        *rbm = 999.99;
        return;
    }

    let dhh = 0.5;
    let mut arlat = [999.99_f32; 182];
    let mut arlon = [999.99_f32; 182];
    let slo = slor;

    let mut jcn = 0;
    let mut jcs = 0;
    let mut rnlat = 999.99_f32;
    let mut rnlon = 999.99_f32;
    let mut rslat = 999.99_f32;
    let mut rslon = 999.99_f32;

    let mut ndir = 0;
    'dir_loop: loop {
        if ndir == 0 {
            for jc in 61..=91 {
                let sla = 90.0 - (jc - 1) as f32;
                let (mut dla, mut dlo) = (0.0, 0.0);
                let (mut cla_val, mut clo_val, mut pmm) = (0.0, 0.0, 0.0);
                geocor(sla, slo, rh, &mut dla, &mut dlo, &mut cla_val, &mut clo_val, &mut pmm, iyr, params);
                if cla_val > 999.0 {
                    ndir = 1;
                    continue 'dir_loop;
                }
                arlat[jc] = cla_val;
                arlon[jc] = clo_val;
            }
            ndir = 1;
        } else {
            for jc in (92..=121).rev() {
                let sla = 90.0 - (jc - 1) as f32;
                let (mut dla, mut dlo) = (0.0, 0.0);
                let (mut cla_val, mut clo_val, mut pmm) = (0.0, 0.0, 0.0);
                geocor(sla, slo, rh, &mut dla, &mut dlo, &mut cla_val, &mut clo_val, &mut pmm, iyr, params);
                if cla_val > 999.0 {
                    ndir = 0;
                    break 'dir_loop;
                }
                arlat[jc] = cla_val;
                arlon[jc] = clo_val;
            }
            ndir = 0;
            break;
        }
    }

    let mut n999 = 0;
    ndir = 0;
    for jc in 61..=121 {
        if arlat[jc] > 999.0 {
            if ndir == 0 {
                jcn = jc - 1;
                rnlat = arlat[jcn];
                rnlon = arlon[jcn];
                ndir = 1;
                n999 = 1;
            }
        }
        if arlat[jc] < 999.0 {
            if ndir == 1 {
                jcs = jc;
                rslat = arlat[jc];
                rslon = arlon[jc];
                ndir = 0;
                break;
            }
        }
    }

    let ih = if n999 == 0 { 3 } else { 1 };

    if n999 != 0 {
        let rdel = (jcs - jcn) as f32;
        let mut delon = 0.0;
        if rdel != 0.0 {
            if rslon > 270.0 && rnlon < 90.0 {
                delon = (rslon - (rnlon + 360.0)) / rdel;
            } else if rslon < 90.0 && rnlon > 270.0 {
                delon = (rslon - (rnlon - 360.0)) / rdel;
            } else {
                delon = (rslon - rnlon) / rdel;
            }
        }
        for jc in jcn + 1..jcs {
            arlon[jc] = rnlon + delon * (jc - jcn) as f32;
            if arlon[jc] < 0.0 {
                arlon[jc] += 360.0;
            }
        }
    }

    let mut rlan = 999.99_f32;
    let mut rlas = 999.99_f32;
    let mut nobm = 0;

    for ihem in ih..=3 {
        let mut rm = rh;
        let mut cla_val = 0.0;
        let mut sla = 0.0;

        if ihem == 1 {
            cla_val = rnlat;
            sla = 90.0 - (jcn - 1) as f32;
        } else if ihem == 2 {
            cla_val = rslat;
            sla = 90.0 - (jcs - 1) as f32;
        } else if ihem == 3 {
            cla_val = 0.0;
            sla = slar;
        }

        let col = (90.0 - cla_val) * 0.017453293;
        let mut slm = (90.0 - sla) * 0.017453293;
        let mut sll = slo * 0.017453293;

        let nm = 10;
        let mut res = crate::igrf::igrf(iyr, nm, rm, slm, sll);
        let mut sz = -res.br;

        let (mut xgeo, mut ygeo, mut zgeo) = (0.0, 0.0, 0.0);
        sphcar(&mut rm, &mut slm, &mut sll, &mut xgeo, &mut ygeo, &mut zgeo, 1);

        let mut bm = (res.br * res.br + res.bt * res.bt + res.bf * res.bf).sqrt();
        let mut xbm = xgeo;
        let mut ybm = ygeo;
        let mut zbm = zgeo;

        let rl = 1.0 / (col.sin()).powi(2);
        let mut frac = 0.03 / (1.0 + 3.0 / (rl - 0.6));
        if sz <= 0.0 {
            frac = -frac;
        }
        let dsd = rl * frac;
        let mut ds = dsd;

        let mut xbm1 = 0.0;
        let mut ybm1 = 0.0;
        let mut zbm1 = 0.0;
        let mut rbm1 = 0.0;

        let mut bmin_count = 0;
        'bmin_loop: loop {
            bmin_count += 1;
            if bmin_count > 500 {
                nobm = 1;
                break 'bmin_loop;
            }
            let mut bc = [0.0_f32; 2];
            for i in 0..2 {
                let mut dd = ds;
                shag(&mut xgeo, &mut ygeo, &mut zgeo, dd, iyr, nm);
                if i == 0 {
                    xbm1 = xgeo;
                    ybm1 = ygeo;
                    zbm1 = zgeo;
                    rbm1 = (xbm1 * xbm1 + ybm1 * ybm1 + zbm1 * zbm1).sqrt();
                }

                sphcar(&mut rm, &mut slm, &mut sll, &mut xgeo, &mut ygeo, &mut zgeo, -1);
                res = crate::igrf::igrf(iyr, nm, rm, slm, sll);

                if rm < rh {
                    nobm = 1;
                    break 'bmin_loop;
                }
                bc[i] = (res.br * res.br + res.bt * res.bt + res.bf * res.bf).sqrt();
            }

            let b2 = bc[0];
            let b3 = bc[1];

            if (bm > b2 && b2 < b3) || (bm >= b2 && b2 < b3) || (bm > b2 && b2 <= b3) {
                let bb3 = (b3 - b2).abs();
                let bb2 = (bm - b2).abs();
                if bb2 < dhh && bb3 < dhh {
                    break 'bmin_loop;
                }
                xgeo = xbm;
                ygeo = ybm;
                zgeo = zbm;
                ds /= 2.0;
                continue 'bmin_loop;
            }

            bm = bc[0];
            xgeo = xbm1;
            ygeo = ybm1;
            zgeo = zbm1;
            xbm = xbm1;
            ybm = ybm1;
            zbm = zbm1;
        }

        if nobm == 0 {
            let mut rla_temp = 0.0;
            let mut rlo_temp = 0.0;
            let mut rbm1_temp = rbm1;
            sphcar(&mut rbm1_temp, &mut rla_temp, &mut rlo_temp, &mut xbm1, &mut ybm1, &mut zbm1, -1);
            let rla = 90.0 - rla_temp * 57.2957751;
            let rlo = rlo_temp * 57.2957751;

            if ihem == 1 {
                rlan = rla;
            }
            if ihem == 2 {
                rlas = rla;
            }

            if ihem == 3 {
                *rbm = rbm1;
                let mut rm_trace = *rbm;
                let mut ds_trace = dsd;
                let mut r1 = rm_trace;
                let mut x1 = xbm1;
                let mut y1 = ybm1;
                let mut z1 = zbm1;

                let mut trace_count = 0;
                loop {
                    trace_count += 1;
                    if trace_count > 500 {
                        *slac = 999.99;
                        *sloc = 999.99;
                        break;
                    }
                    let nm_low = ((1.0 + 9.0 / rm_trace) + 0.5) as i32;
                    shag(&mut xbm1, &mut ybm1, &mut zbm1, ds_trace, iyr, nm_low);
                    let rr = (xbm1 * xbm1 + ybm1 * ybm1 + zbm1 * zbm1).sqrt();

                    if rr > rh {
                        r1 = rr;
                        x1 = xbm1;
                        y1 = ybm1;
                        z1 = zbm1;
                    } else {
                        let dr1 = (rh - r1).abs();
                        let dr0 = (rh - rr).abs();
                        let dr10 = dr1 + dr0;
                        if dr10 != 0.0 {
                            ds_trace = ds_trace * (dr1 / dr10);
                            let nm_low2 = ((1.0 + 9.0 / r1) + 0.5) as i32;
                            shag(&mut x1, &mut y1, &mut z1, ds_trace, iyr, nm_low2);
                        }
                        let mut rr_mut = rr;
                        sphcar(&mut rr_mut, slac, sloc, &mut x1, &mut y1, &mut z1, -1);
                        *slac = 90.0 - *slac * 57.29578;
                        *sloc = *sloc * 57.29578;
                        break;
                    }
                }
            }
        }
    }

    if n999 == 0 {
        return;
    }

    if nobm == 1 {
        let rdel = (jcs - jcn) as f32;
        let delat = if rdel == 0.0 { 0.0 } else { (rslat - rnlat) / rdel };
        let mut jdel = 0;
        for jc in jcn + 1..jcs {
            jdel += 1;
            arlat[jc] = rnlat + delat * jdel as f32;
        }
        *rbm = 999.99;
        *slac = 999.99;
        *sloc = 999.99;
    } else {
        let rla = (rlan + rlas) / 2.0;

        let rdel_n = (90.0 - (jcn - 1) as f32) - rla;
        let delat_n = if rdel_n == 0.0 { 0.0 } else { rnlat / rdel_n };
        let jdn = rdel_n.abs() as usize;
        let mut jdel = 0;
        for jc in jcn + 1..=jcn + jdn {
            jdel += 1;
            arlat[jc] = rnlat - delat_n * jdel as f32;
        }

        let rdel_s = (90.0 - (jcs - 1) as f32) - rla;
        let delat_s = if rdel_s == 0.0 { 0.0 } else { rslat / rdel_s };
        let jds = rdel_s.abs() as usize;
        let mut jdel = 0;
        for jc in (jcs - jds..jcs).rev() {
            jdel += 1;
            arlat[jc] = rslat + delat_s * jdel as f32;
        }
    }

    let l1 = (90.0 - slar + 1.0) as usize;
    let l2 = if slar < 0.0 { l1 - 1 } else { l1 + 1 };
    let dsla = (slar - slar.trunc()).abs();
    let delcla = arlat[l2] - arlat[l1];
    let delclo = arlon[l2] - arlon[l1];
    *clar = arlat[l1] + delcla * dsla;
    *clor = arlon[l1] + delclo * dsla;
}

pub fn dfridr<F>(func: F, x: f32, h: f32, err: &mut f32) -> f32
where
    F: Fn(f32) -> f32,
{
    if h == 0.0 {
        return 0.0;
    }
    const CON: f32 = 1.4;
    const CON2: f32 = CON * CON;
    const BIG: f32 = 1e30;
    const SAFE: f32 = 2.0;

    let mut a = [[0.0_f32; 10]; 10];
    let mut hh = h;
    a[0][0] = (func(x + hh) - func(x - hh)) / (2.0 * hh);
    *err = BIG;
    let mut dfridr_val = a[0][0];

    for i in 1..10 {
        hh /= CON;
        a[0][i] = (func(x + hh) - func(x - hh)) / (2.0 * hh);
        let mut fac = CON2;
        for j in 1..=i {
            a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
            fac *= CON2;
            let errt = (a[j][i] - a[j - 1][i])
                .abs()
                .max((a[j][i] - a[j - 1][i - 1]).abs());
            if errt <= *err {
                *err = errt;
                dfridr_val = a[j][i];
            }
        }
        if (a[i][i] - a[i - 1][i - 1]).abs() >= SAFE * (*err) {
            return dfridr_val;
        }
    }
    dfridr_val
}

pub fn ovl_ang(
    sla: f32,
    slo: f32,
    cla: f32,
    clo: f32,
    rr: f32,
    iyr: i32,
    params: &DipoleParams,
) -> f32 {
    if sla.abs() >= 89.99 || cla.abs() >= 89.99 || sla.abs() < 30.0 {
        return 999.99;
    }

    let cr360 = slo >= 270.0;
    let cr0 = slo <= 90.0;
    let step = 10.0;

    let cgmgla_closure = |clon: f32| {
        let mut clon_val = clon;
        if clon_val > 360.0 {
            clon_val -= 360.0;
        }
        if clon_val < 0.0 {
            clon_val += 360.0;
        }
        let (mut geolat, mut geolon) = (0.0, 0.0);
        let (mut dla, mut dlo, mut pmi) = (0.0, 0.0, 0.0);
        corgeo(
            &mut geolat,
            &mut geolon,
            rr,
            &mut dla,
            &mut dlo,
            cla,
            clon_val,
            &mut pmi,
            iyr,
            params,
        );
        geolat
    };

    let cgmglo_closure = |clon: f32| {
        let mut clon_val = clon;
        loop {
            if clon_val > 360.0 {
                clon_val -= 360.0;
            }
            if clon_val < 0.0 {
                clon_val += 360.0;
            }
            let (mut geolat, mut geolon) = (0.0, 0.0);
            let (mut dla, mut dlo, mut pmi) = (0.0, 0.0, 0.0);
            corgeo(
                &mut geolat,
                &mut geolon,
                rr,
                &mut dla,
                &mut dlo,
                cla,
                clon_val,
                &mut pmi,
                iyr,
                params,
            );

            if geolat.abs() >= 89.99 {
                clon_val -= 0.01;
                continue;
            }

            if cr360 && geolon <= 90.0 {
                return geolon + 360.0;
            } else if cr0 && geolon >= 270.0 {
                return geolon - 360.0;
            } else {
                return geolon;
            }
        }
    };

    let mut err1 = 0.0;
    let hom = dfridr(cgmgla_closure, clo, step, &mut err1);

    let mut err2 = 0.0;
    let mut denom = dfridr(cgmglo_closure, clo, step, &mut err2);

    denom = denom * (sla * 0.017453293).cos();

    let mut ovl = -hom.atan2(denom);
    ovl *= 57.2957751;
    ovl
}

pub fn azm_ang(sla: f32, slo: f32, cla: f32, pla: f32, plo: f32) -> f32 {
    if sla.abs() >= 89.99 || cla.abs() >= 89.99 {
        return 999.99;
    }

    let rad = 0.017453293;
    let am = (90.0 - pla.abs()) * rad;
    let cm = if pla.signum() == sla.signum() {
        (90.0 - sla.abs()) * rad
    } else {
        (90.0 + sla.abs()) * rad
    };

    let bet = if sla >= 0.0 {
        (plo - slo) * rad
    } else {
        (slo - plo) * rad
    };

    let sb = bet.sin();
    let st = cm.sin() / am.tan() - cm.cos() * bet.cos();
    let alfa = sb.atan2(st);
    alfa / rad
}

pub fn mltut(sla: f32, slo: f32, cla: f32, pla: f32, plo: f32, ut: &mut f32) {
    if sla.abs() >= 89.99 || cla.abs() >= 89.99 {
        *ut = 99.99;
        return;
    }

    let tpi = 6.283185307;
    let rad = 0.017453293;

    let qq = plo * rad;
    let mut cff = 90.0 - pla.abs();
    cff *= rad;
    if cff < 0.0000001 {
        cff = 0.0000001;
    }

    let mut cft = if pla.signum() == sla.signum() {
        90.0 - sla.abs()
    } else {
        90.0 + sla.abs()
    };
    cft *= rad;
    if cft < 0.0000001 {
        cft = 0.0000001;
    }

    let qt = slo * rad;
    let a = cff.sin() / cft.sin();
    let y = a * qq.sin() - qt.sin();
    let x = qt.cos() - a * qq.cos();
    let mut ut_val = y.atan2(x);

    if ut_val < 0.0 {
        ut_val += tpi;
    }
    let qqu = qq + ut_val;
    let qtu = qt + ut_val;
    let bp = cff.sin() * qqu.cos();
    let bt = cft.sin() * qtu.cos();
    ut_val = ut_val / rad;
    ut_val = ut_val / 15.0;

    if bp >= bt {
        if ut_val < 12.0 {
            ut_val += 12.0;
        }
        if ut_val > 12.0 {
            ut_val -= 12.0;
        }
    }
    *ut = ut_val;
}

pub fn mfc(sla: f32, slo: f32, r: f32, h: &mut f32, d: &mut f32, z: &mut f32, iyr: i32, nm: i32) {
    if sla >= 999.0 {
        *h = 99999.0;
        *d = 999.99;
        *z = 99999.0;
        return;
    }
    let rla = (90.0 - sla) * 0.017453293;
    let rlo = slo * 0.017453293;

    let res = crate::igrf::igrf(iyr, nm, r, rla, rlo);
    let x = -res.bt;
    let y = res.bf;
    *z = -res.br;
    *h = (x * x + y * y).sqrt();
    *d = 57.2957751 * y.atan2(x);
}

// Flat array index mapping helper
macro_rules! dat {
    ($i:expr, $j:expr) => {
        (($i - 1) + ($j - 1) * 11)
    };
}

pub fn geocgm01(
    icor: i32,
    iyear: i32,
    hi: f32,
    dat: &mut [f32; 44],
    pla: &mut [f32; 4],
    plo: &mut [f32; 4],
) {
    let params = recalc(iyear);
    let re = 6371.2;
    let rh = (re + hi) / re;

    if dat[dat!(1, 1)] > 90.0 {
        dat[dat!(1, 1)] = 180.0 - dat[dat!(1, 1)];
    }
    if dat[dat!(1, 1)] < -90.0 {
        dat[dat!(1, 1)] = -180.0 - dat[dat!(1, 1)];
    }
    if dat[dat!(3, 1)] > 90.0 {
        dat[dat!(3, 1)] = 180.0 - dat[dat!(3, 1)];
    }
    if dat[dat!(3, 1)] < -90.0 {
        dat[dat!(3, 1)] = -180.0 - dat[dat!(3, 1)];
    }

    if dat[dat!(2, 1)] > 360.0 {
        dat[dat!(2, 1)] -= 360.0;
    }
    if dat[dat!(2, 1)] < -360.0 {
        dat[dat!(2, 1)] += 360.0;
    }
    if dat[dat!(4, 1)] > 360.0 {
        dat[dat!(4, 1)] -= 360.0;
    }
    if dat[dat!(4, 1)] < -360.0 {
        dat[dat!(4, 1)] += 360.0;
    }

    let mut slar = dat[dat!(1, 1)];
    let mut slor = dat[dat!(2, 1)];
    let mut clar = dat[dat!(3, 1)];
    let mut clor = dat[dat!(4, 1)];
    let mut pmr = 0.0;
    let (mut dla, mut dlo) = (0.0, 0.0);

    if icor == 1 {
        if slar.abs() == 90.0 {
            slor = 360.0;
        }
        geocor(slar, slor, rh, &mut dla, &mut dlo, &mut clar, &mut clor, &mut pmr, iyear, &params);
        dat[dat!(3, 1)] = clar;
        dat[dat!(4, 1)] = clor;
    } else {
        if clar.abs() == 90.0 {
            clor = 360.0;
        }
        corgeo(&mut slar, &mut slor, rh, &mut dla, &mut dlo, clar, clor, &mut pmr, iyear, &params);
        dat[dat!(1, 1)] = slar;
        dat[dat!(2, 1)] = slor;
    }

    if pmr >= 16.0 {
        pmr = 999.99;
    }
    dat[dat!(5, 1)] = pmr;

    let mut slac = 0.0;
    let mut sloc = 0.0;
    let mut clac = 0.0;
    let mut cloc = 0.0;
    let mut rbm = 0.0;

    if clar > 999.0 {
        geolow(slar, slor, rh, &mut clar, &mut clor, &mut rbm, &mut slac, &mut sloc, iyear, &params);
        dat[dat!(3, 1)] = clar;
        dat[dat!(4, 1)] = clor;
        if rbm >= 16.0 {
            rbm = 999.99;
        }
        dat[dat!(5, 1)] = rbm;

        // Low latitude conjugate points rounding to matching Fortran precision float representation
        let s_slac = format!("{:6.2}", slac);
        let s_sloc = format!("{:6.2}", sloc);
        slac = s_slac.trim().parse::<f32>().unwrap_or(slac);
        sloc = s_sloc.trim().parse::<f32>().unwrap_or(sloc);

        dat[dat!(1, 2)] = slac;
        dat[dat!(2, 2)] = sloc;
        let mut daa = 0.0;
        let mut doo = 0.0;
        geocor(slac, sloc, rh, &mut daa, &mut doo, &mut clac, &mut cloc, &mut rbm, iyear, &params);
        if clac > 999.0 {
            let mut slal = 0.0;
            let mut slol = 0.0;
            geolow(slac, sloc, rh, &mut clac, &mut cloc, &mut rbm, &mut slal, &mut slol, iyear, &params);
        }
        dat[dat!(3, 2)] = clac;
        dat[dat!(4, 2)] = cloc;
        dat[dat!(5, 2)] = rbm;
    } else {
        clac = -clar;
        cloc = clor;
        dat[dat!(3, 2)] = clac;
        dat[dat!(4, 2)] = cloc;
        let mut daa = 0.0;
        let mut doo = 0.0;
        let mut pmc = 0.0;
        corgeo(&mut slac, &mut sloc, rh, &mut daa, &mut doo, clac, cloc, &mut pmc, iyear, &params);
        dat[dat!(1, 2)] = slac;
        dat[dat!(2, 2)] = sloc;
        if pmc >= 16.0 {
            pmc = 999.99;
        }
        dat[dat!(5, 2)] = pmc;
    }

    dat[dat!(5, 3)] = dat[dat!(5, 1)];
    dat[dat!(5, 4)] = dat[dat!(5, 2)];

    let mut aclar = 0.0;
    let mut aclor = 0.0;
    let mut aclac = 0.0;
    let mut acloc = 0.0;

    if rh > 1.0 && clar < 999.0 && clar < 999.0 {
        let mut slarf = 0.0;
        let mut slorf = 0.0;
        ftprnt(
            rh,
            slar,
            slor,
            clar,
            clor,
            &mut aclar,
            &mut aclor,
            &mut slarf,
            &mut slorf,
            1.0,
            iyear,
            &params,
        );
        dat[dat!(1, 3)] = slarf;
        dat[dat!(2, 3)] = slorf;
        dat[dat!(3, 3)] = aclar;
        dat[dat!(4, 3)] = aclor;

        let mut slacf = 0.0;
        let mut slocf = 0.0;
        ftprnt(
            rh,
            slac,
            sloc,
            clac,
            cloc,
            &mut aclac,
            &mut acloc,
            &mut slacf,
            &mut slocf,
            1.0,
            iyear,
            &params,
        );
        dat[dat!(1, 4)] = slacf;
        dat[dat!(2, 4)] = slocf;
        dat[dat!(3, 4)] = aclac;
        dat[dat!(4, 4)] = acloc;
    } else {
        for i in 1..=4 {
            for j in 3..=4 {
                dat[dat!(i, j)] = 999.99;
            }
        }
    }

    let mut plan = 0.0;
    let mut plon = 0.0;
    let mut daa_dummy = 0.0;
    let mut doo_dummy = 0.0;
    let mut pmp = 0.0;
    corgeo(&mut plan, &mut plon, rh, &mut daa_dummy, &mut doo_dummy, 90.0, 360.0, &mut pmp, iyear, &params);
    let mut plan1 = plan;
    let mut plon1 = plon;

    let mut plas = 0.0;
    let mut plos = 0.0;
    corgeo(&mut plas, &mut plos, rh, &mut daa_dummy, &mut doo_dummy, -90.0, 360.0, &mut pmp, iyear, &params);
    let mut plas1 = plas;
    let mut plos1 = plos;

    if rh > 1.0 {
        corgeo(&mut plan1, &mut plon1, 1.0, &mut daa_dummy, &mut doo_dummy, 90.0, 360.0, &mut pmp, iyear, &params);
        let mut pmm = 0.0;
        corgeo(&mut plas1, &mut plos1, 1.0, &mut daa_dummy, &mut doo_dummy, -90.0, 360.0, &mut pmm, iyear, &params);
    }

    if clar < 0.0 {
        pla[0] = plas;
        plo[0] = plos;
    } else {
        pla[0] = plan;
        plo[0] = plon;
    }

    if aclar < 0.0 {
        pla[2] = plas1;
        plo[2] = plos1;
    } else {
        pla[2] = plan1;
        plo[2] = plon1;
    }

    if clac < 0.0 {
        pla[1] = plas;
        plo[1] = plos;
    } else {
        pla[1] = plan;
        plo[1] = plon;
    }

    if aclac < 0.0 {
        pla[3] = plas1;
        plo[3] = plos1;
    } else {
        pla[3] = plan1;
        plo[3] = plon1;
    }

    for j in 1..=4 {
        dat[dat!(6, j)] = 99999.0;
        dat[dat!(7, j)] = 999.99;
        dat[dat!(8, j)] = 99999.0;
        dat[dat!(9, j)] = 999.99;
        dat[dat!(10, j)] = 999.99;
        dat[dat!(11, j)] = 99.99;
    }

    let icount = if rh > 1.0 { 4 } else { 2 };
    let mut rj = rh;
    for j in 1..=icount {
        if j > 2 {
            rj = 1.0;
        }

        let plaj = pla[j - 1];
        let ploj = plo[j - 1];

        let slaj = dat[dat!(1, j)];
        let sloj = dat[dat!(2, j)];
        let claj = dat[dat!(3, j)];
        let cloj = dat[dat!(4, j)];

        let mut btr = 0.0;
        let mut bfr = 0.0;
        let mut brr = 0.0;
        mfc(slaj, sloj, rj, &mut btr, &mut bfr, &mut brr, iyear, 10);
        dat[dat!(6, j)] = btr;
        dat[dat!(7, j)] = bfr;
        dat[dat!(8, j)] = brr;

        let ovl = ovl_ang(slaj, sloj, claj, cloj, rj, iyear, &params);
        dat[dat!(9, j)] = ovl;

        let azm = azm_ang(slaj, sloj, claj, plaj, ploj);
        dat[dat!(10, j)] = azm;

        let mut ut = 0.0;
        mltut(slaj, sloj, claj, plaj, ploj, &mut ut);
        dat[dat!(11, j)] = ut;
    }
}
