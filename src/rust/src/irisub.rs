use crate::irifun_utils::{PI, UMR, ARGMAX};

// --- Invariant Latitude / Invariant Dip Latitude Helpers ---

pub fn invdpc(fl: f32, dimo: f32, b0: f32, dipl: f32) -> f32 {
    let b = [
        1.259921_f64,
        -0.1984259_f64,
        -0.04686632_f64,
        -0.01314096_f64,
        -0.00308824_f64,
        0.00082777_f64,
        -0.00105877_f64,
        0.00183142_f64,
    ];
    let a = ((dimo / b0).powf(1.0 / 3.0) / fl) as f64;
    let mut asa = a * (b[0]
        + b[1] * a
        + b[2] * a.powi(2)
        + b[3] * a.powi(3)
        + b[4] * a.powi(4)
        + b[5] * a.powi(5)
        + b[6] * a.powi(6)
        + b[7] * a.powi(7));
    if asa > 1.0 {
        asa = 1.0;
    }
    if asa < 0.0 {
        asa = 0.0;
    }
    let rinvl = (asa.sqrt()).acos() as f32;
    let invl = rinvl / UMR;
    let alfa = 2.0
        - (1.0 / (((dipl - 25.0) / 2.0).exp() + 1.0)
            + 1.0 / (((-dipl - 25.0) / 2.0).exp() + 1.0));
    let beta = 1.0 / (((invl - 25.0) / 2.0).exp() + 1.0)
        + 1.0 / (((-invl - 25.0) / 2.0).exp() + 1.0)
        - 1.0;
    (alfa * dipl.signum() * invl + beta * dipl) / (alfa + beta)
}

pub fn invdpc_old(fl: f32, dimo: f32, b0: f32, dipl: f32) -> f32 {
    let b = [
        1.259921_f64,
        -0.1984259_f64,
        -0.04686632_f64,
        -0.01314096_f64,
        -0.00308824_f64,
        0.00082777_f64,
        -0.00105877_f64,
        0.00183142_f64,
    ];
    let a = ((dimo / b0).powf(1.0 / 3.0) / fl) as f64;
    let mut asa = a * (b[0]
        + b[1] * a
        + b[2] * a.powi(2)
        + b[3] * a.powi(3)
        + b[4] * a.powi(4)
        + b[5] * a.powi(5)
        + b[6] * a.powi(6)
        + b[7] * a.powi(7));
    if asa > 1.0 {
        asa = 1.0;
    }
    if asa < 0.0 {
        asa = 0.0;
    }
    let rinvl = (asa.sqrt()).acos() as f32;
    let invl = rinvl / UMR;
    let rdipl = dipl * UMR;
    let alfa = rdipl.abs().sin().powi(3);
    let beta = rinvl.cos().powi(3);
    (alfa * dipl.signum() * invl + beta * dipl) / (alfa + beta)
}

// --- SOCO (Solar Zenith Solver) ---

pub fn soco(
    ld: i32,
    t: f32,
    flat: f32,
    elon: f32,
    height: f32,
) -> (f32, f32, f32, f32) {
    let p1 = 0.017203534_f32;
    let p2 = 0.034407068_f32;
    let p3 = 0.051610602_f32;
    let p4 = 0.068814136_f32;
    let p6 = 0.103221204_f32;

    let humr = PI / 12.0;
    let dtr = UMR;

    let wlon = 360.0 - elon;
    let td = (ld as f32) + (t + wlon / 15.0) / 24.0;
    let te = td + 0.9369;

    let dcl = 23.256 * (p1 * (te - 82.242)).sin()
        + 0.381 * (p2 * (te - 44.855)).sin()
        + 0.167 * (p3 * (te - 23.355)).sin()
        - 0.013 * (p4 * (te + 11.97)).sin()
        + 0.011 * (p6 * (te - 10.41)).sin()
        + 0.339137;

    let declin = dcl;
    let dc = dcl * dtr;

    let tf = te - 0.5;
    let eqt = -7.38 * (p1 * (tf - 4.0)).sin()
        - 9.87 * (p2 * (tf + 9.0)).sin()
        + 0.27 * (p3 * (tf - 53.0)).sin()
        - 0.2 * (p4 * (tf - 17.0)).cos();
    let mut et = eqt * dtr / 4.0;

    let fa = flat * dtr;
    let phi = humr * (t - 12.0) + et;

    let a = fa.sin() * dc.sin();
    let b = fa.cos() * dc.cos();
    let mut cosx = a + b * phi.cos();
    if cosx.abs() > 1.0 {
        cosx = cosx.signum();
    }
    let zenith = cosx.acos() / dtr;

    let h = height * 1000.0;
    let chih = 90.83 + 0.0347 * h.sqrt();
    let ch = (chih * dtr).cos();
    let cosphi = (ch - a) / b;
    let mut secphi = 999999.0;
    if cosphi != 0.0 {
        secphi = 1.0 / cosphi;
    }

    let mut sunset = 99.0;
    let mut sunrse = 99.0;

    if secphi > -1.0 && secphi <= 0.0 {
        return (declin, zenith, sunrse, sunset);
    }
    sunset = -99.0;
    sunrse = -99.0;
    if secphi > 0.0 && secphi < 1.0 {
        return (declin, zenith, sunrse, sunset);
    }

    let mut cosx_s = cosphi;
    if cosx_s.abs() > 1.0 {
        cosx_s = cosx_s.signum();
    }
    let mut phi_s = cosx_s.acos();
    et /= humr;
    phi_s /= humr;
    sunrse = 12.0 - phi_s - et;
    sunset = 12.0 + phi_s - et;

    if sunrse < 0.0 {
        sunrse += 24.0;
    }
    if sunset >= 24.0 {
        sunset -= 24.0;
    }

    if sunrse > sunset {
        let sunx = 99.0_f32.copysign(flat);
        if ld > 91 && ld < 273 {
            sunset = sunx;
            sunrse = sunx;
        } else {
            sunset = -sunx;
            sunrse = -sunx;
        }
    }

    (declin, zenith, sunrse, sunset)
}

// --- HPOL (Epstein Interpolation Helper) ---

pub fn hpol(hour: f32, tw: f32, xnw: f32, sa: f32, su: f32, dsa: f32, dsu: f32) -> f32 {
    use crate::irifun_utils::epst;
    if su.abs() > 25.0 {
        if su > 0.0 {
            tw
        } else {
            xnw
        }
    } else {
        xnw + (tw - xnw) * epst(hour, dsa, sa) + (xnw - tw) * epst(hour, dsu, su)
    }
}

// --- tops_cor2 (Topside Solar Activity Correction Coefficients) ---

pub fn tops_cor2(xh: f32, vmod: f32, a01: &mut [[f32; 2]; 2]) {
    use crate::irifun_utils::booker;
    let pa: [f32; 6 * 3 * 2 * 2] = [
        0.0, 0.0, -2.4, -2.4, 0.0, 0.0,
        0.0, 0.0, -1.6, -1.6, 0.0, 0.0,
        0.0, 0.0, -2.2, -2.2, 0.0, 0.0,
        0.0, 0.0, 0.0185, 0.0185, 0.0, 0.0,
        0.0, 0.0, 0.018, 0.018, 0.0, 0.0,
        0.0, 0.0, 0.0175, 0.0175, 0.0, 0.0,
        0.0, 0.0, -1.1, -1.1, 0.0, 0.0,
        0.0, 0.0, -0.7, -0.7, 0.0, 0.0,
        0.0, 0.0, -1.4, -1.4, 0.0, 0.0,
        0.0, 0.0, 0.007, 0.007, 0.0, 0.0,
        0.0, 0.0, 0.005, 0.005, 0.0, 0.0,
        0.0, 0.0, 0.01, 0.01, 0.0, 0.0,
    ];
    
    let ha: [f32; 6 * 3 * 2 * 2] = [
        0.0, 200.0, 600.0, 900.0, 1400.0, 1700.0,
        0.0, 550.0, 700.0, 1100.0, 1400.0, 1700.0,
        0.0, 200.0, 600.0, 950.0, 1600.0, 1700.0,
        0.0, 300.0, 650.0, 750.0, 1300.0, 1700.0,
        0.0, 450.0, 750.0, 850.0, 1400.0, 1700.0,
        0.0, 300.0, 650.0, 750.0, 1500.0, 1700.0,
        0.0, 400.0, 500.0, 900.0, 1200.0, 1700.0,
        0.0, 400.0, 500.0, 900.0, 1200.0, 1700.0,
        0.0, 350.0, 550.0, 800.0, 1200.0, 1700.0,
        0.0, 400.0, 500.0, 750.0, 900.0, 1700.0,
        0.0, 400.0, 550.0, 750.0, 900.0, 1700.0,
        0.0, 400.0, 550.0, 750.0, 900.0, 1700.0,
    ];

    let xmod = [-90.0, -60.0, -25.0, 0.0, 25.0, 60.0, 90.0];
    let thh = [30.0, 30.0, 30.0, 30.0];
    let thhb = [0.1, 0.1, 0.1, 0.1, 0.1];

    let mut ap01 = [[[0.0_f32; 2]; 2]; 3];
    for j2 in 0..3 {
        for k in 0..2 {
            for l3 in 0..2 {
                let mut ah = [0.0_f32; 6];
                let mut av = [0.0_f32; 6];
                let offset = (j2 + 3 * (k + 2 * l3)) * 6;
                for i in 0..6 {
                    ah[i] = ha[offset + i];
                    av[i] = pa[offset + i];
                }
                ap01[j2][k][l3] = booker(xh, 6, &ah, &av, &thh);
            }
        }
    }

    let mut pb = [[[0.0_f32; 2]; 2]; 7];
    for i in 0..2 {
        for k in 0..2 {
            for l in 0..2 {
                pb[l][i][k] = 0.0;
                pb[l + 5][i][k] = 0.0;
                pb[l + 2][i][k] = ap01[l][i][k];
            }
            pb[4][i][k] = ap01[2][i][k];
        }
    }

    for k in 0..2 {
        for l4 in 0..2 {
            let mut bv = [0.0_f32; 7];
            for i in 0..7 {
                bv[i] = pb[i][k][l4];
            }
            a01[k][l4] = booker(vmod, 7, &xmod, &bv, &thhb);
        }
    }
}

// --- tcor2cal (Topside Correction Solver) ---

pub fn tcor2cal(h: f32, hmf2: f32, hour: f32, xmodip: f32, pf107: f32, hcor2: f32, scahei: f32, sax300: f32, sux300: f32) -> f32 {
    let mut a01 = [[0.0_f32; 2]; 2];
    tops_cor2(h, xmodip, &mut a01);
    let tc2d = a01[0][0] + a01[1][0] * pf107;
    let tc2n = a01[0][1] + a01[1][1] * pf107;
    let mut tc2 = hpol(hour, tc2d, tc2n, sax300, sux300, 1.0, 1.0);
    if h < hcor2 {
        tc2 = (((h - hmf2) / scahei).exp() - 1.0) * tc2;
    }
    tc2
}

// --- GEODIP (Spherical Dipole Transform) ---

pub fn geodip(iyr: i32, sla: &mut f32, slo: &mut f32, dla: &mut f32, dlo: &mut f32, j: i32) {
    let params = crate::iriflip::recalc(iyr);
    let mut r = 1.0;
    if j == 0 {
        let mut col = (90.0 - *sla) * UMR;
        let mut rlo = *slo * UMR;
        let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
        crate::iriflip::sphcar(&mut r, &mut col, &mut rlo, &mut x, &mut y, &mut z, 1);
        let (mut xm, mut ym, mut zm) = (0.0, 0.0, 0.0);
        crate::iriflip::geomag(&mut x, &mut y, &mut z, &mut xm, &mut ym, &mut zm, 1, &params);
        let (mut rm, mut th, mut pf) = (0.0, 0.0, 0.0);
        crate::iriflip::sphcar(&mut rm, &mut th, &mut pf, &mut xm, &mut ym, &mut zm, -1);
        *dlo = pf / UMR;
        let dco = th / UMR;
        *dla = 90.0 - dco;
    } else {
        let mut col = (90.0 - *dla) * UMR;
        let mut rlo = *dlo * UMR;
        let (mut xm, mut ym, mut zm) = (0.0, 0.0, 0.0);
        crate::iriflip::sphcar(&mut r, &mut col, &mut rlo, &mut xm, &mut ym, &mut zm, 1);
        let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
        crate::iriflip::geomag(&mut x, &mut y, &mut z, &mut xm, &mut ym, &mut zm, -1, &params);
        let (mut rm, mut th, mut pf) = (0.0, 0.0, 0.0);
        crate::iriflip::sphcar(&mut rm, &mut th, &mut pf, &mut x, &mut y, &mut z, -1);
        *slo = pf / UMR;
        let sco = th / UMR;
        *sla = 90.0 - sco;
    }
}

pub fn ckp(ap: i32) -> f32 {
    let ap_array = [
        0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67,
        80, 94, 111, 132, 154, 179, 207, 236, 300, 400
    ];
    let mut kp_array = [0.0_f32; 29];
    for i in 2..=28 {
        kp_array[i] = ((i - 1) as f32) / 3.0;
    }
    if ap == 0 {
        return 0.0;
    }
    if ap == 1 {
        return kp_array[2] / 2.0;
    }
    if ap < 8 && ap > 1 {
        return kp_array[ap as usize];
    }
    let xl_ap = (ap as f32).ln();
    let mut alap = [0.0_f32; 29];
    let mut i = 8;
    loop {
        alap[i] = (ap_array[i - 1] as f32).ln();
        if xl_ap > alap[i] {
            i += 1;
            if i <= 28 {
                continue;
            }
        }
        break;
    }
    let i_clamped = if i > 28 { 28 } else { i };
    alap[i_clamped - 1] = (ap_array[i_clamped - 2] as f32).ln();
    alap[i_clamped] = (ap_array[i_clamped - 1] as f32).ln();
    let slope = (kp_array[i_clamped] - kp_array[i_clamped - 1]) / (alap[i_clamped] - alap[i_clamped - 1]);
    kp_array[i_clamped] + slope * (xl_ap - alap[i_clamped])
}

fn get_xe1_val(
    h: f32,
    hmf2: f32,
    nmf2s: f32,
    hmf1: f32,
    f1reg: bool,
    b0: f32,
    b1: f32,
    c1: f32,
    hz: f32,
    t: f32,
    hst: f32,
    hme: f32,
    nmes: f32,
    hef: f32,
    enight: bool,
    e: &[f32; 4],
    hmd: f32,
    nmd: f32,
    hdx: f32,
    d1: f32,
    xkk: f32,
    fp30: f32,
    fp3u: f32,
    fp1: f32,
    fp2: f32,
    beta: f32,
    eta: f32,
    delta: f32,
    zeta: f32,
    b2top: f32,
    itopn: i32,
    tcor1: f32,
    tcor2: f32,
) -> f32 {
    let profile = crate::xe_profile::XeProfile {
        hmf2,
        xnmf2: nmf2s,
        hmf1,
        f1reg,
        b0,
        b1,
        c1,
        hz,
        t,
        hst,
        hme,
        xnme: nmes,
        hef,
        night: enight,
        e: *e,
        hmd,
        xnmd: nmd,
        hdx,
        d1,
        xkk,
        fp30,
        fp3u,
        fp1,
        fp2,
        beta,
        eta,
        delta,
        zeta,
        b2top,
        itopn,
        tcor1,
        tcor2,
    };
    profile.xe_1(h)
}

// --- Main IRI_SUB Subroutine ---

pub fn iri_sub(
    jf: &[bool; 50],
    jmag: i32,
    alati: f32,
    along: f32,
    iyyyy: i32,
    mmdd: i32,
    dhour: f32,
    heibeg: f32,
    heiend: f32,
    heistp: f32,
    outf: &mut [f32],
    oarr: &mut [f32; 100],
) {
    let mess = jf[33];
    
    // Clear outf to -1.0
    for val in outf.iter_mut() {
        *val = -1.0;
    }
    
    // Initialize oarr
    oarr[6] = -1.0;
    oarr[7] = -1.0;
    oarr[8] = -1.0;
    for idx in 10..=13 {
        oarr[idx] = -1.0;
    }
    for idx in 16..=21 {
        oarr[idx] = -1.0;
    }
    for idx in 22..=27 {
        oarr[idx] = -99.0;
    }
    for idx in 28..=31 {
        oarr[idx] = -1.0;
    }
    oarr[33] = -1.0;
    oarr[35] = -1.0;
    oarr[36] = -1.0;
    oarr[37] = -1.0;
    oarr[39] = -1.0;
    oarr[41] = -1.0;
    oarr[42] = -1.0;
    oarr[43] = 9999.0;
    for idx in 44..=56 {
        oarr[idx] = -1.0;
    }
    oarr[48] = -99.0;
    oarr[52] = -99.0;
    oarr[54] = -99.0;
    for idx in 57..=81 {
        oarr[idx] = -99.0;
    }
    oarr[83] = -99.0;

    let era = 6371.2_f32;
    let erequ = 6378.16_f32;
    let erpol = 6356.775_f32;
    let aquad = erequ * erequ;
    let bquad = erpol * erpol;
    let eexc = 0.01675_f32;
    let mut dimo = 0.311653_f32;
    
    let nummax = 1000;
    let mut numhei = (((heiend - heibeg).abs() / heistp.abs()) as i32) + 1;
    if numhei > nummax {
        numhei = nummax;
    }
    
    let mut iyear = iyyyy;
    if iyear < 100 {
        iyear += 1900;
    }
    if iyear < 30 {
        iyear += 2000;
    }
    let mut idayy = 365;
    if iyear % 4 == 0 {
        idayy = 366;
    }
    
    let mut month = 0;
    let mut iday = 0;
    let mut daynr = 0;
    let mut nrdaym = 0;
    if mmdd < 0 {
        daynr = -mmdd;
        crate::irifun_utils::moda(1, iyear, &mut month, &mut iday, &mut daynr, &mut nrdaym);
    } else {
        month = mmdd / 100;
        iday = mmdd - month * 100;
        crate::irifun_utils::moda(0, iyear, &mut month, &mut iday, &mut daynr, &mut nrdaym);
    }
    
    let ryear = (iyear as f32) + ((daynr - 1) as f32) / (idayy as f32);
    let iyd = iyear * 1000 + daynr;
    let amx = PI * ((daynr - 3) as f32) / 182.6;
    let radj = 1.0 - eexc * (amx.cos() + eexc * ((2.0 * amx).cos() - 1.0) / 2.0);
    
    let height_center = (heibeg + heiend) / 2.0;
    
    let mut along = along;
    if along < 0.0 {
        along += 360.0;
    }
    let mut lati = alati;
    let mut longi = along;
    let mut mlat = 0.0;
    let mut mlong = 0.0;
    
    if jmag > 0 {
        mlat = alati;
        mlong = along;
        geodip(iyear, &mut lati, &mut longi, &mut mlat, &mut mlong, 1);
    } else {
        lati = alati;
        longi = along;
        geodip(iyear, &mut lati, &mut longi, &mut mlat, &mut mlong, 0);
    }
    
    let data_dir = std::env::var("IRI2020_DATA_DIR").unwrap_or_else(|_| "src/data".to_string());
    let mut igrf_model = crate::igrf::IgrfModel::new();
    igrf_model.feldcof(ryear, &data_dir).unwrap();
    dimo = igrf_model.dimo;
    let apf107_data = crate::data_io::Apf107Data::load(&data_dir).unwrap();
    let mut cira_model = crate::cira::CiraModel::new();
    let rocdrift_model = crate::rocdrift::RocdriftModel::new();
    
    let mut dec = 0.0_f32;
    let mut dip = 0.0_f32;
    let mut magbr = 0.0_f32;
    let mut modip = 0.0_f32;
    
    if jf[17] {
        let res = crate::igrf::igrf_dip(lati, longi, ryear, 300.0, &data_dir).unwrap();
        dec = res.dec;
        dip = res.dip;
        magbr = res.dipl;
        modip = res.ymodip;
    } else {
        let res = igrf_model.feldg(lati, longi, 300.0);
        let hsp = (res.bnorth * res.bnorth + res.beast * res.beast).sqrt();
        let dip_rad = res.bdown.atan2(hsp);
        dip = dip_rad.to_degrees();
        dec = res.beast.atan2(res.bnorth).to_degrees();
        magbr = (0.5 * dip_rad.tan()).atan().to_degrees();
        let cos_lati = (lati * UMR).cos();
        modip = (dip_rad / (cos_lati * cos_lati + dip_rad * dip_rad).sqrt()).asin().to_degrees();
    }
    
    let mut fl = 0.0_f32;
    let mut icode = 0_i32;
    let mut dipl = 0.0_f32;
    let mut babs = 0.0_f32;
    
    let mut invdip = -99.0_f32;
    let mut invdip_old = -99.0_f32;
    
    if (jf[2] && !jf[5]) || ((jf[1] && !jf[22]) || (jf[1] && jf[47])) {
        let res_600 = igrf_model.feldg(lati, longi, 600.0);
        babs = res_600.babs;
        dipl = (res_600.bdown / (2.0 * (res_600.bnorth * res_600.bnorth + res_600.beast * res_600.beast).sqrt())).atan().to_degrees();
        let (fl_val, icode_val, _) = crate::irifun_utils::shellg(&igrf_model, lati, longi, 600.0);
        fl = fl_val;
        icode = icode_val;
        let mut fl_clamped = fl;
        if fl_clamped > 10.0 {
            fl_clamped = 10.0;
        }
        
        if jf[2] && !jf[5] {
            invdip = invdpc(fl_clamped, dimo, babs, dipl);
        }
        if (jf[1] && !jf[22]) || (jf[1] && jf[47]) {
            invdip_old = invdpc_old(fl_clamped, dimo, babs, dipl);
        }
        println!("RUST: 600km: lati={}, longi={}, ryear={}", lati, longi, ryear);
        println!("RUST: 600km: babs={}, dipl={}, fl_clamped={}, dimo={}", babs, dipl, fl_clamped, dimo);
        println!("RUST: 600km results: invdip={}, invdip_old={}", invdip, invdip_old);
    }
    
    let abslat = lati.abs();
    let absmlt = mlat.abs();
    let absmdp = modip.abs();
    let absmbr = magbr.abs();
    
    let mut cgm_lat = -99.0_f32;
    let mut cgm_lon = -99.0_f32;
    let mut cgm_mlt = -1.0_f32;
    
    let mut hourut = 0.0_f32;
    let mut hour = 0.0_f32;
    
    if dhour <= 24.0 {
        hour = dhour;
        hourut = hour - longi / 15.0;
        if hourut < 0.0 {
            hourut += 24.0;
        }
    } else {
        hourut = dhour - 25.0;
        hour = hourut + longi / 15.0;
        if hour > 24.0 {
            hour -= 24.0;
        }
    }
    
    if jf[46] && abslat > 25.0 {
        let mut dat = [0.0_f32; 44];
        dat[0] = lati;
        dat[1] = longi;
        let mut pla = [0.0_f32; 4];
        let mut plo = [0.0_f32; 4];
        crate::iriflip::geocgm01(
            1,
            iyear,
            height_center,
            &mut dat,
            &mut pla,
            &mut plo,
        );
        cgm_lat = dat[24];
        cgm_lon = dat[25];
        let cgm_mlt00_ut = dat[32];
        cgm_mlt = hourut - cgm_mlt00_ut;
        if cgm_mlt < 0.0 {
            cgm_mlt += 24.0;
        }
    }
    
    let xmlt = igrf_model.clcmlt(iyear, daynr, hourut, lati, longi);
    
    let mut season = ((daynr + 45) / 92) as i32;
    if season < 1 {
        season = 4;
    }
    let nseasn = season;
    let mut zmonth = (month as f32) + ((iday - 1) as f32) / (nrdaym as f32);
    let mut sday = (daynr as f32) / (idayy as f32) * 360.0;
    let mut seaday = daynr as f32;
    let mut iseamon = month;
    if lati < 0.0 {
        season -= 2;
        if season < 1 {
            season += 4;
        }
        iseamon = month + 6;
        if iseamon > 12 {
            iseamon -= 12;
        }
        seaday = (daynr as f32) + (idayy as f32) / 2.0;
        if seaday > (idayy as f32) {
            seaday -= idayy as f32;
        }
        sday += 180.0;
        if sday > 360.0 {
            sday -= 360.0;
        }
    }
    
    let ig_rz_data = crate::data_io::IgRzData::load(&data_dir).unwrap();
    
    let mut rzar = [0.0_f32; 3];
    let mut arig = [0.0_f32; 3];
    let mut ttt = 0.0_f32;
    let mut nmonth = 0;
    crate::irifun_utils::tcon(
        iyear,
        month,
        iday,
        daynr,
        &ig_rz_data.aig,
        &ig_rz_data.arz,
        ig_rz_data.iymst,
        ig_rz_data.iymend,
        &mut rzar,
        &mut arig,
        &mut ttt,
        &mut nmonth,
        mess,
    ).unwrap();
    
    let rzin = !jf[16];
    let igin = !jf[26];
    let f107in = !jf[24];
    let f107_81in = !jf[31];
    
    if rzin {
        let rrr = oarr[32];
        oarr[32] = -1.0;
        rzar[0] = rrr;
        rzar[1] = rrr;
        rzar[2] = rrr;
        if !igin {
            let zi = ((-0.0031 * rrr + 1.5332) * rrr) - 11.5634;
            arig[0] = zi;
            arig[1] = zi;
            arig[2] = zi;
        }
    }
    
    if igin {
        let zi = oarr[38];
        oarr[38] = -99.0;
        arig[0] = zi;
        arig[1] = zi;
        arig[2] = zi;
        if !rzin {
            let mut xigin = zi;
            if xigin > 178.0066 {
                xigin = 178.0066;
            }
            let rrr = 247.29 - 17.96 * (178.0066 - xigin).sqrt();
            rzar[0] = rrr;
            rzar[1] = rrr;
            rzar[2] = rrr;
        }
    }
    
    let rssn = rzar[2];
    let gind = arig[2];
    let mut cov = 63.75 + rssn * (0.728 + rssn * 0.00089);
    
    let mut f107_daily = 0.0_f32;
    let mut f107pd = 0.0_f32;
    let mut f107_81 = 0.0_f32;
    let mut f107_365 = 0.0_f32;
    let mut iap_daily = 0;
    let mut isdate = 0;
    
    if let Some((f107_d, f107_pd, f107_81_val, f107_365_val, iap_d, is_val)) = apf107_data.apf_only(iyear, month, iday) {
        f107_daily = f107_d;
        f107pd = f107_pd;
        f107_81 = f107_81_val;
        f107_365 = f107_365_val;
        iap_daily = iap_d;
        isdate = is_val as i32;
    }
    
    let mut f107d = cov;
    let mut f107y = cov;
    let mut f10781 = cov;
    let mut f107365 = cov;
    
    if !f107in || !f107_81in {
        if f107_daily > -11.1 {
            f107d = f107_daily;
            f107y = f107pd;
            f10781 = f107_81;
            f107365 = f107_365;
        }
    }
    
    if f107in {
        let f107din = oarr[40];
        oarr[40] = -1.0;
        f107d = f107din;
        f107y = f107din;
        if !f107_81in {
            f10781 = f107din;
            f107365 = f107din;
        }
    }
    
    if f107_81in {
        let f10781in = oarr[45];
        oarr[45] = -1.0;
        f10781 = f10781in;
        f107365 = f10781in;
        if !f107in {
            f107d = f10781in;
            f107y = f10781in;
        }
    }
    
    let pf107 = (f107d + f10781) / 2.0;
    
    let f_adj = radj * radj;
    let f107yobs = f107y / f_adj;
    let f10781obs = f10781 / f_adj;
    let pf107obs = pf107 / f_adj;
    
    if jf[40] {
        cov = f107365;
    }
    let mut covsat = cov;
    if covsat > 188.0 {
        covsat = 188.0;
    }
    
    let (sundec, xhi1, sax80, sux80) = soco(daynr, hour, lati, longi, 80.0);
    let (_, xhi2, sax110, sux110) = soco(daynr, hour, lati, longi, 110.0);
    let (_, xhi3, sax200, sux200) = soco(daynr, hour, lati, longi, 200.0);
    let (_, sax300, sux300, _) = soco(daynr, hour, lati, longi, 300.0); // Wait, declin, zenith, sunrise, sunset order:
    // soco returns (declin, zenith, sunrse, sunset)
    // So the 3rd and 4th outputs are sunrise and sunset.
    let (_, _, sax300, sux300) = soco(daynr, hour, lati, longi, 300.0);
    
    let (_, xhinon, sax1, sux1) = soco(daynr, 12.0, lati, longi, 110.0);
    let (_, xhinon2, sax2, sux2) = soco(daynr, 12.0, lati, longi, 200.0);
    
    let mut dnight = false;
    if sax80.abs() > 25.0 {
        if sax80 < 0.0 {
            dnight = true;
        }
    } else if sax80 <= sux80 {
        if hour > sux80 || hour < sax80 {
            dnight = true;
        }
    } else {
        if hour > sux80 && hour < sax80 {
            dnight = true;
        }
    }
    
    let mut enight = false;
    if sax110.abs() > 25.0 {
        if sax110 < 0.0 {
            enight = true;
        }
    } else if sax110 <= sux110 {
        if hour > sux110 || hour < sax110 {
            enight = true;
        }
    } else {
        if hour > sux110 && hour < sax110 {
            enight = true;
        }
    }
    
    let mut fnight = false;
    if sax200.abs() > 25.0 {
        if sax200 < 0.0 {
            fnight = true;
        }
    } else if sax200 <= sux200 {
        if hour > sux200 || hour < sax200 {
            fnight = true;
        }
    } else {
        if hour > sux200 && hour < sax200 {
            fnight = true;
        }
    }
    
    let hnee = 2000.0_f32;
    let mut hnea = 65.0_f32;
    if dnight {
        hnea = 80.0_f32;
    }
    let noden = !jf[0];
    
    let mut foes = 0.0_f32;
    let mut nmes = 0.0_f32;
    let mut foes_quiet = 0.0_f32;
    let mut foe = 0.0_f32;
    let mut nme = 0.0_f32;
    let mut hme = 0.0_f32;
    
    let foein = !jf[14];
    let hmein = !jf[15];
    let afoe = oarr[4];
    if foein {
        oarr[4] = -1.0;
    }
    let mut anme = afoe;
    if foein {
        if afoe < 100.0 {
            anme = 1.24e10 * afoe * afoe;
        } else {
            anme = afoe;
        }
    }
    
    if foein {
        foe = if afoe < 100.0 { afoe } else { (anme / 1.24e10).sqrt() };
        nme = if afoe < 100.0 { 1.24e10 * afoe * afoe } else { anme };
    } else {
        foe = crate::e_layer::foeedi(cov, xhi2, xhinon, abslat);
        nme = 1.24e10 * foe * foe;
    }
    
    hme = if hmein {
        let val = oarr[5];
        oarr[5] = -1.0;
        val
    } else {
        110.0
    };
    
    let fof2in = !jf[7];
    let hmf2in = !jf[8];
    let afof2 = oarr[0];
    if fof2in {
        oarr[0] = -1.0;
    }
    let mut anmf2 = afof2;
    if fof2in {
        if afof2 < 100.0 {
            anmf2 = 1.24e10 * afof2 * afof2;
        } else {
            anmf2 = afof2;
        }
    }
    
    let mut fof2 = if fof2in {
        if afof2 < 100.0 { afof2 } else { (anmf2 / 1.24e10).sqrt() }
    } else {
        0.0
    };
    let mut nmf2 = if fof2in {
        if afof2 < 100.0 { 1.24e10 * afof2 * afof2 } else { anmf2 }
    } else {
        0.0
    };
    
    let ursif2 = !jf[5];
    let ahmf2 = oarr[1];
    if hmf2in {
        oarr[1] = -1.0;
    }
    
    let mut yfof2 = 0.0_f32;
    let mut xm3000 = 0.0_f32;
    let mut xm3_ccir = 0.0_f32;
    let mut b2top = 0.0_f32;
    
    let mut itopn = 0;
    if jf[28] {
        if jf[29] {
            itopn = 0;
        } else {
            itopn = 3;
        }
    } else {
        if jf[29] {
            itopn = 1;
        } else {
            itopn = 2;
        }
    }
    
    if !fof2in || !hmf2in {
        let (f2, fm3) = crate::data_io::CcirUrsiData::load(&data_dir, month, !ursif2).unwrap();
        let (f2n, fm3n) = crate::data_io::CcirUrsiData::load(&data_dir, nmonth, !ursif2).unwrap();
        
        let rr2 = arig[0] / 100.0;
        let rr2n = arig[1] / 100.0;
        let rr1 = 1.0 - rr2;
        let rr1n = 1.0 - rr2n;
        
        let mut ff0 = [0.0_f32; 988];
        let mut ff0n = [0.0_f32; 988];
        for i in 0..76 {
            for j in 0..13 {
                let k = j + 13 * i;
                ff0n[k] = f2n[j][i][0] * rr1n + f2n[j][i][1] * rr2n;
                ff0[k] = f2[j][i][0] * rr1 + f2[j][i][1] * rr2;
            }
        }
        
        let rr2_m = rzar[0] / 100.0;
        let rr2n_m = rzar[1] / 100.0;
        let rr1_m = 1.0 - rr2_m;
        let rr1n_m = 1.0 - rr2n_m;
        
        let mut xm0 = [0.0_f32; 441];
        let mut xm0n = [0.0_f32; 441];
        
        let fm3_loaded = if !ursif2 {
            fm3.unwrap()
        } else {
            let (_, fm3_ccir) = crate::data_io::CcirUrsiData::load(&data_dir, month, true).unwrap();
            fm3_ccir.unwrap()
        };
        let fm3n_loaded = if !ursif2 {
            fm3n.unwrap()
        } else {
            let (_, fm3n_ccir) = crate::data_io::CcirUrsiData::load(&data_dir, nmonth, true).unwrap();
            fm3n_ccir.unwrap()
        };
        
        for i in 0..49 {
            for j in 0..9 {
                let k = j + 9 * i;
                xm0n[k] = fm3n_loaded[j][i][0] * rr1n_m + fm3n_loaded[j][i][1] * rr2n_m;
                xm0[k] = fm3_loaded[j][i][0] * rr1_m + fm3_loaded[j][i][1] * rr2_m;
            }
        }
        
        let zfof2 = crate::irifun_utils::fout(modip, lati, longi, hourut, &ff0);
        let fof2n = crate::irifun_utils::fout(modip, lati, longi, hourut, &ff0n);
        let zm3000 = crate::irifun_utils::xmout(modip, lati, longi, hourut, &xm0);
        let xm300n = crate::irifun_utils::xmout(modip, lati, longi, hourut, &xm0n);
        
        let midm = if month == 2 { 14 } else { 15 };
        if iday < midm {
            yfof2 = fof2n + ttt * (zfof2 - fof2n);
            xm3000 = xm300n + ttt * (zm3000 - xm300n);
        } else {
            yfof2 = zfof2 + ttt * (fof2n - zfof2);
            xm3000 = zm3000 + ttt * (xm300n - zm3000);
        }
        xm3_ccir = xm3000;
    }
    
    if !fof2in {
        fof2 = yfof2;
        nmf2 = 1.24e10 * fof2 * fof2;
    }
    
    let mut fof2s = fof2;
    let mut foes = foe;
    let mut nmf2s = nmf2;
    let mut nmes = nme;
    let mut stormcorr = -1.0_f32;
    let mut estormcor = -1.0_f32;
    let fstorm_on = jf[25] && jf[7];
    let estorm_on = jf[34] && jf[14];
    
    let mut indap = [0_i32; 13];
    if let Some(ap_vals) = apf107_data.apf(isdate as usize, hourut) {
        indap = ap_vals;
    }
    let index_3h_ap = indap[12];
    let mut xkp = 3.0_f32;
    if index_3h_ap > -1 {
        xkp = ckp(index_3h_ap);
    }
    
    if fstorm_on && indap[0] > -1 {
        let icoord = 1;
        let kut = hourut as i32;
        let (cf, rgma) = crate::irifun_utils::storm(
            &indap,
            lati,
            longi,
            icoord,
            kut,
            daynr,
        );
        let cglat = rgma;
        stormcorr = cf;
        fof2s = fof2 * stormcorr;
        nmf2s = 1.24e10 * fof2s * fof2s;
    }
    
    if estorm_on && index_3h_ap > -1 {
        estormcor = crate::irifun_utils::storme_ap(daynr, mlat, index_3h_ap as f32);
        if estormcor > -2.0 {
            foes = foe * estormcor;
            nmes = 1.24e10 * foes * foes;
        }
    }
    
    let mut cgmlat = -99.0_f32;
    let mut ab_mlat = [-99.0_f32; 48];
    if jf[32] {
        let mut zmlt = xmlt;
        if zmlt < 0.0 || zmlt > 24.0 {
            zmlt = -1.0;
        }
        crate::irifun_utils::auroral_boundary(
            xkp,
            zmlt,
            &mut cgmlat,
            &mut ab_mlat,
        );
    }
    
    let mut rlat = lati;
    let mut flon = 0.0_f32;
    if (!jf[3] && jf[30]) || (!jf[38] && jf[39]) {
        flon = longi + 15.0 * hourut;
        if flon > 360.0 {
            flon -= 360.0;
        }
        let x11 = -90.0;
        let x22 = 90.0;
        
        let fmodip_closure = |xlat: f32| {
            let res = crate::igrf::igrf_dip(xlat, flon, ryear, 300.0, &data_dir).unwrap();
            res.ymodip
        };
        
        let fx11 = fmodip_closure(x11);
        let fx22 = fmodip_closure(x22);
        
        let (root_err, xrlat) = crate::irifun_utils::regfa1(
            x11,
            x22,
            fx11,
            fx22,
            0.001,
            modip,
            fmodip_closure,
        );
        let mut final_xrlat = xrlat;
        if root_err {
            final_xrlat = lati;
        }
        rlat = final_xrlat;
    }
    
    let mut hmf2 = 0.0_f32;
    if hmf2in {
        if ahmf2 < 50.0 {
            xm3000 = ahmf2;
            let mut ratf = fof2 / foe;
            if !jf[35] {
                ratf = fof2s / foe;
            }
            hmf2 = crate::irifun_utils::hmf2ed(magbr, rssn, ratf, xm3000);
        } else {
            hmf2 = ahmf2;
            xm3000 = -1.0;
        }
    } else if jf[38] {
        let mut ratf = fof2 / foe;
        if !jf[35] {
            ratf = fof2s / foe;
        }
        hmf2 = crate::irifun_utils::hmf2ed(magbr, rssn, ratf, xm3000);
    } else if jf[39] {
        hmf2 = crate::hmf2_coeff::shamdhmf2(rlat, flon, zmonth, rssn);
    } else {
        hmf2 = crate::irifun_utils::model_hmf2(iday, month, hourut, modip, longi, f10781, &data_dir);
    }
    
    let cos2 = (mlat * UMR).cos().powi(2);
    let mut flu = (covsat - 40.0) / 30.0;
    if jf[6] {
        flu = (cov - 40.0) / 30.0;
    }
    let mut fo1 = fof2s;
    if jf[36] {
        fo1 = fof2;
    }
    let ex = (-mlat / 15.0).exp();
    let ex1 = ex + 1.0;
    let epin = 4.0 * ex / (ex1 * ex1);
    let eta1 = -0.02 * epin;
    let eta = 0.058798 + eta1 - flu * (0.014065 - 0.0069724 * cos2) + fo1 * (0.0024287 + 0.004281 * cos2 - 0.0001528 * fo1);
    let zeta = 0.078922 - 0.0046702 * cos2 - flu * (0.019132 - 0.0076545 * cos2) + fo1 * (0.0032513 + 0.006029 * cos2 - 0.00020872 * fo1);
    let beta = -128.03 + 20.253 * cos2 - flu * (8.0755 + 0.65896 * cos2) + fo1 * (0.44041 + 0.71458 * cos2 - 0.042966 * fo1);
    
    let z_val = (94.5 / beta).exp();
    let z1 = z_val + 1.0;
    let z2 = z_val / (beta * z1 * z1);
    let delta = (eta / z1 - zeta / 2.0) / (eta * z2 + zeta / 400.0);
    
    if itopn == 2 {
        let mut fo2 = fof2s;
        if jf[36] {
            fo2 = fof2;
        }
        let mut dndhmx = -3.467 + 1.714 * fo2.ln() + 2.02 * xm3_ccir.ln();
        dndhmx = dndhmx.exp() * 0.01;
        let b2bot = 0.04774 * fo2 * fo2 / dndhmx;
        let mut b2k = 3.22 - 0.0538 * fo2 - 0.00664 * hmf2 + 0.113 * hmf2 / b2bot + 0.00257 * rssn;
        let ee = (2.0 * (b2k - 1.0)).exp();
        b2k = (b2k * ee + 1.0) / (ee + 1.0);
        b2top = b2k * b2bot;
    }
    
    let mut ppmlat = mlat;
    let ppb = 60.0_f32;
    if ppmlat.abs() > ppb {
        ppmlat = ppb.copysign(ppmlat);
    }
    let cosmag = (ppmlat * UMR).cos();
    let cosmag2 = cosmag * cosmag;
    let hpp = 10000.0_f32;
    let xlpp = (1.0 + hpp / era) / cosmag2;
    
    let mut tcor1 = 0.0_f32;
    let mut tcor2 = 0.0_f32;
    
    let mut hcor1 = 0.0_f32;
    let mut hcor2 = 0.0_f32;
    let mut scahei = 0.0_f32;
    
    let mut pah = [0.0_f32; 6];
    let mut palogne = [0.0_f32; 6];
    let dplas = [100.0_f32, 150.0, 10.0, 10.0];
    
    // Default variables for XeProfile initialization (Block variables)
    let mut hmf1 = 0.0_f32;
    let mut f1reg = false;
    let mut b0 = 0.0_f32;
    let mut b1 = 0.0_f32;
    let mut c1 = 0.0_f32;
    let mut hz = 0.0_f32;
    let mut t_val = 0.0_f32;
    let mut hst = 0.0_f32;
    let mut e_arr = [0.0_f32; 4];
    let mut hmd = 0.0_f32;
    let mut nmd = 0.0_f32;
    let mut hdx = 0.0_f32;
    let mut d1 = 0.0_f32;
    let mut xkk = 0.0_f32;
    let mut fp30 = 0.0_f32;
    let mut fp3u = 0.0_f32;
    let mut fp1 = 0.0_f32;
    let mut fp2 = 0.0_f32;
    
    let mut b0_us = 0.0_f32;
    let mut b1_us = 0.0_f32;
    let mut b0in = !jf[42];
    let mut b1in = !jf[43];
    if jf[42] {
        b1in = false;
    }
    if b0in {
        b0_us = oarr[9];
        oarr[9] = -1.0;
    }
    if b1in {
        b1_us = oarr[34];
        oarr[34] = -1.0;
    }
    
    if jf[3] {
        b0 = crate::b0_b1_model::b0_98(hour, sax200, sux200, nseasn, rssn, longi, modip);
        b1 = hpol(hour, 1.9, 2.6, sax200, sux200, 1.0, 1.0);
    } else if jf[30] {
        b0 = crate::b0_b1_model::shamdb0d(rlat, flon, zmonth, rssn);
        b1 = crate::b0_b1_model::shab1d(lati, flon, zmonth, rssn);
    } else {
        let (seax, grat) = crate::xe_profile::rogul(seaday as i32, xhi3);
        b1 = hpol(hour, 1.9, 2.6, sax200, sux200, 1.0, 1.0);
        let bcoef = b1 * (b1 * (0.0046 * b1 - 0.0548) + 0.2546) + 0.3606;
        let b0cnew = hmf2 * (1.0 - grat);
        b0 = b0cnew / bcoef;
    }
    if b0in {
        b0 = b0_us;
    }
    if b1in {
        b1 = b1_us;
    }
    if b1 > 6.0 {
        b1 = 6.0;
    }
    if b1 < 0.6 {
        b1 = 0.6;
    }
    
    let f1_ocpro = jf[18];
    let f1_l_cond = !jf[19];
    let mut fof1 = -1.0_f32;
    let mut nmf1 = -1.0_f32;
    
    if !f1_ocpro && f1_l_cond {
        f1reg = false;
        fof1 = -1.0;
        nmf1 = -1.0;
        hmf1 = 0.0;
        c1 = 0.0;
    } else {
        let fof1in = !jf[12];
        let hmf1in = !jf[13];
        let afoF1 = oarr[2];
        let mut anmf1 = afoF1;
        if fof1in {
            if afoF1 < 100.0 {
                anmf1 = 1.24e10 * afoF1 * afoF1;
            }
        }
        if fof1in {
            fof1 = if afoF1 < 100.0 { afoF1 } else { (anmf1 / 1.24e10).sqrt() };
            nmf1 = if afoF1 < 100.0 { 1.24e10 * afoF1 * afoF1 } else { anmf1 };
        } else {
            oarr[2] = -1.0;
            fof1 = crate::xe_profile::fof1ed(absmbr, rssn, xhi3);
            nmf1 = 1.24e10 * fof1 * fof1;
        }
        if hmf1in {
            let _ahmf1 = oarr[3];
        } else {
            oarr[3] = -1.0;
        }
        c1 = crate::xe_profile::f1_c1(absmdp, hour, sax2, sux2);
        
        let mut f1pb = 0.0_f32;
        if f1_ocpro {
            let (f1pbw, f1pbl) = crate::xe_profile::f1_prob(xhi3, mlat, rssn);
            f1pb = f1pbw;
            if f1_l_cond {
                f1pb = f1pbl;
            }
        } else {
            f1pb = 0.0;
            if !fnight && fof1 > 0.0 {
                f1pb = 1.0;
            }
        }
        f1reg = false;
        if fof1in || f1pb >= 0.5 {
            f1reg = true;
        }
    }
    

    let mut dela = 4.32_f32;
    if absmdp >= 18.0 {
        dela = 1.0 + (-(absmdp - 30.0) / 10.0).exp();
    }
    let dell = 1.0 + (-(abslat - 20.0) / 10.0).exp();
    let xdel = (if season == 1 { 5.0 } else if season == 2 { 5.0 } else if season == 3 { 5.0 } else { 10.0 }) / dela;
    
    let dnds_val = (if season == 1 { 0.016 } else if season == 2 { 0.01 } else if season == 3 { 0.016 } else { 0.016 }) / dela;
    let hdeep = hpol(hour, 10.5 / dela, 28.0, sax110, sux110, 1.0, 1.0);
    let mut width = hpol(hour, 17.8 / dela, 45.0 + 22.0 / dela, sax110, sux110, 1.0, 1.0);
    let mut depth = hpol(hour, xdel, 81.0, sax110, sux110, 1.0, 1.0);
    let dlndh = hpol(hour, dnds_val, 0.06, sax110, sux110, 1.0, 1.0);
    
    let mut hefold = hme;
    let mut hef = hme;
    let mut vner = nmes;
    if depth >= 1.0 {
        let mut depth_val = depth;
        if enight {
            depth_val = -depth_val;
        }
        let mut ext = false;
        crate::xe_profile::tal(hdeep, depth_val, width, dlndh, &mut ext, &mut e_arr);
        if ext {
            width = 0.0;
        }
        hef = hme + width;
        hefold = hef;
        vner = (1.0 - depth_val.abs() / 100.0) * nmes;
    }
    
    let hmex = hme - 9.0;
    nmd = crate::e_layer::xmded(xhi1, rssn, 4.0e8);
    hmd = hpol(hour, 81.0, 88.0, sax80, sux80, 1.0, 1.0);
    let f1_p = hpol(hour, 0.02 + 0.03 / dela, 0.05, sax80, sux80, 1.0, 1.0);
    let f2_p = hpol(hour, 4.6, 4.5, sax80, sux80, 1.0, 1.0);
    let f3_p = hpol(hour, -11.5, -4.0, sax80, sux80, 1.0, 1.0);
    fp1 = f1_p;
    fp2 = -fp1 * fp1 / 2.0;
    fp30 = (-f2_p * fp2 - fp1 + 1.0 / f2_p) / (f2_p * f2_p);
    fp3u = (-f3_p * fp2 - fp1 - 1.0 / f3_p) / (f3_p * f3_p);
    hdx = hmd + f2_p;
    
    let x_xdx = hdx - hmd;
    let xdx = nmd * (x_xdx * (fp1 + x_xdx * (fp2 + x_xdx * fp30))).exp();
    let dxdx = xdx * (fp1 + x_xdx * (2.0 * fp2 + x_xdx * 3.0 * fp30));
    let x_hme = hme - hdx;
    xkk = -dxdx * x_hme / (xdx * (xdx / nmes).ln());
    let xkkmax = 5.0_f32;
    if xkk > xkkmax {
        xkk = xkkmax;
        d1 = -(xdx / nmes).ln() / x_hme.powf(xkk);
    } else {
        d1 = dxdx / (xdx * xkk * x_hme.powf(xkk - 1.0));
    }
    
    let mut ddens = [[-1.0_f32; 11]; 5];
    let dreg = jf[23];
    if !dreg {
        let vkp = 1.0_f32;
        let mut elg = [0.0_f32; 7];
        
        let f5sw_vals = [0.0_f32, 0.5, 1.0, 0.0, 0.0];
        let f6wa_vals = [0.0_f32, 0.0, 0.0, 0.5, 1.0];
        for step_idx in 0..5 {
                crate::d_region::dregion(
                    xhi1,
                    month,
                    f107d,
                    vkp,
                    f5sw_vals[step_idx],
                    f6wa_vals[step_idx],
                    &mut elg,
                );
            for ii in 0..11 {
                if ii < 7 {
                    ddens[step_idx][ii] = 10.0_f32.powf(elg[ii] + 6.0);
                } else {
                    ddens[step_idx][ii] = -1.0;
                }
            }
        }
    }
    
    let layver = !jf[10];
    if !layver {
        let hmf2_copy = hmf2;
        let nmf2s_copy = nmf2s;
        let b0_copy = b0;
        let b1_copy = b1;
        let xe2_eval = move |h: f32| {
            let p = crate::xe_profile::XeProfile {
                hmf2: hmf2_copy,
                xnmf2: nmf2s_copy,
                b0: b0_copy,
                b1: b1_copy,
                hmf1: 0.0,
                f1reg: false,
                c1: 0.0,
                hz: 0.0,
                t: 0.0,
                hst: 0.0,
                hme: 0.0,
                xnme: 0.0,
                hef: 0.0,
                night: false,
                e: [0.0; 4],
                hmd: 0.0,
                xnmd: 0.0,
                hdx: 0.0,
                d1: 0.0,
                xkk: 0.0,
                fp30: 0.0,
                fp3u: 0.0,
                fp1: 0.0,
                fp2: 0.0,
                beta: 0.0,
                eta: 0.0,
                delta: 0.0,
                zeta: 0.0,
                b2top: 0.0,
                itopn: 0,
                tcor1: 0.0,
                tcor2: 0.0,
            };
            p.xe2(h)
        };

        hmf1 = 0.0;
        let mut loop_f1 = f1reg;
        if loop_f1 {
            let bnmf1 = 0.9 * nmf1;
            if nmes >= bnmf1 {
                if mess {
                    println!(" *Ne* hmf1 is not evaluated by the function xe2\n corr.: no f1 region");
                }
                hmf1 = 0.0;
                f1reg = false;
                c1 = 0.0;
                loop_f1 = false;
            }
            
            if loop_f1 {
                'search_loop: loop {
                    let mut xe2h = xe2_eval(hef);
                    let mut omit_f1 = false;
                    while xe2h > bnmf1 {
                        hef -= 1.0;
                        if hef <= hme {
                            hef = hme;
                            width = 0.0;
                            hefold = hef;
                            omit_f1 = true;
                            break;
                        }
                        xe2h = xe2_eval(hef);
                    }
                    
                    if omit_f1 {
                        if mess {
                            println!(" *Ne* hmf1 is not evaluated by the function xe2\n corr.: no f1 region");
                        }
                        hmf1 = 0.0;
                        f1reg = false;
                        c1 = 0.0;
                    } else {
                        let (root_err, xhmf1) = crate::irifun_utils::regfa1(
                            hef,
                            hmf2,
                            xe2h,
                            nmf2s,
                            0.001,
                            nmf1,
                            xe2_eval,
                        );
                        if root_err {
                            if mess {
                                println!(" *Ne* hmf1 is not evaluated by the function xe2\n corr.: no f1 region");
                            }
                            hmf1 = 0.0;
                            f1reg = false;
                            c1 = 0.0;
                        } else {
                            hmf1 = xhmf1;
                        }
                    }
                    
                    if hef != hefold {
                        width = hef - hme;
                        let mut dep = depth;
                        if enight {
                            dep = -depth;
                        }
                        let mut ext_val = false;
                            crate::xe_profile::tal(
                                hdeep,
                                dep,
                                width,
                                dlndh,
                                &mut ext_val,
                                &mut e_arr,
                            );
                        if !ext_val {
                            break 'search_loop;
                        } else {
                            if mess {
                                println!(" *Ne* E-region valley can not be modelled");
                            }
                            width = 0.0;
                            hef = hme;
                            hefold = hef;
                            if !f1reg {
                                break 'search_loop;
                            }
                            continue 'search_loop;
                        }
                    } else {
                        break 'search_loop;
                    }
                }
            }
        }

        let mut hf1 = hmf1;
        let mut xf1 = nmf1;
        if !f1reg {
            hf1 = (hmf2 + hef) / 2.0;
            xf1 = xe2_eval(hf1);
        }
        
        let mut hf2 = 100.0_f32;
        
        let hz_cell = std::cell::Cell::new(hz);
        let hst_cell = std::cell::Cell::new(hst);

        let xe3_1_eval = |h: f32| {
            let p = crate::xe_profile::XeProfile {
                hmf2,
                xnmf2: nmf2s,
                hmf1,
                f1reg,
                b0,
                b1,
                c1,
                hz: hz_cell.get(),
                t: t_val,
                hst: hst_cell.get(),
                hme,
                xnme: nmes,
                hef,
                night: enight,
                e: e_arr,
                hmd,
                xnmd: nmd,
                hdx,
                d1,
                xkk,
                fp30,
                fp3u,
                fp1,
                fp2,
                beta,
                eta,
                delta,
                zeta,
                b2top,
                itopn,
                tcor1,
                tcor2,
            };
            p.xe3_1(h)
        };
        
        let xf2 = xe3_1_eval(hf2);
        if xf2 <= nmes {
            let mut schalt = false;
            let (root_err, xhst) = crate::irifun_utils::regfa1(
                hf1,
                hf2,
                xf1,
                xf2,
                0.001,
                nmes,
                xe3_1_eval,
            );
            schalt = root_err;
            if !schalt {
                hst_cell.set(xhst);
                hz_cell.set((hst_cell.get() + hf1) / 2.0);
            } else {
                hz_cell.set((hef + hf1) / 2.0);
                let xnehz = xe3_1_eval(hz_cell.get());
                t_val = (xnehz - nmes) / (hz_cell.get() - hef);
                hst_cell.set(-333.0);
            }
        } else {
            hz_cell.set((hef + hf1) / 2.0);
            let xnehz = xe3_1_eval(hz_cell.get());
            t_val = (xnehz - nmes) / (hz_cell.get() - hef);
            hst_cell.set(-333.0);
        }
        hz = hz_cell.get();
        hst = hst_cell.get();
    }
    
    let mut hxl = [0.0_f32; 4];
    let mut scl = [0.0_f32; 4];
    let mut amp = [0.0_f32; 4];
    let mut hhalf = 0.0_f32;
    let mut vner_lay = vner;
    let mut hxl1_choice = 0;
    
    if layver {
        let mut hmf1m = 0.0_f32;
        if !jf[13] {
            hmf1m = oarr[3];
            oarr[3] = -1.0;
        } else {
            hmf1m = 165.0 + 0.6428 * xhi3;
        }
        hhalf = grat_from_rogul(seaday, xhi3) * hmf2; // wait, grat is computed from rogul
        // In Fortran:
        // HHALF = GRAT * HMF2
        // where GRAT is computed inside ROGUL or ROGUL call. Let's make sure we have a helper for it.
        let hv1r = hme + width;
        let hv2r = hme + hdeep;
        let mut iqu_c = 0;
        crate::xe_profile::inilay(
            fnight,
            f1reg,
            nmf2s,
            nmf1,
            nmes,
            vner,
            hmf2,
            hmf1m,
            hme,
            hv1r,
            hv2r,
            hhalf,
            &mut hxl,
            &mut scl,
            &mut amp,
            &mut iqu_c,
        );
        hxl1_choice = iqu_c;
    }
    
    if itopn == 1 || itopn == 3 {
        let mut hppo = 20000.0_f32;
        let mut hpt = 30000.0_f32;

        let zmp111 = crate::irifun_utils::epla(modip, 10.0, 0.0);
        let zmp222 = crate::irifun_utils::epla(modip, 19.0, 0.0);
        let r2n = -0.84 - 1.6 * zmp111;
        let r2d = -0.84 - 0.64 * zmp111;
        let x1n = 230.0 - 700.0 * zmp222;
        let x1d = 550.0 - 1900.0 * zmp222;
        let r2 = hpol(hour, r2d, r2n, sax300, sux300, 1.0, 1.0);
        let x1 = hpol(hour, x1d, x1n, sax300, sux300, 1.0, 1.0);
        hcor1 = hmf2 + x1;
        let x12 = 1500.0 - x1;
        let tc3 = r2 / x12;
        hcor2 = (hcor1 + hmf2) / 2.0;
        scahei = (hcor2 - hmf2) / 2.0_f32.ln();

        let hmid = 5000.0_f32;
        let xlmid = (1.0 + hmid / era) / cosmag2;
        let xlpt = (1.0 + hpt / era) / cosmag2;
        let xlppo = (1.0 + hppo / era) / cosmag2;

        let (xnepp, xnemid, xneppo, xnept) = if jf[48] {
            (
                ohzden(xlpp, ppmlat),
                ohzden(xlmid, ppmlat),
                ohzden(xlppo, ppmlat),
                ohzden(xlpt, ppmlat),
            )
        } else {
            (
                gallden(xlpp, daynr as f32, rssn),
                gallden(xlmid, daynr as f32, rssn),
                gallden(xlppo, daynr as f32, rssn),
                gallden(xlpt, daynr as f32, rssn),
            )
        };

        pah[0] = hcor1;
        pah[1] = hmf2 + 1500.0;
        pah[2] = hmid;
        pah[3] = hpp;
        pah[4] = hppo;
        pah[5] = hpt;
        palogne[0] = 0.0;
        palogne[1] = r2 * 10.0_f32.ln();
        palogne[2] = (xnemid / nmf2s).ln();
        palogne[3] = (xnepp / nmf2s).ln();
        palogne[4] = (xneppo / nmf2s).ln();
        palogne[5] = (xnept / nmf2s).ln();

        tcor1 = 0.0;
        let mut a01_mid = [[0.0_f32; 2]; 2];
        tops_cor2(hmid_val(), modip, &mut a01_mid);
        let tc2d = a01_mid[0][0] + a01_mid[1][0] * pf107;
        let tc2n = a01_mid[0][1] + a01_mid[1][1] * pf107;
        let (_, _, sap_hmid, sup_hmid) = soco(daynr, hour, lati, longi, hmid_val());
        tcor2 = hpol(hour, tc2d, tc2n, sap_hmid, sup_hmid, 1.0, 1.0);
        let xe1_hmid = get_xe1_val(hmid_val(), hmf2, nmf2s, hmf1, f1reg, b0, b1, c1, hz, t_val, hst, hme, nmes, hef, enight, &e_arr, hmd, nmd, hdx, d1, xkk, fp30, fp3u, fp1, fp2, beta, eta, delta, zeta, b2top, itopn, tcor1, tcor2);
        let znemid = (xe1_hmid / nmf2s).ln();
        palogne[2] -= znemid;
        
        let (_, _, sap_hpp, sup_hpp) = soco(daynr, hour, lati, longi, hpp);
        tcor2 = tcor2cal(hpp, hmf2, hour, modip, pf107, hcor2, scahei, sap_hpp, sup_hpp);
        tcor1 = crate::irifun_utils::booker(hpp, 6, &pah, &palogne, &dplas);
        let xe1_hpp = get_xe1_val(hpp, hmf2, nmf2s, hmf1, f1reg, b0, b1, c1, hz, t_val, hst, hme, nmes, hef, enight, &e_arr, hmd, nmd, hdx, d1, xkk, fp30, fp3u, fp1, fp2, beta, eta, delta, zeta, b2top, itopn, tcor1, tcor2);
        let znepp = (xe1_hpp / nmf2s).ln();
        palogne[3] -= znepp;
        
        let (_, _, sap_hppo, sup_hppo) = soco(daynr, hour, lati, longi, hppo);
        tcor2 = tcor2cal(hppo, hmf2, hour, modip, pf107, hcor2, scahei, sap_hppo, sup_hppo);
        tcor1 = crate::irifun_utils::booker(hppo, 6, &pah, &palogne, &dplas);
        let xe1_hppo = get_xe1_val(hppo, hmf2, nmf2s, hmf1, f1reg, b0, b1, c1, hz, t_val, hst, hme, nmes, hef, enight, &e_arr, hmd, nmd, hdx, d1, xkk, fp30, fp3u, fp1, fp2, beta, eta, delta, zeta, b2top, itopn, tcor1, tcor2);
        let zneppo = (xe1_hppo / nmf2s).ln();
        palogne[4] -= zneppo;
        
        let (_, _, sap_hpt, sup_hpt) = soco(daynr, hour, lati, longi, hpt);
        tcor2 = tcor2cal(hpt, hmf2, hour, modip, pf107, hcor2, scahei, sap_hpt, sup_hpt);
        tcor1 = crate::irifun_utils::booker(hpt, 6, &pah, &palogne, &dplas);
        let xe1_hpt = get_xe1_val(hpt, hmf2, nmf2s, hmf1, f1reg, b0, b1, c1, hz, t_val, hst, hme, nmes, hef, enight, &e_arr, hmd, nmd, hdx, d1, xkk, fp30, fp3u, fp1, fp2, beta, eta, delta, zeta, b2top, itopn, tcor1, tcor2);
        let znept = (xe1_hpt / nmf2s).ln();
        palogne[5] -= znept;
    }
    
    let hta = 60.0_f32;
    let hequi = 120.0_f32;
    let hte = 3000.0_f32;
    
    let mut tn120 = 0.0_f32;
    let mut iapo = [0.0_f32; 7];
    let mut d_msis = [0.0_f32; 9];
    let mut t_msis = [0.0_f32; 2];
    
    let notem = !jf[1];
    let noion = !jf[2];
    let rbtt = !jf[5];
    
    let mut sec = hourut * 3600.0;
    
    if !notem || (!noion && rbtt) {
        if let Some(apf_msis_vals) = apf107_data.apfmsis(isdate as usize, hourut) {
            iapo = apf_msis_vals;
        }
        println!("RUST iapo: {:?}", iapo);
        let mut swmi = [1.0_f32; 25];
        if iapo[1] < 0.0 {
            swmi[8] = 0.0;
            iapo[0] = 0.0;
        } else {
            swmi[8] = -1.0;
        }
        cira_model.tselec(&swmi);
        cira_model.gtd7(
            iyd,
            sec,
            hequi,
            lati,
            longi,
            hour,
            f10781obs,
            f107yobs,
            &iapo,
            0,
            &mut d_msis,
            &mut t_msis,
        );
        tn120 = t_msis[1];
    }
    
    let mut ahh = [0.0_f32; 8];
    let mut ate = [0.0_f32; 8];
    let mut stte = [0.0_f32; 7];
    let mut dte = [5.0_f32, 5.0, 10.0, 20.0, 20.0, 20.0]; // size 6
    let mut dti = [5.0_f32, 5.0, 10.0, 20.0, 20.0, 20.0];
    
    let mut ti1 = 0.0_f32;
    let mut xteti = 0.0_f32;
    let mut xsm = [0.0_f32; 7];
    let mut mm_ti = [0.0_f32; 6];
    let mut mxsm = 0;
    let mut tnhs = 0.0_f32;
    let mut hs = 200.0_f32;
    let mut ate1 = 0.0_f32;
    
    let teneop = !jf[9];
    
    if !notem {
        ahh[1] = 120.0;
        ate[1] = tn120;
        
        let hmaxd = 60.0 * (-((mlat / 22.41).powi(2))).exp() + 210.0;
        let hmaxn = 150.0;
        ahh[2] = hpol(hour, hmaxd, hmaxn, sax200, sux200, 1.0, 1.0);
        let tmaxd = 800.0 * (-((mlat / 33.0).powi(2))).exp() + 1500.0;
        let secni = (24.0 - longi / 15.0) * 3600.0;
        cira_model.gtd7(
            iyd,
            secni,
            hmaxn,
            lati,
            longi,
            0.0,
            f10781obs,
            f107yobs,
            &iapo,
            0,
            &mut d_msis,
            &mut t_msis,
        );
        let tmaxn = t_msis[1];
        ate[2] = hpol(hour, tmaxd, tmaxn, sax200, sux200, 1.0, 1.0);
        
        if jf[22] {
            let mut tea = [0.0_f32; 4];
            crate::irifun_utils::teba(magbr, hour, nseasn, &mut tea);
            ahh[3] = 300.0;
            ahh[4] = 400.0;
            ahh[5] = 600.0;
            ahh[6] = 1400.0;
            ahh[7] = 3000.0;
            ate[3] = tea[0];
            ate[4] = tea[1];
            ate[6] = tea[2];
            ate[7] = tea[3];
            
            let ett = (-mlat / 11.35).exp();
            let tet = 2900.0 - 5600.0 * ett / ((ett + 1.0).powi(2));
            let ten = 839.0 + 1161.0 / (1.0 + (-(absmlt - 45.0) / 5.0).exp());
            ate[5] = hpol(hour, tet, ten, sax300, sux300, 1.5, 1.5);
        } else {
            ahh[3] = 350.0;
            ahh[4] = 550.0;
            ahh[5] = 850.0;
            ahh[6] = 1400.0;
            ahh[7] = 2000.0;
            
            let isa = if jf[41] { 1 } else { 0 };
            let mut teva = [0.0_f32; 5];
            let mut sdteva = [0.0_f32; 5];
            println!("RUST: elteik inputs: isa={}, invdip_old={}, xmlt={}, daynr={}, pf107obs={}", isa, invdip_old, xmlt, daynr, pf107obs);
            crate::irifun_utils::elteik(
                isa,
                invdip_old,
                xmlt,
                daynr as f32,
                pf107obs,
                &mut teva,
                &mut sdteva,
            );
            println!("RUST: elteik_c outputs: teva={:?}, sdteva={:?}", teva, sdteva);
            for ijk in 3..=7 {
                ate[ijk] = teva[ijk - 3];
            }
        }
        
        if teneop {
            for i in 1..=2 {
                let xnar_i = oarr[i + 13];
                if xnar_i > 0.0 {
                    ate[i + 2] = crate::irifun_utils::tede(ahh[i + 2], xnar_i, -cov);
                }
            }
            oarr[14] = -1.0;
            oarr[15] = -1.0;
        }
        
        // Enforce Te > Tn at ahh[2]
        cira_model.gtd7(
            iyd,
            sec,
            ahh[2],
            lati,
            longi,
            hour,
            f10781obs,
            f107yobs,
            &iapo,
            0,
            &mut d_msis,
            &mut t_msis,
        );
        let tnahh2 = t_msis[1];
        println!("RUST: tnahh2={}, ate[2] (before clamp)={}, ate[1]={}", tnahh2, ate[2], ate[1]);
        if ate[2] < tnahh2 {
            ate[2] = tnahh2;
        }
        let mut stte1 = (ate[2] - ate[1]) / (ahh[2] - ahh[1]);
        for i in 2..=6 {
            let orig_ate_i_plus_1 = ate[i + 1];
            cira_model.gtd7(
                iyd,
                sec,
                ahh[i + 1],
                lati,
                longi,
                hour,
                f10781obs,
                f107yobs,
                &iapo,
                0,
                &mut d_msis,
                &mut t_msis,
            );
            let tnahhi = t_msis[1];
            if ate[i + 1] < tnahhi {
                ate[i + 1] = tnahhi;
            }
            println!("RUST: loop i={}: ahh[i+1]={}, tnahhi={}, ate[i+1] (before clamp)={}, ate[i+1] (after clamp)={}", i, ahh[i + 1], tnahhi, orig_ate_i_plus_1, ate[i + 1]);
            let stte2 = (ate[i + 1] - ate[i]) / (ahh[i + 1] - ahh[i]);
            let orig_ate_i = ate[i];
            ate[i] = ate[i] - (stte2 - stte1) * dte[i - 2] * 2.0_f32.ln();
            println!("RUST: loop i={}: stte2={}, stte1={}, dte[i-2]={}, ate[i] (before correction)={}, ate[i] (after correction)={}", i, stte2, stte1, dte[i - 2], orig_ate_i, ate[i]);
            stte1 = stte2;
        }
        
        for i in 1..=6 {
            stte[i] = (ate[i + 1] - ate[i]) / (ahh[i + 1] - ahh[i]);
        }
        ate1 = ate[1];
        
        hs = 200.0_f32;
        xsm[1] = hs;
        cira_model.gtd7(
            iyd,
            sec,
            hs,
            lati,
            longi,
            hour,
            f10781obs,
            f107yobs,
            &iapo,
            0,
            &mut d_msis,
            &mut t_msis,
        );
        tnhs = t_msis[1];
        
        if jf[47] {
            xsm[1] = 200.0;
            xsm[2] = 350.0;
            xsm[3] = 430.0;
            xsm[4] = 600.0;
            xsm[5] = 850.0;
            let mut tiv = [0.0_f32; 4];
            let mut sigtv = [0.0_f32; 4];
            crate::irifun_utils::iontif(
                1,
                invdip_old,
                xmlt,
                daynr as f32,
                pf107obs,
                &mut tiv,
                &mut sigtv,
            );
            mm_ti[1] = (tiv[0] - tnhs) / (xsm[2] - xsm[1]);
            mm_ti[2] = (tiv[1] - tiv[0]) / (xsm[3] - xsm[2]);
            mm_ti[3] = (tiv[2] - tiv[1]) / (xsm[4] - xsm[3]);
            mm_ti[4] = (tiv[3] - tiv[2]) / (xsm[5] - xsm[4]);
            mxsm = 3;
        } else {
            let xsm1 = 430.0_f32;
            xsm[2] = xsm1;
            let z1_ti = (-0.09 * mlat).exp();
            let z2_ti = z1_ti + 1.0;
            let tid1 = 1240.0 - 1400.0 * z1_ti / (z2_ti * z2_ti);
            
            let z1_n = absmlt;
            let z2_n = z1_n * (0.47 + z1_n * 0.024) * UMR;
            let z3_n = z2_n.cos();
            let tin1 = 1200.0 - 300.0 * z3_n.signum() * z3_n.abs().sqrt();
            
            let mut ti1_val = tin1;
            if tid1 > tin1 {
                ti1_val = hpol(hour, tid1, tin1, sax300, sux300, 1.0, 1.0);
            }
            
            let ten1 = crate::irifun_utils::booker1(xsm1, 5, ate1, &ahh[1..], &stte[1..], &dte);
            cira_model.gtd7(
                iyd,
                secni,
                xsm1,
                lati,
                longi,
                0.0,
                f10781obs,
                f107yobs,
                &iapo,
                0,
                &mut d_msis,
                &mut t_msis,
            );
            let tnn1 = t_msis[1];
            
            let mut final_ten1 = ten1;
            if final_ten1 < tnn1 {
                final_ten1 = tnn1;
            }
            if ti1_val > final_ten1 {
                ti1_val = final_ten1;
            }
            if ti1_val < tnn1 {
                ti1_val = tnn1;
            }
            ti1 = ti1_val;
            
            xsm[1] = hs;
            mm_ti[1] = (ti1 - tnhs) / (xsm1 - hs);
            mxsm = 1;
            
            mm_ti[2] = hpol(hour, 3.0, 0.0, sax300, sux300, 1.0, 1.0);
            xsm[3] = hte;
            
            let mut xtts = 500.0_f32;
            let mut x_ti = 500.0_f32;
            loop {
                x_ti += xtts;
                if x_ti >= ahh[7] {
                    break;
                }
                let tex = crate::irifun_utils::booker1(x_ti, 5, ate1, &ahh[1..], &stte[1..], &dte);
                let tix = crate::irifun_utils::booker1(x_ti, mxsm, tnhs, &xsm[1..], &mm_ti[1..], &dti);
                if tix >= tex {
                    x_ti -= xtts;
                    xtts /= 10.0;
                    if xtts <= 0.1 {
                        xteti = x_ti + xtts * 5.0;
                        break;
                    }
                }
            }
            
            mxsm = 2;
            mm_ti[3] = stte[6];
            xsm[3] = xteti;
            if xteti <= ahh[6] {
                mxsm = 3;
                mm_ti[3] = stte[5];
                mm_ti[4] = stte[6];
                xsm[4] = ahh[6];
                if xteti <= ahh[5] {
                    mxsm = 4;
                    mm_ti[3] = stte[4];
                    mm_ti[4] = stte[5];
                    mm_ti[5] = stte[6];
                    xsm[4] = ahh[5];
                    xsm[5] = ahh[6];
                }
            }
        }
    }
    
    let mut hnia = 75.0_f32;
    if rbtt {
        hnia = 80.0;
    }
    let hnie = 2000.0_f32;
    
    let mut xhmf1 = hmf1;
    let mut actual_hmf1 = hmf1;
    if hmf1 <= 0.0 {
        actual_hmf1 = hz;
    }
    
    // Core Step Loop
    let mut height = heibeg;
    let mut kk = 1;
    let mut jfirsta = 1;
    let mut jfirste = 1;
    
    loop {
        let (_, xhi_step, sax_step, sux_step) = soco(daynr, hour, lati, longi, height);
        if !noden {
            if height < hnea {
                jfirsta = 1;
            } else {
                if height > hnee {
                    jfirste = 1;
                }
                
                let elede = if layver {
                    let mut elede_val = -9.0_f32;
                    if hxl1_choice < 2 {
                        // Lay-function evaluation
                        let profile = crate::xe_profile::XeProfile {
                            hmf2,
                            xnmf2: nmf2s,
                            hmf1: actual_hmf1,
                            f1reg,
                            b0,
                            b1,
                            c1,
                            hz,
                            t: t_val,
                            hst,
                            hme,
                            xnme: nmes,
                            hef,
                            night: enight,
                            e: e_arr,
                            hmd,
                            xnmd: nmd,
                            hdx,
                            d1,
                            xkk,
                            fp30,
                            fp3u,
                            fp1,
                            fp2,
                            beta,
                            eta,
                            delta,
                            zeta,
                            b2top,
                            itopn,
                            tcor1,
                            tcor2,
                        };
                        elede_val = crate::xe_profile::xen(
                            height,
                            hmf2,
                            nmf2s,
                            hme,
                            4,
                            &hxl,
                            &scl,
                            &amp,
                            &profile,
                        );
                    }
                    elede_val
                } else {
                    if itopn == 1 || itopn == 3 {
                        tcor1 = 0.0;
                        if height >= hcor1 {
                            tcor1 = crate::irifun_utils::booker(height, 6, &pah, &palogne, &dplas);
                        }
                        tcor2 = 0.0;
                        if itopn == 3 && height > hmf2 {
                            let (_, _, sap, sup) = soco(daynr, hour, lati, longi, height);
                            tcor2 = tcor2cal(height, hmf2, hour, modip, pf107, hcor2, scahei, sap, sup);
                        }
                    }
                    
                    let mut density = get_xe1_val(height, hmf2, nmf2s, actual_hmf1, f1reg, b0, b1, c1, hz, t_val, hst, hme, nmes, hef, enight, &e_arr, hmd, nmd, hdx, d1, xkk, fp30, fp3u, fp1, fp2, beta, eta, delta, zeta, b2top, itopn, tcor1, tcor2);
                    
                    if !dreg && height <= 140.0 {
                        let mut edens = 0.0_f32;
                        let mut ierror = 0_i32;
                        if let Ok((edens_val, ierror_val)) = crate::iridreg::f00(height, lati, daynr, xhi_step, f107d) {
                            edens = edens_val;
                            ierror = ierror_val;
                        }
                        if ierror == 0 || ierror == 2 {
                            density = edens;
                        }
                    }
                    
                    density
                };
                
                outf[((kk - 1) * 20) as usize] = elede;
                
                let pf_gf = 3.2045e-3 * (elede / 1e6).sqrt() / babs;
                outf[((kk - 1) * 20 + 14) as usize] = pf_gf;
                
                if !notem {
                    if height <= hte && height >= hta {
                        let mut d_msis_t = [0.0_f32; 9];
                        let mut t_msis_t = [0.0_f32; 2];
                        cira_model.gtd7(
                            iyd,
                            sec,
                            height,
                            lati,
                            longi,
                            hour,
                            f10781obs,
                            f107yobs,
                            &iapo,
                            0,
                            &mut d_msis_t,
                            &mut t_msis_t,
                        );
                        let tnh = t_msis_t[1];
                        let mut tih = tnh;
                        if height > hs {
                            tih = crate::irifun_utils::booker1(height, mxsm, tnhs, &xsm[1..], &mm_ti[1..], &dti);
                        }
                        let mut teh = tnh;
                        if height > hequi {
                            teh = crate::irifun_utils::booker1(height, 5, ate1, &ahh[1..], &stte[1..], &dte);
                        }
                        
                        if tih < tnh {
                            tih = tnh;
                        }
                        if teh < tnh {
                            teh = tnh;
                        }
                        if tih > teh {
                            tih = teh;
                        }
                        
                        outf[((kk - 1) * 20 + 1) as usize] = tnh;
                        outf[((kk - 1) * 20 + 2) as usize] = tih;
                        outf[((kk - 1) * 20 + 3) as usize] = teh;
                        
                        if !noion {
                            if height <= hnie && height >= hnia {
                                let mut rox = -1.0_f32;
                                let mut rhx = -1.0_f32;
                                let mut rnx = -1.0_f32;
                                let mut rhex = -1.0_f32;
                                let mut rnox = -1.0_f32;
                                let mut ro2x = -1.0_f32;
                                let mut rclust = -1.0_f32;
                                
                                 if rbtt {
                                    if height >= 300.0 {
                                        let res_h = igrf_model.feldg(lati, longi, height);
                                        babs = res_h.babs;
                                        dipl = (res_h.bdown / (2.0 * (res_h.bnorth * res_h.bnorth + res_h.beast * res_h.beast).sqrt())).atan().to_degrees();
                                        let (fl_val, icode_val, _) = crate::irifun_utils::shellg(&igrf_model, lati, longi, height);
                                        fl = fl_val;
                                        icode = icode_val;
                                        let mut fl_clamped = fl;
                                        if fl_clamped > 10.0 {
                                            fl_clamped = 10.0;
                                        }
                                        invdip = invdpc(fl_clamped, dimo, babs, dipl);
                                        
                                        let mut xic_o = 0.0_f32;
                                        let mut xic_h = 0.0_f32;
                                        let mut xic_he = 0.0_f32;
                                        let mut xic_n = 0.0_f32;
                                        crate::ioncom::calion(
                                            invdip,
                                            xmlt,
                                            height,
                                            daynr,
                                            pf107obs,
                                            &mut xic_o,
                                            &mut xic_h,
                                            &mut xic_he,
                                            &mut xic_n,
                                        );
                                        rox = xic_o * 100.0;
                                        rhx = xic_h * 100.0;
                                        rnx = xic_n * 100.0;
                                        rhex = xic_he * 100.0;
                                        rnox = 0.0;
                                        ro2x = 0.0;
                                    } else {
                                        let mut d_msis_chem = [0.0_f32; 9];
                                        let mut t_msis_chem = [0.0_f32; 2];
                                         cira_model.gtd7(
                                             iyd,
                                             sec,
                                             height,
                                             lati,
                                             longi,
                                             hour,
                                             f10781obs,
                                             f107yobs,
                                             &iapo,
                                             48,
                                             &mut d_msis_chem,
                                             &mut t_msis_chem,
                                         );
                                        let xn4s = 0.5 * d_msis_chem[7];
                                        let edens = elede / 1e6;
                                        let jprint = if jf[37] { 0 } else { 1 };
                                        let mut ro = 0.0_f32;
                                        let mut ro2 = 0.0_f32;
                                        let mut rno = 0.0_f32;
                                        let mut rn2 = 0.0_f32;
                                        let mut rn = 0.0_f32;
                                        let mut den_no = 0.0_f32;
                                        let mut den_n2d = 0.0_f32;
                                        let mut inewt = 0_i32;
                                        
                                        crate::irifun_utils::chemion(
                                             jprint,
                                             height,
                                             f107yobs,
                                             f10781obs,
                                             teh,
                                             tih,
                                             tnh,
                                             d_msis_chem[1],
                                             d_msis_chem[3],
                                             d_msis_chem[2],
                                             d_msis_chem[0],
                                             d_msis_chem[6],
                                             -1.0,
                                             xn4s,
                                             edens,
                                             -1.0,
                                             xhi_step,
                                             &mut ro,
                                             &mut ro2,
                                             &mut rno,
                                             &mut rn2,
                                             &mut rn,
                                             &mut den_no,
                                             &mut den_n2d,
                                             &mut inewt,
                                         );
                                        if inewt > 0 {
                                            let sumion = edens / 100.0;
                                            rox = ro / sumion;
                                            rhx = 0.0;
                                            rhex = 0.0;
                                            rnx = rn / sumion;
                                            rnox = rno / sumion;
                                            ro2x = ro2 / sumion;
                                        }
                                    }
                                } else {
                                    let dion = crate::ioncom::iondani(
                                        daynr,
                                        iseamon,
                                        height,
                                        xhi_step,
                                        lati,
                                        f107365,
                                    );
                                    rox = dion[0];
                                    rhx = dion[1];
                                    rnx = dion[2];
                                    rhex = dion[3];
                                    rnox = dion[4];
                                    ro2x = dion[5];
                                    rclust = dion[6];
                                }
                                
                                let xnorm = if jf[21] { 1.0 } else { elede / 100.0 };
                                outf[((kk - 1) * 20 + 4) as usize] = rox * xnorm;
                                outf[((kk - 1) * 20 + 5) as usize] = rhx * xnorm;
                                outf[((kk - 1) * 20 + 6) as usize] = rhex * xnorm;
                                outf[((kk - 1) * 20 + 7) as usize] = ro2x * xnorm;
                                outf[((kk - 1) * 20 + 8) as usize] = rnox * xnorm;
                                outf[((kk - 1) * 20 + 9) as usize] = rclust * xnorm;
                                outf[((kk - 1) * 20 + 10) as usize] = rnx * xnorm;
                            }
                        }
                    }
                }
            }
        }
        
        height += heistp;
        kk += 1;
        if kk > numhei {
            break;
        }
    }
    
    if !dreg {
        for ii in 1..=11 {
            let htemp = 55.0 + (ii as f32) * 5.0;
            outf[((ii - 1) * 20 + 13) as usize] = -1.0;
            if htemp >= 65.0 {
                outf[((ii - 1) * 20 + 13) as usize] = get_xe1_val(htemp, hmf2, nmf2s, actual_hmf1, f1reg, b0, b1, c1, hz, t_val, hst, hme, nmes, hef, enight, &e_arr, hmd, nmd, hdx, d1, xkk, fp30, fp3u, fp1, fp2, beta, eta, delta, zeta, b2top, itopn, tcor1, tcor2);
            }
            outf[((10 + ii) * 20 + 13) as usize] = -1.0;
            let mut edens = 0.0_f32;
            let mut ierror = 0_i32;
            if let Ok((edens_val, ierror_val)) = crate::iridreg::f00(htemp, lati, daynr, xhi1, f107d) {
                edens = edens_val;
                ierror = ierror_val;
            }
            if ierror == 0 || ierror == 2 {
                outf[((10 + ii) * 20 + 13) as usize] = edens;
            }
            outf[((21 + ii) * 20 + 13) as usize] = ddens[0][(ii - 1) as usize];
            outf[((32 + ii) * 20 + 13) as usize] = ddens[1][(ii - 1) as usize];
            outf[((43 + ii) * 20 + 13) as usize] = ddens[2][(ii - 1) as usize];
            outf[((54 + ii) * 20 + 13) as usize] = ddens[3][(ii - 1) as usize];
            outf[((65 + ii) * 20 + 13) as usize] = ddens[4][(ii - 1) as usize];
        }
    }
    
    let mut drift = -1.0_f32;
    if jf[20] && magbr.abs() < 25.0 {
        let (jsea, jf107) = rocdrift_model.vfjmodelrocinit(f107d, daynr);
        drift = rocdrift_model.vfjmodelroc(hour, longi, jsea, jf107);
    }
    
    let mut spreadf = -1.0_f32;
    if jf[27] {
        if !(hour > 7.25 && hour < 17.75) && lati.abs() <= 25.0 {
            let mut spfhour = hour;
            let mut daynr1 = daynr;
            if hour < 12.0 {
                spfhour = hour + 24.0;
                daynr1 = daynr - 1;
                if daynr1 < 1 {
                    daynr1 = idayy;
                }
            }
            let mut osfbr = [0.0_f32; 25];
            crate::spreadf_brazil::spreadf_brazil(daynr1, idayy, f107d, lati, &mut osfbr);
            let ispf = ((spfhour - 17.75) / 0.5) as i32 + 1;
            if ispf > 0 && ispf < 26 {
                spreadf = osfbr[(ispf - 1) as usize];
            }
        }
    }
    
    if !noden {
        oarr[0] = nmf2s;
        oarr[1] = hmf2;
        if f1reg {
            oarr[2] = nmf1;
            oarr[3] = xhmf1;
        }
        oarr[4] = nmes;
        oarr[5] = hme;
        oarr[6] = nmd;
        oarr[7] = hmd;
        oarr[8] = hhalf;
        oarr[9] = b0;
        oarr[10] = vner;
        oarr[11] = hef;
    }
    
    if !notem && (!noion || rbtt) {
        oarr[12] = ate[2];
        oarr[13] = ahh[2];
        oarr[14] = ate[3];
        oarr[15] = ate[4];
        oarr[16] = ate[5];
        oarr[17] = ate[6];
        oarr[18] = ate[7];
        oarr[19] = ate[1];
        oarr[20] = ti1;
        oarr[21] = xteti;
    }
    
    oarr[22] = xhi3;
    oarr[23] = sundec;
    oarr[24] = dip;
    oarr[25] = magbr;
    oarr[26] = modip;
    oarr[27] = lati;
    oarr[28] = sax200;
    oarr[29] = sux200;
    oarr[30] = season as f32;
    oarr[31] = longi;
    oarr[32] = rssn;
    oarr[33] = cov;
    oarr[34] = b1;
    oarr[35] = xm3000;
    oarr[38] = gind;
    // oarr(40) = f1pb in Fortran.
    // wait! how is f1pb defined in Fortran?
    // 1700: f1pb = 0.0
    // so oarr[39] = f1pb!
    let mut f1pb_val = 0.0_f32;
    if !f1_ocpro && f1_l_cond {
        f1pb_val = 0.0;
    } else {
        if f1_ocpro {
            let (f1pbw, f1pbl) = crate::xe_profile::f1_prob(xhi3, mlat, rssn);
            f1pb_val = f1pbw;
            if f1_l_cond {
                f1pb_val = f1pbl;
            }
        } else {
            f1pb_val = 0.0;
            if !fnight && fof1 > 0.0 {
                f1pb_val = 1.0;
            }
        }
    }
    oarr[39] = f1pb_val;
    oarr[40] = f107d;
    oarr[41] = c1;
    oarr[42] = daynr as f32;
    oarr[43] = drift;
    oarr[44] = stormcorr;
    oarr[45] = f10781;
    oarr[46] = estormcor;
    oarr[47] = spreadf;
    oarr[48] = mlat;
    oarr[49] = mlong;
    oarr[50] = index_3h_ap as f32;
    oarr[51] = iap_daily as f32;
    oarr[52] = invdip;
    oarr[53] = xmlt;
    oarr[54] = cgm_lat;
    oarr[55] = cgm_lon;
    oarr[56] = cgm_mlt;
    oarr[57] = cgmlat;
    
    let mut jjj = 57;
    for iii in (0..47).step_by(2) {
        jjj += 1;
        oarr[jjj] = ab_mlat[iii];
    }
    
    oarr[82] = xkp;
    oarr[83] = dec;
    oarr[84] = fl;
    oarr[85] = dimo;
    oarr[86] = sax300;
    oarr[87] = sux300;
    oarr[88] = hnea;
    oarr[89] = hnee;
    
    oarr[99] = fof2;
}

// Helpers for other models

fn hmid_val() -> f32 {
    5000.0
}

fn grat_from_rogul(seaday: f32, xhi3: f32) -> f32 {
    let (_, grat) = crate::xe_profile::rogul(seaday as i32, xhi3);
    grat
}

fn ohzden(l: f32, lat: f32) -> f32 {
    let y1 = 4.4693 - 0.4903 * l;
    let y1_clamped = if y1.abs() > 38.0 { y1.signum() * 38.0 } else { y1 };
    let xneq = 10.0_f32.powf(y1_clamped);
    let xinv = (1.0 / l).sqrt().acos() / UMR;
    let y2 = 1.01 * lat / xinv;
    let y3 = (PI * y2 / 2.0).cos();
    let y4 = y3.powf(-0.75);
    1.0e6 * xneq * y4
}

fn gallden(l: f32, day: f32, rz12: f32) -> f32 {
    let dumr = PI / 182.5;
    let y1 = -0.79 * l + 5.3;
    let y2 = dumr * (day + 9.0);
    let y5 = 0.15 * (y2.cos() - 0.5 * (2.0 * y2).cos());
    let y6 = y5 + 0.00127 * rz12 - 0.0635;
    let y7 = y6 * (-(l - 2.0) / 1.5).exp();
    let mut xlog_ne = y1 + y7;
    if xlog_ne.abs() > 38.0 {
        xlog_ne = xlog_ne.signum() * 38.0;
    }
    10.0_f32.powf(xlog_ne + 6.0)
}
