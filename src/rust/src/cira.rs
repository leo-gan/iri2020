use crate::cira_coeff::*;

pub struct CiraModel {
    // Switches
    pub sw: [f32; 26],
    pub swc: [f32; 26],
    pub sav: [f32; 26],
    pub imr: i32,
    pub mess: bool,

    // Caching state for VTST7
    iydl: [i32; 3],
    secl: [f32; 3],
    glatl: [f32; 3],
    gll: [f32; 3],
    stll: [f32; 3],
    fal: [f32; 3],
    fl: [f32; 3],
    apl: [[f32; 8]; 3],
    swl: [[f32; 26]; 3],
    swcl: [[f32; 26]; 3],

    // Caching state for GTD7
    gtd7_alast: f32,
    gtd7_mssl: i32,
    gtd7_ds: [f32; 10],
    gtd7_ts: [f32; 3],

    // Caching state for GTS7
    gts7_alast: f32,
    gts7_tinf: f32,
    gts7_g0: f32,
    gts7_tlb: f32,
    gts7_s: f32,
    gts7_tn1: [f32; 6],
    gts7_tgn1: [f32; 3],

    // Common /GTS3C/
    tlb: f32,
    s: f32,
    db04: f32,
    db16: f32,
    db28: f32,
    db32: f32,
    db40: f32,
    db48: f32,
    db01: f32,
    za: f32,
    t0: f32,
    z0: f32,
    g0: f32,
    rl: f32,
    dd: f32,
    db14: f32,
    tr12: f32,

    // Common /MESO7/
    tn1: [f32; 6],
    tn2: [f32; 5],
    tn3: [f32; 6],
    tgn1: [f32; 3],
    tgn2: [f32; 3],
    tgn3: [f32; 3],

    // Constant node arrays (ZN1 is modified, ZN2/3 are read-only)
    zn1: [f32; 6],
    zn2: [f32; 5],
    zn3: [f32; 6],

    // Common /DMIX/
    dm04: f32,
    dm16: f32,
    dm28: f32,
    dm32: f32,
    dm40: f32,
    dm01: f32,
    dm14: f32,

    // Common /PARMB/
    gsurf: f32,
    re: f32,

    // Common /LPOLY/
    plg: [[f32; 5]; 10], // 1-based indexing plg[1..9][1..4]
    ctloc: f32,
    stloc: f32,
    c2tloc: f32,
    s2tloc: f32,
    c3tloc: f32,
    s3tloc: f32,
    iyr: i32,
    day: f32,
    df: f32,
    dfa: f32,
    apd: f32,
    apdf: f32,
    apt: [f32; 5],
    xlong: f32,

    // Common /TTEST/
    ttest_tinf: f32,
    ttest_gb: f32,
    ttest_rout: f32,
    ttest_t: [f32; 16],

    // Cache for GLOBE7
    globe7_xl: f32,
    globe7_tll: f32,
    globe7_sw9: f32,
    globe7_dayl: f32,
    globe7_p14: f32,
    globe7_p18: f32,
    globe7_p32: f32,
    globe7_p39: f32,
    globe7_cd14: f32,
    globe7_cd18: f32,
    globe7_cd32: f32,
    globe7_cd39: f32,

    // Cache for GLOB7S
    glob7s_dayl: f32,
    glob7s_p32: f32,
    glob7s_p18: f32,
    glob7s_p14: f32,
    glob7s_p39: f32,
    glob7s_cd32: f32,
    glob7s_cd18: f32,
    glob7s_cd14: f32,
    glob7s_cd39: f32,
}

#[inline(always)]
fn zeta(zz: f32, zl: f32, re: f32) -> f32 {
    (zz - zl) * (re + zl) / (re + zz)
}

impl CiraModel {
    pub fn new() -> Self {
        let mut sw = [0.0; 26];
        let mut swc = [0.0; 26];
        let mut sav = [1.0; 26];
        for i in 1..=25 {
            sw[i] = 1.0;
            swc[i] = 1.0;
        }

        let mut zn1 = [0.0; 6];
        zn1[1] = 120.0;
        zn1[2] = 110.0;
        zn1[3] = 100.0;
        zn1[4] = 90.0;
        zn1[5] = 72.5;

        let mut zn2 = [0.0; 5];
        zn2[1] = 72.5;
        zn2[2] = 55.0;
        zn2[3] = 45.0;
        zn2[4] = 32.5;

        let mut zn3 = [0.0; 6];
        zn3[1] = 32.5;
        zn3[2] = 20.0;
        zn3[3] = 15.0;
        zn3[4] = 10.0;
        zn3[5] = 0.0;

        CiraModel {
            sw,
            swc,
            sav,
            imr: 0,
            mess: true,

            iydl: [-999; 3],
            secl: [-999.0; 3],
            glatl: [-999.0; 3],
            gll: [-999.0; 3],
            stll: [-999.0; 3],
            fal: [-999.0; 3],
            fl: [-999.0; 3],
            apl: [[-999.0; 8]; 3],
            swl: [[-999.0; 26]; 3],
            swcl: [[-999.0; 26]; 3],

            gtd7_alast: 99999.0,
            gtd7_mssl: -999,
            gtd7_ds: [0.0; 10],
            gtd7_ts: [0.0; 3],

            gts7_alast: -999.0,
            gts7_tinf: 0.0,
            gts7_g0: 0.0,
            gts7_tlb: 0.0,
            gts7_s: 0.0,
            gts7_tn1: [0.0; 6],
            gts7_tgn1: [0.0; 3],

            tlb: 0.0,
            s: 0.0,
            db04: 0.0,
            db16: 0.0,
            db28: 0.0,
            db32: 0.0,
            db40: 0.0,
            db48: 0.0,
            db01: 0.0,
            za: 0.0,
            t0: 0.0,
            z0: 0.0,
            g0: 0.0,
            rl: 0.0,
            dd: 0.0,
            db14: 0.0,
            tr12: 0.0,

            tn1: [0.0; 6],
            tn2: [0.0; 5],
            tn3: [0.0; 6],
            tgn1: [0.0; 3],
            tgn2: [0.0; 3],
            tgn3: [0.0; 3],

            zn1,
            zn2,
            zn3,

            dm04: 0.0,
            dm16: 0.0,
            dm28: 0.0,
            dm32: 0.0,
            dm40: 0.0,
            dm01: 0.0,
            dm14: 0.0,

            gsurf: 0.0,
            re: 0.0,

            plg: [[0.0; 5]; 10],
            ctloc: 0.0,
            stloc: 0.0,
            c2tloc: 0.0,
            s2tloc: 0.0,
            c3tloc: 0.0,
            s3tloc: 0.0,
            iyr: 0,
            day: 0.0,
            df: 0.0,
            dfa: 0.0,
            apd: 0.0,
            apdf: 0.0,
            apt: [0.0; 5],
            xlong: 0.0,

            ttest_tinf: 0.0,
            ttest_gb: 0.0,
            ttest_rout: 0.0,
            ttest_t: [0.0; 16],

            globe7_xl: 1000.0,
            globe7_tll: 1000.0,
            globe7_sw9: 1.0,
            globe7_dayl: -1.0,
            globe7_p14: -1000.0,
            globe7_p18: -1000.0,
            globe7_p32: -1000.0,
            globe7_p39: -1000.0,
            globe7_cd14: 0.0,
            globe7_cd18: 0.0,
            globe7_cd32: 0.0,
            globe7_cd39: 0.0,

            glob7s_dayl: -1.0,
            glob7s_p32: -1000.0,
            glob7s_p18: -1000.0,
            glob7s_p14: -1000.0,
            glob7s_p39: -1000.0,
            glob7s_cd32: 0.0,
            glob7s_cd18: 0.0,
            glob7s_cd14: 0.0,
            glob7s_cd39: 0.0,
        }
    }

    #[inline(always)]
    fn pdm(&self, i: usize, j: usize) -> f32 {
        PDM[(j - 1) * 10 + (i - 1)]
    }

    #[inline(always)]
    fn pdl(&self, i: usize, j: usize) -> f32 {
        PDL[(j - 1) * 25 + (i - 1)]
    }

    #[inline(always)]
    fn ptl(&self, i: usize, j: usize) -> f32 {
        PTL[(j - 1) * 100 + (i - 1)]
    }

    #[inline(always)]
    fn pma(&self, i: usize, j: usize) -> f32 {
        PMA[(j - 1) * 100 + (i - 1)]
    }

    pub fn tselec(&mut self, sv: &[f32]) {
        for i in 1..=25 {
            let val = sv[i - 1];
            self.sav[i] = val;
            self.sw[i] = val % 2.0;
            if val.abs() == 1.0 || val.abs() == 2.0 {
                self.swc[i] = 1.0;
            } else {
                self.swc[i] = 0.0;
            }
        }
    }

    pub fn tretrv(&self, svv: &mut [f32]) {
        for i in 1..=25 {
            svv[i - 1] = self.sav[i];
        }
    }

    pub fn meters(&mut self, meter: bool) {
        if meter {
            self.imr = 1;
        } else {
            self.imr = 0;
        }
    }

    pub fn glatf(&mut self, lat: f32) {
        let dgtr = 1.74533e-2;
        let c2 = (2.0 * dgtr * lat).cos();
        self.gsurf = 980.616 * (1.0 - 0.0026373 * c2);
        self.re = 2.0 * self.gsurf / (3.085462e-6 + 2.27e-9 * c2) * 1e-5;
    }

    pub fn scalh(&self, alt: f32, xm: f32, temp: f32) -> f32 {
        let rgas = 831.4;
        let g = self.gsurf / (1.0 + alt / self.re).powi(2);
        rgas * temp / (g * xm)
    }

    pub fn dnet(&self, dd: f32, dm: f32, zhm: f32, xmm: f32, xm: f32) -> f32 {
        let a = zhm / (xmm - xm);
        if dm > 0.0 && dd > 0.0 {
            let ylog = a * (dm / dd).ln();
            if ylog < -10.0 {
                dd
            } else if ylog > 10.0 {
                dm
            } else {
                dd * (1.0 + ylog.exp()).powf(1.0 / a)
            }
        } else {
            if dd == 0.0 && dm == 0.0 {
                1.0
            } else if dm == 0.0 {
                dd
            } else {
                dm
            }
        }
    }

    pub fn ccor(&self, alt: f32, r: f32, h1: f32, zh: f32) -> f32 {
        let e = (alt - zh) / h1;
        if e > 70.0 {
            1.0
        } else if e < -70.0 {
            r.exp()
        } else {
            let ex = e.exp();
            (r / (1.0 + ex)).exp()
        }
    }

    pub fn ccor2(&self, alt: f32, r: f32, h1: f32, zh: f32, h2: f32) -> f32 {
        let e1 = (alt - zh) / h1;
        let e2 = (alt - zh) / h2;
        if e1 > 70.0 || e2 > 70.0 {
            1.0
        } else if e1 < -70.0 && e2 < -70.0 {
            r.exp()
        } else {
            let ex1 = e1.exp();
            let ex2 = e2.exp();
            (r / (1.0 + 0.5 * (ex1 + ex2))).exp()
        }
    }

    pub fn vtst7(
        &mut self,
        iyd: i32,
        sec: f32,
        glat: f32,
        glong: f32,
        stl: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        ic: usize,
    ) -> f32 {
        if iyd != self.iydl[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        if sec != self.secl[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        if glat != self.glatl[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        if glong != self.gll[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        if stl != self.stll[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        if f107a != self.fal[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        if f107 != self.fl[ic] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        for i in 1..=7 {
            if ap[i - 1] != self.apl[ic][i] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        }
        for i in 1..=25 {
            if self.sw[i] != self.swl[ic][i] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
            if self.swc[i] != self.swcl[ic][i] { return self.vtst7_update(iyd, sec, glat, glong, stl, f107a, f107, ap, ic); }
        }
        0.0
    }

    fn vtst7_update(
        &mut self,
        iyd: i32,
        sec: f32,
        glat: f32,
        glong: f32,
        stl: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        ic: usize,
    ) -> f32 {
        self.iydl[ic] = iyd;
        self.secl[ic] = sec;
        self.glatl[ic] = glat;
        self.gll[ic] = glong;
        self.stll[ic] = stl;
        self.fal[ic] = f107a;
        self.fl[ic] = f107;
        for i in 1..=7 {
            self.apl[ic][i] = ap[i - 1];
        }
        for i in 1..=25 {
            self.swl[ic][i] = self.sw[i];
            self.swcl[ic][i] = self.swc[i];
        }
        1.0
    }

    pub fn splinem(&self, x: &[f32], y: &[f32], n: usize, yp1: f32, ypn: f32, y2: &mut [f32]) {
        let mut u = [0.0; 101];
        if yp1 > 0.99e30 {
            y2[1] = 0.0;
            u[1] = 0.0;
        } else {
            y2[1] = -0.5;
            u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
        }
        for i in 2..n {
            let sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
            let p = sig * y2[i - 1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1])
                / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
        }
        let (qn, un) = if ypn > 0.99e30 {
            (0.0, 0.0)
        } else {
            (0.5, (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1])))
        };
        y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);
        for k in (1..=n - 1).rev() {
            y2[k] = y2[k] * y2[k + 1] + u[k];
        }
    }

    pub fn splintm(&self, xa: &[f32], ya: &[f32], y2a: &[f32], n: usize, x: f32) -> f32 {
        let mut klo = 1;
        let mut khi = n;
        while khi - klo > 1 {
            let k = (khi + klo) / 2;
            if xa[k] > x {
                khi = k;
            } else {
                klo = k;
            }
        }
        let h = xa[khi] - xa[klo];
        if h == 0.0 {
            if self.mess {
                eprintln!("BAD XA INPUT TO SPLINT");
            }
        }
        let a = (xa[khi] - x) / h;
        let b = (x - xa[klo]) / h;
        a * ya[klo] + b * ya[khi] +
            ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0
    }

    pub fn splini(&self, xa: &[f32], ya: &[f32], y2a: &[f32], n: usize, x: f32) -> f32 {
        let mut yi = 0.0;
        let mut klo = 1;
        let mut khi = 2;
        while klo < n && khi <= n && x > xa[klo] {
            let xx = if khi < n { x.min(xa[khi]) } else { x };
            let h = xa[khi] - xa[klo];
            let a = (xa[khi] - xx) / h;
            let b = (xx - xa[klo]) / h;
            let a2 = a * a;
            let b2 = b * b;
            yi += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 +
                ((-(1.0 + a2 * a2) / 4.0 + a2 / 2.0) * y2a[klo] +
                (b2 * b2 / 4.0 - b2 / 2.0) * y2a[khi]) * h * h / 6.0) * h;
            klo += 1;
            khi += 1;
        }
        yi
    }

    pub fn densu(
        &mut self,
        alt: f32,
        dlb: f32,
        tinf: f32,
        tlb: f32,
        xm: f32,
        alpha: f32,
        tz: &mut f32,
        zlb: f32,
        s2: f32,
        mn1: usize,
    ) -> f32 {
        let re = self.re;
        let gsurf = self.gsurf;
        let rgas = 831.4;
        
        let za = self.zn1[1];
        let z = alt.max(za);
        let zg2 = zeta(z, zlb, re);
        let tt = tinf - (tinf - tlb) * (-s2 * zg2).exp();
        let ta = tt;
        *tz = tt;
        let mut densu_val = *tz;
        
        if alt < za {
            let dta = (tinf - ta) * s2 * ((re + zlb) / (re + za)).powi(2);
            self.tgn1[1] = dta;
            self.tn1[1] = ta;
            
            let z_low = alt.max(self.zn1[mn1]);
            let z1 = self.zn1[1];
            let z2 = self.zn1[mn1];
            let t1 = self.tn1[1];
            let t2 = self.tn1[mn1];
            let zg = zeta(z_low, z1, re);
            let zgdif = zeta(z2, z1, re);
            
            let mut xs = [0.0; 6];
            let mut ys = [0.0; 6];
            for k in 1..=mn1 {
                xs[k] = zeta(self.zn1[k], z1, re) / zgdif;
                ys[k] = 1.0 / self.tn1[k];
            }
            let yd1 = -self.tgn1[1] / (t1 * t1) * zgdif;
            let yd2 = -self.tgn1[2] / (t2 * t2) * zgdif * ((re + z2) / (re + z1)).powi(2);
            
            let mut y2out = [0.0; 6];
            self.splinem(&xs, &ys, mn1, yd1, yd2, &mut y2out);
            let x = zg / zgdif;
            let y = self.splintm(&xs, &ys, &y2out, mn1, x);
            *tz = 1.0 / y;
            densu_val = *tz;

            
            if xm != 0.0 {
                let glb = gsurf / (1.0 + zlb / re).powi(2);
                let gamma = xm * glb / (s2 * rgas * tinf);
                let mut expl = (-s2 * gamma * zg2).exp();
                if expl > 50.0 || tt <= 0.0 {
                    expl = 50.0;
                }
                let densa = dlb * (tlb / tt).powf(1.0 + alpha + gamma) * expl;
                densu_val = densa;
                
                let glb2 = gsurf / (1.0 + z1 / re).powi(2);
                let gamm = xm * glb2 * zgdif / rgas;
                let yi = self.splini(&xs, &ys, &y2out, mn1, x);
                let mut expl2 = gamm * yi;
                if expl2 > 50.0 || *tz <= 0.0 {
                    expl2 = 50.0;
                }
                densu_val = densu_val * (t1 / *tz).powf(1.0 + alpha) * (-expl2).exp();
            }
        } else {
            if xm != 0.0 {
                let glb = gsurf / (1.0 + zlb / re).powi(2);
                let gamma = xm * glb / (s2 * rgas * tinf);
                let mut expl = (-s2 * gamma * zg2).exp();
                if expl > 50.0 || tt <= 0.0 {
                    expl = 50.0;
                }
                let densa = dlb * (tlb / tt).powf(1.0 + alpha + gamma) * expl;
                densu_val = densa;
            }
        }
        
        densu_val
    }

    pub fn densm(
        &mut self,
        alt: f32,
        d0: f32,
        xm: f32,
        tz: &mut f32,
        mn3: usize,
        mn2: usize,
    ) -> f32 {
        let re = self.re;
        let gsurf = self.gsurf;
        let rgas = 831.4;
        
        let mut densm_val = d0;
        if alt > self.zn2[1] {
            if xm == 0.0 {
                densm_val = *tz;
            }
            return densm_val;
        }
        
        // STRATOSPHERE/MESOSPHERE TEMPERATURE
        let z = alt.max(self.zn2[mn2]);
        let z1 = self.zn2[1];
        let z2 = self.zn2[mn2];
        let t1 = self.tn2[1];
        let t2 = self.tn2[mn2];
        let zg = zeta(z, z1, re);
        let zgdif = zeta(z2, z1, re);
        
        let mut xs = [0.0; 11];
        let mut ys = [0.0; 11];
        for k in 1..=mn2 {
            xs[k] = zeta(self.zn2[k], z1, re) / zgdif;
            ys[k] = 1.0 / self.tn2[k];
        }
        let yd1 = -self.tgn2[1] / (t1 * t1) * zgdif;
        let yd2 = -self.tgn2[2] / (t2 * t2) * zgdif * ((re + z2) / (re + z1)).powi(2);
        
        let mut y2out = [0.0; 11];
        self.splinem(&xs, &ys, mn2, yd1, yd2, &mut y2out);
        let x = zg / zgdif;
        let y = self.splintm(&xs, &ys, &y2out, mn2, x);
        *tz = 1.0 / y;
        
        if xm != 0.0 {
            let glb = gsurf / (1.0 + z1 / re).powi(2);
            let gamm = xm * glb * zgdif / rgas;
            let yi = self.splini(&xs, &ys, &y2out, mn2, x);
            let mut expl = gamm * yi;
            if expl > 50.0 {
                expl = 50.0;
            }
            densm_val = densm_val * (t1 / *tz) * (-expl).exp();
        }
        
        if alt > self.zn3[1] {
            if xm == 0.0 {
                densm_val = *tz;
            }
            return densm_val;
        }
        
        // TROPOSPHERE/STRATOSPHERE TEMPERATURE
        let z_trop = alt;
        let z1_trop = self.zn3[1];
        let z2_trop = self.zn3[mn3];
        let t1_trop = self.tn3[1];
        let t2_trop = self.tn3[mn3];
        let zg_trop = zeta(z_trop, z1_trop, re);
        let zgdif_trop = zeta(z2_trop, z1_trop, re);
        
        for k in 1..=mn3 {
            xs[k] = zeta(self.zn3[k], z1_trop, re) / zgdif_trop;
            ys[k] = 1.0 / self.tn3[k];
        }
        let yd1_trop = -self.tgn3[1] / (t1_trop * t1_trop) * zgdif_trop;
        let yd2_trop = -self.tgn3[2] / (t2_trop * t2_trop) * zgdif_trop * ((re + z2_trop) / (re + z1_trop)).powi(2);
        
        self.splinem(&xs, &ys, mn3, yd1_trop, yd2_trop, &mut y2out);
        let x_trop = zg_trop / zgdif_trop;
        let y_trop = self.splintm(&xs, &ys, &y2out, mn3, x_trop);
        *tz = 1.0 / y_trop;
        
        if xm != 0.0 {
            let glb = gsurf / (1.0 + z1_trop / re).powi(2);
            let gamm = xm * glb * zgdif_trop / rgas;
            let yi = self.splini(&xs, &ys, &y2out, mn3, x_trop);
            let mut expl = gamm * yi;
            if expl > 50.0 {
                expl = 50.0;
            }
            densm_val = densm_val * (t1_trop / *tz) * (-expl).exp();
        }
        
        if xm == 0.0 {
            densm_val = *tz;
        }
        
        densm_val
    }

    pub fn globe7(
        &mut self,
        yrd: i32,
        sec: f32,
        lat: f32,
        long: f32,
        tloc: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        p: &[f32],
    ) -> f32 {
        let dgtr = 1.74533e-2;
        let dr = 1.72142e-2;
        let hr = 0.2618;
        let sr = 7.2722e-5;
        
        let g0 = |a: f32| {
            let p25_abs = p[25 - 1].abs();
            let a_minus_4 = a - 4.0;
            a_minus_4 + (p[26 - 1] - 1.0) * (a_minus_4 + ((-p25_abs * a_minus_4).exp() - 1.0) / p25_abs)
        };
        
        let sumex = |ex: f32| {
            1.0 + (1.0 - ex.powi(19)) / (1.0 - ex) * ex.sqrt()
        };
        
        let sg0 = |ex: f32| {
            (g0(ap[2 - 1]) + (g0(ap[3 - 1]) * ex + g0(ap[4 - 1]) * ex * ex + g0(ap[5 - 1]) * ex.powi(3)
                + (g0(ap[6 - 1]) * ex.powi(4) + g0(ap[7 - 1]) * ex.powi(12)) * (1.0 - ex.powi(8)) / (1.0 - ex)
            )) / sumex(ex)
        };
        
        let iyr = yrd / 1000;
        let day = (yrd - iyr * 1000) as f32;
        self.iyr = iyr;
        self.day = day;
        self.xlong = long;
        
        if self.globe7_xl != lat {
            let c = (lat * dgtr).sin();
            let s = (lat * dgtr).cos();
            let c2 = c * c;
            let c4 = c2 * c2;
            let s2 = s * s;
            self.plg[2][1] = c;
            self.plg[3][1] = 0.5 * (3.0 * c2 - 1.0);
            self.plg[4][1] = 0.5 * (5.0 * c * c2 - 3.0 * c);
            self.plg[5][1] = (35.0 * c4 - 30.0 * c2 + 3.0) / 8.0;
            self.plg[6][1] = (63.0 * c2 * c2 * c - 70.0 * c2 * c + 15.0 * c) / 8.0;
            self.plg[7][1] = (11.0 * c * self.plg[6][1] - 5.0 * self.plg[5][1]) / 6.0;
            
            self.plg[2][2] = s;
            self.plg[3][2] = 3.0 * c * s;
            self.plg[4][2] = 1.5 * (5.0 * c2 - 1.0) * s;
            self.plg[5][2] = 2.5 * (7.0 * c2 * c - 3.0 * c) * s;
            self.plg[6][2] = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s;
            self.plg[7][2] = (11.0 * c * self.plg[6][2] - 6.0 * self.plg[5][2]) / 5.0;
            
            self.plg[3][3] = 3.0 * s2;
            self.plg[4][3] = 15.0 * s2 * c;
            self.plg[5][3] = 7.5 * (7.0 * c2 - 1.0) * s2;
            self.plg[6][3] = 3.0 * c * self.plg[5][3] - 2.0 * self.plg[4][3];
            self.plg[7][3] = (11.0 * c * self.plg[6][3] - 7.0 * self.plg[5][3]) / 4.0;
            self.plg[8][3] = (13.0 * c * self.plg[7][3] - 8.0 * self.plg[6][3]) / 5.0;
            
            self.plg[4][4] = 15.0 * s2 * s;
            self.plg[5][4] = 105.0 * s2 * s * c;
            self.plg[6][4] = (9.0 * c * self.plg[5][4] - 7.0 * self.plg[4][4]) / 2.0;
            self.plg[7][4] = (11.0 * c * self.plg[6][4] - 8.0 * self.plg[5][4]) / 3.0;
            
            self.globe7_xl = lat;
        }
        
        if self.globe7_tll != tloc {
            if self.sw[7] != 0.0 || self.sw[8] != 0.0 || self.sw[14] != 0.0 {
                self.stloc = (hr * tloc).sin();
                self.ctloc = (hr * tloc).cos();
                self.s2tloc = (2.0 * hr * tloc).sin();
                self.c2tloc = (2.0 * hr * tloc).cos();
                self.s3tloc = (3.0 * hr * tloc).sin();
                self.c3tloc = (3.0 * hr * tloc).cos();
                self.globe7_tll = tloc;
            }
        }
        
        if day != self.globe7_dayl || p[14 - 1] != self.globe7_p14 {
            self.globe7_cd14 = (dr * (day - p[14 - 1])).cos();
        }
        if day != self.globe7_dayl || p[18 - 1] != self.globe7_p18 {
            self.globe7_cd18 = (2.0 * dr * (day - p[18 - 1])).cos();
        }
        if day != self.globe7_dayl || p[32 - 1] != self.globe7_p32 {
            self.globe7_cd32 = (dr * (day - p[32 - 1])).cos();
        }
        if day != self.globe7_dayl || p[39 - 1] != self.globe7_p39 {
            self.globe7_cd39 = (2.0 * dr * (day - p[39 - 1])).cos();
        }
        
        self.globe7_dayl = day;
        self.globe7_p14 = p[14 - 1];
        self.globe7_p18 = p[18 - 1];
        self.globe7_p32 = p[32 - 1];
        self.globe7_p39 = p[39 - 1];
        
        let df = f107 - f107a;
        let dfa = f107a - 150.0;
        self.df = df;
        self.dfa = dfa;
        
        let mut t_local = [0.0; 16];
        t_local[1] = p[20 - 1] * df * (1.0 + p[60 - 1] * dfa) + p[21 - 1] * df * df + p[22 - 1] * dfa + p[30 - 1] * dfa * dfa;
        let f1 = 1.0 + (p[48 - 1] * dfa + p[20 - 1] * df + p[21 - 1] * df * df) * self.swc[1];
        let f2 = 1.0 + (p[50 - 1] * dfa + p[20 - 1] * df + p[21 - 1] * df * df) * self.swc[1];
        
        t_local[2] = p[2 - 1] * self.plg[3][1] + p[3 - 1] * self.plg[5][1] + p[23 - 1] * self.plg[7][1]
            + (p[15 - 1] * self.plg[3][1]) * dfa * self.swc[1]
            + p[27 - 1] * self.plg[2][1];
            
        t_local[3] = p[19 - 1] * self.globe7_cd32;
        
        t_local[4] = (p[16 - 1] + p[17 - 1] * self.plg[3][1]) * self.globe7_cd18;
        
        t_local[5] = f1 * (p[10 - 1] * self.plg[2][1] + p[11 - 1] * self.plg[4][1]) * self.globe7_cd14;
        
        t_local[6] = p[38 - 1] * self.plg[2][1] * self.globe7_cd39;
        
        if self.sw[7] != 0.0 {
            let t71 = (p[12 - 1] * self.plg[3][2]) * self.globe7_cd14 * self.swc[5];
            let t72 = (p[13 - 1] * self.plg[3][2]) * self.globe7_cd14 * self.swc[5];
            t_local[7] = f2 * ((p[4 - 1] * self.plg[2][2] + p[5 - 1] * self.plg[4][2] + p[28 - 1] * self.plg[6][2] + t71) * self.ctloc
                + (p[7 - 1] * self.plg[2][2] + p[8 - 1] * self.plg[4][2] + p[29 - 1] * self.plg[6][2] + t72) * self.stloc);
        }
        
        if self.sw[8] != 0.0 {
            let t81 = (p[24 - 1] * self.plg[4][3] + p[36 - 1] * self.plg[6][3]) * self.globe7_cd14 * self.swc[5];
            let t82 = (p[34 - 1] * self.plg[4][3] + p[37 - 1] * self.plg[6][3]) * self.globe7_cd14 * self.swc[5];
            t_local[8] = f2 * ((p[6 - 1] * self.plg[3][3] + p[42 - 1] * self.plg[5][3] + t81) * self.c2tloc
                + (p[9 - 1] * self.plg[3][3] + p[43 - 1] * self.plg[5][3] + t82) * self.s2tloc);
        }
        
        if self.sw[14] != 0.0 {
            t_local[14] = f2 * ((p[40 - 1] * self.plg[4][4] + (p[94 - 1] * self.plg[5][4] + p[47 - 1] * self.plg[7][4]) * self.globe7_cd14 * self.swc[5]) * self.s3tloc
                + (p[41 - 1] * self.plg[4][4] + (p[95 - 1] * self.plg[5][4] + p[49 - 1] * self.plg[7][4]) * self.globe7_cd14 * self.swc[5]) * self.c3tloc);
        }
        
        let sw9 = if self.sw[9] > 0.0 { 1.0 } else if self.sw[9] < 0.0 { -1.0 } else { 0.0 };
        if sw9 != -1.0 {
            let apd = ap[0] - 4.0;
            let mut p44 = p[44 - 1];
            let p45 = p[45 - 1];
            if p44 < 0.0 { p44 = 1e-5; }
            let apdf = apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44);
            self.apd = apd;
            self.apdf = apdf;
            if self.sw[9] != 0.0 {
                t_local[9] = apdf * (p[33 - 1] + p[46 - 1] * self.plg[3][1] + p[35 - 1] * self.plg[5][1]
                    + (p[101 - 1] * self.plg[2][1] + p[102 - 1] * self.plg[4][1] + p[103 - 1] * self.plg[6][1]) * self.globe7_cd14 * self.swc[5]
                    + (p[122 - 1] * self.plg[2][2] + p[123 - 1] * self.plg[4][2] + p[124 - 1] * self.plg[6][2]) * self.swc[7]
                        * (hr * (tloc - p[125 - 1])).cos());
            }
        } else {
            if p[52 - 1] != 0.0 {
                let mut exp1 = (-10800.0 * p[52 - 1].abs() / (1.0 + p[139 - 1] * (45.0 - lat.abs()))).exp();
                if exp1 > 0.99999 { exp1 = 0.99999; }
                let apt1 = sg0(exp1);
                self.apt[1] = apt1;
                if self.sw[9] != 0.0 {
                    t_local[9] = apt1 * (p[51 - 1] + p[97 - 1] * self.plg[3][1] + p[55 - 1] * self.plg[5][1]
                        + (p[126 - 1] * self.plg[2][1] + p[127 - 1] * self.plg[4][1] + p[128 - 1] * self.plg[6][1]) * self.globe7_cd14 * self.swc[5]
                        + (p[129 - 1] * self.plg[2][2] + p[130 - 1] * self.plg[4][2] + p[131 - 1] * self.plg[6][2]) * self.swc[7]
                            * (hr * (tloc - p[132 - 1])).cos());
                }
            }
        }
        
        if self.sw[10] != 0.0 && long > -1000.0 {
            if self.sw[11] != 0.0 {
                t_local[11] = (1.0 + p[81 - 1] * dfa * self.swc[1])
                    * (((p[65 - 1] * self.plg[3][2] + p[66 - 1] * self.plg[5][2] + p[67 - 1] * self.plg[7][2]
                        + p[104 - 1] * self.plg[2][2] + p[105 - 1] * self.plg[4][2] + p[106 - 1] * self.plg[6][2]
                        + self.swc[5] * (p[110 - 1] * self.plg[2][2] + p[111 - 1] * self.plg[4][2] + p[112 - 1] * self.plg[6][2]) * self.globe7_cd14)
                        * (dgtr * long).cos())
                    + ((p[91 - 1] * self.plg[3][2] + p[92 - 1] * self.plg[5][2] + p[93 - 1] * self.plg[7][2]
                        + p[107 - 1] * self.plg[2][2] + p[108 - 1] * self.plg[4][2] + p[109 - 1] * self.plg[6][2]
                        + self.swc[5] * (p[113 - 1] * self.plg[2][2] + p[114 - 1] * self.plg[4][2] + p[115 - 1] * self.plg[6][2]) * self.globe7_cd14)
                        * (dgtr * long).sin()));
            }
            
            if self.sw[12] != 0.0 {
                t_local[12] = (1.0 + p[96 - 1] * self.plg[2][1]) * (1.0 + p[82 - 1] * dfa * self.swc[1])
                    * (1.0 + p[120 - 1] * self.plg[2][1] * self.swc[5] * self.globe7_cd14)
                    * ((p[69 - 1] * self.plg[2][1] + p[70 - 1] * self.plg[4][1] + p[71 - 1] * self.plg[6][1])
                        * (sr * (sec - p[72 - 1])).cos())
                    + self.swc[11] * (p[77 - 1] * self.plg[4][3] + p[78 - 1] * self.plg[6][3] + p[79 - 1] * self.plg[8][3])
                        * (sr * (sec - p[80 - 1]) + 2.0 * dgtr * long).cos() * (1.0 + p[138 - 1] * dfa * self.swc[1]);
            }
            
            if self.sw[13] != 0.0 {
                if sw9 != -1.0 {
                    let apdf = self.apdf;
                    t_local[13] = apdf * self.swc[11] * (1.0 + p[121 - 1] * self.plg[2][1])
                        * ((p[61 - 1] * self.plg[3][2] + p[62 - 1] * self.plg[5][2] + p[63 - 1] * self.plg[7][2])
                            * (dgtr * (long - p[64 - 1])).cos())
                        + apdf * self.swc[11] * self.swc[5]
                            * (p[116 - 1] * self.plg[2][2] + p[117 - 1] * self.plg[4][2] + p[118 - 1] * self.plg[6][2])
                            * self.globe7_cd14 * (dgtr * (long - p[119 - 1])).cos()
                        + apdf * self.swc[12]
                            * (p[84 - 1] * self.plg[2][1] + p[85 - 1] * self.plg[4][1] + p[86 - 1] * self.plg[6][1])
                            * (sr * (sec - p[76 - 1])).cos();
                } else {
                    if p[52 - 1] != 0.0 {
                        let apt1 = self.apt[1];
                        t_local[13] = apt1 * self.swc[11] * (1.0 + p[133 - 1] * self.plg[2][1])
                            * ((p[53 - 1] * self.plg[3][2] + p[99 - 1] * self.plg[5][2] + p[68 - 1] * self.plg[7][2])
                                * (dgtr * (long - p[98 - 1])).cos())
                            + apt1 * self.swc[11] * self.swc[5]
                                * (p[134 - 1] * self.plg[2][2] + p[135 - 1] * self.plg[4][2] + p[136 - 1] * self.plg[6][2])
                                * self.globe7_cd14 * (dgtr * (long - p[137 - 1])).cos()
                            + apt1 * self.swc[12]
                                * (p[56 - 1] * self.plg[2][1] + p[57 - 1] * self.plg[4][1] + p[58 - 1] * self.plg[6][1])
                                * (sr * (sec - p[59 - 1])).cos();
                    }
                }
            }
        }
        
        let mut tinf_val = p[31 - 1];
        for i in 1..=14 {
            tinf_val += self.sw[i].abs() * t_local[i];
            self.ttest_t[i] = t_local[i];
        }
        self.ttest_tinf = tinf_val;
        
        tinf_val
    }

    pub fn glob7s(&mut self, p: &[f32]) -> f32 {
        let dr = 1.72142e-2;
        let dgtr = 1.74533e-2;
        let pset = 2.0;
        let day = self.day;
        
        let mut p_val = p[100 - 1];
        if p_val == 0.0 {
            p_val = pset;
        }
        if p_val != pset {
            if self.mess {
                eprintln!("WRONG PARAMETER SET FOR GLOB7S");
            }
            panic!("WRONG PARAMETER SET FOR GLOB7S");
        }
        
        if day != self.glob7s_dayl || p[32 - 1] != self.glob7s_p32 {
            self.glob7s_cd32 = (dr * (day - p[32 - 1])).cos();
        }
        if day != self.glob7s_dayl || p[18 - 1] != self.glob7s_p18 {
            self.glob7s_cd18 = (2.0 * dr * (day - p[18 - 1])).cos();
        }
        if day != self.glob7s_dayl || p[14 - 1] != self.glob7s_p14 {
            self.glob7s_cd14 = (dr * (day - p[14 - 1])).cos();
        }
        if day != self.glob7s_dayl || p[39 - 1] != self.glob7s_p39 {
            self.glob7s_cd39 = (2.0 * dr * (day - p[39 - 1])).cos();
        }
        
        self.glob7s_dayl = day;
        self.glob7s_p32 = p[32 - 1];
        self.glob7s_p18 = p[18 - 1];
        self.glob7s_p14 = p[14 - 1];
        self.glob7s_p39 = p[39 - 1];
        
        let mut t_local = [0.0; 15];
        let dfa = self.dfa;
        let apdf = self.apdf;
        let apt1 = self.apt[1];
        
        t_local[1] = p[22 - 1] * dfa;
        
        t_local[2] = p[2 - 1] * self.plg[3][1] + p[3 - 1] * self.plg[5][1] + p[23 - 1] * self.plg[7][1]
            + p[27 - 1] * self.plg[2][1] + p[15 - 1] * self.plg[4][1] + p[60 - 1] * self.plg[6][1];
            
        t_local[3] = (p[19 - 1] + p[48 - 1] * self.plg[3][1] + p[30 - 1] * self.plg[5][1]) * self.glob7s_cd32;
        
        t_local[4] = (p[16 - 1] + p[17 - 1] * self.plg[3][1] + p[31 - 1] * self.plg[5][1]) * self.glob7s_cd18;
        
        t_local[5] = (p[10 - 1] * self.plg[2][1] + p[11 - 1] * self.plg[4][1] + p[21 - 1] * self.plg[6][1]) * self.glob7s_cd14;
        
        t_local[6] = (p[38 - 1] * self.plg[2][1]) * self.glob7s_cd39;
        
        if self.sw[7] != 0.0 {
            let t71 = p[12 - 1] * self.plg[3][2] * self.glob7s_cd14 * self.swc[5];
            let t72 = p[13 - 1] * self.plg[3][2] * self.glob7s_cd14 * self.swc[5];
            t_local[7] = (p[4 - 1] * self.plg[2][2] + p[5 - 1] * self.plg[4][2] + t71) * self.ctloc
                + (p[7 - 1] * self.plg[2][2] + p[8 - 1] * self.plg[4][2] + t72) * self.stloc;
        }
        
        if self.sw[8] != 0.0 {
            let t81 = (p[24 - 1] * self.plg[4][3] + p[36 - 1] * self.plg[6][3]) * self.glob7s_cd14 * self.swc[5];
            let t82 = (p[34 - 1] * self.plg[4][3] + p[37 - 1] * self.plg[6][3]) * self.glob7s_cd14 * self.swc[5];
            t_local[8] = (p[6 - 1] * self.plg[3][3] + p[42 - 1] * self.plg[5][3] + t81) * self.c2tloc
                + (p[9 - 1] * self.plg[3][3] + p[43 - 1] * self.plg[5][3] + t82) * self.s2tloc;
        }
        
        if self.sw[14] != 0.0 {
            t_local[14] = p[40 - 1] * self.plg[4][4] * self.s3tloc + p[41 - 1] * self.plg[4][4] * self.c3tloc;
        }
        
        if self.sw[9] != 0.0 {
            if self.sw[9] == 1.0 {
                t_local[9] = apdf * (p[33 - 1] + p[46 - 1] * self.plg[3][1] * self.swc[2]);
            } else if self.sw[9] == -1.0 {
                t_local[9] = p[51 - 1] * apt1 + p[97 - 1] * self.plg[3][1] * apt1 * self.swc[2];
            }
        }
        
        if self.sw[10] != 0.0 && self.sw[11] != 0.0 && self.xlong > -1000.0 {
            t_local[11] = (1.0 + self.plg[2][1] * (p[81 - 1] * self.swc[5] * (dr * (day - p[82 - 1])).cos()
                + p[86 - 1] * self.swc[6] * (2.0 * dr * (day - p[87 - 1])).cos())
                + p[84 - 1] * self.swc[3] * (dr * (day - p[85 - 1])).cos()
                + p[88 - 1] * self.swc[4] * (2.0 * dr * (day - p[89 - 1])).cos())
                * ((p[65 - 1] * self.plg[3][2] + p[66 - 1] * self.plg[5][2] + p[67 - 1] * self.plg[7][2]
                    + p[75 - 1] * self.plg[2][2] + p[76 - 1] * self.plg[4][2] + p[77 - 1] * self.plg[6][2]) * (dgtr * self.xlong).cos()
                + (p[91 - 1] * self.plg[3][2] + p[92 - 1] * self.plg[5][2] + p[93 - 1] * self.plg[7][2]
                    + p[78 - 1] * self.plg[2][2] + p[79 - 1] * self.plg[4][2] + p[80 - 1] * self.plg[6][2]) * (dgtr * self.xlong).sin());
        }
        
        let mut tt = 0.0;
        for i in 1..=14 {
            tt += self.sw[i].abs() * t_local[i];
        }
        
        tt
    }

    pub fn gts7(
        &mut self,
        iyd: i32,
        sec: f32,
        alt: f32,
        glat: f32,
        glong: f32,
        stl: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        mass: i32,
        d: &mut [f32],
        t: &mut [f32],
    ) {
        let alpha = [-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0];
        let mut tnmod = 0.0;
        let mut b28 = 0.0;
        if d[0] < 0.0 {
            tnmod = -d[0];
        }
        for j in 1..=9 {
            d[j - 1] = 0.0;
        }
        
        let v2 = self.vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 2);
        
        let yrd = iyd;
        self.za = self.pdl(16, 2);
        self.zn1[1] = self.za;
        
        // VARIATIONS NOT IMPORTANT BELOW ZA
        let mut tinf = self.gts7_tinf;
        if alt > self.zn1[1] {
            if v2 == 1.0 || self.gts7_alast <= self.zn1[1] {
                tinf = PTM[0] * PT[0] * (1.0 + self.sw[16] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PT));
                self.gts7_tinf = tinf;
            }
        } else {
            tinf = PTM[0] * PT[0];
            self.gts7_tinf = tinf;
        }
        
        if tnmod > 600.0 {
            tinf = tnmod;
            self.gts7_tinf = tinf;
        }
        
        t[0] = tinf;
        
        // GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
        if alt > self.zn1[5] {
            if v2 == 1.0 || self.gts7_alast <= self.zn1[5] {
                self.g0 = PTM[3] * PS[0] * (1.0 + self.sw[19] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PS));
            }
        } else {
            self.g0 = PTM[3] * PS[0];
        }
        
        // Calculate these temperatures only if input changed
        if v2 == 1.0 || alt < 300.0 {
            self.tlb = PTM[1] * (1.0 + self.sw[17] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA4)) * PDA4[0];
        }
        self.s = self.g0 / (tinf - self.tlb);
        
        // Lower thermosphere temp variations not significant for density above 300 km
        if alt < 300.0 {
            if v2 == 1.0 || self.gts7_alast >= 300.0 {
                self.tn1[2] = PTM[6] * PTL[0] / (1.0 - self.sw[18] * self.glob7s(&PTL[0..100]));
                self.tn1[3] = PTM[2] * PTL[100] / (1.0 - self.sw[18] * self.glob7s(&PTL[100..200]));
                self.tn1[4] = PTM[7] * PTL[200] / (1.0 - self.sw[18] * self.glob7s(&PTL[200..300]));
                self.tn1[5] = PTM[4] * PTL[300] / (1.0 - self.sw[18] * self.sw[20] * self.glob7s(&PTL[300..400]));
                self.tgn1[2] = PTM[8] * self.pma(1, 9) * (1.0 + self.sw[18] * self.sw[20] * self.glob7s(&PMA[800..900]))
                    * self.tn1[5] * self.tn1[5] / (PTM[4] * PTL[300]).powi(2);
            }
        } else {
            self.tn1[2] = PTM[6] * PTL[0];
            self.tn1[3] = PTM[2] * PTL[100];
            self.tn1[4] = PTM[7] * PTL[200];
            self.tn1[5] = PTM[4] * PTL[300];
            self.tgn1[2] = PTM[8] * self.pma(1, 9) * self.tn1[5] * self.tn1[5] / (PTM[4] * PTL[300]).powi(2);
        }
        
        self.z0 = self.zn1[4];
        self.t0 = self.tn1[4];
        self.tr12 = 1.0;
        
        let mut tz = 0.0;
        let _ddum = self.densu(alt.abs(), 1.0, tinf, self.tlb, 0.0, 0.0, &mut tz, PTM[5], self.s, 5);
        t[1] = tz;
        
        if mass == 0 {
            self.gts7_alast = alt;
            return;
        }
        
        let g28 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA3);
        let day = (yrd % 1000) as f32;
        let dgtr = 1.74533e-2;
        let dr = 1.72142e-2;
        
        let zhf = self.pdl(25, 2) * (1.0 + self.sw[5] * self.pdl(25, 1) * (dgtr * glat).sin() * (dr * (day - PT[14 - 1])).cos());
        self.gts7_alast = alt;
        let xmm = self.pdm(5, 3);
        let z = alt;
        
        let mt = [0, 48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17];
        let mut j_opt = None;
        for j in 1..=11 {
            if mass == mt[j] {
                j_opt = Some(j);
                break;
            }
        }
        if j_opt.is_none() {
            if self.mess {
                eprintln!("MASS {} NOT VALID", mass);
            }
            return;
        }
        let j = j_opt.unwrap();
        
        let altl = [0.0, 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0];
        let mut dd = 0.0;
        
        if z <= altl[6] || mass == 28 || mass == 48 {
            self.db28 = self.pdm(1, 3) * g28.exp() * PDA3[0];
            let mut tz = 0.0;
            d[2] = self.densu(z, self.db28, tinf, self.tlb, 28.0, alpha[3 - 1], &mut tz, PTM[6 - 1], self.s, 5);
            t[1] = tz;
            dd = d[2];
            
            let zh28 = self.pdm(3, 3) * zhf;
            let zhm28 = self.pdm(4, 3) * self.pdl(6, 2);
            let xmd = 28.0 - xmm;
            
            b28 = self.densu(zh28, self.db28, tinf, self.tlb, xmd, alpha[3 - 1] - 1.0, &mut tz, PTM[6 - 1], self.s, 5);
            if z <= altl[3] && self.sw[15] != 0.0 {
                self.dm28 = self.densu(z, b28, tinf, self.tlb, xmm, alpha[3 - 1], &mut tz, PTM[6 - 1], self.s, 5);
                t[1] = tz;
                d[2] = self.dnet(d[2], self.dm28, zhm28, xmm, 28.0);
            }
        }

        
        if j == 1 || j == 3 {
            let g4 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA1);
            self.db04 = self.pdm(1, 1) * g4.exp() * PDA1[0];
            let mut t2_val = t[1];
            d[0] = self.densu(z, self.db04, tinf, self.tlb, 4.0, alpha[1 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            t[1] = t2_val;
            dd = d[0];
            
            if z <= altl[1] && self.sw[15] != 0.0 {
                let zh04 = self.pdm(3, 1);
                let b04 = self.densu(zh04, self.db04, tinf, self.tlb, 4.0 - xmm, alpha[1 - 1] - 1.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                self.dm04 = self.densu(z, b04, tinf, self.tlb, xmm, 0.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                let zhm04 = self.pdm(4, 3) * self.pdl(6, 2);
                d[0] = self.dnet(d[0], self.dm04, zhm04, xmm, 4.0);
                
                let rl = (b28 * self.pdm(2, 1) / b04).ln();
                let zc04 = self.pdm(5, 1) * self.pdl(1, 2);
                let hc04 = self.pdm(6, 1) * self.pdl(2, 2);
                d[0] = d[0] * self.ccor(z, rl, hc04, zc04);
            }
            if mass != 48 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }

        
        if j == 1 || j == 4 || j == 9 {
            let g16 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA2);
            self.db16 = self.pdm(1, 2) * g16.exp() * PDA2[0];
            let mut t2_val = t[1];
            d[1] = self.densu(z, self.db16, tinf, self.tlb, 16.0, alpha[2 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            t[1] = t2_val;
            dd = d[1];
            
            if z <= altl[2] && self.sw[15] != 0.0 {
                let zh16 = self.pdm(3, 2);
                let b16 = self.densu(zh16, self.db16, tinf, self.tlb, 16.0 - xmm, alpha[2 - 1] - 1.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                self.dm16 = self.densu(z, b16, tinf, self.tlb, xmm, 0.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                let zhm16 = self.pdm(4, 3) * self.pdl(6, 2); // ZHM16 = ZHM28
                d[1] = self.dnet(d[1], self.dm16, zhm16, xmm, 16.0);
                
                let rl = self.pdm(2, 2) * self.pdl(17, 2) * (1.0 + self.sw[1] * self.pdl(24, 1) * (f107a - 150.0));
                let hc16 = self.pdm(6, 2) * self.pdl(4, 2);
                let zc16 = self.pdm(5, 2) * self.pdl(3, 2);
                let hc216 = self.pdm(6, 2) * self.pdl(5, 2);
                d[1] = d[1] * self.ccor2(z, rl, hc16, zc16, hc216);
                
                let hcc16 = self.pdm(8, 2) * self.pdl(14, 2);
                let zcc16 = self.pdm(7, 2) * self.pdl(13, 2);
                let rc16 = self.pdm(4, 2) * self.pdl(15, 2);
                d[1] = d[1] * self.ccor(z, rc16, hcc16, zcc16);
            }
            if mass != 48 && mass != 49 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }

        
        if j == 1 || j == 6 {
            let g32 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA5);
            self.db32 = self.pdm(1, 4) * g32.exp() * PDA5[0];
            let mut t2_val = t[1];
            d[3] = self.densu(z, self.db32, tinf, self.tlb, 32.0, alpha[4 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            t[1] = t2_val;
            dd = d[3];
            
            if self.sw[15] != 0.0 {
                if z <= altl[4] {
                    let zh32 = self.pdm(3, 4);
                    let b32 = self.densu(zh32, self.db32, tinf, self.tlb, 32.0 - xmm, alpha[4 - 1] - 1.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                    t[1] = t2_val;
                    self.dm32 = self.densu(z, b32, tinf, self.tlb, xmm, 0.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                    t[1] = t2_val;
                    let zhm32 = self.pdm(4, 3) * self.pdl(6, 2); // ZHM32 = ZHM28
                    d[3] = self.dnet(d[3], self.dm32, zhm32, xmm, 32.0);
                    
                    let rl = (b28 * self.pdm(2, 4) / b32).ln();
                    let hc32 = self.pdm(6, 4) * self.pdl(8, 2);
                    let zc32 = self.pdm(5, 4) * self.pdl(7, 2);
                    d[3] = d[3] * self.ccor(z, rl, hc32, zc32);
                }
                
                let hcc32 = self.pdm(8, 4) * self.pdl(23, 2);
                let hcc232 = self.pdm(8, 4) * self.pdl(23, 1);
                let zcc32 = self.pdm(7, 4) * self.pdl(22, 2);
                let rc32 = self.pdm(4, 4) * self.pdl(24, 2) * (1.0 + self.sw[1] * self.pdl(24, 1) * (f107a - 150.0));
                d[3] = d[3] * self.ccor2(z, rc32, hcc32, zcc32, hcc232);
            }
            if mass != 48 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }

        
        if j == 1 || j == 7 {
            let g40 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA6);
            self.db40 = self.pdm(1, 5) * g40.exp() * PDA6[0];
            let mut t2_val = t[1];
            d[4] = self.densu(z, self.db40, tinf, self.tlb, 40.0, alpha[5 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            t[1] = t2_val;
            dd = d[4];
            
            if z <= altl[5] && self.sw[15] != 0.0 {
                let zh40 = self.pdm(3, 5);
                let b40 = self.densu(zh40, self.db40, tinf, self.tlb, 40.0 - xmm, alpha[5 - 1] - 1.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                self.dm40 = self.densu(z, b40, tinf, self.tlb, xmm, 0.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                let zhm40 = self.pdm(4, 3) * self.pdl(6, 2);
                d[4] = self.dnet(d[4], self.dm40, zhm40, xmm, 40.0);
                
                let rl = (b28 * self.pdm(2, 5) / b40).ln();
                let hc40 = self.pdm(6, 5) * self.pdl(10, 2);
                let zc40 = self.pdm(5, 5) * self.pdl(9, 2);
                d[4] = d[4] * self.ccor(z, rl, hc40, zc40);
            }
            if mass != 48 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }

        
        if j == 1 || j == 8 {
            let g1 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA7);
            self.db01 = self.pdm(1, 6) * g1.exp() * PDA7[0];
            let mut t2_val = t[1];
            d[6] = self.densu(z, self.db01, tinf, self.tlb, 1.0, alpha[7 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            t[1] = t2_val;
            dd = d[6];
            
            if z <= altl[7] && self.sw[15] != 0.0 {
                let zh01 = self.pdm(3, 6);
                let b01 = self.densu(zh01, self.db01, tinf, self.tlb, 1.0 - xmm, alpha[7 - 1] - 1.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                self.dm01 = self.densu(z, b01, tinf, self.tlb, xmm, 0.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                let zhm01 = self.pdm(4, 3) * self.pdl(6, 2);
                d[6] = self.dnet(d[6], self.dm01, zhm01, xmm, 1.0);
                
                let rl = (b28 * self.pdm(2, 6) * self.pdl(18, 2).abs() / b01).ln();
                let hc01 = self.pdm(6, 6) * self.pdl(12, 2);
                let zc01 = self.pdm(5, 6) * self.pdl(11, 2);
                d[6] = d[6] * self.ccor(z, rl, hc01, zc01);
                
                let hcc01 = self.pdm(8, 6) * self.pdl(20, 2);
                let zcc01 = self.pdm(7, 6) * self.pdl(19, 2);
                let rc01 = self.pdm(4, 6) * self.pdl(21, 2);
                d[6] = d[6] * self.ccor(z, rc01, hcc01, zcc01);
            }
            if mass != 48 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }

        
        if j == 1 || j == 10 {
            let g14 = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA8);
            self.db14 = self.pdm(1, 7) * g14.exp() * PDA8[0];
            let mut t2_val = t[1];
            d[7] = self.densu(z, self.db14, tinf, self.tlb, 14.0, alpha[8 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            t[1] = t2_val;
            dd = d[7];
            
            if z <= altl[8] && self.sw[15] != 0.0 {
                let zh14 = self.pdm(3, 7);
                let b14 = self.densu(zh14, self.db14, tinf, self.tlb, 14.0 - xmm, alpha[8 - 1] - 1.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                self.dm14 = self.densu(z, b14, tinf, self.tlb, xmm, 0.0, &mut t2_val, PTM[6 - 1], self.s, 5);
                t[1] = t2_val;
                let zhm14 = self.pdm(4, 3) * self.pdl(6, 2);
                d[7] = self.dnet(d[7], self.dm14, zhm14, xmm, 14.0);
                
                let rl = (b28 * self.pdm(2, 7) * self.pdl(3, 1).abs() / b14).ln();
                let hc14 = self.pdm(6, 7) * self.pdl(2, 1);
                let zc14 = self.pdm(5, 7) * self.pdl(1, 1);
                d[7] = d[7] * self.ccor(z, rl, hc14, zc14);
                
                let hcc14 = self.pdm(8, 7) * self.pdl(5, 1);
                let zcc14 = self.pdm(7, 7) * self.pdl(4, 1);
                let rc14 = self.pdm(4, 7) * self.pdl(6, 1);
                d[7] = d[7] * self.ccor(z, rc14, hcc14, zcc14);
            }
            if mass != 48 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }

        
        if j == 1 || j == 11 {
            let g16h = self.sw[21] * self.globe7(iyd, sec, glat, glong, stl, f107a, f107, ap, &PDA9);
            let db16h = self.pdm(1, 8) * g16h.exp() * PDA9[0];
            let tho = self.pdm(10, 8) * self.pdl(7, 1);
            let mut t2_val = 0.0;
            let dd_anom = self.densu(z, db16h, tho, tho, 16.0, alpha[9 - 1], &mut t2_val, PTM[6 - 1], self.s, 5);
            let zsht = self.pdm(6, 8);
            let zmho = self.pdm(5, 8);
            let zsho = self.scalh(zmho, 16.0, tho);
            d[8] = dd_anom * (-zsht / zsho * ((-(z - zmho) / zsht).exp() - 1.0)).exp();
            
            if mass != 48 {
                self.gts7_adjust_and_exit(mass, d, t);
                return;
            }
        }
        
        if mass == 48 {
            d[5] = 1.66e-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7]);
            self.db48 = 1.66e-24 * (4.0 * self.db04 + 16.0 * self.db16 + 28.0 * self.db28 + 32.0 * self.db32 + 40.0 * self.db40 + self.db01 + 14.0 * self.db14);
        }
        
        self.gts7_adjust_and_exit(mass, d, t);
    }

    pub fn gtd7(
        &mut self,
        iyd: i32,
        sec: f32,
        alt: f32,
        glat: f32,
        glong: f32,
        stl: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        mass: i32,
        d: &mut [f32],
        t: &mut [f32],
    ) {
        let v1 = self.vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 1);
        let mut xlat = glat;
        if self.sw[2] == 0.0 {
            xlat = 45.0;
        }
        self.glatf(xlat);
        
        let xmm = self.pdm(5, 3);
        let altt = alt.max(self.zn2[1]);
        let mut mss = mass;
        let zmix = 62.5;
        if alt < zmix && mass > 0 {
            mss = 28;
        }
        
        if v1 == 1.0 || alt > self.zn2[1] || self.gtd7_alast > self.zn2[1] || mss != self.gtd7_mssl {
            self.gtd7_ds[1] = d[0];
            let mut ds_tmp = [0.0; 10];
            let mut ts_tmp = [0.0; 3];
            self.gts7(iyd, sec, altt, glat, glong, stl, f107a, f107, ap, mss, &mut ds_tmp, &mut ts_tmp);
            for i in 1..=9 {
                self.gtd7_ds[i] = ds_tmp[i - 1];
            }
            self.gtd7_ts[1] = ts_tmp[0];
            self.gtd7_ts[2] = ts_tmp[1];
            self.gtd7_mssl = mss;
        }
        
        t[0] = self.gtd7_ts[1];
        t[1] = self.gtd7_ts[2];
        
        if alt >= self.zn2[1] {
            for j in 1..=9 {
                d[j - 1] = self.gtd7_ds[j];
            }
            self.gtd7_alast = alt;
            return;
        }
        
        if v1 == 1.0 || self.gtd7_alast >= self.zn2[1] {
            self.tgn2[1] = self.tgn1[2];
            self.tn2[1] = self.tn1[5];
            self.tn2[2] = self.pma(1, 1) * PAVGM[1 - 1] / (1.0 - self.sw[20] * self.glob7s(&PMA[0..100]));
            self.tn2[3] = self.pma(1, 2) * PAVGM[2 - 1] / (1.0 - self.sw[20] * self.glob7s(&PMA[100..200]));
            self.tn2[4] = self.pma(1, 3) * PAVGM[3 - 1] / (1.0 - self.sw[20] * self.sw[22] * self.glob7s(&PMA[200..300]));
            self.tgn2[2] = PAVGM[9 - 1] * self.pma(1, 10) * (1.0 + self.sw[20] * self.sw[22] * self.glob7s(&PMA[900..1000]))
                * self.tn2[4] * self.tn2[4] / (self.pma(1, 3) * PAVGM[3 - 1]).powi(2);
            self.tn3[1] = self.tn2[4];
        }
        
        if alt < self.zn3[1] {
            if v1 == 1.0 || self.gtd7_alast >= self.zn3[1] {
                self.tgn3[1] = self.tgn2[2];
                self.tn3[2] = self.pma(1, 4) * PAVGM[4 - 1] / (1.0 - self.sw[22] * self.glob7s(&PMA[300..400]));
                self.tn3[3] = self.pma(1, 5) * PAVGM[5 - 1] / (1.0 - self.sw[22] * self.glob7s(&PMA[400..500]));
                self.tn3[4] = self.pma(1, 6) * PAVGM[6 - 1] / (1.0 - self.sw[22] * self.glob7s(&PMA[500..600]));
                self.tn3[5] = self.pma(1, 7) * PAVGM[7 - 1] / (1.0 - self.sw[22] * self.glob7s(&PMA[600..700]));
                self.tgn3[2] = self.pma(1, 8) * PAVGM[8 - 1] * (1.0 + self.sw[22] * self.glob7s(&PMA[700..800]))
                    * self.tn3[5] * self.tn3[5] / (self.pma(1, 7) * PAVGM[7 - 1]).powi(2);
            }
        }
        
        if mass == 0 {
            let mut tz = 0.0;
            let _dd = self.densm(alt, 1.0, 0.0, &mut tz, 5, 4);
            t[1] = tz;
            self.gtd7_alast = alt;
            return;
        }
        
        let mut dmc = 0.0;
        if alt > zmix {
            dmc = 1.0 - (self.zn2[1] - alt) / (self.zn2[1] - zmix);
        }
        let dz28 = self.gtd7_ds[3];
        let mut dm28m = self.dm28;
        if self.imr == 1 {
            dm28m = self.dm28 * 1e6;
        }
        
        let mut tz = 0.0;
        let dmr = self.gtd7_ds[3] / dm28m - 1.0;
        d[2] = self.densm(alt, dm28m, xmm, &mut tz, 5, 4);
        d[2] = d[2] * (1.0 + dmr * dmc);
        
        d[0] = 0.0;
        if mass == 4 || mass == 48 {
            let dmr = self.gtd7_ds[1] / (dz28 * self.pdm(2, 1)) - 1.0;
            d[0] = d[2] * self.pdm(2, 1) * (1.0 + dmr * dmc);
        }
        
        d[1] = 0.0;
        d[8] = 0.0;
        
        d[3] = 0.0;
        if mass == 32 || mass == 48 {
            let dmr = self.gtd7_ds[4] / (dz28 * self.pdm(2, 4)) - 1.0;
            d[3] = d[2] * self.pdm(2, 4) * (1.0 + dmr * dmc);
        }
        
        d[4] = 0.0;
        if mass == 40 || mass == 48 {
            let dmr = self.gtd7_ds[5] / (dz28 * self.pdm(2, 5)) - 1.0;
            d[4] = d[2] * self.pdm(2, 5) * (1.0 + dmr * dmc);
        }
        
        d[6] = 0.0;
        d[7] = 0.0;
        
        if mass == 48 {
            d[5] = 1.66e-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7]);
            if self.imr == 1 {
                d[5] = d[5] / 1000.0;
            }
        }
        t[1] = tz;
        self.gtd7_alast = alt;
    }

    pub fn gtd7d(
        &mut self,
        iyd: i32,
        sec: f32,
        alt: f32,
        glat: f32,
        glong: f32,
        stl: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        mass: i32,
        d: &mut [f32],
        t: &mut [f32],
    ) {
        self.gtd7(iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, d, t);
        if mass == 48 {
            d[5] = 1.66e-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7] + 16.0 * d[8]);
            if self.imr == 1 {
                d[5] = d[5] / 1000.0;
            }
        }
    }

    pub fn ghp7(
        &mut self,
        iyd: i32,
        sec: f32,
        alt: &mut f32,
        glat: f32,
        glong: f32,
        stl: f32,
        f107a: f32,
        f107: f32,
        ap: &[f32],
        d: &mut [f32],
        t: &mut [f32],
        press: f32,
    ) {
        let mut d_tmp = [0.0; 9];
        let mut t_tmp = [0.0; 2];
        let bm = 1.3806e-19;
        let rgas = 831.4;
        let test = 0.00043;
        let ltest = 12;
        let pl = press.log10();
        let mut z = 0.0;
        if pl >= -5.0 {
            let mut zi = 0.0;
            if pl > 2.5 { zi = 18.06 * (3.00 - pl); }
            else if pl > 0.75 && pl <= 2.5 { zi = 14.98 * (3.08 - pl); }
            else if pl > -1.0 && pl <= 0.75 { zi = 17.8 * (2.72 - pl); }
            else if pl > -2.0 && pl <= -1.0 { zi = 14.28 * (3.64 - pl); }
            else if pl > -4.0 && pl <= -2.0 { zi = 12.72 * (4.32 - pl); }
            else { zi = 25.3 * (0.11 - pl); }
            let iday = iyd % 1000;
            let cl = glat / 90.0;
            let cl2 = cl * cl;
            let mut cd = 0.0;
            if iday < 182 { cd = 1.0 - (iday as f32) / 91.25; }
            else { cd = (iday as f32) / 91.25 - 3.0; }
            let mut ca = 0.0;
            if pl > -1.11 && pl <= -0.23 { ca = 1.0; }
            else if pl > -0.23 { ca = (2.79 - pl) / (2.79 + 0.23); }
            else if pl <= -1.11 && pl > -3.0 { ca = (-2.93 - pl) / (-2.93 + 1.11); }
            z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl;
        } else {
            z = 22.0 * (pl + 4.0).powi(2) + 110.0;
        }
        
        let mut l = 0;
        loop {
            l += 1;
            self.gtd7(iyd, sec, z, glat, glong, stl, f107a, f107, ap, 48, &mut d_tmp, &mut t_tmp);
            let xn = d_tmp[0] + d_tmp[1] + d_tmp[2] + d_tmp[3] + d_tmp[4] + d_tmp[6] + d_tmp[7];
            let mut p = bm * xn * t_tmp[1];
            if self.imr == 1 {
                p = p * 1e-6;
            }
            let diff = pl - p.log10();
            if diff.abs() < test || l == ltest {
                if self.mess && l == ltest {
                    eprintln!("GHP7 NOT CONVERGING FOR PRESS {:.2e} {:.2e}", press, diff);
                }
                for i in 0..9 { d[i] = d_tmp[i]; }
                for i in 0..2 { t[i] = t_tmp[i]; }
                break;
            }
            let mut xm = d_tmp[5] / xn / 1.66e-24;
            if self.imr == 1 {
                xm = xm * 1e3;
            }
            let g = self.gsurf / (1.0 + z / self.re).powi(2);
            let sh = rgas * t_tmp[1] / (xm * g);
            if l < 6 {
                z = z - sh * diff * 2.302;
            } else {
                z = z - sh * diff;
            }
        }
        *alt = z;
    }

    fn gts7_adjust_and_exit(&self, _mass: i32, d: &mut [f32], _t: &mut [f32]) {
        if self.imr == 1 {
            for i in 0..9 {
                d[i] = d[i] * 1e6;
            }
            d[5] = d[5] / 1000.0;
        }
    }


}

pub fn cira_gtd7(
    iyd: i32,
    sec: f32,
    alt: f32,
    glat: f32,
    glong: f32,
    stl: f32,
    f107a: f32,
    f107: f32,
    ap: &[f32],
    mass: i32,
    d: &mut [f32],
    t: &mut [f32],
) {
    let mut model = CiraModel::new();
    model.gtd7(iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, d, t);
}

pub fn cira_gtd7d(
    iyd: i32,
    sec: f32,
    alt: f32,
    glat: f32,
    glong: f32,
    stl: f32,
    f107a: f32,
    f107: f32,
    ap: &[f32],
    mass: i32,
    d: &mut [f32],
    t: &mut [f32],
) {
    let mut model = CiraModel::new();
    model.gtd7d(iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, d, t);
}
