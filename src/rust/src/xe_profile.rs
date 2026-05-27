use crate::irifun_utils::ARGMAX;
use std::f32::consts::PI;

const UMR: f32 = PI / 180.0;

// --- F1 Layer Probability ---

pub fn f1_prob(sza: f32, glat: f32, rz12: f32) -> (f32, f32) {
    let xarg = 0.5 + 0.5 * (sza * UMR).cos();
    let a = 2.98 + 0.0854 * rz12;
    let b = 0.0107 - 0.0022 * rz12;
    let c = -0.000256 + 0.0000147 * rz12;
    let gamma = a + (b + c * glat) * glat;
    
    let mut f1prob = if xarg > 0.0 { xarg.powf(gamma) } else { 0.0 };
    if f1prob < 1e-3 {
        f1prob = 0.0;
    }
    
    let mut f1probl = if xarg > 0.0 { xarg.powf(2.36) } else { 0.0 };
    if f1probl < 1e-3 {
        f1probl = 0.0;
    }
    
    (f1prob, f1probl)
}

// --- Electron Density Profile Models ---

pub struct XeProfile {
    // Block 1
    pub hmf2: f32,
    pub xnmf2: f32,
    pub hmf1: f32,
    pub f1reg: bool,
    // Block 2
    pub b0: f32,
    pub b1: f32,
    pub c1: f32,
    // Block 3
    pub hz: f32,
    pub t: f32,
    pub hst: f32,
    // Block 4
    pub hme: f32,
    pub xnme: f32,
    pub hef: f32,
    // Block 5
    pub night: bool,
    pub e: [f32; 4],
    // Block 6
    pub hmd: f32,
    pub xnmd: f32,
    pub hdx: f32,
    // Block 7
    pub d1: f32,
    pub xkk: f32,
    pub fp30: f32,
    pub fp3u: f32,
    pub fp1: f32,
    pub fp2: f32,
    // Block 10
    pub beta: f32,
    pub eta: f32,
    pub delta: f32,
    pub zeta: f32,
    // Block 11
    pub b2top: f32,
    pub itopn: i32,
    pub tcor1: f32,
    pub tcor2: f32,
}

impl XeProfile {
    pub fn topq(&self, h: f32, no: f32, hmax: f32, ho: f32) -> f32 {
        let g = 0.125;
        let rfac = 100.0;
        let dh = h - hmax;
        let g1 = g * dh;
        let denom = rfac * ho + g1;
        let z = if denom != 0.0 {
            dh / (ho * (1.0 + rfac * g1 / denom))
        } else {
            0.0
        };
        if z > 40.0 {
            return 0.0;
        }
        let ee = z.exp();
        let ep = if ee > 1.0e7 {
            4.0 / ee
        } else {
            4.0 * ee / (1.0 + ee).powi(2)
        };
        no * ep
    }

    pub fn zero(&self, delta: f32) -> f32 {
        let arg1 = delta / 100.0;
        let z1 = if arg1.abs() < ARGMAX {
            1.0 / (1.0 + arg1.exp())
        } else if arg1 < 0.0 {
            1.0
        } else {
            0.0
        };

        let arg2 = (delta + 94.5) / self.beta;
        let z2 = if arg2.abs() < ARGMAX {
            1.0 / (1.0 + arg2.exp())
        } else if arg2 < 0.0 {
            1.0
        } else {
            0.0
        };

        self.zeta * (1.0 - z1) - self.eta * z2
    }

    pub fn dxe1n(&self, h: f32) -> f32 {
        use crate::irifun_utils::epst;
        let x0 = 300.0 - self.delta;
        let x = (h - self.hmf2) / (1000.0 - self.hmf2) * 700.0 + x0;
        let epst2 = epst(x, 100.0, 300.0);
        let epst1 = epst(x, self.beta, 394.5);
        -self.eta * epst1 + self.zeta * (1.0 - epst2)
    }

    pub fn xe1(&self, h: f32) -> f32 {
        use crate::irifun_utils::eptr;
        if self.itopn == 2 {
            return self.topq(h, self.xnmf2, self.hmf2, self.b2top);
        }

        let dxdh = (1000.0 - self.hmf2) / 700.0;
        let x0 = 300.0 - self.delta;
        let xmx0 = (h - self.hmf2) / dxdh;
        let x = xmx0 + x0;
        let eptr1 = eptr(x, self.beta, 394.5) - eptr(x0, self.beta, 394.5);
        let eptr2 = eptr(x, 100.0, 300.0) - eptr(x0, 100.0, 300.0);
        let y = self.beta * self.eta * eptr1 + self.zeta * (100.0 * eptr2 - xmx0);
        let y = y * dxdh;
        let mut yc = -y + self.tcor1 + self.tcor2;
        if yc.abs() > ARGMAX {
            yc = yc.signum() * ARGMAX;
        }
        self.xnmf2 * yc.exp()
    }

    pub fn xe2(&self, h: f32) -> f32 {
        let mut x = (self.hmf2 - h) / self.b0;
        if x <= 0.0 {
            x = 0.0;
        }
        let mut z = x.powf(self.b1);
        if z > ARGMAX {
            z = ARGMAX;
        }
        self.xnmf2 * (-z).exp() / x.cosh()
    }

    pub fn xe3_1(&self, h: f32) -> f32 {
        let mut h1bar = h;
        if self.f1reg {
            let ratio = (self.hmf1 - h) / self.hmf1;
            let ratio_clamped = if ratio < 0.0 { 0.0 } else { ratio };
            h1bar = self.hmf1 * (1.0 - ratio_clamped.powf(1.0 + self.c1));
        }
        self.xe2(h1bar)
    }

    pub fn xe4_1(&self, h: f32) -> f32 {
        if self.hst < 0.0 {
            return self.xnme + self.t * (h - self.hef);
        }
        let haha = if self.hst == self.hef {
            h
        } else {
            let d = self.hz - self.hst;
            let t_local = d * d / (self.hst - self.hef);
            let term = t_local * (0.25 * t_local + self.hz - h);
            let term_clamped = if term < 0.0 { 0.0 } else { term };
            if self.hst > self.hef {
                self.hz + 0.5 * t_local - term_clamped.sqrt()
            } else {
                self.hz + 0.5 * t_local + term_clamped.sqrt()
            }
        };
        let mut h2bar = haha;
        if self.f1reg {
            let ratio = (self.hmf1 - haha) / self.hmf1;
            let ratio_clamped = if ratio < 0.0 { 0.0 } else { ratio };
            h2bar = self.hmf1 * (1.0 - ratio_clamped.powf(1.0 + self.c1));
        }
        self.xe2(h2bar)
    }

    pub fn xe5(&self, h: f32) -> f32 {
        let t3 = h - self.hme;
        let t1 = t3 * t3 * (self.e[0] + t3 * (self.e[1] + t3 * (self.e[2] + t3 * self.e[3])));
        if self.night {
            self.xnme * t1.exp()
        } else {
            self.xnme * (1.0 + t1)
        }
    }

    pub fn xe6(&self, h: f32) -> f32 {
        if h > self.hdx {
            let z = self.hme - h;
            let base = if z < 0.0 { 0.0 } else { z };
            self.xnme * (-self.d1 * base.powf(self.xkk)).exp()
        } else {
            let z = h - self.hmd;
            let fp3 = if z > 0.0 { self.fp30 } else { self.fp3u };
            self.xnmd * (z * (self.fp1 + z * (self.fp2 + z * fp3))).exp()
        }
    }

    pub fn xe_1(&self, h: f32) -> f32 {
        let hmf1 = if self.f1reg { self.hmf1 } else { self.hmf2 };
        if h >= self.hmf2 {
            self.xe1(h)
        } else if h >= hmf1 {
            self.xe2(h)
        } else if h >= self.hz {
            self.xe3_1(h)
        } else if h >= self.hef {
            self.xe4_1(h)
        } else if h >= self.hme {
            self.xe5(h)
        } else {
            self.xe6(h)
        }
    }
}

pub fn xe2to5(h: f32, hmf2: f32, nl: usize, hx: &[f32], sc: &[f32], amp: &[f32]) -> f32 {
    let mut sum = 1.0_f32;
    for i in 0..nl {
        let ylay = amp[i] * crate::irifun_utils::rlay(h, hmf2, sc[i], hx[i]);
        let zlay = 10.0_f32.powf(ylay);
        sum *= zlay;
    }
    sum
}

pub fn xen(
    h: f32,
    hmf2: f32,
    xnmf2: f32,
    hme: f32,
    nl: usize,
    hx: &[f32],
    sc: &[f32],
    amp: &[f32],
    profile: &XeProfile,
) -> f32 {
    if h >= hmf2 {
        profile.xe1(h)
    } else if h >= hme {
        xnmf2 * xe2to5(h, hmf2, nl, hx, sc, amp)
    } else {
        profile.xe6(h)
    }
}

pub fn fof1ed(ylati: f32, r: f32, chi: f32) -> f32 {
    if chi > 90.0 {
        return 0.0;
    }
    let umr = std::f32::consts::PI / 180.0;
    let dla = ylati;
    let f0 = 4.35 + dla * (0.0058 - 1.2e-4 * dla);
    let f100 = 5.348 + dla * (0.011 - 2.3e-4 * dla);
    let fs = f0 + (f100 - f0) * r / 100.0;
    let xmue = 0.093 + dla * (0.0046 - 5.4e-5 * dla) + 3.0e-4 * r;
    let mut fof1 = fs * (chi * umr).cos().powf(xmue);
    let chi0 = 49.84733 + 0.349504 * dla;
    let chi100 = 38.96113 + 0.509932 * dla;
    let chim = chi0 + (chi100 - chi0) * r / 100.0;
    if chi > chim {
        fof1 = -fof1;
    }
    fof1
}

pub fn f1_c1(absmdp: f32, hour: f32, saxnon: f32, suxnon: f32) -> f32 {
    let pi = std::f32::consts::PI;
    let mut dela = 4.32_f32;
    if absmdp >= 18.0 {
        dela = 1.0 + (-(absmdp - 30.0) / 10.0).exp();
    }
    let c1old = 0.09 + 0.11 / dela;
    let mut c1 = if suxnon == saxnon {
        2.5 * c1old
    } else {
        2.5 * c1old * ((hour - 12.0) / (suxnon - saxnon) * pi).cos()
    };
    if c1 < 0.0 {
        c1 = 0.0;
    }
    c1
}

pub fn rogul(iday: i32, xhi: f32) -> (f32, f32) {
    let dumr = std::f32::consts::PI / 182.5;
    let sx = 2.0 - (iday as f32 * dumr).cos();
    let xs = (xhi - 20.0 * sx) / 15.0;
    let gro = 0.8 - 0.2 / (1.0 + xs.exp());
    (sx, gro)
}

pub fn tal(
    shabr: f32,
    sdelta: f32,
    shbr: f32,
    sdtdh0: f32,
    aus6: &mut bool,
    spt: &mut [f32; 4],
) {
    *aus6 = false;
    if shbr <= 0.0 {
        *aus6 = true;
        return;
    }
    let mut z1 = -sdelta / (100.0 * shabr * shabr);
    let mut sdelta_abs = sdelta;
    if sdelta <= 0.0 {
        sdelta_abs = -sdelta;
        z1 = (1.0 - sdelta_abs / 100.0).ln() / (shabr * shabr);
    }
    let mut z3 = sdtdh0 / (2.0 * shbr);
    let z4 = shabr - shbr;
    spt[3] = 2.0 * (z1 * (shbr - 2.0 * shabr) * shbr + z3 * z4 * shabr)
        / (shabr * shabr * z4 * z4 * z4);
    spt[2] = z1 * (2.0 * shbr - 3.0 * shabr) / (shabr * z4 * z4) - (2.0 * shabr + shbr) * spt[3];
    spt[1] = -2.0 * z1 / shabr - 2.0 * shabr * spt[2] - 3.0 * shabr * shabr * spt[3];
    spt[0] = z1 - shabr * (spt[1] + shabr * (spt[2] + shabr * spt[3]));
    let b = 4.0 * spt[2] / (5.0 * spt[3]) + shabr;
    let c = -2.0 * spt[0] / (5.0 * spt[3] * shabr);
    let z2_sq = b * b / 4.0 - c;
    if z2_sq < 0.0 {
        return;
    }
    z3 = z2_sq.sqrt();
    let z1_val = b / 2.0;
    let mut z2 = -z1_val + z3;
    if z2 > 0.0 && z2 < shbr {
        *aus6 = true;
    }
    if z3.abs() > 1.0e-15 {
        z2 = -z1_val - z3;
        if z2 > 0.0 && z2 < shbr {
            *aus6 = true;
        }
        return;
    }
    z2 = c / z2;
    if z2 > 0.0 && z2 < shbr {
        *aus6 = true;
    }
}

pub fn lsknm(
    n: usize,
    m: usize,
    m0: usize,
    m1: usize,
    hm: f32,
    sc: &[f32],
    hx: &[f32],
    w: &[f32],
    x: &[f32],
    y: &[f32],
    var: &mut [f32],
    sing: &mut bool,
) {
    let m01 = m0 + m1;
    let mut bli = [0.0_f32; 5];
    let mut ali = [[0.0_f32; 5]; 5];
    let mut xli = [[0.0_f32; 10]; 5];
    
    for i in 0..n {
        for k in 0..m0 {
            xli[i][k] = crate::irifun_utils::rlay(x[k], hm, sc[i], hx[i]);
        }
        for k in m0..m01 {
            xli[i][k] = crate::irifun_utils::d1lay(x[k], hm, sc[i], hx[i]);
        }
        for k in m01..m {
            xli[i][k] = crate::irifun_utils::d2lay(x[k], hm, sc[i], hx[i]);
        }
    }
    
    for j in 0..n {
        for k in 0..m {
            bli[j] += w[k] * y[k] * xli[j][k];
            for i in 0..n {
                ali[j][i] += w[k] * xli[i][k] * xli[j][k];
            }
        }
    }
    
    *sing = crate::irifun_utils::lnglsn(n, &mut ali, &mut bli);
    if !*sing {
        for i in 0..n {
            var[i] = ali[n - 1][i];
        }
    }
}

pub fn inilay(
    night: bool,
    f1reg: bool,
    xnmf2: f32,
    xnmf1: f32,
    xnme: f32,
    vne: f32,
    hmf2: f32,
    hmf1: f32,
    hme: f32,
    hv1: f32,
    hv2: f32,
    hhalf: f32,
    hxl: &mut [f32; 4],
    scl: &mut [f32; 4],
    amp: &mut [f32; 4],
    iqual: &mut i32,
) {
    let numlay = 4;
    let nc1 = 2;
    let alg102 = 2.0_f32.log10();
    
    let alogf = xnmf2.log10();
    let alogef = xnme.log10() - alogf;
    let xhalf = xnmf2 / 2.0;
    
    let mut xx = [0.0_f32; 8];
    let mut yy = [0.0_f32; 8];
    let mut ww = [0.0_f32; 8];
    
    xx[0] = hhalf;
    xx[1] = hv1;
    xx[2] = hv2;
    xx[3] = hme;
    xx[4] = hme - (hv2 - hme);
    
    yy[0] = -alg102;
    yy[1] = alogef;
    yy[2] = vne.log10() - alogf;
    yy[3] = alogef;
    yy[4] = yy[2];
    yy[6] = 0.0;
    
    ww[1] = 1.0;
    ww[2] = 2.0;
    ww[3] = 5.0;
    
    let scl0 = 0.7 * (0.216 * (hmf2 - hhalf) + 56.8);
    scl[0] = 0.8 * scl0;
    scl[1] = 10.0;
    scl[2] = 9.0;
    scl[3] = 6.0;
    hxl[2] = hv2;
    
    let mut hfff = hhalf;
    let mut xfff = xhalf;
    
    let mut numcon = 0;
    let mut hxl1t = 0.0;
    
    if night {
        numcon = 7;
        hxl[0] = hhalf;
        hxl1t = 0.4 * hmf2 + 30.0;
        hxl[1] = (hmf2 + hv1) / 2.0;
        hxl[3] = hme;
        xx[5] = hv2;
        xx[6] = hme;
        yy[5] = 0.0;
        ww[0] = 1.0;
        ww[2] = 3.0;
        ww[4] = 0.5;
        ww[5] = 50.0;
        ww[6] = 500.0;
    } else {
        numcon = 8;
        hxl[0] = 0.9 * hmf2;
        hxl1t = hhalf;
        hxl[1] = hmf1;
        hxl[3] = hme - scl[3];
        xx[5] = hmf1;
        xx[6] = hv2;
        xx[7] = hme;
        yy[7] = 0.0;
        ww[4] = 1.0;
        ww[6] = 50.0;
        ww[7] = 500.0;
        
        if f1reg {
            yy[5] = xnmf1.log10() - alogf;
            ww[5] = 3.0;
            if (xnmf1 - xhalf) * (hmf1 - hhalf) < 0.0 {
                ww[0] = 0.5;
            } else {
                let zet = yy[0] - yy[5];
                ww[0] = crate::irifun_utils::epst(zet, 0.1, 0.15);
            }
            if hhalf > hmf1 {
                hfff = hmf1;
                xfff = xnmf1;
            } else {
                hfff = hhalf;
                xfff = xhalf;
            }
        } else {
            hxl[1] = (hmf2 + hhalf) / 2.0;
            yy[5] = 0.0;
            ww[5] = 0.0;
            ww[0] = 1.0;
        }
    }
    
    if (hv1 - hfff) * (xnme - xfff) < 0.0 {
        ww[1] = 0.5;
    }
    if hv1 <= hv2 + 5.0 {
        ww[1] = 0.5;
    }
    
    let nc0 = numcon - nc1;
    *iqual = 0;
    
    loop {
        let mut ssin = false;
        lsknm(numlay, numcon, nc0, nc1, hmf2, scl, hxl, &ww[0..numcon], &xx[0..numcon], &yy[0..numcon], amp, &mut ssin);
        
        if *iqual > 0 {
            if ssin {
                *iqual = 2;
            }
            break;
        }
        
        if amp[0].abs() > 10.0 || ssin {
            *iqual = 1;
            hxl[0] = hxl1t;
            continue;
        }
        break;
    }
}
