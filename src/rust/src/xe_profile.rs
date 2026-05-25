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
