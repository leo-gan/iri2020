use std::path::Path;

pub const PI: f32 = 3.141592653589793;
pub const UMR: f32 = PI / 180.0;
pub const ERA: f32 = 6371.2;
pub const EREQU: f32 = 6378.16;
pub const ERPOL: f32 = 6356.775;
pub const AQUAD: f32 = EREQU * EREQU;
pub const BQUAD: f32 = ERPOL * ERPOL;

const FILMOD: [&str; 17] = [
    "dgrf1945.dat", "dgrf1950.dat", "dgrf1955.dat",
    "dgrf1960.dat", "dgrf1965.dat", "dgrf1970.dat", "dgrf1975.dat",
    "dgrf1980.dat", "dgrf1985.dat", "dgrf1990.dat", "dgrf1995.dat",
    "dgrf2000.dat", "dgrf2005.dat", "dgrf2010.dat", "dgrf2015.dat",
    "igrf2020.dat", "igrf2020s.dat"
];

const DTEMOD: [f32; 17] = [
    1945.0, 1950.0, 1955.0, 1960.0, 1965.0,
    1970.0, 1975.0, 1980.0, 1985.0, 1990.0, 1995.0, 2000.0, 2005.0,
    2010.0, 2015.0, 2020.0, 2025.0
];

const NUMYE: usize = 16;

#[derive(Debug, Clone)]
pub struct FeldgResult {
    pub bnorth: f32,
    pub beast: f32,
    pub bdown: f32,
    pub babs: f32,
}

#[derive(Debug, Clone)]
pub struct IgrfModel {
    pub nmax: i32,
    pub time: f32,
    pub g: [f32; 197],
    pub era: f32,
    pub dimo: f32,
    pub ghi1: f32,
    pub ghi2: f32,
    pub ghi3: f32,
}

impl IgrfModel {
    pub fn new() -> Self {
        IgrfModel {
            nmax: 0,
            time: 0.0,
            g: [0.0; 197],
            era: ERA,
            dimo: 0.311653,
            ghi1: 0.0,
            ghi2: 0.0,
            ghi3: 0.0,
        }
    }

    pub fn feldcof(&mut self, year: f32, data_dir: &str) -> Result<(), String> {
        self.time = year;
        let iyea = ((year / 5.0) as i32) * 5;
        let mut l = ((iyea - 1945) / 5) + 1;
        if l < 1 {
            l = 1;
        }
        if l > (NUMYE as i32) {
            l = NUMYE as i32;
        }

        let idx1 = (l - 1) as usize;
        let dte1 = DTEMOD[idx1];
        let fil1 = FILMOD[idx1];
        let dte2 = DTEMOD[idx1 + 1];
        let fil2 = FILMOD[idx1 + 1];

        let mut gh1 = [0.0_f32; 197];
        let mut gh2 = [0.0_f32; 197];
        let mut gha = [0.0_f32; 197];

        let (nmax1, erad1, _) = getshc(data_dir, fil1, &mut gh1)?;
        let (nmax2, _, _) = getshc(data_dir, fil2, &mut gh2)?;

        self.era = erad1;

        if l <= (NUMYE as i32) - 1 {
            intershc(year, dte1, nmax1, &gh1, dte2, nmax2, &gh2, &mut self.nmax, &mut gha);
        } else {
            extrashc(year, dte1, nmax1, &gh1, nmax2, &gh2, &mut self.nmax, &mut gha);
        }

        let mut f0 = 0.0_f64;
        for j in 1..=3 {
            let f = (gha[j] as f64) * 1e-5;
            f0 += f * f;
        }
        self.dimo = f0.sqrt() as f32;
        self.ghi1 = gha[1];
        self.ghi2 = gha[2];
        self.ghi3 = gha[3];

        self.g[1] = 0.0;
        let mut i = 2;
        let mut f0_norm = -1e-5_f64; // IS = 0, so negative
        let sqrt2 = 2.0_f64.sqrt();

        for n in 1..=(self.nmax as usize) {
            let x = n as f64;
            f0_norm = f0_norm * x * x / (4.0 * x - 2.0);
            f0_norm = f0_norm * (2.0 * x - 1.0) / x;
            let mut f = f0_norm * 0.5 * sqrt2;
            self.g[i] = gha[i - 1] * (f0_norm as f32);
            i += 1;

            for m in 1..=n {
                f = f * (x + m as f64) / (x - m as f64 + 1.0);
                f = f * ((x - m as f64 + 1.0) / (x + m as f64)).sqrt();
                self.g[i] = gha[i - 1] * (f as f32);
                self.g[i + 1] = gha[i] * (f as f32);
                i += 2;
            }
        }

        Ok(())
    }

    pub fn feldg(&self, glat: f32, glon: f32, alt: f32) -> FeldgResult {
        let rlat = glat * UMR;
        let ct = rlat.sin();
        let st = rlat.cos();
        let d = (AQUAD - (AQUAD - BQUAD) * ct * ct).sqrt();
        let rlon = glon * UMR;
        let cp = rlon.cos();
        let sp = rlon.sin();
        let zzz = (alt + BQUAD / d) * ct / ERA;
        let rho = (alt + AQUAD / d) * st / ERA;
        let xxx = rho * cp;
        let yyy = rho * sp;

        let rq = 1.0 / (xxx * xxx + yyy * yyy + zzz * zzz);
        let mut xi = [0.0_f32; 4];
        xi[1] = xxx * rq;
        xi[2] = yyy * rq;
        xi[3] = zzz * rq;

        let mut h = [0.0_f32; 197];
        let ihmax = (self.nmax * self.nmax + 1) as isize;
        let last = ihmax + (self.nmax + self.nmax) as isize;
        let imax = (self.nmax + self.nmax - 1) as isize;

        for idx in ihmax..=last {
            h[idx as usize] = self.g[idx as usize];
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
                    h[il_u] = self.g[il_u] + z * h[ih_u] + 2.0 * (x * h[ih_u + 1] + y * h[ih_u + 2]);
                } else if i == 1 {
                    h[il_u + 2] = self.g[il_u + 2] + z * h[ih_u + 2] + x * h[ih_u + 4] - y * (h[ih_u + 3] + h[ih_u]);
                    h[il_u + 1] = self.g[il_u + 1] + z * h[ih_u + 1] + y * h[ih_u + 4] + x * (h[ih_u + 3] - h[ih_u]);
                    h[il_u] = self.g[il_u] + z * h[ih_u] + 2.0 * (x * h[ih_u + 1] + y * h[ih_u + 2]);
                } else {
                    let mut m = 3_isize;
                    while m <= i {
                        let m_u = m as usize;
                        h[il_u + m_u + 1] = self.g[il_u + m_u + 1] + z * h[ih_u + m_u + 1] + x * (h[ih_u + m_u + 3] - h[ih_u + m_u - 1]) - y * (h[ih_u + m_u + 2] + h[ih_u + m_u - 2]);
                        h[il_u + m_u] = self.g[il_u + m_u] + z * h[ih_u + m_u] + x * (h[ih_u + m_u + 2] - h[ih_u + m_u - 2]) + y * (h[ih_u + m_u + 3] + h[ih_u + m_u - 1]);
                        m += 2;
                    }
                    h[il_u + 2] = self.g[il_u + 2] + z * h[ih_u + 2] + x * h[ih_u + 4] - y * (h[ih_u + 3] + h[ih_u]);
                    h[il_u + 1] = self.g[il_u + 1] + z * h[ih_u + 1] + y * h[ih_u + 4] + x * (h[ih_u + 3] - h[ih_u]);
                    h[il_u] = self.g[il_u] + z * h[ih_u] + 2.0 * (x * h[ih_u + 1] + y * h[ih_u + 2]);
                }

                ih = il;
                if i < k {
                    break;
                }
            }
        }

        let s = 0.5 * h[1] + 2.0 * (h[2] * xi[3] + h[3] * xi[1] + h[4] * xi[2]);
        let t = (rq + rq) * rq.sqrt();
        let bxxx = t * (h[3] - s * xxx);
        let byyy = t * (h[4] - s * yyy);
        let bzzz = t * (h[2] - s * zzz);

        let babs = (bxxx * bxxx + byyy * byyy + bzzz * bzzz).sqrt();
        let beast = byyy * cp - bxxx * sp;
        let brho = byyy * sp + bxxx * cp;
        let bnorth = bzzz * st - brho * ct;
        let bdown = -bzzz * ct - brho * st;

        FeldgResult {
            bnorth,
            beast,
            bdown,
            babs,
        }
    }
}

fn getshc(data_dir: &str, fspec: &str, gh: &mut [f32; 197]) -> Result<(i32, f32, f32), String> {
    let path = Path::new(data_dir).join(fspec);
    let content = std::fs::read_to_string(&path)
        .map_err(|e| format!("Failed to read {}: {}", path.display(), e))?;

    let mut tokens = content.split_whitespace();

    let _name = tokens.next().ok_or_else(|| format!("Empty file {}", fspec))?;

    let nmax_str = tokens.next().ok_or_else(|| format!("Missing NMAX in {}", fspec))?;
    let erad_str = tokens.next().ok_or_else(|| format!("Missing ERAD in {}", fspec))?;
    let xmyear_str = tokens.next().ok_or_else(|| format!("Missing XMYEAR in {}", fspec))?;

    let nmax: i32 = nmax_str.parse().map_err(|e| format!("Invalid NMAX in {}: {}", fspec, e))?;
    let erad: f32 = erad_str.parse().map_err(|e| format!("Invalid ERAD in {}: {}", fspec, e))?;
    let xmyear: f32 = xmyear_str.parse().map_err(|e| format!("Invalid XMYEAR in {}: {}", fspec, e))?;

    let nm = nmax * (nmax + 2);
    gh.fill(0.0);

    for i in 1..=(nm as usize) {
        let val_str = tokens.next().ok_or_else(|| format!("Expected {} coefficients, found only {} in {}", nm, i - 1, fspec))?;
        let val: f32 = val_str.parse().map_err(|e| format!("Invalid coefficient at index {} in {}: {}", i, fspec, e))?;
        gh[i] = val;
    }

    Ok((nmax, erad, xmyear))
}

fn intershc(date: f32, dte1: f32, nmax1: i32, gh1: &[f32; 197], dte2: f32, nmax2: i32, gh2: &[f32; 197], nmax: &mut i32, gh: &mut [f32; 197]) {
    let factor = (date - dte1) / (dte2 - dte1);
    let k: usize;
    if nmax1 == nmax2 {
        k = (nmax1 * (nmax1 + 2)) as usize;
        *nmax = nmax1;
    } else if nmax1 > nmax2 {
        k = (nmax2 * (nmax2 + 2)) as usize;
        let l = (nmax1 * (nmax1 + 2)) as usize;
        for i in (k + 1)..=l {
            gh[i] = gh1[i] + factor * (-gh1[i]);
        }
        *nmax = nmax1;
    } else {
        k = (nmax1 * (nmax1 + 2)) as usize;
        let l = (nmax2 * (nmax2 + 2)) as usize;
        for i in (k + 1)..=l {
            gh[i] = factor * gh2[i];
        }
        *nmax = nmax2;
    }

    for i in 1..=k {
        gh[i] = gh1[i] + factor * (gh2[i] - gh1[i]);
    }
}

fn extrashc(date: f32, dte1: f32, nmax1: i32, gh1: &[f32; 197], nmax2: i32, gh2: &[f32; 197], nmax: &mut i32, gh: &mut [f32; 197]) {
    let factor = date - dte1;
    let k: usize;
    if nmax1 == nmax2 {
        k = (nmax1 * (nmax1 + 2)) as usize;
        *nmax = nmax1;
    } else if nmax1 > nmax2 {
        k = (nmax2 * (nmax2 + 2)) as usize;
        let l = (nmax1 * (nmax1 + 2)) as usize;
        for i in (k + 1)..=l {
            gh[i] = gh1[i];
        }
        *nmax = nmax1;
    } else {
        k = (nmax1 * (nmax1 + 2)) as usize;
        let l = (nmax2 * (nmax2 + 2)) as usize;
        for i in (k + 1)..=l {
            gh[i] = factor * gh2[i];
        }
        *nmax = nmax2;
    }

    for i in 1..=k {
        gh[i] = gh1[i] + factor * gh2[i];
    }
}

fn get_coeff_array(year: i32) -> (&'static [f32; 66], &'static [f32; 66]) {
    match year {
        1900 => (&crate::igrf_coeff::G1900, &crate::igrf_coeff::H1900),
        1905 => (&crate::igrf_coeff::G1905, &crate::igrf_coeff::H1905),
        1910 => (&crate::igrf_coeff::G1910, &crate::igrf_coeff::H1910),
        1915 => (&crate::igrf_coeff::G1915, &crate::igrf_coeff::H1915),
        1920 => (&crate::igrf_coeff::G1920, &crate::igrf_coeff::H1920),
        1925 => (&crate::igrf_coeff::G1925, &crate::igrf_coeff::H1925),
        1930 => (&crate::igrf_coeff::G1930, &crate::igrf_coeff::H1930),
        1935 => (&crate::igrf_coeff::G1935, &crate::igrf_coeff::H1935),
        1940 => (&crate::igrf_coeff::G1940, &crate::igrf_coeff::H1940),
        1945 => (&crate::igrf_coeff::G1945, &crate::igrf_coeff::H1945),
        1950 => (&crate::igrf_coeff::G1950, &crate::igrf_coeff::H1950),
        1955 => (&crate::igrf_coeff::G1955, &crate::igrf_coeff::H1955),
        1960 => (&crate::igrf_coeff::G1960, &crate::igrf_coeff::H1960),
        1965 => (&crate::igrf_coeff::G1965, &crate::igrf_coeff::H1965),
        1970 => (&crate::igrf_coeff::G1970, &crate::igrf_coeff::H1970),
        1975 => (&crate::igrf_coeff::G1975, &crate::igrf_coeff::H1975),
        1980 => (&crate::igrf_coeff::G1980, &crate::igrf_coeff::H1980),
        1985 => (&crate::igrf_coeff::G1985, &crate::igrf_coeff::H1985),
        1990 => (&crate::igrf_coeff::G1990, &crate::igrf_coeff::H1990),
        1995 => (&crate::igrf_coeff::G1995, &crate::igrf_coeff::H1995),
        2000 => (&crate::igrf_coeff::G2000, &crate::igrf_coeff::H2000),
        2005 => (&crate::igrf_coeff::G2005, &crate::igrf_coeff::H2005),
        2010 => (&crate::igrf_coeff::G2010, &crate::igrf_coeff::H2010),
        _ => (&crate::igrf_coeff::G2015, &crate::igrf_coeff::H2015),
    }
}

fn get_standalone_coeffs(iy: i32) -> ([f32; 67], [f32; 67]) {
    let mut iyr = iy;
    if iyr < 1900 {
        iyr = 1900;
    }
    if iyr > 2020 {
        iyr = 2020;
    }

    let mut g = [0.0_f32; 67];
    let mut h = [0.0_f32; 67];

    if iyr >= 2015 {
        let dt = (iyr - 2015) as f32;
        for n in 1..=66 {
            g[n] = crate::igrf_coeff::G2015[n - 1];
            h[n] = crate::igrf_coeff::H2015[n - 1];
            if n <= 45 {
                g[n] += crate::igrf_coeff::DG[n - 1] * dt;
                h[n] += crate::igrf_coeff::DH[n - 1] * dt;
            }
        }
    } else {
        let idx = ((iyr - 1900) / 5) as usize;
        let year_base = 1900 + (idx * 5) as i32;
        let f2 = (iyr - year_base) as f32 / 5.0;
        let f1 = 1.0 - f2;

        let (g1, h1) = get_coeff_array(year_base);
        let (g2, h2) = get_coeff_array(year_base + 5);

        for n in 1..=66 {
            g[n] = g1[n - 1] * f1 + g2[n - 1] * f2;
            h[n] = h1[n - 1] * f1 + h2[n - 1] * f2;
        }
    }

    (g, h)
}

#[derive(Debug, Clone)]
pub struct IgrfResult {
    pub br: f32,
    pub bt: f32,
    pub bf: f32,
}

pub fn igrf(iy: i32, nm: i32, r: f32, t: f32, f: f32) -> IgrfResult {
    let mut rec = [0.0_f32; 67];
    for n in 1..=11 {
        let n = n as i32;
        let mut n2 = 2 * n - 1;
        n2 = n2 * (n2 - 2);
        for m in 1..=n {
            let mn = (n * (n - 1) / 2 + m) as usize;
            rec[mn] = ((n - m) * (n + m - 2)) as f32 / n2 as f32;
        }
    }

    let (mut g, mut h) = get_standalone_coeffs(iy);

    // Schmidt normalization
    let mut s_val = 1.0_f32;
    for n in 2..=11 {
        let mn = n * (n - 1) / 2 + 1;
        s_val = s_val * (2 * n - 3) as f32 / (n - 1) as f32;
        g[mn] = g[mn] * s_val;
        h[mn] = h[mn] * s_val;
        let mut p = s_val;
        for m in 2..=n {
            let aa = if m == 2 { 2.0_f32 } else { 1.0_f32 };
            p = p * (aa * (n - m + 1) as f32 / (n + m - 2) as f32).sqrt();
            let mnn = mn + m - 1;
            g[mnn] = g[mnn] * p;
            h[mnn] = h[mnn] * p;
        }
    }

    let pp = 1.0 / r;
    let mut p_pow = pp;
    let k = nm + 1;

    let mut a = vec![0.0_f32; (k + 1) as usize];
    let mut b = vec![0.0_f32; (k + 1) as usize];

    for n in 1..=(k as usize) {
        p_pow = p_pow * pp;
        a[n] = p_pow;
        b[n] = p_pow * n as f32;
    }

    let mut p = 1.0_f32;
    let mut d = 0.0_f32;
    let mut bbr = 0.0_f32;
    let mut bbt = 0.0_f32;
    let mut bbf = 0.0_f32;

    let u = t;
    let cf = f.cos();
    let sf = f.sin();
    let c = u.cos();
    let s = u.sin();

    let mut x = 0.0_f32;
    let mut y = 0.0_f32;

    for m in 1..=(k as usize) {
        if m == 1 {
            x = 0.0;
            y = 1.0;
        } else {
            let w = x;
            x = w * cf + y * sf;
            y = y * cf - w * sf;
        }

        let mut q = p;
        let mut z = d;
        let mut bi = 0.0_f32;
        let mut p2 = 0.0_f32;
        let mut d2 = 0.0_f32;

        for n in m..=(k as usize) {
            let an = a[n];
            let mn = n * (n - 1) / 2 + m;
            let e = g[mn];
            let hh = h[mn];
            let w = e * y + hh * x;
            bbr = bbr + b[n] * w * q;
            bbt = bbt - an * w * z;
            if m > 1 {
                let mut qq = q;
                if s < 1e-5 {
                    qq = z;
                }
                bi = bi + an * (e * x - hh * y) * qq;
            }
            let xk = rec[mn];
            let dp = c * z - s * q - xk * d2;
            let pm = c * q - xk * p2;
            d2 = z;
            p2 = q;
            z = dp;
            q = pm;
        }

        d = s * d + c * p;
        p = s * p;
        if m > 1 {
            let mm = (m - 1) as f32;
            bi = bi * mm;
            bbf = bbf + bi;
        }
    }

    let br = bbr;
    let bt = bbt;
    let bf;
    if s < 1e-5 {
        if c < 0.0 {
            bbf = -bbf;
        }
        bf = bbf;
    } else {
        bf = bbf / s;
    }

    IgrfResult { br, bt, bf }
}

#[derive(Debug, Clone)]
pub struct IgrfDipResult {
    pub dec: f32,
    pub dip: f32,
    pub dipl: f32,
    pub ymodip: f32,
}

pub fn igrf_dip(xlat: f32, xlong: f32, year: f32, height: f32, data_dir: &str) -> Result<IgrfDipResult, String> {
    let mut model = IgrfModel::new();
    model.feldcof(year, data_dir)?;

    let res = model.feldg(xlat, xlong, height);

    let dec_arg = res.beast / (res.beast * res.beast + res.bnorth * res.bnorth).sqrt();
    let mut dec_arg_clamped = dec_arg;
    if dec_arg_clamped.abs() > 1.0 {
        dec_arg_clamped = dec_arg_clamped.signum();
    }
    let dec = dec_arg_clamped.asin();

    let bd_ba = res.bdown / res.babs;
    let mut bd_ba_clamped = bd_ba;
    if bd_ba_clamped.abs() > 1.0 {
        bd_ba_clamped = bd_ba_clamped.signum();
    }
    let dip = bd_ba_clamped.asin();

    let dip_div = dip / (dip * dip + (xlat * UMR).cos()).sqrt();
    let mut dip_div_clamped = dip_div;
    if dip_div_clamped.abs() > 1.0 {
        dip_div_clamped = dip_div_clamped.signum();
    }
    let smodip = dip_div_clamped.asin();

    let dipl = (res.bdown / (2.0 * (res.bnorth * res.bnorth + res.beast * res.beast).sqrt())).atan() / UMR;
    let ymodip = smodip / UMR;

    Ok(IgrfDipResult {
        dec: dec / UMR,
        dip: dip / UMR,
        dipl,
        ymodip,
    })
}
