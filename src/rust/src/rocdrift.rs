//! ROCSAT-1 Vertical Plasma Drift Model
//! Quiet time equatorial F-region vertical plasma drift model derived from ROCSAT-1 observations.

use crate::rocdrift_coeff::FJROCVZ;

pub struct RocdriftModel {
    pub vzm: [[[[f32; 11]; 4]; 25]; 59],
}

impl RocdriftModel {
    /// Creates a new `RocdriftModel` instance and resolves the drift arrays (fills in missing values).
    /// Mimics the Fortran `vfjmodelrocstart` subroutine.
    pub fn new() -> Self {
        let mut vzm = [[[[0.0; 11]; 4]; 25]; 59];
        let mlt = 58;

        for is in 0..11 {
            for isn in 0..4 {
                for il in 0..25 {
                    for it in 0..59 {
                        let mut vzmx = FJROCVZ[it][il][isn][is];
                        if vzmx < -900.0 {
                            let mut itm = it as i32 - 1;
                            if itm < 0 {
                                itm += mlt;
                            }
                            vzmx = FJROCVZ[itm as usize][il][isn][is];
                        }
                        vzm[it][il][isn][is] = vzmx;
                    }
                }
            }
        }

        RocdriftModel { vzm }
    }

    /// Solves the season and solar activity indices based on day of year and F10.7 flux.
    /// Mimics `vfjmodelrocinit` and returns 0-based indices.
    pub fn vfjmodelrocinit(&self, f107: f32, idoy: i32) -> (usize, usize) {
        let seas = [59.0, 120.0, 243.0, 304.0];
        let sfl = [100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0];

        let doy = idoy as f32;
        let mut jseas = 1;
        if doy > seas[0] && doy <= seas[1] {
            jseas = 2;
        } else if doy > seas[1] && doy <= seas[2] {
            jseas = 3;
        } else if doy > seas[2] && doy <= seas[3] {
            jseas = 4;
        }

        let msol = 11;
        let mut jsfl = 1;
        if f107 < sfl[1] {
            jsfl = 1;
        } else if f107 >= sfl[msol - 1] {
            jsfl = msol;
        } else {
            for i in 0..(msol - 1) {
                if sfl[i] <= f107 && f107 < sfl[i + 1] {
                    jsfl = i + 1;
                    break;
                }
            }
        }

        (jseas - 1, jsfl - 1)
    }

    /// Evaluates the vertical plasma drift for a given local time and geographic longitude.
    /// Mimics `vfjmodelroc`.
    pub fn vfjmodelroc(&self, ttl: f32, gglon: f32, jseas: usize, jsfl: usize) -> f32 {
        let glon = [
            -180.0, -165.0, -150.0, -135.0, -120.0, -105.0, -90.0, -75.0, -60.0,
            -45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0, 60.0, 75.0,
            90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0,
        ];
        let tl = [
            0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00,
            4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50,
            9.00, 9.50, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00,
            13.50, 14.00, 14.50, 15.00, 15.50, 16.00, 16.50, 17.00, 17.25,
            17.50, 17.75, 18.00, 18.25, 18.50, 18.75, 19.00, 19.25, 19.50,
            19.75, 20.00, 20.25, 20.50, 20.75, 21.00, 21.25, 21.50, 21.75,
            22.00, 22.50, 23.00, 23.50, 24.00,
        ];

        let mut xgglon = gglon;
        if xgglon > 180.0 {
            xgglon -= 360.0;
        }

        let mut xtl = ttl;
        if xtl < 0.0 {
            xtl += 24.0;
        }
        if xtl > 24.0 {
            xtl -= 24.0;
        }

        // Copy subset slice for interpolation
        let mut aa = [[0.0; 25]; 59];
        for j in 0..59 {
            for jj in 0..25 {
                aa[j][jj] = self.vzm[j][jj][jseas][jsfl];
            }
        }

        fjlin22dex(xtl, xgglon, &tl, &glon, &aa)
    }
}

fn fjlocate(x: f32, xa: &[f32]) -> usize {
    let na = xa.len();
    if na <= 1 {
        return 0;
    }

    let mut klo = 0;
    if xa[1] > xa[0] {
        // Ascending
        if x < xa[0] {
            return 0;
        }
        if x >= xa[na - 1] {
            return na - 1;
        }
        if xa[klo] <= x && x < xa[klo + 1] {
            return klo;
        }
        let mut khi = na - 1;
        while khi - klo > 1 {
            let k = (khi + klo) / 2;
            if xa[k] > x {
                khi = k;
            } else {
                klo = k;
            }
        }
        klo
    } else {
        // Descending
        if x > xa[0] {
            return 0;
        }
        if x < xa[na - 1] {
            return na - 1;
        }
        if xa[klo] >= x && x > xa[klo + 1] {
            return klo;
        }
        let mut khi = na - 1;
        while khi - klo > 1 {
            let k = (khi + klo) / 2;
            if xa[k] >= x {
                klo = k;
            } else {
                khi = k;
            }
        }
        klo
    }
}

fn fjlin22dex(
    x: f32,
    y: f32,
    xd: &[f32],
    yd: &[f32],
    d: &[[f32; 25]; 59],
) -> f32 {
    let nx = xd.len();
    let ny = yd.len();

    let kx = fjlocate(x, xd).min(nx - 2);
    let hx = xd[kx + 1] - xd[kx];
    let ax = (xd[kx + 1] - x) / hx;
    let bx = (x - xd[kx]) / hx;

    let ky = fjlocate(y, yd).min(ny - 2);
    let hy = yd[ky + 1] - yd[ky];
    let ay = (yd[ky + 1] - y) / hy;
    let by = (y - yd[ky]) / hy;

    ax * ay * d[kx][ky]
        + bx * ay * d[kx + 1][ky]
        + ax * by * d[kx][ky + 1]
        + bx * by * d[kx + 1][ky + 1]
}
