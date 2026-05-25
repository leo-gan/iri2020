use iri2020::e_layer::*;
use iri2020::spharm::*;
use iri2020::d_region::*;
use iri2020::b0_b1_model::*;
use iri2020::ffi::*;

#[test]
fn test_e_layer_equivalence() {
    unsafe {
        init_igrf_c();
    }

    // Grid of inputs for foeedi
    let covs = [70.0, 100.0, 150.0];
    let xhis = [10.0, 45.0, 75.0, 89.0, 95.0, 120.0];
    let xhims = [20.0, 50.0, 80.0];
    let lats = [-60.0_f32, -30.0_f32, 0.0_f32, 30.0_f32, 60.0_f32];

    for &cov in &covs {
        for &xhi in &xhis {
            for &xhim in &xhims {
                for &lat in &lats {
                    let rust_val = foeedi(cov, xhi, xhim, lat.abs());
                    let fortran_val = unsafe { foeedi_c(cov, xhi, xhim, lat.abs()) };
                    let diff = (rust_val - fortran_val).abs();
                    assert!(
                        diff < 1e-4,
                        "foeedi mismatch for cov={}, xhi={}, xhim={}, lat={}: rust = {}, fortran = {}",
                        cov, xhi, xhim, lat, rust_val, fortran_val
                    );
                }
            }
        }
    }

    // Grid of inputs for xmded
    let rs = [10.0, 100.0, 200.0];
    let yws = [1.0e7, 4.0e8];
    for &r in &rs {
        for &xhi in &xhis {
            for &yw in &yws {
                let rust_val = xmded(xhi, r, yw);
                let fortran_val = unsafe { xmded_c(xhi, r, yw) };
                let diff = (rust_val - fortran_val).abs() / fortran_val;
                assert!(
                    diff < 1e-4,
                    "xmded mismatch for xhi={}, r={}, yw={}: rust = {}, fortran = {}",
                    xhi, r, yw, rust_val, fortran_val
                );
            }
        }
    }

    // Grid of inputs for valgul
    for &xhi in &xhis {
        let (rust_hvb, rust_vwu, rust_vwa, rust_vdp) = valgul(xhi);
        let (mut fort_hvb, mut fort_vwu, mut fort_vwa, mut fort_vdp) = (0.0, 0.0, 0.0, 0.0);
        unsafe {
            valgul_c(xhi, &mut fort_hvb, &mut fort_vwu, &mut fort_vwa, &mut fort_vdp);
        }
        assert!((rust_hvb - fort_hvb).abs() < 1e-3, "valgul hvb mismatch for xhi={}: rust={}, fort={}", xhi, rust_hvb, fort_hvb);
        assert!((rust_vwu - fort_vwu).abs() < 1e-3, "valgul vwu mismatch for xhi={}: rust={}, fort={}", xhi, rust_vwu, fort_vwu);
        assert!((rust_vwa - fort_vwa).abs() < 1e-3, "valgul vwa mismatch for xhi={}: rust={}, fort={}", xhi, rust_vwa, fort_vwa);
        assert!((rust_vdp - fort_vdp).abs() < 1e-4, "valgul vdp mismatch for xhi={}: rust={}, fort={}", xhi, rust_vdp, fort_vdp);
    }
}

#[test]
fn test_spharm_equivalence() {
    let colats = [0.1, 0.5, 1.0, 1.5, 2.0];
    let azs = [0.0, 0.5, 1.5, 3.0, 5.0];
    
    for &colat in &colats {
        for &az in &azs {
            // Test SPHARM
            let mut rust_c = [0.0_f32; 82];
            let mut fort_c = [0.0_f32; 82];
            spharm(&mut rust_c, 8, 8, colat, az);
            unsafe {
                spharm_c(fort_c.as_mut_ptr(), 8, 8, colat, az);
            }
            for i in 0..82 {
                let diff = (rust_c[i] - fort_c[i]).abs();
                assert!(diff < 1e-5, "spharm mismatch at idx {}: rust = {}, fortran = {}, colat = {}, az = {}", i, rust_c[i], fort_c[i], colat, az);
            }
            
            // Test SPHARM_IK
            let mut rust_c_ik = [0.0_f32; 82];
            let mut fort_c_ik = [0.0_f32; 82];
            spharm_ik(&mut rust_c_ik, 8, 8, colat, az);
            unsafe {
                spharm_ik_c(fort_c_ik.as_mut_ptr(), 8, 8, colat, az);
            }
            for i in 0..82 {
                let diff = (rust_c_ik[i] - fort_c_ik[i]).abs();
                assert!(diff < 1e-5, "spharm_ik mismatch at idx {}: rust = {}, fortran = {}, colat = {}, az = {}", i, rust_c_ik[i], fort_c_ik[i], colat, az);
            }
        }
    }
}

#[test]
fn test_d_region_equivalence() {
    // Test dregion
    let zs = [10.0, 45.0, 75.0, 85.0, 95.0, 120.0];
    let its = [1, 3, 6, 10, 12];
    let fs = [70.0, 150.0, 220.0];
    let vkps = [0.0, 1.5, 3.0];
    let f5sws = [0.0, 0.5, 1.0];
    let f6was = [0.0, 0.5, 1.0];
    
    for &z in &zs {
        for &it in &its {
            for &f in &fs {
                for &vkp in &vkps {
                    for &f5sw in &f5sws {
                        for &f6wa in &f6was {
                            let mut rust_elg = [0.0_f32; 7];
                            let mut fort_elg = [0.0_f32; 7];
                            dregion(z, it, f, vkp, f5sw, f6wa, &mut rust_elg);
                            unsafe {
                                dregion_c(z, it, f, vkp, f5sw, f6wa, fort_elg.as_mut_ptr());
                            }
                            for i in 0..7 {
                                let diff = (rust_elg[i] - fort_elg[i]).abs();
                                assert!(diff < 1e-5, "dregion mismatch at idx {}: rust = {}, fortran = {}, z = {}, it = {}, f = {}, vkp = {}", i, rust_elg[i], fort_elg[i], z, it, f, vkp);
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Test f00 (FIRI model)
    let hgts = [60.5, 75.0, 90.3, 110.0, 139.5];
    let glats = [-60.0, -30.0, 0.0, 30.0, 60.0];
    let idays = [10, 80, 172, 266, 355];
    let zangs = [10.0, 45.0, 75.0, 85.0, 95.0, 120.0];
    let f107ts = [70.0, 130.0, 220.0];
    
    for &hgt in &hgts {
        for &glat in &glats {
            for &iday in &idays {
                for &zang in &zangs {
                    for &f107t in &f107ts {
                        let rust_res = f00(hgt, glat, iday, zang, f107t);
                        
                        let mut fort_edens = 0.0_f32;
                        let mut fort_ierror = 0_i32;
                        unsafe {
                            f00_c(hgt, glat, iday, zang, f107t, &mut fort_edens, &mut fort_ierror);
                        }
                        
                        match rust_res {
                            Ok((rust_edens, rust_ierror)) => {
                                assert_eq!(rust_ierror, fort_ierror, "f00 ierror mismatch: rust = {}, fortran = {}", rust_ierror, fort_ierror);
                                let rel_diff = if fort_edens > 0.0 {
                                    (rust_edens - fort_edens).abs() / fort_edens
                                } else {
                                    (rust_edens - fort_edens).abs()
                                };
                                assert!(rel_diff < 1e-4, "f00 edens mismatch for hgt={}, glat={}, iday={}, zang={}, f107t={}: rust = {}, fortran = {}", hgt, glat, iday, zang, f107t, rust_edens, fort_edens);
                            }
                            Err(rust_ierror) => {
                                assert_eq!(rust_ierror, fort_ierror, "f00 ierror mismatch (err path): rust = {}, fortran = {}", rust_ierror, fort_ierror);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn test_b0_b1_equivalence() {
    unsafe {
        init_igrf_c();
    }
    let lats = [-80.0, -40.0, 0.0, 40.0, 80.0];
    let lons = [-160.0, -80.0, 0.0, 80.0, 160.0];
    let ts = [1.0, 4.0, 7.0, 10.0];
    let rzs = [10.0, 80.0, 150.0];
    
    for &lat in &lats {
        for &lon in &lons {
            for &t in &ts {
                for &rz in &rzs {
                    // Test shamdb0d
                    let rust_b0 = shamdb0d(lat, lon, t, rz);
                    let fort_b0 = unsafe { shamdb0d_c(lat, lon, t, rz) };
                    let diff_b0 = (rust_b0 - fort_b0).abs();
                    assert!(diff_b0 < 1e-4, "shamdb0d mismatch for lat={}, lon={}, t={}, rz={}: rust = {}, fortran = {}", lat, lon, t, rz, rust_b0, fort_b0);
                    
                    // Test shab1d
                    let rust_b1 = shab1d(lat, lon, t, rz);
                    let fort_b1 = unsafe { shab1d_c(lat, lon, t, rz) };
                    let diff_b1 = (rust_b1 - fort_b1).abs();
                    assert!(diff_b1 < 1e-4, "shab1d mismatch for lat={}, lon={}, t={}, rz={}: rust = {}, fortran = {}", lat, lon, t, rz, rust_b1, fort_b1);
                }
            }
        }
    }
}
