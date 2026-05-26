use iri2020::iriflip::*;
use iri2020::ffi::*;

#[test]
fn test_geocgm01_equivalence() {
    unsafe {
        init_igrf_c();
    }

    // Grid of inputs for GEOCGM01
    let icors = [1, -1];
    let iyears = [1990, 2010, 2025];
    let heights = [0.0, 300.0, 800.0];
    
    // Grid of latitudes and longitudes
    // Avoid near-equator region (abs < 30 deg where geolow is used and some ranges can be undefined or highly sensitive)
    // Avoid near-pole region (abs >= 89.99 deg)
    let latitudes = [-60.0, -35.0, 45.0, 80.0];
    let longitudes = [-120.0, 0.0, 120.0];

    for &icor in &icors {
        for &iyear in &iyears {
            for &hi in &heights {
                for &lat in &latitudes {
                    for &lon in &longitudes {
                        let mut rust_dat = [0.0_f32; 44];
                        let mut fort_dat = [0.0_f32; 44];
                        
                        let mut rust_pla = [0.0_f32; 4];
                        let mut fort_pla = [0.0_f32; 4];
                        
                        let mut rust_plo = [0.0_f32; 4];
                        let mut fort_plo = [0.0_f32; 4];

                        if icor == 1 {
                            rust_dat[0] = lat;
                            rust_dat[1] = lon;
                            fort_dat[0] = lat;
                            fort_dat[1] = lon;
                        } else {
                            rust_dat[2] = lat;
                            rust_dat[3] = lon;
                            fort_dat[2] = lat;
                            fort_dat[3] = lon;
                        }

                        // Call Rust version
                        geocgm01(icor, iyear, hi, &mut rust_dat, &mut rust_pla, &mut rust_plo);

                        // Call Fortran version via FFI
                        unsafe {
                            geocgm01_c(
                                icor,
                                iyear,
                                hi,
                                fort_dat.as_mut_ptr(),
                                fort_pla.as_mut_ptr(),
                                fort_plo.as_mut_ptr(),
                            );
                        }

                        // Compare start points and conjugate points (location indices 1 and 2, which are the first 2 columns of 11 elements each)
                        for col in 0..2 {
                            // Elements 0: lat, 1: lon
                            let r_lat = rust_dat[0 + col * 11];
                            let f_lat = fort_dat[0 + col * 11];
                            let r_lon = rust_dat[1 + col * 11];
                            let f_lon = fort_dat[1 + col * 11];

                            if f_lat < 900.0 && r_lat < 900.0 {
                                let diff_lat = (r_lat - f_lat).abs();
                                assert!(
                                    diff_lat < 0.01,
                                    "Latitude mismatch for icor={}, year={}, hi={}, lat={}, lon={}, col={}: rust={}, fort={}",
                                    icor, iyear, hi, lat, lon, col, r_lat, f_lat
                                );

                                let mut diff_lon = (r_lon - f_lon).abs();
                                if diff_lon > 180.0 {
                                    diff_lon = 360.0 - diff_lon;
                                }
                                assert!(
                                    diff_lon < 0.01,
                                    "Longitude mismatch for icor={}, year={}, hi={}, lat={}, lon={}, col={}: rust={}, fort={}",
                                    icor, iyear, hi, lat, lon, col, r_lon, f_lon
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
