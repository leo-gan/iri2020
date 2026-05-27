use iri2020::iritec::*;
use iri2020::ffi::*;

fn get_data_dir() -> String {
    if let Ok(dir) = std::env::var("IRI2020_DATA_DIR") {
        return dir;
    }
    let paths = ["src/data", "../data", "../../src/data", "data"];
    for path in &paths {
        if std::path::Path::new(path).exists() {
            if let Ok(abs_path) = std::fs::canonicalize(path) {
                if let Some(path_str) = abs_path.to_str() {
                    std::env::set_var("IRI2020_DATA_DIR", path_str);
                    return path_str.to_string();
                }
            }
        }
    }
    panic!("Could not find data directory");
}

#[test]
fn test_iritec_equivalence() {
    let _data_dir = get_data_dir();
    unsafe {
        init_igrf_c();
        read_ig_rz_c();
        readapf107_c();
    }

    // Default jf switches
    let mut jf = [true; 50];
    // turn off some options as per standard driver
    jf[3] = false;
    jf[4] = false;
    jf[5] = false;
    jf[21] = false;
    jf[22] = false;
    jf[29] = false;
    jf[32] = false;
    jf[33] = false;
    jf[34] = false;
    jf[38] = false;
    jf[39] = false;
    jf[46] = false;

    // Grid of scenarios
    let latitudes = [-45.0, 0.0, 45.0, 75.0];
    let longitudes = [-90.0, 0.0, 90.0];
    let years = [2000, 2015, 2020];
    let mmdds = [315, 615];
    let hours = [12.0, 0.0, 18.5];
    let altitude_ends = [1000.0, 2000.0];
    let steps = [0.1, 0.2];

    for &lat in &latitudes {
        for &lon in &longitudes {
            for &year in &years {
                for &mmdd in &mmdds {
                    for &hour in &hours {
                        for &hend in &altitude_ends {
                            for &hstep in &steps {
                                let mut rust_oarr = [0.0_f32; 100];
                                let mut fort_oarr = [0.0_f32; 100];

                                let mut rust_tecbo = 0.0_f32;
                                let mut rust_tecto = 0.0_f32;

                                let mut fort_tecbo = 0.0_f32;
                                let mut fort_tecto = 0.0_f32;

                                // Call FFI version
                                unsafe {
                                    iritec_c(
                                        lat,
                                        lon,
                                        0, // geographic
                                        jf.as_ptr(),
                                        year,
                                        mmdd,
                                        hour + 25.0,
                                        0.0,
                                        hend,
                                        hstep,
                                        fort_oarr.as_mut_ptr(),
                                        &mut fort_tecbo,
                                        &mut fort_tecto,
                                    );
                                }

                                // Call Rust version
                                iritec(
                                    lat,
                                    lon,
                                    0, // geographic
                                    &jf,
                                    year,
                                    mmdd,
                                    hour + 25.0,
                                    0.0,
                                    hend,
                                    hstep,
                                    &mut rust_oarr,
                                    &mut rust_tecbo,
                                    &mut rust_tecto,
                                );

                                // 0.1 TECU = 1e15 m^-2
                                let tolerance = 1e15_f32;

                                let diff_bo = (rust_tecbo - fort_tecbo).abs();
                                let diff_to = (rust_tecto - fort_tecto).abs();

                                if diff_bo >= tolerance || diff_to >= tolerance {
                                    println!("lat={}, lon={}, year={}, hour={}", lat, lon, year, hour);
                                    for idx in 0..100 {
                                        println!("oarr[{}] -> Rust={}, Fortran={}", idx, rust_oarr[idx], fort_oarr[idx]);
                                    }
                                }
                                assert!(
                                    diff_bo < tolerance,
                                    "Bottomside TEC mismatch for lat={}, lon={}, year={}, hour={}, hend={}, hstep={}: Rust={}, Fortran={}, diff={}",
                                    lat, lon, year, hour, hend, hstep, rust_tecbo, fort_tecbo, diff_bo
                                );

                                assert!(
                                    diff_to < tolerance,
                                    "Topside TEC mismatch for lat={}, lon={}, year={}, hour={}, hend={}, hstep={}: Rust={}, Fortran={}, diff={}",
                                    lat, lon, year, hour, hend, hstep, rust_tecto, fort_tecto, diff_to
                                );

                                // Check that oarr outputs (like hmF2, xnmF2) also match closely
                                // especially because oarr[36] (index 36) contains the computed TEC
                                // but wait, the driver overwrites oarr[36] with tecbo + tecto
                                let diff_hm = (rust_oarr[1] - fort_oarr[1]).abs();
                                let diff_xnm = (rust_oarr[0] - fort_oarr[0]).abs();

                                assert!(
                                    diff_hm < 0.1,
                                    "hmF2 mismatch: Rust={}, Fortran={}",
                                    rust_oarr[1], fort_oarr[1]
                                );
                                assert!(
                                    diff_xnm < 1e-4 * fort_oarr[0] || (diff_xnm < 1.0 && fort_oarr[0] == 0.0),
                                    "xnmF2 mismatch: Rust={}, Fortran={}",
                                    rust_oarr[0], fort_oarr[0]
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn test_ioncorr() {
    let tec = 1e17_f32;
    let freq = 1.5e9_f32; // 1.5 GHz
    let corr = ioncorr(tec, freq);
    let expected = 40.3 * tec / (freq * freq);
    assert!((corr - expected).abs() < 1e-5);
}


