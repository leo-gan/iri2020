use iri2020::ffi::*;
use iri2020::igrf::{IgrfModel, igrf as rust_igrf, igrf_dip as rust_igrf_dip, UMR};
use std::sync::Mutex;
use std::sync::OnceLock;

static FFI_MUTEX: OnceLock<Mutex<()>> = OnceLock::new();

fn get_ffi_mutex() -> &'static Mutex<()> {
    FFI_MUTEX.get_or_init(|| Mutex::new(()))
}

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

fn assert_close(computed: f32, expected: f32, tol: f32) {
    if computed.is_nan() && expected.is_nan() {
        return;
    }
    let abs_diff = (computed - expected).abs();
    if abs_diff <= tol {
        return;
    }
    if expected.abs() > 1e-5 {
        let rel_diff = abs_diff / expected.abs();
        if rel_diff <= tol {
            return;
        }
        panic!(
            "Values differ: computed {}, expected {} (abs_diff = {}, rel_diff = {}, tol = {})",
            computed, expected, abs_diff, rel_diff, tol
        );
    } else {
        panic!(
            "Values differ: computed {}, expected {} (abs_diff = {}, tol = {})",
            computed, expected, abs_diff, tol
        );
    }
}

fn run_side_by_side(glat: f32, glon: f32, alt: f32, year: f32, data_dir: &str) {
    let _lock = get_ffi_mutex().lock().unwrap();
    
    // Set data dir for Fortran FFI (via env var check)
    std::env::set_var("IRI2020_DATA_DIR", data_dir);

    // 1. Call Fortran FELDCOF and FELDG via FFI
    let mut bnorth_f = 0.0_f32;
    let mut beast_f = 0.0_f32;
    let mut bdown_f = 0.0_f32;
    let mut babs_f = 0.0_f32;
    
    unsafe {
        feldcof_c(year);
        feldg_c(glat, glon, alt, &mut bnorth_f, &mut beast_f, &mut bdown_f, &mut babs_f);
    }

    // 2. Call Rust FELDCOF and FELDG
    let mut model = IgrfModel::new();
    model.feldcof(year, data_dir).unwrap();
    let res = model.feldg(glat, glon, alt);

    // Assert FELDCOF / FELDG Equivalence
    assert_close(res.bnorth, bnorth_f, 1e-4);
    assert_close(res.beast, beast_f, 1e-4);
    assert_close(res.bdown, bdown_f, 1e-4);
    assert_close(res.babs, babs_f, 1e-4);

    // 3. Call Fortran standalone IGRF via FFI
    // Subroutine IGRF inputs: IY, NM, R, T, F
    // T (geographic colatitude in radians) = PI/2 - lat_radians
    // F (longitude in radians)
    // R (radius in units of RE=6371.2)
    let lat_rad = glat * UMR;
    let lon_rad = glon * UMR;
    let colat_rad = (std::f32::consts::PI / 2.0) - lat_rad;
    let re = 6371.2_f32;
    let r_re = (alt + re) / re;
    let iy = year as i32;
    let nm = 10; // order <= 10

    let mut br_f = 0.0_f32;
    let mut bt_f = 0.0_f32;
    let mut bf_f = 0.0_f32;

    unsafe {
        igrf_c(iy, nm, r_re, colat_rad, lon_rad, &mut br_f, &mut bt_f, &mut bf_f);
    }

    // Call Rust standalone igrf
    let res_standalone = rust_igrf(iy, nm, r_re, colat_rad, lon_rad);

    // Assert standalone IGRF Equivalence
    assert_close(res_standalone.br, br_f, 1e-4);
    assert_close(res_standalone.bt, bt_f, 1e-4);
    assert_close(res_standalone.bf, bf_f, 1e-4);

    // 4. Call Fortran igrf_dip via FFI
    let mut dec_f = 0.0_f32;
    let mut dip_f = 0.0_f32;
    let mut dipl_f = 0.0_f32;
    let mut ymodip_f = 0.0_f32;

    unsafe {
        igrf_dip_c(glat, glon, year, alt, &mut dec_f, &mut dip_f, &mut dipl_f, &mut ymodip_f);
    }

    // Call Rust igrf_dip
    let res_dip = rust_igrf_dip(glat, glon, year, alt, data_dir).unwrap();

    // Assert igrf_dip Equivalence
    assert_close(res_dip.dec, dec_f, 1e-4);
    assert_close(res_dip.dip, dip_f, 1e-4);
    assert_close(res_dip.dipl, dipl_f, 1e-4);
    assert_close(res_dip.ymodip, ymodip_f, 1e-4);
}

#[test]
fn test_igrf_equivalence() {
    let data_dir = get_data_dir();
    unsafe {
        init_igrf_c();
    }
    
    // Grid of coordinates, altitudes, and years to cover interpolation/extrapolation limits
    let years = [1947.5, 1968.2, 1985.0, 2003.1, 2018.7, 2024.0];
    let locations = [
        (-11.95, -76.77, 300.0), // Scenario 1 location
        (65.1, -147.5, 100.0),   // Scenario 2 location
        (0.0, 0.0, 600.0),       // Equator
        (80.0, 120.0, 80.0),     // Polar region
    ];

    for &year in &years {
        for &loc in &locations {
            run_side_by_side(loc.0, loc.1, loc.2, year, &data_dir);
        }
    }
}
