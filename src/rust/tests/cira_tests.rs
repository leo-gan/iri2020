use iri2020::ffi::*;
use iri2020::cira::{cira_gtd7, cira_gtd7d, CiraModel};
use std::sync::Mutex;
use std::sync::OnceLock;

static FFI_MUTEX: OnceLock<Mutex<()>> = OnceLock::new();

fn get_ffi_mutex() -> &'static Mutex<()> {
    FFI_MUTEX.get_or_init(|| Mutex::new(()))
}

fn assert_close(computed: f32, expected: f32, tol: f32, label: &str, inputs: &str) {
    if computed.is_nan() && expected.is_nan() {
        return;
    }
    let abs_diff = (computed - expected).abs();
    if abs_diff <= 1e-7 {
        return;
    }
    if expected.abs() > 1e-5 {
        let rel_diff = abs_diff / expected.abs();
        if rel_diff <= tol {
            return;
        }
        panic!(
            "[{}] Values differ: computed {}, expected {} (abs_diff = {}, rel_diff = {}, tol = {}) for {}",
            label, computed, expected, abs_diff, rel_diff, tol, inputs
        );
    } else {
        if abs_diff <= tol {
            return;
        }
        panic!(
            "[{}] Values differ: computed {}, expected {} (abs_diff = {}, tol = {}) for {}",
            label, computed, expected, abs_diff, tol, inputs
        );
    }
}

fn run_side_by_side(
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
) {
    let _lock = get_ffi_mutex().lock().unwrap();

    let mut d_f = [0.0_f32; 9];
    let mut t_f = [0.0_f32; 2];
    let mut d_r = [0.0_f32; 9];
    let mut t_r = [0.0_f32; 2];

    // Reset switches to default (1.0) and meters to false (cgs)
    let sv_default = [1.0_f32; 25];
    unsafe {
        tselec_c(sv_default.as_ptr());
        meters_c(false);
    }

    // Call Fortran GTD7 via FFI
    unsafe {
        gtd7_c(
            iyd,
            sec,
            alt,
            glat,
            glong,
            stl,
            f107a,
            f107,
            ap.as_ptr(),
            mass,
            d_f.as_mut_ptr(),
            t_f.as_mut_ptr(),
        );
    }

    // Call Rust GTD7
    let mut model = CiraModel::new();
    model.gtd7(
        iyd,
        sec,
        alt,
        glat,
        glong,
        stl,
        f107a,
        f107,
        ap,
        mass,
        &mut d_r,
        &mut t_r,
    );



    // Verify GTD7 equivalence under 1e-5 relative tolerance
    let inputs_str = format!(
        "alt={}, lat={}, lon={}, stl={}, f107a={}, f107={}, mass={}",
        alt, glat, glong, stl, f107a, f107, mass
    );
    for i in 0..9 {
        assert_close(d_r[i], d_f[i], 1e-5, &format!("GTD7 D[{}]", i), &inputs_str);
    }
    for i in 0..2 {
        assert_close(t_r[i], t_f[i], 1e-5, &format!("GTD7 T[{}]", i), &inputs_str);
    }

    // Call Fortran GTD7D via FFI
    let mut d_fd = [0.0_f32; 9];
    let mut t_fd = [0.0_f32; 2];
    let mut d_rd = [0.0_f32; 9];
    let mut t_rd = [0.0_f32; 2];

    unsafe {
        gtd7d_c(
            iyd,
            sec,
            alt,
            glat,
            glong,
            stl,
            f107a,
            f107,
            ap.as_ptr(),
            mass,
            d_fd.as_mut_ptr(),
            t_fd.as_mut_ptr(),
        );
    }

    // Call Rust GTD7D
    let mut model_d = CiraModel::new();
    model_d.gtd7d(
        iyd,
        sec,
        alt,
        glat,
        glong,
        stl,
        f107a,
        f107,
        ap,
        mass,
        &mut d_rd,
        &mut t_rd,
    );

    // Verify GTD7D equivalence under 1e-5 relative tolerance
    for i in 0..9 {
        assert_close(d_rd[i], d_fd[i], 1e-5, &format!("GTD7D D[{}]", i), &inputs_str);
    }
    for i in 0..2 {
        assert_close(t_rd[i], t_fd[i], 1e-5, &format!("GTD7D T[{}]", i), &inputs_str);
    }
}

#[test]
fn test_cira_equivalence() {
    let altitudes: [f32; 11] = [0.0, 50.0, 75.0, 80.0, 100.0, 150.0, 200.0, 300.0, 500.0, 800.0, 1000.0];
    let latitudes: [f32; 5] = [-90.0, -45.0, 0.0, 45.0, 90.0];
    let longitudes: [f32; 5] = [-180.0, -90.0, 0.0, 90.0, 180.0];
    let f107_vals: [f32; 3] = [150.0, 80.0, 220.0];
    let ap_vals = [
        [4.0_f32; 7],
        [10.0, 20.0, 5.0, 8.0, 9.0, 12.0, 15.0],
    ];
    let masses = [48, 0, 4, 16, 28, 32, 40, 1, 14, 17];

    // Run a variety of coordinate, date, and activity configurations
    for &alt in &altitudes {
        for &lat in &latitudes {
            for &lon in &longitudes {
                let stl = (12.0 + lon / 15.0).rem_euclid(24.0_f32);
                for &f107 in &f107_vals {
                    for ap in &ap_vals {
                        for &mass in &masses {
                            run_side_by_side(
                                2003120, // YYDDD
                                43200.0, // UT sec
                                alt,
                                lat,
                                lon,
                                stl,
                                f107,     // F10.7A
                                f107,     // F10.7 daily
                                ap,
                                mass,
                            );
                        }
                    }
                }
            }
        }
    }
}
