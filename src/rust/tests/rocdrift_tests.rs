use iri2020::ffi::*;
use iri2020::rocdrift::RocdriftModel;
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

#[test]
fn test_rocdrift_equivalence() {
    let _lock = get_ffi_mutex().lock().unwrap();

    // 1. Verify initialization (vfjmodelrocstart)
    let mut vzm_f = [0.0_f32; 64900];
    unsafe {
        vfjmodelrocstart_c(vzm_f.as_mut_ptr());
    }

    let model = RocdriftModel::new();

    // Flatten model.vzm in column-major order to compare against vzm_f
    // Fortran: vzm(59,25,4,11) is column-major:
    // vzm_f[it + 59 * (il + 25 * (isn + 4 * is))]
    for is in 0..11 {
        for isn in 0..4 {
            for il in 0..25 {
                for it in 0..59 {
                    let idx = it + 59 * (il + 25 * (isn + 4 * is));
                    let val_f = vzm_f[idx];
                    let val_r = model.vzm[it][il][isn][is];
                    assert_close(
                        val_r,
                        val_f,
                        1e-5,
                        "VZM initialization",
                        &format!("it={}, il={}, isn={}, is={}", it, il, isn, is),
                    );
                }
            }
        }
    }

    // 2. Grid test over multiple configurations
    let f107_vals = [80.0, 150.0, 220.0];
    let doy_vals = [1, 80, 180, 280, 355];
    let ttl_vals = [0.0, 6.0, 12.0, 18.0, 23.5];
    let gglon_vals = [-180.0, -90.0, 0.0, 90.0, 180.0];

    for &f107 in &f107_vals {
        for &doy in &doy_vals {
            // Solve season and solar indices
            let mut jseas_f: i32 = 0;
            let mut jsfl_f: i32 = 0;
            unsafe {
                vfjmodelrocinit_c(f107, doy, &mut jseas_f, &mut jsfl_f);
            }

            let (jseas_r, jsfl_r) = model.vfjmodelrocinit(f107, doy);

            assert_eq!(jseas_r + 1, jseas_f as usize, "jseas mismatch for f107={}, doy={}", f107, doy);
            assert_eq!(jsfl_r + 1, jsfl_f as usize, "jsfl mismatch for f107={}, doy={}", f107, doy);

            for &ttl in &ttl_vals {
                for &gglon in &gglon_vals {
                    let mut viv_f: f32 = 0.0;
                    unsafe {
                        vfjmodelroc_c(
                            vzm_f.as_ptr(),
                            ttl,
                            gglon,
                            jseas_f,
                            jsfl_f,
                            &mut viv_f,
                        );
                    }

                    let viv_r = model.vfjmodelroc(ttl, gglon, jseas_r, jsfl_r);

                    let inputs_str = format!(
                        "f107={}, doy={}, ttl={}, gglon={}",
                        f107, doy, ttl, gglon
                    );

                    assert_close(
                        viv_r,
                        viv_f,
                        1e-5,
                        "Vertical Drift (viv)",
                        &inputs_str,
                    );
                }
            }
        }
    }
}
