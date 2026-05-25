use iri2020::ffi::*;
use iri2020::data_io::{IgRzData, Apf107Data, McsatData, CcirUrsiData};
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

#[test]
fn test_ig_rz_equivalence() {
    let _lock = get_ffi_mutex().lock().unwrap();
    let data_dir = get_data_dir();
    std::env::set_var("IRI2020_DATA_DIR", &data_dir);
    
    // 1. Load via Fortran FFI
    unsafe {
        read_ig_rz_c();
    }
    
    let mut aig_f = [0.0_f32; 806];
    let mut arz_f = [0.0_f32; 806];
    let mut iymst_f = 0;
    let mut iymend_f = 0;
    unsafe {
        get_igrz_c(aig_f.as_mut_ptr(), arz_f.as_mut_ptr(), &mut iymst_f, &mut iymend_f);
    }
    
    // 2. Load via Rust
    let rust_data = IgRzData::load(&data_dir).unwrap();
    
    assert_eq!(rust_data.iymst, iymst_f, "iymst mismatch");
    assert_eq!(rust_data.iymend, iymend_f, "iymend mismatch");
    
    for i in 0..806 {
        assert_eq!(
            rust_data.aig[i],
            aig_f[i],
            "aig mismatch at index {}",
            i
        );
        assert_eq!(
            rust_data.arz[i],
            arz_f[i],
            "arz mismatch at index {}",
            i
        );
    }
}

#[test]
fn test_apf107_equivalence() {
    let _lock = get_ffi_mutex().lock().unwrap();
    let data_dir = get_data_dir();
    std::env::set_var("IRI2020_DATA_DIR", &data_dir);
    
    // 1. Load via Fortran FFI
    unsafe {
        readapf107_c();
    }
    
    let mut aap_f_flat = vec![0_i32; 27000 * 9];
    let mut af107_f_flat = vec![0.0_f32; 27000 * 3];
    let mut n_f = 0;
    unsafe {
        get_apfa_c(aap_f_flat.as_mut_ptr(), af107_f_flat.as_mut_ptr(), &mut n_f);
    }
    
    // Transpose column-major Fortran arrays to row-major
    let mut aap_f = vec![[0_i32; 9]; 27000];
    let mut af107_f = vec![[0.0_f32; 3]; 27000];
    for i in 0..27000 {
        for j in 0..9 {
            aap_f[i][j] = aap_f_flat[j * 27000 + i];
        }
        for j in 0..3 {
            af107_f[i][j] = af107_f_flat[j * 27000 + i];
        }
    }
    
    // 2. Load via Rust
    let rust_data = Apf107Data::load(&data_dir).unwrap();
    
    assert_eq!(rust_data.n, n_f, "n mismatch");
    
    for i in 0..27000 {
        assert_eq!(
            rust_data.aap[i],
            aap_f[i],
            "aap mismatch at index {}",
            i
        );
        assert_eq!(
            rust_data.af107[i],
            af107_f[i],
            "af107 mismatch at index {}",
            i
        );
    }
}

#[test]
fn test_mcsat_equivalence() {
    let _lock = get_ffi_mutex().lock().unwrap();
    let data_dir = get_data_dir();
    std::env::set_var("IRI2020_DATA_DIR", &data_dir);
    
    for month in 1..=12 {
        let mut coeff_f_flat = vec![0.0_f64; 149 * 48];
        unsafe {
            read_data_sd_c(month, coeff_f_flat.as_mut_ptr());
        }
        
        // Transpose column-major Fortran arrays to row-major
        let mut coeff_f = [[0.0_f64; 48]; 149];
        for i in 0..149 {
            for j in 0..48 {
                coeff_f[i][j] = coeff_f_flat[j * 149 + i];
            }
        }
        
        let coeff_r = McsatData::load(&data_dir, month).unwrap();
        
        for i in 0..149 {
            for j in 0..48 {
                assert_eq!(
                    coeff_r[i][j],
                    coeff_f[i][j],
                    "mcsat mismatch for month {} at i={}, j={}",
                    month,
                    i,
                    j
                );
            }
        }
    }
}

#[test]
fn test_ccir_ursi_equivalence() {
    let _lock = get_ffi_mutex().lock().unwrap();
    let data_dir = get_data_dir();
    std::env::set_var("IRI2020_DATA_DIR", &data_dir);
    
    for month in 1..=12 {
        for &is_ccir in &[true, false] {
            let mut f2_f_flat = vec![0.0_f32; 13 * 76 * 2];
            let mut fm3_f_flat = vec![0.0_f32; 9 * 49 * 2];
            unsafe {
                read_coeff_c(month, is_ccir, f2_f_flat.as_mut_ptr(), fm3_f_flat.as_mut_ptr());
            }
            
            // Transpose column-major Fortran arrays to row-major
            let mut f2_f = [[[0.0_f32; 2]; 76]; 13];
            for i in 0..13 {
                for j in 0..76 {
                    for k in 0..2 {
                        let idx = i + 13 * (j + 76 * k);
                        f2_f[i][j][k] = f2_f_flat[idx];
                    }
                }
            }
            
            let mut fm3_f = [[[0.0_f32; 2]; 49]; 9];
            for i in 0..9 {
                for j in 0..49 {
                    for k in 0..2 {
                        let idx = i + 9 * (j + 49 * k);
                        fm3_f[i][j][k] = fm3_f_flat[idx];
                    }
                }
            }
            
            let (f2_r, fm3_r) = CcirUrsiData::load(&data_dir, month, is_ccir).unwrap();
            
            // Check F2
            for i in 0..13 {
                for j in 0..76 {
                    for k in 0..2 {
                        assert_eq!(
                            f2_r[i][j][k],
                            f2_f[i][j][k],
                            "F2 mismatch for month {}, ccir={} at i={}, j={}, k={}",
                            month,
                            is_ccir,
                            i,
                            j,
                            k
                        );
                    }
                }
            }
            
            // Check FM3
            if is_ccir {
                let fm3_r_val = fm3_r.expect("Expected FM3 for CCIR");
                for i in 0..9 {
                    for j in 0..49 {
                        for k in 0..2 {
                            assert_eq!(
                                fm3_r_val[i][j][k],
                                fm3_f[i][j][k],
                                "FM3 mismatch for month {} at i={}, j={}, k={}",
                                month,
                                i,
                                j,
                                k
                            );
                        }
                    }
                }
            } else {
                assert!(fm3_r.is_none(), "Expected no FM3 for URSI");
            }
        }
    }
}
