use iri2020::irifun_utils::*;
use iri2020::igrf::IgrfModel;
use std::path::Path;

fn get_data_dir() -> String {
    if let Ok(dir) = std::env::var("IRI2020_DATA_DIR") {
        return dir;
    }
    let paths = ["src/data", "../data", "../../src/data", "data", "src/rust/src/data"];
    for path in &paths {
        if Path::new(path).exists() {
            if let Ok(abs_path) = std::fs::canonicalize(path) {
                if let Some(path_str) = abs_path.to_str() {
                    std::env::set_var("IRI2020_DATA_DIR", path_str);
                    return path_str.to_string();
                }
            }
        }
    }
    // Try tracing up if nested in test run
    let mut current = std::env::current_dir().unwrap();
    loop {
        let candidate = current.join("src/data");
        if candidate.exists() {
            let path_str = candidate.to_str().unwrap().to_string();
            std::env::set_var("IRI2020_DATA_DIR", &path_str);
            return path_str;
        }
        if !current.pop() {
            break;
        }
    }
    panic!("Could not find data directory");
}

#[test]
fn test_trig_wrappers() {
    // acos clamping
    assert_eq!(acos(1.5), 0.0);
    assert_eq!(acos(-1.5), PI);
    assert!((acos(0.5) - PI / 3.0).abs() < 1e-5);

    // asin clamping
    assert_eq!(asin(1.5), PI / 2.0);
    assert_eq!(asin(-1.5), -PI / 2.0);
    assert!((asin(0.5) - PI / 6.0).abs() < 1e-5);

    // atan2d
    assert!((atan2d(1.0, 1.0) - 45.0).abs() < 1e-5);
    assert!((atan2d(1.0, 0.0) - 90.0).abs() < 1e-5);
}

#[test]
fn test_epstein_transitions() {
    // eptr
    assert!((eptr(5.0, 2.0, 5.0) - 2.0_f32.ln()).abs() < 1e-5);
    // extreme values
    assert_eq!(eptr(5.0 + 90.0 * 2.0, 2.0, 5.0), 90.0);
    assert_eq!(eptr(5.0 - 90.0 * 2.0, 2.0, 5.0), 0.0);

    // epst
    assert!((epst(5.0, 2.0, 5.0) - 0.5).abs() < 1e-5);
    assert_eq!(epst(5.0 + 90.0 * 2.0, 2.0, 5.0), 1.0);
    assert_eq!(epst(5.0 - 90.0 * 2.0, 2.0, 5.0), 0.0);

    // epstep
    assert!((epstep(10.0, 2.0, 2.0, 5.0, 5.0) - 6.0).abs() < 1e-5);

    // epla
    assert!((epla(5.0, 2.0, 5.0) - 0.25).abs() < 1e-5);
    assert_eq!(epla(5.0 + 90.0 * 2.0, 2.0, 5.0), 0.0);
}

#[test]
fn test_lay_functions() {
    let x = 6.0;
    let xm = 5.0;
    let sc = 2.0;
    let hx = 5.0;

    let y1 = eptr(x, sc, hx);
    let y1m = eptr(xm, sc, hx);
    let y2m = epst(xm, sc, hx);
    let expected_rlay = y1 - y1m - (x - xm) * y2m / sc;
    assert!((rlay(x, xm, sc, hx) - expected_rlay).abs() < 1e-5);

    let expected_d1lay = (epst(x, sc, hx) - epst(xm, sc, hx)) / sc;
    assert!((d1lay(x, xm, sc, hx) - expected_d1lay).abs() < 1e-5);

    let expected_d2lay = epla(x, sc, hx) / (sc * sc);
    assert!((d2lay(x, xm, sc, hx) - expected_d2lay).abs() < 1e-5);
}

#[test]
fn test_booker() {
    let ah = [100.0, 110.0, 120.0];
    let av = [1.0, 2.0, 3.0];
    let d = [2.0];
    let val = booker(105.0, 3, &ah, &av, &d);
    assert!(val > 0.0);

    let st = [0.1, 0.1];
    let val1 = booker1(105.0, 1, 1.0, &ah, &st, &d);
    assert!(val1 > 0.0);
}

#[test]
fn test_regfa1() {
    // Find where x^2 - 4.0 = 0
    let f = |x: f32| x * x;
    let (schalt, x_root) = regfa1(1.0, 3.0, 1.0, 9.0, 0.001, 4.0, f);
    assert!(!schalt);
    assert!((x_root - 2.0).abs() < 0.01);
}

#[test]
fn test_lnglsn() {
    // Solve:
    //  x +  y = 3
    // 2x +  y = 4
    // Solution is x=1, y=2.
    let mut a = [[0.0_f32; 5]; 5];
    let mut b = [0.0_f32; 5];
    
    a[0][0] = 1.0; a[0][1] = 1.0;
    a[1][0] = 2.0; a[1][1] = 1.0;
    
    b[0] = 3.0;
    b[1] = 4.0;
    
    let n = 2;
    let aus = lnglsn(n, &mut a, &mut b);
    assert!(!aus);
    
    // Solution is stored in A(N, J) for J=1, N
    // In Rust, that is row index n-1, columns 0..n
    let x = a[1][0];
    let y = a[1][1];
    assert!((x - 1.0).abs() < 1e-4);
    assert!((y - 2.0).abs() < 1e-4);
}

#[test]
fn test_moda() {
    // Leap year (2000)
    let mut month = 3;
    let mut iday = 1;
    let mut idoy = 0;
    let mut nrdaymo = 0;
    moda(0, 2000, &mut month, &mut iday, &mut idoy, &mut nrdaymo);
    assert_eq!(idoy, 61);
    assert_eq!(nrdaymo, 31);

    let mut month_inv = 0;
    let mut iday_inv = 0;
    let mut idoy_mut = 61;
    moda(1, 2000, &mut month_inv, &mut iday_inv, &mut idoy_mut, &mut nrdaymo);
    assert_eq!(month_inv, 3);
    assert_eq!(iday_inv, 1);

    // Non-leap year (2001)
    month = 3;
    iday = 1;
    moda(0, 2001, &mut month, &mut iday, &mut idoy, &mut nrdaymo);
    assert_eq!(idoy, 60);

    month_inv = 0;
    iday_inv = 0;
    idoy_mut = 60;
    moda(1, 2001, &mut month_inv, &mut iday_inv, &mut idoy_mut, &mut nrdaymo);
    assert_eq!(month_inv, 3);
    assert_eq!(iday_inv, 1);
}

#[test]
fn test_tcon() {
    let mut ionoindx = [0.0_f32; 806];
    let mut indrz = [0.0_f32; 806];
    // Fill indices for year 2020 month 1 and 2
    // Let's set start to 202001
    let iymst = 202001;
    let iymend = 202012;
    // index 0 -> 201912 (month before start month)
    // index 1 -> 202001 (start month)
    // index 2 -> 202002 (month after start month)
    indrz[0] = 5.0;
    indrz[1] = 10.0;
    indrz[2] = 20.0;
    ionoindx[0] = 90.0;
    ionoindx[1] = 100.0;
    ionoindx[2] = 110.0;

    let mut rz = [0.0_f32; 3];
    let mut ig = [0.0_f32; 3];
    let mut rsn = 0.0_f32;
    let mut nmonth = 0;

    // Test on 2020-01-15 (mid month)
    let res = tcon(
        2020, 1, 15, 15,
        &ionoindx, &indrz,
        iymst, iymend,
        &mut rz, &mut ig,
        &mut rsn, &mut nmonth,
        false
    );
    assert!(res.is_ok());
    assert_eq!(rz[0], 10.0);
    assert_eq!(ig[0], 100.0);
}

#[test]
fn test_tbfit() {
    let (thint, tzero) = tbfit(1.0, 5.0, 0);
    assert_eq!(tzero, 3.0);
    assert_eq!(thint, 2.0);

    let (thint2, tzero2) = tbfit(1.0, 5.0, 2);
    assert_eq!(tzero2, 1.0);
    assert_eq!(thint2, 2.0);
}

#[test]
fn test_shellg_integration() {
    use iri2020::ffi::shellg_c;
    use std::sync::Mutex;
    use std::sync::OnceLock;

    static FFI_MUTEX: OnceLock<Mutex<()>> = OnceLock::new();
    let _lock = FFI_MUTEX.get_or_init(|| Mutex::new(())).lock().unwrap();

    let data_dir = get_data_dir();
    std::env::set_var("IRI2020_DATA_DIR", &data_dir);

    // Call Fortran FELDCOF to initialize the global variables for shellg
    unsafe {
        iri2020::ffi::init_igrf_c();
        iri2020::ffi::feldcof_c(2020.0);
    }

    // Call Fortran shellg
    let mut f_fl = 0.0_f32;
    let mut f_icode = 0_i32;
    let mut f_b0 = 0.0_f32;
    unsafe {
        shellg_c(-11.95, -76.77, 300.0, &mut f_fl, &mut f_icode, &mut f_b0);
    }

    // Call Rust shellg
    let mut model = IgrfModel::new();
    model.feldcof(2020.0, &data_dir).unwrap();
    let (r_fl, r_icode, r_b0) = shellg(&model, -11.95, -76.77, 300.0);

    // Assert exact numerical match!
    assert_eq!(r_icode, f_icode, "ICODE mismatch");
    assert!((r_fl - f_fl).abs() < 1e-4, "FL mismatch: rust = {}, fortran = {}", r_fl, f_fl);
    assert!((r_b0 - f_b0).abs() < 1e-5, "B0 mismatch: rust = {}, fortran = {}", r_b0, f_b0);
}
