use crate::ffi::*;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::sync::{Mutex, OnceLock};

static IRI_MUTEX: OnceLock<Mutex<()>> = OnceLock::new();

fn get_iri_mutex() -> &'static Mutex<()> {
    IRI_MUTEX.get_or_init(|| Mutex::new(()))
}

pub struct IriResult {
    pub altkm: Vec<f32>,
    pub outf: Vec<f32>,
    pub oarr: Vec<f32>,
}

pub fn run_iri(
    year: i32,
    month: i32,
    day: i32,
    hour: f32,
    glat: f32,
    glon: f32,
    alt_range: [f32; 3],
) -> Result<IriResult, PyErr> {
    let _lock = get_iri_mutex().lock().unwrap();
    if std::env::var("IRI2020_DATA_DIR").is_err() {
        let paths = ["src/data", "../data", "../../src/data", "data"];
        for path in &paths {
            if std::path::Path::new(path).exists() {
                if let Ok(abs_path) = std::fs::canonicalize(path) {
                    if let Some(path_str) = abs_path.to_str() {
                        std::env::set_var("IRI2020_DATA_DIR", path_str);
                        break;
                    }
                }
            }
        }
    }

    let mut jf = [true; 50];
    for item in jf.iter_mut().take(6).skip(3) { *item = false; }
    jf[21] = false;
    jf[22] = false;
    jf[29] = false;
    jf[32] = false;
    jf[33] = false;
    jf[34] = false;
    for item in jf.iter_mut().take(40).skip(38) { *item = false; }
    jf[46] = false;

    let jmag = 0;
    let mmdd = month * 100 + day;
    let dhour = hour;
    let dhour_plus_25 = dhour + 25.0;

    unsafe {
        init_igrf_c();
        read_ig_rz_c();
        readapf107_c();
    }

    let heibeg = alt_range[0];
    let heiend = alt_range[1];
    let heistp = alt_range[2];

    if heistp <= 0.0 {
        return Err(PyValueError::new_err("heistp must be greater than 0.0"));
    }
    if heiend < heibeg {
        return Err(PyValueError::new_err("heiend must be greater than or equal to heibeg"));
    }

    let num_alt = ((heiend - heibeg) / heistp) as usize + 1;
    if num_alt > 1000 {
        return Err(PyValueError::new_err("num_alt must be less than or equal to 1000"));
    }

    let mut outf = vec![0.0_f32; 20 * 1000];
    let mut oarr = [0.0_f32; 100];

    crate::irisub::iri_sub(
        &jf,
        jmag,
        glat,
        glon,
        year,
        mmdd,
        dhour_plus_25,
        heibeg,
        heiend,
        heistp,
        &mut outf,
        &mut oarr,
    );

    let mut tecbo = 0.0_f32;
    let mut tecto = 0.0_f32;

    crate::iritec::iritec(
        glat,
        glon,
        jmag,
        &jf,
        year,
        mmdd,
        dhour_plus_25,
        0.0,
        heiend,
        0.1,
        &mut oarr,
        &mut tecbo,
        &mut tecto,
    );

    oarr[36] = tecbo + tecto;
    if oarr[36] != 0.0 {
        oarr[37] = tecto / oarr[36] * 100.0;
    }

    let mut altkm = Vec::with_capacity(num_alt);
    for i in 0..num_alt {
        altkm.push(heibeg + (i as f32) * heistp);
    }

    Ok(IriResult {
        altkm,
        outf,
        oarr: oarr.to_vec(),
    })
}

