pub mod ffi;
pub mod driver;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use numpy::{IntoPyArray, PyArray2};

#[pyfunction]
#[pyo3(signature = (year, month, day, hour, glat, glon, altkmrange))]
fn run_iri_py<'py>(
    py: Python<'py>,
    year: i32,
    month: i32,
    day: i32,
    hour: f32,
    glat: f32,
    glon: f32,
    altkmrange: [f32; 3],
) -> PyResult<&'py PyDict> {

    let result = driver::run_iri(year, month, day, hour, glat, glon, altkmrange);

    let dict = PyDict::new(py);

    let altkm_py = result.altkm.into_pyarray(py);
    dict.set_item("altkm", altkm_py)?;

    let oarr_py = result.oarr.into_pyarray(py);
    dict.set_item("oarr", oarr_py)?;

    let outf_py = PyArray2::from_vec2(py, &result.outf).unwrap();
    dict.set_item("outf", outf_py)?;

    Ok(dict)
}

#[pymodule]
fn iri2020(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_iri_py, m)?)?;
    Ok(())
}
