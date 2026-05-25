pub mod ffi;
pub mod driver;
pub mod igrf_coeff;
pub mod igrf;
pub mod cira_coeff;
pub mod cira;
pub mod rocdrift_coeff;
pub mod rocdrift;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use numpy::{IntoPyArray, ToPyArray};
use numpy::ndarray::{Array2, s, ShapeBuilder};

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

    let result = driver::run_iri(year, month, day, hour, glat, glon, altkmrange)?;

    let dict = PyDict::new(py);

    let num_alt = result.altkm.len();

    let altkm_py = result.altkm.into_pyarray(py);
    dict.set_item("altkm", altkm_py)?;

    let oarr_py = result.oarr.into_pyarray(py);
    dict.set_item("oarr", oarr_py)?;

    let outf_arr = Array2::from_shape_vec((20, 1000).f(), result.outf)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("Reshape error: {}", e)))?;

    let outf_sliced = outf_arr.slice(s![.., 0..num_alt]);
    let outf_py = outf_sliced.to_pyarray(py);
    dict.set_item("outf", outf_py)?;

    Ok(dict)
}

#[pymodule]
fn iri2020(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_iri_py, m)?)?;
    Ok(())
}

