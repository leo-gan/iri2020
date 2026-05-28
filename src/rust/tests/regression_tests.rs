use serde::Deserialize;
use iri2020::driver::run_iri;

#[derive(Deserialize)]
struct GoldenData {
    altkm: Vec<Option<f32>>,
    oarr: Vec<Option<f32>>,
    outf: Vec<Vec<Option<f32>>>,
}

fn assert_approx_eq(a: &[f32], b: &[Option<f32>], tol: f32) {
    assert_eq!(a.len(), b.len(), "Vectors have different lengths");
    for (i, (&x, &y_opt)) in a.iter().zip(b.iter()).enumerate() {
        let y = y_opt.unwrap_or(f32::NAN);
        if x.is_nan() && y.is_nan() {
            continue;
        }
        let diff = (x - y).abs();
        let max_val = x.abs().max(y.abs());
        let allowed_diff = if max_val > 1.0 { tol * max_val } else { tol };
        assert!(
            diff <= allowed_diff,
            "Mismatch at index {}: {} vs {} (diff = {}, tol = {})",
            i, x, y, diff, allowed_diff
        );
    }
}

fn verify_scenario(
    golden_str: &str,
    year: i32,
    month: i32,
    day: i32,
    hour: f32,
    glat: f32,
    glon: f32,
    alt_range: [f32; 3]
) {
    let sanitized = golden_str.replace("NaN", "null");
    let golden: GoldenData = serde_json::from_str(&sanitized).unwrap();
    let result = run_iri(year, month, day, hour, glat, glon, alt_range).unwrap();

    // Verify altkm
    assert_approx_eq(&result.altkm, &golden.altkm, 1e-4);

    // Verify oarr
    assert_approx_eq(&result.oarr, &golden.oarr, 1e-4);

    // Verify outf (columns matching the altitude range)
    let num_alt = result.altkm.len();
    assert_eq!(golden.outf.len(), 20, "Golden outf must contain 20 variables");
    for row in 0..20 {
        assert_eq!(
            golden.outf[row].len(),
            num_alt,
            "Golden outf row {} length does not match num_alt ({})",
            row,
            num_alt
        );
        for col in 0..num_alt {
            let computed_val = result.outf[col * 20 + row];
            let expected_val = golden.outf[row][col].unwrap_or(f32::NAN);
            if computed_val.is_nan() && expected_val.is_nan() {
                continue;
            }
            let diff = (computed_val - expected_val).abs();
            let max_val = computed_val.abs().max(expected_val.abs());
            let allowed = if max_val > 1.0 { 1e-4 * max_val } else { 1e-4 };
            assert!(
                diff <= allowed,
                "outf mismatch at row {}, col {}: computed {}, expected {} (diff = {}, allowed = {})",
                row,
                col,
                computed_val,
                expected_val,
                diff,
                allowed
            );
        }
    }
}

#[test]
fn test_all_scenarios() {
    // Scenario 1: lat=−11.95, lon=−76.77, date=2003-11-21, alt=100–300 km
    let golden_str1 = include_str!("fixtures/scenario1.json");
    verify_scenario(
        golden_str1,
        2003,
        11,
        21,
        12.0,
        -11.95,
        -76.77,
        [100.0, 300.0, 20.0]
    );

    // Scenario 2: lat=65.1, lon=-147.5, date=2015-12-13, alt=100-1000 km
    let golden_str2 = include_str!("fixtures/scenario2.json");
    verify_scenario(
        golden_str2,
        2015,
        12,
        13,
        10.0,
        65.1,
        -147.5,
        [100.0, 1000.0, 10.0]
    );

    // Scenario 3: lat=0.0, lon=0.0, date=2020-06-21, alt=100-200 km
    let golden_str3 = include_str!("fixtures/scenario3.json");
    verify_scenario(
        golden_str3,
        2020,
        6,
        21,
        12.0,
        0.0,
        0.0,
        [100.0, 200.0, 50.0]
    );
}
