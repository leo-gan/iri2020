# IRI2020 Neural Surrogates

This document describes the neural network and gradient-boosting **surrogate models**
that approximate the IRI2020 empirical ionosphere at pointwise
`(time, latitude, longitude, altitude)` queries, plus the benchmark that compares
them to the original empirical model (Rust runtime; Fortran noted but not wired
into the current Python package).

## Motivation

IRI2020 is accurate as an empirical standard but relatively expensive for Monte Carlo,
real-time assimilation, or large grids. A learned surrogate can:

1. Predict key ionospheric quantities orders of magnitude faster than running IRI.
2. Provide **ensemble uncertainty** as a diagnostic for where the surrogate is untrustworthy.
3. Serve as a differentiable (NN) or fast tabular (XGBoost) stand-in for optimization loops.

**Important:** surrogates approximate the **model**, not nature. Error vs IRI is not
error vs measurements.

## Architecture overview

```
IRI2020 (Rust) ──generate──► training batches (.npz)
                                │
                                ▼
                     IRIPreprocessor (fit on train)
                       • sin/cos periodic encodings
                       • Gaussian Fourier features
                       • log10 for multi-decade outputs
                       • z-score standardization
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
 ResidualFourierMLP      FiLMConditionedMLP      XGBoost baseline
        │                       │
        └──────── Deep ensemble (FiLM members) ─────────┘
                                │
                                ▼
                     benchmark_report.json
              (error by regime, speed vs IRI Rust)
```

### 1. Preprocessing (`iri2020.surrogate.preprocessing`)

| Concern | Approach | Rationale |
|--------|----------|-----------|
| Periodic inputs (DOY, UT hour, longitude) | `(sin, cos)` pairs | Avoids 0/360 and 24 h discontinuities |
| Latitude / altitude | Scaled linear features | Bounded geometry; alt normalized over broad physical range |
| Densities, TEC, Nm* | `log10(max(x, floor))` then z-score | Outputs span 6+ orders of magnitude; MSE without log is dominated by peaks |
| Temperatures, heights, foF2, drifts | z-score only | Log would distort near-zero/negative drifts; heights are not multiplicative |
| High-frequency spatial/temporal structure | Fixed Gaussian Fourier features (Tancik et al.) | MLPs are spectrally biased toward low frequencies |
| Solar / geomagnetic drivers | Separate **condition** vector `(f107, ap, sin/cos doy)` | Enables FiLM modulation without polluting geometry features |

**Trade-offs**

- Global (not altitude-conditional) normalization is simpler to invert but ignores
  heteroscedasticity with height.
- `log_floor` flattens the lowest densities (often IRI zero-pads absent species).
- Fourier feature count/scale is a bias–variance knob; too large overfits small datasets.

### 2. Deep Residual MLP + Fourier features (`ResidualFourierMLP`)

- Stem linear → LayerNorm → GELU
- Stacked **residual blocks** (two linear layers + skip)
- Linear multi-target head in **normalized** space

Residual connections allow modest depth without vanishing gradients; Fourier features
are applied in preprocessing so the module stays a plain MLP.

### 3. FiLM-conditioned MLP (`FiLMConditionedMLP`)

Feature-wise Linear Modulation predicts per-channel `γ, β` from drivers and applies
`h ← γ ⊙ h + β` inside each block. This targets **driver-dependent** amplitude/shape
changes (solar max / storm) that a monolithic MLP tends to average away.

Conditioning uses IRI-returned `f107` and `ap` (from index files), not user-specified
overrides — matching the default Python `IRI()` API.

### 4. Loss

**Huber loss** on normalized multi-target vectors (after log for density-like targets).

- Equalizes targets regardless of physical units.
- Huber is more robust than MSE to occasional IRI outliers / numerical spikes.
- Not physics-informed; no explicit peak-height constraints.

### 5. Deep ensemble UQ (`DeepEnsemble`)

Train `M` FiLM networks with different seeds; report mean and member std at inference.

| Strength | Weakness |
|----------|----------|
| Simple epistemic uncertainty proxy | Cost ×M at train/infer |
| No variational machinery | Std is **not** calibrated coverage |
| Short training ⇒ underfit members | UQ quality poor until epochs/data scale up |

Use ensemble variance as a **warning**, not a confidence interval for science products.

### 6. XGBoost baseline (`XGBoostBaseline`)

`MultiOutputRegressor(XGBRegressor)` on concatenated `[X_features, cond]`.

- Strong tabular baseline, fast to train, few dependencies beyond `xgboost`.
- No shared multi-task representation; no first-class FiLM.
- No native ensemble UQ in this implementation.

## Data generation

`iri2020.surrogate.data.generate_samples` queries the **Rust** IRI through
`iri2020.base.IRI` at random times/locations/altitudes.

- **Nominal**: uniform lat/lon/alt, years 2000–2019 (index file coverage).
- **Extreme drivers**: solar-max years, equinox/storm months, polar + equatorial latitudes.

Labels are stored in `surrogate_artifacts/*_batch.npz`. Regenerate by deleting those files.

## Training & benchmark (short runs)

Default pipeline uses **small N and few epochs** so CI / smoke tests finish quickly and
you can verify **convergence direction** (train/val loss decreasing), not final accuracy.

```bash
source .venv/bin/activate
export IRI2020_DATA_DIR=$PWD/src/data   # required for reliable bulk generation
python -m iri2020.surrogate.scripts.run_pipeline \
  --n-train 600 --n-val 120 --n-test 120 --n-extreme 80 \
  --epochs 6 --ensemble-size 3 \
  --artifact-dir surrogate_artifacts
```

### Example smoke-run results (illustrative, not publication-grade)

Run on ~215 train / 40 val / 55 eval samples, 5 epochs, ensemble size 3 (CPU):

| Model | Val loss trend | All mean MAPE | Extreme mean MAPE | Speed (samples/s) |
|-------|----------------|---------------|-------------------|-------------------|
| IRI Rust (reference) | — | 0 (labels) | 0 | ~21 |
| XGBoost baseline | norm MSE ~0.52 | ~0.65 | ~0.24 | ~12k |
| Residual Fourier MLP | 1.26 → 0.69 train / 1.03 → 0.84 val | high (underfit) | high | ~33k |
| FiLM MLP | 1.17 → 0.62 / 0.89 → 0.76 | high (underfit) | better than all-regime | ~29k |
| FiLM ensemble | members converge | moderate | ~0.81 | ~13k |

**Takeaways from the smoke run**

1. **Losses decrease** — optimizers and preprocessing are wired correctly.
2. **XGBoost wins accuracy** at this data scale; NNs need more data/epochs to compete.
3. **Surrogates are ~500–1500× faster** than live IRI at batch inference.
4. **FiLM / ensemble help extreme regimes** relative to single residual MLP in this run,
   but absolute NN errors remain large — do not deploy without scaling data.
5. **Fortran** not timed (not in Python path); Rust is the empirical reference.

Data generation uses **batched subprocess workers** (`iri_worker`) because occasional
Rust panics (missing coeff edge cases) poison an in-process mutex; subprocesses
isolate failures. Set `IRI2020_DATA_DIR` explicitly.

Artifacts written under `surrogate_artifacts/`:

| File | Contents |
|------|----------|
| `preprocessor.json` | Fitted scalers, Fourier `B`, target list |
| `config.json` | Hyperparameters |
| `residual_mlp.pt` / `film_mlp.pt` | Single-network weights |
| `film_ensemble.pt` | Ensemble member state dicts |
| `xgboost.joblib` | Baseline |
| `train_summary.json` | Loss histories |
| `benchmark_report.json` | Error by regime + speed |

### Metrics

Per target: MAE, RMSE, MAPE (relative), median AE; for log targets also `log10_mae/rmse`.

Reported for regimes: **all**, **nominal**, **extreme**.

Speed: mean wall time over 3 timed runs after 1 warmup; `samples/s` and `ms/sample`.

### Fortran vs Rust

The Python package currently runs the **Rust** port only (`run_iri_py`). The original
Fortran sources remain under `src/fortran/` for reference and historical builds, but
are **not** invoked from the surrogate benchmark. The report records
`fortran_available: false` and explains this explicitly.

If you restore an f2py/CMake Fortran extension, extend `benchmark.py` with a
`_predict_iri_fortran_batch` parallel to the Rust path and add a speed entry.

## Critique & known limitations (honest)

1. **Sample efficiency**: hundreds of points is nowhere near enough for global IRI fidelity;
   expect large MAPE especially on ion species and TEC. Treat numbers as **pipeline sanity**,
   not model quality.
2. **Pointwise not profile-aware**: we predict at a single altitude per sample; vertical
   structure consistency (monotonic topside, peak relationships) is not enforced.
3. **Drivers are endogenous**: f107/ap come from IRI's index lookup at the query date,
   so FiLM sees drivers correlated with season/time — not fully independent interventions.
4. **No negative density constraints** in physical space beyond post-hoc log inversion.
5. **Ensemble UQ is undertrained** in short runs; do not publish σ as calibrated uncertainty.
6. **Fortran not benchmarked** in this package state — document gap, don't pretend parity.
7. **XGBoost** trains one model per target via `MultiOutputRegressor`; memory/time scales
   with target count.

## Proposed improvements (prioritized)

| Priority | Improvement | Why |
|----------|-------------|-----|
| P0 | Scale dataset to 10⁴–10⁵+ samples; 50–200 epochs | Everything else is noise without data |
| P0 | Altitude-profile multitask: predict full height vector per geo/time | Matches IRI use cases; enforces structure |
| P1 | Physics / monotonicity penalties; non-negativity in output head | Reduces implausible profiles |
| P1 | Calibrate ensemble via temperature scaling or CRPS on val | Makes UQ usable |
| P1 | Wire Fortran driver for true backend speed/accuracy matrix | Completes user-requested Fortran/Rust comparison |
| P2 | Conditional normalization by altitude bins | Better multi-scale fit |
| P2 | Learned Fourier/positional encodings | Less hand-tuned `B` |
| P2 | Quantile regression or NGBoost baseline | UQ without deep ensembles |
| P3 | ONNX / TorchScript export for deployment | Production latency |

Several **P0/P1 polish fixes** are implemented in code (Huber loss, regime splits,
sanity check of labels vs live Rust, explicit Fortran unavailable handling, FiLM γ
initialized near identity, joblib/xgb paths, artifact config reload in benchmark).

## API sketch

```python
from iri2020.surrogate.config import SurrogateConfig
from iri2020.surrogate.preprocessing import IRIPreprocessor
from iri2020.surrogate.train import run_full_training_pipeline
from iri2020.surrogate.benchmark import run_benchmark

cfg = SurrogateConfig(epochs=6, n_train=600)
# ... generate / load batches ...
run_full_training_pipeline(train_batch, val_batch, cfg)
run_benchmark("surrogate_artifacts")
```

## Dependencies

Optional ML stack (not required for core IRI):

- `torch`
- `xgboost`
- `scikit-learn`
- `joblib`

Install into the project venv:

```bash
source .venv/bin/activate
uv pip install torch xgboost scikit-learn joblib tqdm pyyaml
```

## References

- Bilitza et al., IRI-2020, JSWSC (2022)
- Tancik et al., Fourier Features Let Networks Learn High Frequency Functions, NeurIPS 2020
- Perez et al., FiLM, AAAI 2018
- Lakshminarayanan et al., Deep Ensembles, NeurIPS 2017
- Chen & Guestrin, XGBoost, KDD 2016
