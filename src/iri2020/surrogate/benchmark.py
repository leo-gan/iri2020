"""Benchmark harness: surrogates vs IRI backends (Rust primary, Fortran optional).

Compares accuracy on held-out / extreme samples and throughput (samples/sec).
Runs are kept short by default — enough to produce comparable metrics, not
publication-grade statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable
import json

import numpy as np
import torch

from iri2020.base import IRI

from .config import SurrogateConfig
from .data import IRISampleBatch, generate_samples, load_batch
from .preprocessing import IRIPreprocessor
from .metrics import regime_split_metrics, timed_call, per_target_metrics, aggregate_metrics
from .models.residual_mlp import ResidualFourierMLP
from .models.film_mlp import FiLMConditionedMLP
from .models.ensemble import DeepEnsemble, load_ensemble_into, build_ensemble
from .models.xgboost_baseline import XGBoostBaseline
from .train import batch_to_tensors


@dataclass
class Predictor:
    name: str
    predict_fn: Callable[[IRISampleBatch], np.ndarray]
    supports_batch: bool = True


def _predict_iri_rust_batch(batch: IRISampleBatch, targets: list[str]) -> np.ndarray:
    """Re-run IRI (Rust via Python) at stored coordinates — ground truth / reference speed."""
    y = np.full((len(batch.doy), len(targets)), np.nan, dtype=np.float64)
    for i in range(len(batch.doy)):
        t = datetime(int(batch.year[i]), int(batch.month[i]), int(batch.day[i]),
                     int(batch.hour[i]) % 24, int((batch.hour[i] % 1) * 60), 0)
        a = float(batch.alt_km[i])
        try:
            ds = IRI(t, [a, a, 1.0], float(batch.glat[i]), float(batch.glon[i]))
            for j, name in enumerate(targets):
                if name in ds and "alt_km" in ds[name].dims:
                    y[i, j] = float(ds[name].values.ravel()[0])
                else:
                    y[i, j] = float(np.asarray(ds[name].values).ravel()[0])
        except BaseException:
            pass
    return y


def try_fortran_available() -> bool:
    """Fortran driver is not wired into the Python package anymore (Rust port).

    We still report it as unavailable unless a future binding appears.
    Benchmark records this explicitly rather than silently skipping.
    """
    return False


def load_models(artifact_dir: Path, config: SurrogateConfig, pre: IRIPreprocessor, device: str = "cpu"):
    models = {}
    res_path = artifact_dir / "residual_mlp.pt"
    if res_path.exists():
        m = ResidualFourierMLP(pre.input_dim(), pre.output_dim(), config.res_hidden, config.res_blocks, config.res_dropout)
        m.load_state_dict(torch.load(res_path, map_location=device, weights_only=True))
        m.eval()
        models["residual_mlp"] = m

    film_path = artifact_dir / "film_mlp.pt"
    if film_path.exists():
        m = FiLMConditionedMLP(pre.input_dim(), pre.cond_dim(), pre.output_dim(),
                               config.film_hidden, config.film_blocks, config.film_dropout)
        m.load_state_dict(torch.load(film_path, map_location=device, weights_only=True))
        m.eval()
        models["film_mlp"] = m

    ens_path = artifact_dir / "film_ensemble.pt"
    if ens_path.exists():
        def factory():
            return FiLMConditionedMLP(pre.input_dim(), pre.cond_dim(), pre.output_dim(),
                                      config.film_hidden, config.film_blocks, config.film_dropout)
        ens = build_ensemble(factory, config.ensemble_size)
        load_ensemble_into(ens, ens_path, map_location=device)
        ens.eval()
        models["film_ensemble"] = ens

    xgb_path = artifact_dir / "xgboost.joblib"
    if xgb_path.exists():
        xgb = XGBoostBaseline()
        xgb.load(xgb_path)
        models["xgboost"] = xgb

    return models


def nn_predict_physical(model, pre: IRIPreprocessor, batch: IRISampleBatch, use_cond: bool, device: str = "cpu") -> np.ndarray:
    X, C, _ = batch_to_tensors(batch, pre)
    xt = torch.from_numpy(X).to(device)
    ct = torch.from_numpy(C).to(device)
    model.eval()
    with torch.no_grad():
        if isinstance(model, DeepEnsemble):
            y_norm, _ = model.predict_with_uncertainty(xt, ct if use_cond else None)
            y_norm = y_norm.cpu().numpy()
        else:
            y_norm = model(xt, ct if use_cond else None).cpu().numpy()
    return pre.inverse_y(y_norm)


def run_benchmark(
    artifact_dir: Path | str,
    test_batch: IRISampleBatch | None = None,
    extreme_batch: IRISampleBatch | None = None,
    config: SurrogateConfig | None = None,
    n_speed_samples: int = 32,
    out_path: Path | str | None = None,
) -> dict[str, Any]:
    artifact_dir = Path(artifact_dir)
    config = config or SurrogateConfig()
    if (artifact_dir / "config.json").exists():
        with open(artifact_dir / "config.json") as f:
            config = SurrogateConfig.from_dict(json.load(f))

    pre = IRIPreprocessor.load(artifact_dir / "preprocessor.json")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    models = load_models(artifact_dir, config, pre, device=device)

    # Build evaluation set
    if test_batch is None:
        test_path = artifact_dir / "test_batch.npz"
        if test_path.exists():
            test_batch = load_batch(test_path)
        else:
            print("generating small test batch ...", flush=True)
            test_batch = generate_samples(min(200, config.n_test), config, extreme=False)

    if extreme_batch is None:
        ext_path = artifact_dir / "extreme_batch.npz"
        if ext_path.exists():
            extreme_batch = load_batch(ext_path)
        else:
            print("generating small extreme batch ...", flush=True)
            extreme_batch = generate_samples(min(100, config.n_extreme), config, extreme=True)

    eval_batch = IRISampleBatch(
        doy=np.concatenate([test_batch.doy, extreme_batch.doy]),
        hour=np.concatenate([test_batch.hour, extreme_batch.hour]),
        year=np.concatenate([test_batch.year, extreme_batch.year]),
        month=np.concatenate([test_batch.month, extreme_batch.month]),
        day=np.concatenate([test_batch.day, extreme_batch.day]),
        glat=np.concatenate([test_batch.glat, extreme_batch.glat]),
        glon=np.concatenate([test_batch.glon, extreme_batch.glon]),
        alt_km=np.concatenate([test_batch.alt_km, extreme_batch.alt_km]),
        f107=np.concatenate([test_batch.f107, extreme_batch.f107]),
        ap=np.concatenate([test_batch.ap, extreme_batch.ap]),
        y=np.concatenate([test_batch.y, extreme_batch.y], axis=0),
        targets=test_batch.targets,
        regime=np.concatenate([
            np.zeros(len(test_batch.doy), dtype=np.int8),
            np.ones(len(extreme_batch.doy), dtype=np.int8),
        ]),
    )

    targets = eval_batch.targets
    y_true = eval_batch.y
    report: dict[str, Any] = {
        "n_eval": int(len(eval_batch.doy)),
        "n_nominal": int((eval_batch.regime == 0).sum()),
        "n_extreme": int((eval_batch.regime == 1).sum()),
        "targets": targets,
        "fortran_available": try_fortran_available(),
        "models": {},
        "speed": {},
        "notes": [],
    }

    if not report["fortran_available"]:
        report["notes"].append(
            "Fortran backend is not exposed through the current Python package "
            "(Rust port is the active runtime). Speed/accuracy vs Fortran would "
            "require restoring the f2py/CMake driver; Rust IRI is the reference."
        )

    # --- Accuracy: each surrogate vs stored IRI labels (Rust-generated) ---
    predictors: list[tuple[str, Callable[[], np.ndarray]]] = []

    if "xgboost" in models:
        xgb = models["xgboost"]

        def _xgb():
            X, C, _ = batch_to_tensors(eval_batch, pre)
            return pre.inverse_y(xgb.predict(X, cond=C))

        predictors.append(("xgboost", _xgb))

    if "residual_mlp" in models:
        m = models["residual_mlp"].to(device)
        predictors.append(("residual_mlp", lambda m=m: nn_predict_physical(m, pre, eval_batch, use_cond=False, device=device)))

    if "film_mlp" in models:
        m = models["film_mlp"].to(device)
        predictors.append(("film_mlp", lambda m=m: nn_predict_physical(m, pre, eval_batch, use_cond=True, device=device)))

    if "film_ensemble" in models:
        m = models["film_ensemble"].to(device)
        predictors.append(("film_ensemble", lambda m=m: nn_predict_physical(m, pre, eval_batch, use_cond=True, device=device)))

    # Self-consistency of stored labels vs live Rust (sanity, small subsample)
    n_check = min(16, len(eval_batch.doy))
    sub = IRISampleBatch(
        doy=eval_batch.doy[:n_check], hour=eval_batch.hour[:n_check],
        year=eval_batch.year[:n_check], month=eval_batch.month[:n_check], day=eval_batch.day[:n_check],
        glat=eval_batch.glat[:n_check], glon=eval_batch.glon[:n_check], alt_km=eval_batch.alt_km[:n_check],
        f107=eval_batch.f107[:n_check], ap=eval_batch.ap[:n_check],
        y=eval_batch.y[:n_check], targets=targets, regime=eval_batch.regime[:n_check],
    )
    live = _predict_iri_rust_batch(sub, targets)
    live_err = per_target_metrics(sub.y, live, targets)
    report["rust_label_sanity"] = {
        "aggregate": aggregate_metrics(live_err),
        "note": "stored labels vs re-run Rust on subset; should be ~0",
    }

    if len(eval_batch.doy) == 0:
        report["notes"].append("Empty evaluation batch — skip accuracy metrics.")
        predictors = []

    for name, fn in predictors:
        y_pred = fn()
        split = regime_split_metrics(y_true, y_pred, targets, eval_batch.regime)
        report["models"][name] = split
        all_m = split.get("all", {}).get("aggregate", {})
        ext_m = split.get("extreme", {}).get("aggregate", {})
        print(
            f"  {name}: all mean_mape={all_m.get('mean_mape', float('nan')):.4f}  "
            f"extreme mean_mape={ext_m.get('mean_mape', float('nan')):.4f}",
            flush=True,
        )

    # --- Speed ---
    speed_n = min(n_speed_samples, len(eval_batch.doy))
    speed_batch = IRISampleBatch(
        doy=eval_batch.doy[:speed_n], hour=eval_batch.hour[:speed_n],
        year=eval_batch.year[:speed_n], month=eval_batch.month[:speed_n], day=eval_batch.day[:speed_n],
        glat=eval_batch.glat[:speed_n], glon=eval_batch.glon[:speed_n], alt_km=eval_batch.alt_km[:speed_n],
        f107=eval_batch.f107[:speed_n], ap=eval_batch.ap[:speed_n],
        y=eval_batch.y[:speed_n], targets=targets, regime=eval_batch.regime[:speed_n],
    )

    def bench_speed(label: str, fn: Callable[[], Any]):
        _, secs = timed_call(fn, n_warmup=1, n_runs=3)
        sps = speed_n / secs if secs > 0 else float("inf")
        report["speed"][label] = {
            "n_samples": speed_n,
            "mean_seconds": secs,
            "samples_per_sec": sps,
            "ms_per_sample": 1000.0 * secs / speed_n,
        }
        print(f"  speed {label}: {sps:.1f} samples/s ({1000*secs/speed_n:.2f} ms/sample)", flush=True)

    bench_speed("iri_rust", lambda: _predict_iri_rust_batch(speed_batch, targets))

    if "xgboost" in models:
        xgb = models["xgboost"]
        X, C, _ = batch_to_tensors(speed_batch, pre)
        bench_speed("xgboost", lambda: xgb.predict(X, cond=C))

    if "residual_mlp" in models:
        m = models["residual_mlp"].to(device)
        bench_speed("residual_mlp", lambda: nn_predict_physical(m, pre, speed_batch, False, device))

    if "film_mlp" in models:
        m = models["film_mlp"].to(device)
        bench_speed("film_mlp", lambda: nn_predict_physical(m, pre, speed_batch, True, device))

    if "film_ensemble" in models:
        m = models["film_ensemble"].to(device)
        bench_speed("film_ensemble", lambda: nn_predict_physical(m, pre, speed_batch, True, device))

    report["speed"]["fortran"] = {
        "available": False,
        "note": "Not wired in Python package; see docs/surrogate.md",
    }

    out_path = Path(out_path) if out_path else artifact_dir / "benchmark_report.json"
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"wrote {out_path}", flush=True)
    return report
