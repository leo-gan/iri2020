"""Error and speed metrics for IRI surrogate evaluation."""

from __future__ import annotations

from typing import Any
import time

import numpy as np

from .config import LOG_TARGETS


def per_target_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    targets: list[str],
) -> dict[str, dict[str, float]]:
    """Compute MAE, RMSE, MAPE-like, and log10-MAE for log targets."""
    out: dict[str, dict[str, float]] = {}
    y_true = np.asarray(y_true, dtype=np.float64)
    y_pred = np.asarray(y_pred, dtype=np.float64)
    for j, t in enumerate(targets):
        yt = y_true[:, j]
        yp = y_pred[:, j]
        err = yp - yt
        mae = float(np.mean(np.abs(err)))
        rmse = float(np.sqrt(np.mean(err**2)))
        # Symmetric relative error clipped to avoid /0
        denom = np.maximum(np.abs(yt), 1e-12)
        mape = float(np.mean(np.abs(err) / denom))
        med_ae = float(np.median(np.abs(err)))
        m: dict[str, float] = {"mae": mae, "rmse": rmse, "mape": mape, "med_ae": med_ae}
        if t in LOG_TARGETS:
            yt_l = np.log10(np.maximum(yt, 1e-3))
            yp_l = np.log10(np.maximum(yp, 1e-3))
            m["log10_mae"] = float(np.mean(np.abs(yp_l - yt_l)))
            m["log10_rmse"] = float(np.sqrt(np.mean((yp_l - yt_l) ** 2)))
        out[t] = m
    return out


def aggregate_metrics(per_t: dict[str, dict[str, float]]) -> dict[str, float]:
    maes = [v["mae"] for v in per_t.values()]
    rmses = [v["rmse"] for v in per_t.values()]
    mapes = [v["mape"] for v in per_t.values()]
    log_maes = [v["log10_mae"] for v in per_t.values() if "log10_mae" in v]
    agg = {
        "mean_mae": float(np.mean(maes)),
        "mean_rmse": float(np.mean(rmses)),
        "mean_mape": float(np.mean(mapes)),
    }
    if log_maes:
        agg["mean_log10_mae"] = float(np.mean(log_maes))
    return agg


def timed_call(fn, n_warmup: int = 1, n_runs: int = 3) -> tuple[Any, float]:
    """Return (last_result, mean_seconds_per_call)."""
    for _ in range(n_warmup):
        fn()
    times = []
    result = None
    for _ in range(n_runs):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return result, float(np.mean(times))


def regime_split_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    targets: list[str],
    regime: np.ndarray,
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name, mask in [
        ("all", np.ones(len(regime), dtype=bool)),
        ("nominal", regime == 0),
        ("extreme", regime == 1),
    ]:
        if mask.sum() == 0:
            continue
        pt = per_target_metrics(y_true[mask], y_pred[mask], targets)
        out[name] = {
            "per_target": pt,
            "aggregate": aggregate_metrics(pt),
            "n": int(mask.sum()),
        }
    return out
