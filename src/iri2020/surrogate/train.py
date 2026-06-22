"""Training loop for IRI neural surrogates.

Loss
----
We train on **normalized multi-target MSE** after log-transform for density-like
outputs. This equalizes contribution across targets spanning many orders of
magnitude. Optionally we add a small **relative L1 in normalized space** for
robustness to outliers (extreme regimes).

Short training
--------------
Default epochs are intentionally low (smoke / convergence check). For production
surrogates, increase `config.epochs`, data volume, and ensemble size.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import time

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from .config import SurrogateConfig
from .data import IRISampleBatch
from .preprocessing import IRIPreprocessor
from .models.residual_mlp import ResidualFourierMLP
from .models.film_mlp import FiLMConditionedMLP
from .models.ensemble import DeepEnsemble, build_ensemble, save_ensemble
from .models.xgboost_baseline import XGBoostBaseline


def batch_to_tensors(
    batch: IRISampleBatch, pre: IRIPreprocessor
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = pre.transform_X(batch.doy, batch.hour, batch.glat, batch.glon, batch.alt_km)
    C = pre.transform_cond(batch.f107, batch.ap, batch.doy)
    Y = pre.transform_y(batch.y)
    return X, C, Y


def make_loader(
    X: np.ndarray, C: np.ndarray, Y: np.ndarray, batch_size: int, shuffle: bool
) -> DataLoader:
    ds = TensorDataset(
        torch.from_numpy(X),
        torch.from_numpy(C),
        torch.from_numpy(Y),
    )
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, drop_last=False)


def multi_target_loss(
    pred: torch.Tensor, target: torch.Tensor, huber_delta: float = 1.0
) -> torch.Tensor:
    """Huber on normalized targets — robust to occasional extreme IRI outliers."""
    return nn.functional.huber_loss(pred, target, delta=huber_delta, reduction="mean")


def train_one_model(
    model: nn.Module,
    train_loader: DataLoader,
    val_loader: DataLoader,
    config: SurrogateConfig,
    device: str,
    use_cond: bool,
) -> dict[str, Any]:
    model = model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(1, config.epochs))

    history: dict[str, list[float]] = {"train_loss": [], "val_loss": []}
    best_val = float("inf")
    best_state = None
    patience_left = config.patience

    for epoch in range(config.epochs):
        model.train()
        tr_losses = []
        for xb, cb, yb in train_loader:
            xb, cb, yb = xb.to(device), cb.to(device), yb.to(device)
            opt.zero_grad(set_to_none=True)
            pred = model(xb, cb if use_cond else None)
            loss = multi_target_loss(pred, yb)
            loss.backward()
            if config.grad_clip > 0:
                nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
            opt.step()
            tr_losses.append(float(loss.item()))
        sched.step()

        model.eval()
        va_losses = []
        with torch.no_grad():
            for xb, cb, yb in val_loader:
                xb, cb, yb = xb.to(device), cb.to(device), yb.to(device)
                pred = model(xb, cb if use_cond else None)
                va_losses.append(float(multi_target_loss(pred, yb).item()))

        tr_m = float(np.mean(tr_losses)) if tr_losses else float("nan")
        va_m = float(np.mean(va_losses)) if va_losses else float("nan")
        history["train_loss"].append(tr_m)
        history["val_loss"].append(va_m)
        print(f"    epoch {epoch+1}/{config.epochs}  train={tr_m:.5f}  val={va_m:.5f}", flush=True)

        if va_m < best_val:
            best_val = va_m
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            patience_left = config.patience
        else:
            patience_left -= 1
            if patience_left <= 0:
                print("    early stop", flush=True)
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return {"history": history, "best_val": best_val}


def fit_preprocessor(train: IRISampleBatch, config: SurrogateConfig) -> IRIPreprocessor:
    pre = IRIPreprocessor(config)
    pre.fit(
        train.doy, train.hour, train.glat, train.glon, train.alt_km,
        train.f107, train.ap, train.y, targets=train.targets,
    )
    return pre


def train_residual_mlp(
    train: IRISampleBatch,
    val: IRISampleBatch,
    config: SurrogateConfig,
    pre: IRIPreprocessor | None = None,
    device: str | None = None,
) -> tuple[ResidualFourierMLP, IRIPreprocessor, dict]:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    pre = pre or fit_preprocessor(train, config)
    Xtr, Ctr, Ytr = batch_to_tensors(train, pre)
    Xva, Cva, Yva = batch_to_tensors(val, pre)
    tr_loader = make_loader(Xtr, Ctr, Ytr, config.batch_size, shuffle=True)
    va_loader = make_loader(Xva, Cva, Yva, config.batch_size, shuffle=False)

    model = ResidualFourierMLP(
        in_dim=pre.input_dim(),
        out_dim=pre.output_dim(),
        hidden=config.res_hidden,
        n_blocks=config.res_blocks,
        dropout=config.res_dropout,
    )
    print("  training ResidualFourierMLP ...", flush=True)
    info = train_one_model(model, tr_loader, va_loader, config, device, use_cond=False)
    return model, pre, info


def train_film_mlp(
    train: IRISampleBatch,
    val: IRISampleBatch,
    config: SurrogateConfig,
    pre: IRIPreprocessor | None = None,
    device: str | None = None,
) -> tuple[FiLMConditionedMLP, IRIPreprocessor, dict]:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    pre = pre or fit_preprocessor(train, config)
    Xtr, Ctr, Ytr = batch_to_tensors(train, pre)
    Xva, Cva, Yva = batch_to_tensors(val, pre)
    tr_loader = make_loader(Xtr, Ctr, Ytr, config.batch_size, shuffle=True)
    va_loader = make_loader(Xva, Cva, Yva, config.batch_size, shuffle=False)

    model = FiLMConditionedMLP(
        in_dim=pre.input_dim(),
        cond_dim=pre.cond_dim(),
        out_dim=pre.output_dim(),
        hidden=config.film_hidden,
        n_blocks=config.film_blocks,
        dropout=config.film_dropout,
    )
    print("  training FiLMConditionedMLP ...", flush=True)
    info = train_one_model(model, tr_loader, va_loader, config, device, use_cond=True)
    return model, pre, info


def train_ensemble(
    kind: str,
    train: IRISampleBatch,
    val: IRISampleBatch,
    config: SurrogateConfig,
    pre: IRIPreprocessor,
    device: str | None = None,
) -> tuple[DeepEnsemble, dict]:
    """Train a deep ensemble of residual or film models."""
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    Xtr, Ctr, Ytr = batch_to_tensors(train, pre)
    Xva, Cva, Yva = batch_to_tensors(val, pre)
    tr_loader = make_loader(Xtr, Ctr, Ytr, config.batch_size, shuffle=True)
    va_loader = make_loader(Xva, Cva, Yva, config.batch_size, shuffle=False)
    use_cond = kind == "film"

    def factory():
        if kind == "film":
            return FiLMConditionedMLP(
                in_dim=pre.input_dim(), cond_dim=pre.cond_dim(), out_dim=pre.output_dim(),
                hidden=config.film_hidden, n_blocks=config.film_blocks, dropout=config.film_dropout,
            )
        return ResidualFourierMLP(
            in_dim=pre.input_dim(), out_dim=pre.output_dim(),
            hidden=config.res_hidden, n_blocks=config.res_blocks, dropout=config.res_dropout,
        )

    ensemble = build_ensemble(factory, config.ensemble_size)
    histories = []
    print(f"  training {kind} ensemble (M={config.ensemble_size}) ...", flush=True)
    for i, member in enumerate(ensemble.members):
        print(f"  -- member {i+1}/{config.ensemble_size}", flush=True)
        info = train_one_model(member, tr_loader, va_loader, config, device, use_cond=use_cond)
        histories.append(info)
    return ensemble, {"members": histories}


def train_xgboost(
    train: IRISampleBatch,
    val: IRISampleBatch,
    config: SurrogateConfig,
    pre: IRIPreprocessor,
) -> tuple[XGBoostBaseline, dict]:
    Xtr, Ctr, Ytr = batch_to_tensors(train, pre)
    Xva, Cva, Yva = batch_to_tensors(val, pre)
    model = XGBoostBaseline(
        n_estimators=config.xgb_n_estimators,
        max_depth=config.xgb_max_depth,
        learning_rate=config.xgb_learning_rate,
        subsample=config.xgb_subsample,
        random_state=config.seed,
    )
    t0 = time.perf_counter()
    print("  training XGBoost baseline ...", flush=True)
    model.fit(Xtr, Ytr, cond=Ctr)
    dt = time.perf_counter() - t0
    pred = model.predict(Xva, cond=Cva)
    val_mse = float(np.mean((pred - Yva) ** 2))
    print(f"    val MSE (norm)={val_mse:.5f}  time={dt:.2f}s", flush=True)
    return model, {"val_mse_norm": val_mse, "train_seconds": dt}


def save_torch_model(model: nn.Module, path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), path)


def run_full_training_pipeline(
    train: IRISampleBatch,
    val: IRISampleBatch,
    config: SurrogateConfig,
    artifact_dir: Path | None = None,
) -> dict[str, Any]:
    """Train all models, write artifacts, return summary metadata."""
    artifact_dir = artifact_dir or config.artifact_path()
    artifact_dir = Path(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}  artifacts={artifact_dir}", flush=True)

    if len(train.doy) < 2:
        raise RuntimeError(
            f"Need at least 2 training samples, got {len(train.doy)}. "
            "IRI data generation likely failed (check IRI2020_DATA_DIR / subprocess worker)."
        )
    # If val is empty, carve a tiny holdout from train so metrics are defined.
    if len(val.doy) < 1:
        print("  warning: empty val set — holding out last train sample", flush=True)
        from .data import IRISampleBatch as _B

        def _take(b, sl):
            return _B(
                doy=b.doy[sl], hour=b.hour[sl], year=b.year[sl], month=b.month[sl], day=b.day[sl],
                glat=b.glat[sl], glon=b.glon[sl], alt_km=b.alt_km[sl], f107=b.f107[sl], ap=b.ap[sl],
                y=b.y[sl], targets=b.targets, regime=b.regime[sl],
            )

        val = _take(train, slice(-1, None))
        train = _take(train, slice(None, -1))

    pre = fit_preprocessor(train, config)
    pre.save(artifact_dir / "preprocessor.json")
    config.save(artifact_dir / "config.json")

    results: dict[str, Any] = {"device": device}

    res_model, _, res_info = train_residual_mlp(train, val, config, pre=pre, device=device)
    save_torch_model(res_model, artifact_dir / "residual_mlp.pt")
    results["residual_mlp"] = res_info

    film_model, _, film_info = train_film_mlp(train, val, config, pre=pre, device=device)
    save_torch_model(film_model, artifact_dir / "film_mlp.pt")
    results["film_mlp"] = film_info

    ens, ens_info = train_ensemble("film", train, val, config, pre=pre, device=device)
    save_ensemble(ens, artifact_dir / "film_ensemble.pt")
    results["film_ensemble"] = ens_info

    xgb, xgb_info = train_xgboost(train, val, config, pre=pre)
    xgb.save(artifact_dir / "xgboost.joblib")
    results["xgboost"] = xgb_info

    with open(artifact_dir / "train_summary.json", "w") as f:
        # histories only; avoid non-serializable tensors
        json.dump(_jsonable(results), f, indent=2)

    print("training pipeline complete.", flush=True)
    return results


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (float, int, str, bool)) or obj is None:
        return obj
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    return str(obj)
