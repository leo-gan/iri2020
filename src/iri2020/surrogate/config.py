"""Configuration for IRI2020 neural surrogates."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import json


# Outputs that span many orders of magnitude -> log10 transform before scaling.
LOG_TARGETS: tuple[str, ...] = (
    "ne",
    "nO+",
    "nH+",
    "nHe+",
    "nO2+",
    "nNO+",
    "nCI",
    "nN+",
    "NmF2",
    "NmF1",
    "NmE",
    "TEC",
)

# Linear targets (temperatures, heights, drifts, frequencies).
LINEAR_TARGETS: tuple[str, ...] = (
    "Tn",
    "Ti",
    "Te",
    "hmF2",
    "hmF1",
    "hmE",
    "EqVertIonDrift",
    "foF2",
)

# Profile quantities vary with altitude; scalars are altitude-independent
# diagnostics (still predicted pointwise for a consistent interface).
PROFILE_TARGETS: tuple[str, ...] = (
    "ne",
    "Tn",
    "Ti",
    "Te",
    "nO+",
    "nH+",
    "nHe+",
    "nO2+",
    "nNO+",
    "nCI",
    "nN+",
)

SCALAR_TARGETS: tuple[str, ...] = (
    "NmF2",
    "hmF2",
    "NmF1",
    "hmF1",
    "NmE",
    "hmE",
    "TEC",
    "EqVertIonDrift",
    "foF2",
)

ALL_TARGETS: tuple[str, ...] = PROFILE_TARGETS + SCALAR_TARGETS

# Default subset for fast training / benchmarking (still multi-scale).
DEFAULT_TARGETS: tuple[str, ...] = ("ne", "Tn", "Ti", "Te", "NmF2", "hmF2", "TEC", "foF2")

# Periodic input dimensions (radians after encoding).
PERIODIC_INPUTS: dict[str, float] = {
    "doy": 365.25,  # day of year period
    "hour": 24.0,  # UT hour period
    "glon": 360.0,  # longitude degrees
}


@dataclass
class SurrogateConfig:
    """Central hyperparameter / path configuration."""

    # Data
    targets: tuple[str, ...] = DEFAULT_TARGETS
    alt_km_min: float = 100.0
    alt_km_max: float = 600.0
    n_train: int = 4000
    n_val: int = 800
    n_test: int = 800
    n_extreme: int = 400
    seed: int = 42

    # Preprocessing
    fourier_features: int = 16  # Gaussian Fourier feature count (per spatial/time group)
    fourier_scale: float = 1.0
    log_floor: float = 1e-3  # floor for log targets (avoids -inf)

    # Residual Fourier MLP
    res_hidden: int = 128
    res_blocks: int = 3
    res_dropout: float = 0.05

    # FiLM MLP
    film_hidden: int = 128
    film_blocks: int = 3
    film_cond_dim: int = 4  # f107, ap, sin(doy), cos(doy) style cond vector size after prep
    film_dropout: float = 0.05

    # Training (intentionally short by default for smoke / convergence checks)
    batch_size: int = 256
    lr: float = 1e-3
    weight_decay: float = 1e-5
    epochs: int = 8
    patience: int = 4
    grad_clip: float = 1.0
    num_workers: int = 0

    # Ensemble
    ensemble_size: int = 3

    # XGBoost baseline
    xgb_n_estimators: int = 80
    xgb_max_depth: int = 6
    xgb_learning_rate: float = 0.1
    xgb_subsample: float = 0.9

    # Paths
    artifact_dir: str = "surrogate_artifacts"

    def artifact_path(self) -> Path:
        p = Path(self.artifact_dir)
        p.mkdir(parents=True, exist_ok=True)
        return p

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["targets"] = list(self.targets)
        return d

    def save(self, path: Path | str) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "SurrogateConfig":
        if "targets" in d:
            d = {**d, "targets": tuple(d["targets"])}
        known = {f.name for f in cls.__dataclass_fields__.values()}  # type: ignore[attr-defined]
        return cls(**{k: v for k, v in d.items() if k in known})


def default_config(**overrides: Any) -> SurrogateConfig:
    cfg = SurrogateConfig()
    for k, v in overrides.items():
        if not hasattr(cfg, k):
            raise ValueError(f"Unknown config key: {k}")
        setattr(cfg, k, v)
    return cfg
