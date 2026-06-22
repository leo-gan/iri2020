"""Preprocessing pipeline for IRI2020 surrogate models.

Design goals
------------
1. **Periodic structure**: day-of-year, UT hour, and longitude are circular.
   We encode them as (sin, cos) pairs so the model never sees a discontinuity
   at the 0/360 or 24h boundary.
2. **Orders-of-magnitude outputs**: electron/ion densities and TEC vary by
   6+ orders of magnitude. We apply log10(x + floor) then standardize so
   MSE on normalized targets is not dominated by the largest absolute values.
3. **Linear outputs**: temperatures, heights, foF2, drifts are standardized
   only (z-score), preserving sign/scale where log is inappropriate.
4. **Fourier features**: optional Gaussian Fourier feature lift (Tancik et al.)
   applied to continuous coordinates (lat, alt, time-of-year embedding) to
   help an MLP resolve high-frequency spatial/temporal structure that pure
   MLPs are biased against (spectral bias).
5. **FiLM conditioning**: solar/geomagnetic drivers (f107, ap) and seasonal
   context are split into a separate "condition" vector for FiLM models,
   while geometry/time go into the main feature stream.

Trade-offs
----------
- Log + floor clips densities below `log_floor`, which is physically fine for
  IRI (zero-padded species regions) but slightly distorts the lowest densities.
- Global z-score ignores heteroscedasticity across altitude; altitude-conditional
  normalization would be more accurate but harder to invert and document.
- Fourier features add parameters and can overfit with very few samples; scale
  is modest by default.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import json

import numpy as np

from .config import (
    LOG_TARGETS,
    PERIODIC_INPUTS,
    SurrogateConfig,
)


def _sin_cos(values: np.ndarray, period: float) -> tuple[np.ndarray, np.ndarray]:
    ang = 2.0 * np.pi * values / period
    return np.sin(ang), np.cos(ang)


def encode_periodic_row(
    doy: float,
    hour: float,
    glat: float,
    glon: float,
    alt_km: float,
) -> np.ndarray:
    """Encode a single sample's geometry/time into a fixed feature vector.

    Order: sin_doy, cos_doy, sin_hour, cos_hour, glat/90, sin_glon, cos_glon,
    alt_norm (0-1 over [100,1000] km for stability), year_frac (unused here).
    """
    s_doy, c_doy = _sin_cos(np.asarray(doy, dtype=np.float64), PERIODIC_INPUTS["doy"])
    s_hr, c_hr = _sin_cos(np.asarray(hour, dtype=np.float64), PERIODIC_INPUTS["hour"])
    s_lon, c_lon = _sin_cos(np.asarray(glon, dtype=np.float64), PERIODIC_INPUTS["glon"])
    lat_n = np.clip(glat / 90.0, -1.2, 1.2)
    alt_n = (alt_km - 60.0) / (2000.0 - 60.0)  # broad physical range
    alt_n = np.clip(alt_n, -0.1, 1.2)
    return np.array(
        [
            float(s_doy),
            float(c_doy),
            float(s_hr),
            float(c_hr),
            float(lat_n),
            float(s_lon),
            float(c_lon),
            float(alt_n),
        ],
        dtype=np.float64,
    )


def encode_condition_row(f107: float, ap: float, doy: float) -> np.ndarray:
    """Conditioning vector for FiLM: solar/geomag + coarse season."""
    s_doy, c_doy = _sin_cos(np.asarray(doy, dtype=np.float64), PERIODIC_INPUTS["doy"])
    # Rough normalizations around climatological means
    f107_n = (f107 - 100.0) / 80.0
    ap_n = (ap - 10.0) / 30.0
    return np.array(
        [float(f107_n), float(ap_n), float(s_doy), float(c_doy)], dtype=np.float64
    )


def gaussian_fourier_features(
    x: np.ndarray,
    B: np.ndarray,
) -> np.ndarray:
    """Map x (N, D) through fixed random Fourier features: [cos(2π x B), sin(2π x B)].

    B has shape (D, M). Output shape (N, 2M).
    """
    proj = 2.0 * np.pi * x @ B  # (N, M)
    return np.concatenate([np.cos(proj), np.sin(proj)], axis=-1)


@dataclass
class StandardizerState:
    mean: np.ndarray
    std: np.ndarray

    def transform(self, x: np.ndarray) -> np.ndarray:
        return (x - self.mean) / self.std

    def inverse(self, z: np.ndarray) -> np.ndarray:
        return z * self.std + self.mean


class IRIPreprocessor:
    """Fit/transform pipeline for features and multi-target outputs."""

    def __init__(self, config: SurrogateConfig | None = None):
        self.config = config or SurrogateConfig()
        self.targets: list[str] = list(self.config.targets)
        self.log_mask: np.ndarray | None = None
        self.y_state: StandardizerState | None = None
        self.x_state: StandardizerState | None = None
        self.cond_state: StandardizerState | None = None
        self.fourier_B: np.ndarray | None = None
        self._fitted = False

    # ------------------------------------------------------------------
    # Target transforms
    # ------------------------------------------------------------------
    def _apply_log_targets(self, y: np.ndarray) -> np.ndarray:
        """y: (N, T) in physical units -> log10 where applicable."""
        assert self.log_mask is not None
        out = y.copy().astype(np.float64)
        floor = self.config.log_floor
        for j, is_log in enumerate(self.log_mask):
            if is_log:
                out[:, j] = np.log10(np.maximum(out[:, j], floor))
        return out

    def _invert_log_targets(self, y_loglike: np.ndarray) -> np.ndarray:
        assert self.log_mask is not None
        out = y_loglike.copy().astype(np.float64)
        for j, is_log in enumerate(self.log_mask):
            if is_log:
                out[:, j] = 10.0 ** out[:, j]
        return out

    def build_raw_features(
        self,
        doy: np.ndarray,
        hour: np.ndarray,
        glat: np.ndarray,
        glon: np.ndarray,
        alt_km: np.ndarray,
    ) -> np.ndarray:
        """Stack periodic encodings; shape (N, F0)."""
        n = len(doy)
        feats = np.zeros((n, 8), dtype=np.float64)
        for i in range(n):
            feats[i] = encode_periodic_row(
                float(doy[i]),
                float(hour[i]),
                float(glat[i]),
                float(glon[i]),
                float(alt_km[i]),
            )
        return feats

    def build_raw_conditions(
        self, f107: np.ndarray, ap: np.ndarray, doy: np.ndarray
    ) -> np.ndarray:
        n = len(doy)
        cond = np.zeros((n, 4), dtype=np.float64)
        for i in range(n):
            cond[i] = encode_condition_row(float(f107[i]), float(ap[i]), float(doy[i]))
        return cond

    def fit(
        self,
        doy: np.ndarray,
        hour: np.ndarray,
        glat: np.ndarray,
        glon: np.ndarray,
        alt_km: np.ndarray,
        f107: np.ndarray,
        ap: np.ndarray,
        y: np.ndarray,
        targets: Sequence[str] | None = None,
    ) -> "IRIPreprocessor":
        if targets is not None:
            self.targets = list(targets)
        self.log_mask = np.array([t in LOG_TARGETS for t in self.targets], dtype=bool)

        x_raw = self.build_raw_features(doy, hour, glat, glon, alt_km)
        cond_raw = self.build_raw_conditions(f107, ap, doy)

        rng = np.random.default_rng(self.config.seed)
        d_in = x_raw.shape[1]
        m = self.config.fourier_features
        self.fourier_B = rng.normal(
            0.0, self.config.fourier_scale, size=(d_in, m)
        ).astype(np.float64)

        x_ff = gaussian_fourier_features(x_raw, self.fourier_B)
        x_full = np.concatenate([x_raw, x_ff], axis=-1)

        x_mean = x_full.mean(axis=0)
        x_std = x_full.std(axis=0)
        x_std = np.where(x_std < 1e-8, 1.0, x_std)
        self.x_state = StandardizerState(x_mean, x_std)

        c_mean = cond_raw.mean(axis=0)
        c_std = cond_raw.std(axis=0)
        c_std = np.where(c_std < 1e-8, 1.0, c_std)
        self.cond_state = StandardizerState(c_mean, c_std)

        y_t = self._apply_log_targets(y)
        y_mean = y_t.mean(axis=0)
        y_std = y_t.std(axis=0)
        y_std = np.where(y_std < 1e-8, 1.0, y_std)
        self.y_state = StandardizerState(y_mean, y_std)

        self._fitted = True
        return self

    def transform_X(
        self,
        doy: np.ndarray,
        hour: np.ndarray,
        glat: np.ndarray,
        glon: np.ndarray,
        alt_km: np.ndarray,
    ) -> np.ndarray:
        self._check_fitted()
        assert self.fourier_B is not None and self.x_state is not None
        x_raw = self.build_raw_features(doy, hour, glat, glon, alt_km)
        x_ff = gaussian_fourier_features(x_raw, self.fourier_B)
        x_full = np.concatenate([x_raw, x_ff], axis=-1)
        return self.x_state.transform(x_full).astype(np.float32)

    def transform_cond(
        self, f107: np.ndarray, ap: np.ndarray, doy: np.ndarray
    ) -> np.ndarray:
        self._check_fitted()
        assert self.cond_state is not None
        cond_raw = self.build_raw_conditions(f107, ap, doy)
        return self.cond_state.transform(cond_raw).astype(np.float32)

    def transform_y(self, y: np.ndarray) -> np.ndarray:
        self._check_fitted()
        assert self.y_state is not None
        y_t = self._apply_log_targets(y)
        return self.y_state.transform(y_t).astype(np.float32)

    def inverse_y(self, y_norm: np.ndarray) -> np.ndarray:
        """Normalized predictions -> physical units."""
        self._check_fitted()
        assert self.y_state is not None
        y_t = self.y_state.inverse(np.asarray(y_norm, dtype=np.float64))
        return self._invert_log_targets(y_t)

    def input_dim(self) -> int:
        self._check_fitted()
        assert self.x_state is not None
        return int(self.x_state.mean.shape[0])

    def cond_dim(self) -> int:
        self._check_fitted()
        assert self.cond_state is not None
        return int(self.cond_state.mean.shape[0])

    def output_dim(self) -> int:
        return len(self.targets)

    def save(self, path: Path | str) -> None:
        self._check_fitted()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "targets": self.targets,
            "log_mask": self.log_mask.tolist() if self.log_mask is not None else None,
            "x_mean": self.x_state.mean.tolist(),  # type: ignore[union-attr]
            "x_std": self.x_state.std.tolist(),  # type: ignore[union-attr]
            "cond_mean": self.cond_state.mean.tolist(),  # type: ignore[union-attr]
            "cond_std": self.cond_state.std.tolist(),  # type: ignore[union-attr]
            "y_mean": self.y_state.mean.tolist(),  # type: ignore[union-attr]
            "y_std": self.y_state.std.tolist(),  # type: ignore[union-attr]
            "fourier_B": self.fourier_B.tolist()
            if self.fourier_B is not None
            else None,
            "config": self.config.to_dict(),
        }
        with open(path, "w") as f:
            json.dump(payload, f)

    @classmethod
    def load(cls, path: Path | str) -> "IRIPreprocessor":
        with open(path) as f:
            payload = json.load(f)
        cfg = SurrogateConfig.from_dict(payload.get("config", {}))
        obj = cls(cfg)
        obj.targets = list(payload["targets"])
        obj.log_mask = np.array(payload["log_mask"], dtype=bool)
        obj.x_state = StandardizerState(
            np.array(payload["x_mean"], dtype=np.float64),
            np.array(payload["x_std"], dtype=np.float64),
        )
        obj.cond_state = StandardizerState(
            np.array(payload["cond_mean"], dtype=np.float64),
            np.array(payload["cond_std"], dtype=np.float64),
        )
        obj.y_state = StandardizerState(
            np.array(payload["y_mean"], dtype=np.float64),
            np.array(payload["y_std"], dtype=np.float64),
        )
        obj.fourier_B = np.array(payload["fourier_B"], dtype=np.float64)
        obj._fitted = True
        return obj

    def _check_fitted(self) -> None:
        if not self._fitted:
            raise RuntimeError("IRIPreprocessor is not fitted")
