"""XGBoost multi-output baseline for IRI surrogates.

XGBoost is a strong tabular baseline: it handles nonlinear interactions well,
trains fast, and needs less tuning than deep nets on moderate datasets.

Limitations vs neural surrogates
--------------------------------
- No natural uncertainty quantification without quantile/extra tricks.
- Separate model per target (or MultiOutputRegressor wrapping) — no shared
  representation across ionospheric variables.
- Does not exploit FiLM-style driver conditioning as a first-class mechanism
  (drivers are just concatenated features).
"""

from __future__ import annotations

from pathlib import Path
import json

import numpy as np
from sklearn.multioutput import MultiOutputRegressor

try:
    from xgboost import XGBRegressor
except ImportError as e:  # pragma: no cover
    raise ImportError("xgboost is required for the baseline: pip install xgboost") from e


class XGBoostBaseline:
    def __init__(
        self,
        n_estimators: int = 80,
        max_depth: int = 6,
        learning_rate: float = 0.1,
        subsample: float = 0.9,
        random_state: int = 42,
        n_jobs: int = -1,
    ):
        base = XGBRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            subsample=subsample,
            random_state=random_state,
            n_jobs=n_jobs,
            objective="reg:squarederror",
            tree_method="hist",
        )
        self.model = MultiOutputRegressor(base)
        self._fitted = False

    def fit(self, X: np.ndarray, y: np.ndarray, cond: np.ndarray | None = None) -> "XGBoostBaseline":
        if cond is not None:
            X = np.concatenate([X, cond], axis=-1)
        self.model.fit(X, y)
        self._fitted = True
        return self

    def predict(self, X: np.ndarray, cond: np.ndarray | None = None) -> np.ndarray:
        if not self._fitted:
            raise RuntimeError("XGBoostBaseline not fitted")
        if cond is not None:
            X = np.concatenate([X, cond], axis=-1)
        return self.model.predict(X)

    def save(self, path: Path | str) -> None:
        import joblib

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self.model, path)

    def load(self, path: Path | str) -> "XGBoostBaseline":
        import joblib

        self.model = joblib.load(path)
        self._fitted = True
        return self
