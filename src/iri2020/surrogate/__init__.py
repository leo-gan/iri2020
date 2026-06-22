"""Neural network surrogates for the IRI2020 empirical ionosphere model.

This package provides:
- Preprocessing with normalization, periodic encodings, and log-scale outputs
- Deep Residual MLP with Fourier features
- FiLM-conditioned MLP for solar/geomagnetic driver conditioning
- Deep ensemble for epistemic uncertainty quantification
- XGBoost baseline
- Training, evaluation, and benchmark harnesses

Primary prediction targets are pointwise ionospheric quantities at
(time, lat, lon, alt) with optional solar/geomagnetic drivers derived
from the IRI run (f107, ap) or treated as external conditioning inputs.
"""

from .config import SurrogateConfig, default_config
from .preprocessing import IRIPreprocessor
from .models.residual_mlp import ResidualFourierMLP
from .models.film_mlp import FiLMConditionedMLP
from .models.ensemble import DeepEnsemble
from .models.xgboost_baseline import XGBoostBaseline

__all__ = [
    "SurrogateConfig",
    "default_config",
    "IRIPreprocessor",
    "ResidualFourierMLP",
    "FiLMConditionedMLP",
    "DeepEnsemble",
    "XGBoostBaseline",
]
