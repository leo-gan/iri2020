"""Deep ensemble for epistemic uncertainty quantification.

Lakshminarayanan et al. (2017): train M networks with different random seeds /
initializations; predictive mean and variance approximate posterior uncertainty.

Caveats for this application
----------------------------
- Ensemble variance captures **model** uncertainty, not IRI's inherent empirical
  uncertainty or index-file noise.
- With short training (few epochs), members may be underfit and variances
  poorly calibrated — treat UQ as directional, not calibrated σ.
- Computational cost scales linearly with M at train and inference time.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import torch
import torch.nn as nn


class DeepEnsemble(nn.Module):
    def __init__(self, members: list[nn.Module]):
        super().__init__()
        self.members = nn.ModuleList(members)

    @property
    def n_members(self) -> int:
        return len(self.members)

    def forward(
        self, x: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        """Return mean prediction over ensemble members."""
        preds = [m(x, cond) for m in self.members]
        return torch.stack(preds, dim=0).mean(dim=0)

    def predict_with_uncertainty(
        self, x: torch.Tensor, cond: torch.Tensor | None = None
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Return (mean, std) over members; std is epistemic uncertainty proxy."""
        preds = torch.stack([m(x, cond) for m in self.members], dim=0)  # (M, B, T)
        mean = preds.mean(dim=0)
        std = preds.std(dim=0, unbiased=False)
        return mean, std

    def predict_numpy(
        self, x: np.ndarray, cond: np.ndarray | None = None, device: str = "cpu"
    ) -> tuple[np.ndarray, np.ndarray]:
        self.eval()
        xt = torch.from_numpy(np.asarray(x, dtype=np.float32)).to(device)
        ct = None
        if cond is not None:
            ct = torch.from_numpy(np.asarray(cond, dtype=np.float32)).to(device)
        with torch.no_grad():
            mean, std = self.predict_with_uncertainty(xt, ct)
        return mean.cpu().numpy(), std.cpu().numpy()


def build_ensemble(
    factory: Callable[[], nn.Module],
    n_members: int,
    seeds: list[int] | None = None,
) -> DeepEnsemble:
    if seeds is None:
        seeds = list(range(n_members))
    members = []
    for s in seeds[:n_members]:
        torch.manual_seed(s)
        np.random.seed(s)
        members.append(factory())
    return DeepEnsemble(members)


def save_ensemble(ensemble: DeepEnsemble, path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dicts": [m.state_dict() for m in ensemble.members]}, path)


def load_ensemble_into(
    ensemble: DeepEnsemble, path: Path | str, map_location: str = "cpu"
) -> DeepEnsemble:
    payload = torch.load(path, map_location=map_location, weights_only=True)
    state_dicts = payload["state_dicts"]
    if len(state_dicts) != ensemble.n_members:
        raise ValueError(
            f"ensemble size mismatch: checkpoint has {len(state_dicts)} members, "
            f"but ensemble was built with {ensemble.n_members}"
        )
    for m, sd in zip(ensemble.members, state_dicts, strict=True):
        m.load_state_dict(sd)
    return ensemble
