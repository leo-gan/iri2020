"""Deep Residual MLP operating on Fourier-augmented, normalized features.

Why residual blocks
-------------------
Plain deep MLPs suffer from vanishing gradients and optimization difficulty.
Residual connections (He et al., 2016) let us stack more layers while keeping
gradients flowing — useful because ionospheric fields have multi-scale structure
(smooth climatology + sharper gradients near terminators / peaks).

Why Fourier features live in preprocessing
------------------------------------------
We lift coordinates *before* the network (fixed B matrix) rather than learning
positional encodings end-to-end. This is simpler, matches the Tancik et al.
recipe for spectral bias mitigation, and keeps the PyTorch module a pure
function of already-featurized inputs.
"""

from __future__ import annotations

import torch
import torch.nn as nn


class ResidualBlock(nn.Module):
    def __init__(self, dim: int, dropout: float = 0.0):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, dim),
            nn.LayerNorm(dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim, dim),
            nn.LayerNorm(dim),
        )
        self.act = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.act(x + self.net(x))


class ResidualFourierMLP(nn.Module):
    """Residual MLP: stem -> N residual blocks -> linear head."""

    def __init__(
        self,
        in_dim: int,
        out_dim: int,
        hidden: int = 128,
        n_blocks: int = 3,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.stem = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.LayerNorm(hidden),
            nn.GELU(),
        )
        self.blocks = nn.ModuleList(
            [ResidualBlock(hidden, dropout=dropout) for _ in range(n_blocks)]
        )
        self.head = nn.Linear(hidden, out_dim)
        self._init_weights()

    def _init_weights(self) -> None:
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(
        self, x: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        # cond ignored (interface compatibility with FiLM models)
        h = self.stem(x)
        for blk in self.blocks:
            h = blk(h)
        return self.head(h)
