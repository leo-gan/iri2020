"""FiLM-conditioned MLP for driver-aware ionosphere prediction.

FiLM (Feature-wise Linear Modulation; Perez et al., 2018) predicts per-channel
scale (γ) and shift (β) from a conditioning vector and applies:
    h' = γ ⊙ h + β

Why FiLM here
-------------
Solar flux (F10.7) and geomagnetic activity (ap) modulate the *entire* ionosphere
profile amplitude/shape, not just a local additive effect. Conditioning the main
geometry/time stream via FiLM lets the network keep a shared spatial backbone
while specializing behavior under different driver regimes — important for
extreme-driver evaluation where a monolithic MLP tends to regress to the mean.

Conditioning inputs (after preprocessing): normalized f107, ap, sin/cos(doy).
Main inputs: periodic geometry/time + Fourier features.
"""

from __future__ import annotations

import torch
import torch.nn as nn


class FiLMGenerator(nn.Module):
    def __init__(self, cond_dim: int, hidden: int, n_layers: int):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(cond_dim, hidden),
            nn.GELU(),
            nn.Linear(hidden, hidden),
            nn.GELU(),
        )
        # Predict (gamma, beta) for each block's channel dim
        self.to_gb = nn.ModuleList(
            [nn.Linear(hidden, 2 * hidden) for _ in range(n_layers)]
        )

    def forward(self, cond: torch.Tensor) -> list[tuple[torch.Tensor, torch.Tensor]]:
        h = self.mlp(cond)
        out = []
        for layer in self.to_gb:
            gb = layer(h)
            gamma, beta = gb.chunk(2, dim=-1)
            # gamma centered at 1 for identity-at-init friendliness
            out.append((gamma + 1.0, beta))
        return out


class FiLMBlock(nn.Module):
    def __init__(self, dim: int, dropout: float = 0.0):
        super().__init__()
        self.fc1 = nn.Linear(dim, dim)
        self.fc2 = nn.Linear(dim, dim)
        self.norm = nn.LayerNorm(dim)
        self.drop = nn.Dropout(dropout)
        self.act = nn.GELU()

    def forward(
        self, x: torch.Tensor, gamma: torch.Tensor, beta: torch.Tensor
    ) -> torch.Tensor:
        h = self.fc1(x)
        h = self.act(h)
        h = self.drop(h)
        h = self.fc2(h)
        h = self.norm(h)
        h = gamma * h + beta
        return self.act(x + h)


class FiLMConditionedMLP(nn.Module):
    def __init__(
        self,
        in_dim: int,
        cond_dim: int,
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
        self.film_gen = FiLMGenerator(cond_dim, hidden, n_blocks)
        self.blocks = nn.ModuleList(
            [FiLMBlock(hidden, dropout=dropout) for _ in range(n_blocks)]
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
        if cond is None:
            raise ValueError("FiLMConditionedMLP requires cond tensor")
        h = self.stem(x)
        films = self.film_gen(cond)
        for blk, (gamma, beta) in zip(self.blocks, films):
            h = blk(h, gamma, beta)
        return self.head(h)
