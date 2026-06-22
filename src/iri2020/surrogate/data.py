"""Dataset generation by querying the IRI2020 (Rust) backend.

Sampling strategy
-----------------
- **Nominal domain**: Latin-hypercube-like independent uniforms over lat/lon/alt,
  random times in a multi-year window, reflecting climatological coverage.
- **Extreme-driver regimes**: oversample high/low F10.7 proxies via dates known
  for high solar activity / storm seasons and polar latitudes. Since F10.7/ap
  are read from IRI internal index files (not user inputs in the default
  Python API), we approximate extremes by sampling dates/locations where IRI
  returns elevated f107/ap attributes, plus polar/equatorial extremes.

We store both inputs and IRI outputs so training is fully offline after generation.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterator
import json

import numpy as np

from iri2020.base import IRI

from .config import SurrogateConfig, DEFAULT_TARGETS

import os
import json
import subprocess
import sys


def ensure_iri_data_dir() -> str | None:
    """Point Rust backend at bundled coefficient/index files."""
    if os.environ.get("IRI2020_DATA_DIR"):
        return os.environ["IRI2020_DATA_DIR"]
    here = Path(__file__).resolve()
    candidates = [
        here.parents[2] / "data",
        here.parents[3] / "src" / "data",
        Path.cwd() / "src" / "data",
        Path.cwd() / "data",
    ]
    for c in candidates:
        if c.is_dir() and (c / "ig_rz.dat").exists():
            os.environ["IRI2020_DATA_DIR"] = str(c.resolve())
            return os.environ["IRI2020_DATA_DIR"]
    return None


@dataclass
class IRISampleBatch:
    doy: np.ndarray
    hour: np.ndarray
    year: np.ndarray
    month: np.ndarray
    day: np.ndarray
    glat: np.ndarray
    glon: np.ndarray
    alt_km: np.ndarray
    f107: np.ndarray
    ap: np.ndarray
    y: np.ndarray  # (N, T)
    targets: list[str]
    regime: np.ndarray  # 0=nominal, 1=extreme


def _sample_times(rng: np.random.Generator, n: int, extreme: bool = False) -> list[datetime]:
    times: list[datetime] = []
    # Years with index data coverage in bundled apf107/ig_rz files.
    years = list(range(2000, 2020))
    for _ in range(n):
        if extreme:
            # Favor solar max years and storm-prone seasons
            year = int(rng.choice([2001, 2002, 2003, 2012, 2014, 2015]))
            month = int(rng.choice([3, 4, 9, 10, 11]))  # equinox/storm seasons
        else:
            year = int(rng.choice(years))
            month = int(rng.integers(1, 13))
        day = int(rng.integers(1, 28))
        hour = int(rng.integers(0, 24))
        minute = int(rng.integers(0, 60))
        times.append(datetime(year, month, day, hour, minute, 0))
    return times


def _sample_geo(rng: np.random.Generator, n: int, extreme: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if extreme:
        # Polar + equatorial extremes, full altitude band
        lat_mode = rng.integers(0, 3, size=n)
        glat = np.where(
            lat_mode == 0,
            rng.uniform(-90, -60, size=n),
            np.where(lat_mode == 1, rng.uniform(60, 90, size=n), rng.uniform(-15, 15, size=n)),
        )
        glon = rng.uniform(-180, 180, size=n)
        alt = rng.uniform(100, 800, size=n)
    else:
        glat = rng.uniform(-80, 80, size=n)
        glon = rng.uniform(-180, 180, size=n)
        alt = rng.uniform(100, 600, size=n)
    return glat.astype(np.float64), glon.astype(np.float64), alt.astype(np.float64)


def _extract_targets(ds, alt_km: float, targets: list[str]) -> np.ndarray:
    """Extract target vector from an IRI xarray Dataset at a single altitude."""
    row = []
    for t in targets:
        if t in ds and "alt_km" in ds[t].dims:
            # nearest altitude in the returned profile (we query single alt)
            vals = ds[t].values
            row.append(float(vals.ravel()[0]))
        elif t in ds:
            row.append(float(np.asarray(ds[t].values).ravel()[0]))
        else:
            row.append(np.nan)
    return np.array(row, dtype=np.float64)


def _safe_iri_call_inline(t: datetime, a: float, glat: float, glon: float, targets: list[str]):
    """In-process IRI call (fast, but a Rust panic can poison the extension mutex)."""
    ensure_iri_data_dir()
    try:
        ds = IRI(t, [a, a, 1.0], float(glat), float(glon))
        f107 = float(ds.attrs.get("f107", 100.0))
        ap = float(ds.attrs.get("ap", 10.0))
        y_row = _extract_targets(ds, a, targets)
        if not (np.isfinite(f107) and np.isfinite(ap) and np.isfinite(y_row).all()):
            return None
        return f107, ap, y_row
    except BaseException:
        return None


def _safe_iri_call_subprocess(t: datetime, a: float, glat: float, glon: float, targets: list[str]):
    """Subprocess IRI call — isolates panics; preferred for bulk dataset generation."""
    import tempfile

    ensure_iri_data_dir()
    req = {
        "year": t.year,
        "month": t.month,
        "day": t.day,
        "hour": t.hour,
        "minute": t.minute,
        "alt_km": a,
        "glat": glat,
        "glon": glon,
        "targets": targets,
    }
    try:
        with tempfile.TemporaryDirectory() as td:
            req_p = Path(td) / "req.json"
            out_p = Path(td) / "out.json"
            req_p.write_text(json.dumps(req))
            proc = subprocess.run(
                [sys.executable, "-m", "iri2020.surrogate.iri_worker",
                 "--req", str(req_p), "--out", str(out_p)],
                capture_output=True,
                text=True,
                timeout=90,
                env=os.environ.copy(),
            )
            if not out_p.exists():
                return None
            payload = json.loads(out_p.read_text())
        if not payload.get("ok"):
            return None
        y_row = np.asarray(payload["y"], dtype=np.float64)
        return float(payload["f107"]), float(payload["ap"]), y_row
    except (subprocess.SubprocessError, json.JSONDecodeError, KeyError, ValueError, OSError):
        return None


def _safe_iri_call(
    t: datetime,
    a: float,
    glat: float,
    glon: float,
    targets: list[str],
    *,
    use_subprocess: bool = True,
):
    """Call IRI; return (f107, ap, y_row) or None on failure."""
    if use_subprocess:
        return _safe_iri_call_subprocess(t, a, glat, glon, targets)
    return _safe_iri_call_inline(t, a, glat, glon, targets)


def generate_samples(
    n: int,
    config: SurrogateConfig | None = None,
    extreme: bool = False,
    progress: bool = True,
    max_attempts_factor: int = 4,
    batch_size: int = 40,
) -> IRISampleBatch:
    """Query IRI n successful times via batched subprocess workers."""
    import tempfile

    config = config or SurrogateConfig()
    targets = list(config.targets)
    rng = np.random.default_rng(config.seed + (1000 if extreme else 0))
    ensure_iri_data_dir()

    rows: list[dict] = []
    failed = 0
    attempts = 0
    max_attempts = max(n * max_attempts_factor, n + 50)

    while len(rows) < n and attempts < max_attempts:
        need = min(batch_size, (n - len(rows)) * 2 + 4)
        batch_reqs = []
        for _ in range(need):
            attempts += 1
            if attempts > max_attempts:
                break
            times = _sample_times(rng, 1, extreme=extreme)
            glat_a, glon_a, alt_a = _sample_geo(rng, 1, extreme=extreme)
            t = times[0]
            batch_reqs.append({
                "year": t.year, "month": t.month, "day": t.day,
                "hour": t.hour, "minute": t.minute,
                "alt_km": float(alt_a[0]), "glat": float(glat_a[0]), "glon": float(glon_a[0]),
                "targets": targets,
                "meta": {
                    "doy": t.timetuple().tm_yday + t.hour / 24.0,
                    "hour": t.hour + t.minute / 60.0,
                    "year": t.year, "month": t.month, "day": t.day,
                    "glat": float(glat_a[0]), "glon": float(glon_a[0]),
                    "alt_km": float(alt_a[0]),
                },
            })

        try:
            with tempfile.TemporaryDirectory() as td:
                req_p = Path(td) / "req.json"
                out_p = Path(td) / "out.json"
                req_p.write_text(json.dumps({"batch": batch_reqs}))
                subprocess.run(
                    [sys.executable, "-m", "iri2020.surrogate.iri_worker",
                     "--req", str(req_p), "--out", str(out_p)],
                    capture_output=True, text=True, timeout=600, env=os.environ.copy(),
                )
                if not out_p.exists():
                    failed += len(batch_reqs)
                    continue
                payload = json.loads(out_p.read_text())
            results = payload.get("results", [])
            for item, res in zip(batch_reqs, results):
                if res.get("ok"):
                    meta = item["meta"]
                    y_row = np.asarray(res["y"], dtype=np.float64)
                    if not np.isfinite(y_row).all():
                        failed += 1
                        continue
                    rows.append({**meta, "f107": float(res["f107"]), "ap": float(res["ap"]), "y": y_row})
                    if len(rows) >= n:
                        break
                else:
                    failed += 1
        except (subprocess.SubprocessError, json.JSONDecodeError, OSError):
            failed += len(batch_reqs)

        if progress:
            print(f"  generated {min(len(rows), n)}/{n} ({'extreme' if extreme else 'nominal'}; failed={failed})", flush=True)

    rows = rows[:n]
    if len(rows) < n:
        print(f"  warning: only collected {len(rows)}/{n} samples (failed={failed})", flush=True)
    if not rows:
        return IRISampleBatch(
            doy=np.zeros(0), hour=np.zeros(0), year=np.zeros(0, dtype=np.int32),
            month=np.zeros(0, dtype=np.int32), day=np.zeros(0, dtype=np.int32),
            glat=np.zeros(0), glon=np.zeros(0), alt_km=np.zeros(0),
            f107=np.zeros(0), ap=np.zeros(0), y=np.zeros((0, len(targets))),
            targets=targets, regime=np.zeros(0, dtype=np.int8),
        )
    regime = np.full(len(rows), 1 if extreme else 0, dtype=np.int8)
    return IRISampleBatch(
        doy=np.array([r["doy"] for r in rows], dtype=np.float64),
        hour=np.array([r["hour"] for r in rows], dtype=np.float64),
        year=np.array([r["year"] for r in rows], dtype=np.int32),
        month=np.array([r["month"] for r in rows], dtype=np.int32),
        day=np.array([r["day"] for r in rows], dtype=np.int32),
        glat=np.array([r["glat"] for r in rows], dtype=np.float64),
        glon=np.array([r["glon"] for r in rows], dtype=np.float64),
        alt_km=np.array([r["alt_km"] for r in rows], dtype=np.float64),
        f107=np.array([r["f107"] for r in rows], dtype=np.float64),
        ap=np.array([r["ap"] for r in rows], dtype=np.float64),
        y=np.stack([r["y"] for r in rows], axis=0),
        targets=targets, regime=regime,
    )



def concat_batches(*batches: IRISampleBatch) -> IRISampleBatch:
    assert batches
    targets = batches[0].targets
    return IRISampleBatch(
        doy=np.concatenate([b.doy for b in batches]),
        hour=np.concatenate([b.hour for b in batches]),
        year=np.concatenate([b.year for b in batches]),
        month=np.concatenate([b.month for b in batches]),
        day=np.concatenate([b.day for b in batches]),
        glat=np.concatenate([b.glat for b in batches]),
        glon=np.concatenate([b.glon for b in batches]),
        alt_km=np.concatenate([b.alt_km for b in batches]),
        f107=np.concatenate([b.f107 for b in batches]),
        ap=np.concatenate([b.ap for b in batches]),
        y=np.concatenate([b.y for b in batches], axis=0),
        targets=targets,
        regime=np.concatenate([b.regime for b in batches]),
    )


def save_batch(batch: IRISampleBatch, path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        doy=batch.doy,
        hour=batch.hour,
        year=batch.year,
        month=batch.month,
        day=batch.day,
        glat=batch.glat,
        glon=batch.glon,
        alt_km=batch.alt_km,
        f107=batch.f107,
        ap=batch.ap,
        y=batch.y,
        regime=batch.regime,
        targets=np.array(batch.targets),
    )


def load_batch(path: Path | str) -> IRISampleBatch:
    z = np.load(path, allow_pickle=False)
    targets = [str(t) for t in z["targets"].tolist()]
    return IRISampleBatch(
        doy=z["doy"],
        hour=z["hour"],
        year=z["year"],
        month=z["month"],
        day=z["day"],
        glat=z["glat"],
        glon=z["glon"],
        alt_km=z["alt_km"],
        f107=z["f107"],
        ap=z["ap"],
        y=z["y"],
        targets=targets,
        regime=z["regime"],
    )


def train_val_test_split(
    batch: IRISampleBatch,
    n_val: int,
    n_test: int,
    seed: int = 42,
) -> tuple[IRISampleBatch, IRISampleBatch, IRISampleBatch]:
    n = len(batch.doy)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    n_test = min(n_test, n // 5)
    n_val = min(n_val, n // 5)
    te = idx[:n_test]
    va = idx[n_test : n_test + n_val]
    tr = idx[n_test + n_val :]

    def take(ii: np.ndarray) -> IRISampleBatch:
        return IRISampleBatch(
            doy=batch.doy[ii],
            hour=batch.hour[ii],
            year=batch.year[ii],
            month=batch.month[ii],
            day=batch.day[ii],
            glat=batch.glat[ii],
            glon=batch.glon[ii],
            alt_km=batch.alt_km[ii],
            f107=batch.f107[ii],
            ap=batch.ap[ii],
            y=batch.y[ii],
            targets=batch.targets,
            regime=batch.regime[ii],
        )

    return take(tr), take(va), take(te)
