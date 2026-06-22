#!/usr/bin/env python3
"""End-to-end short pipeline: generate data, train, benchmark.

Usage (from repo root, venv active):
    python -m iri2020.surrogate.scripts.run_pipeline
    python -m iri2020.surrogate.scripts.run_pipeline --n-train 500 --epochs 5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure src layout works when run as script
_REPO = Path(__file__).resolve().parents[4]
if str(_REPO / "src") not in sys.path:
    sys.path.insert(0, str(_REPO / "src"))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="IRI2020 surrogate short training + benchmark"
    )
    p.add_argument("--artifact-dir", default="surrogate_artifacts")
    p.add_argument("--n-train", type=int, default=600)
    p.add_argument("--n-val", type=int, default=120)
    p.add_argument("--n-test", type=int, default=120)
    p.add_argument("--n-extreme", type=int, default=80)
    p.add_argument("--epochs", type=int, default=6)
    p.add_argument("--ensemble-size", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--skip-generate", action="store_true")
    p.add_argument("--skip-train", action="store_true")
    p.add_argument("--skip-benchmark", action="store_true")
    args = p.parse_args(argv)

    from iri2020.surrogate.config import SurrogateConfig
    from iri2020.surrogate.data import (
        generate_samples,
        concat_batches,
        save_batch,
        load_batch,
        train_val_test_split,
        ensure_iri_data_dir,
    )
    from iri2020.surrogate.train import run_full_training_pipeline
    from iri2020.surrogate.benchmark import run_benchmark

    data_dir = ensure_iri_data_dir()
    print(f"IRI2020_DATA_DIR={data_dir}", flush=True)

    cfg = SurrogateConfig(
        n_train=args.n_train,
        n_val=args.n_val,
        n_test=args.n_test,
        n_extreme=args.n_extreme,
        epochs=args.epochs,
        ensemble_size=args.ensemble_size,
        seed=args.seed,
        artifact_dir=args.artifact_dir,
    )
    art = Path(args.artifact_dir)
    art.mkdir(parents=True, exist_ok=True)

    train_p = art / "train_batch.npz"
    val_p = art / "val_batch.npz"
    test_p = art / "test_batch.npz"
    ext_p = art / "extreme_batch.npz"

    if not args.skip_generate and not train_p.exists():
        print("=== generating training data (Rust IRI) ===", flush=True)
        # Slightly oversample then split
        n_nom = args.n_train + args.n_val + args.n_test
        nom = generate_samples(n_nom, cfg, extreme=False)
        ext = generate_samples(args.n_extreme, cfg, extreme=True)
        # Mix some extreme into training distribution
        n_ext_tr = min(len(ext.doy) // 2, max(20, args.n_extreme // 3))
        from iri2020.surrogate.data import IRISampleBatch
        import numpy as np

        ext_tr_idx = np.arange(n_ext_tr)
        ext_ev_idx = np.arange(n_ext_tr, len(ext.doy))

        def take(b, ii):
            return IRISampleBatch(
                doy=b.doy[ii],
                hour=b.hour[ii],
                year=b.year[ii],
                month=b.month[ii],
                day=b.day[ii],
                glat=b.glat[ii],
                glon=b.glon[ii],
                alt_km=b.alt_km[ii],
                f107=b.f107[ii],
                ap=b.ap[ii],
                y=b.y[ii],
                targets=b.targets,
                regime=b.regime[ii],
            )

        full_nom = nom
        tr, va, te = train_val_test_split(
            full_nom, args.n_val, args.n_test, seed=args.seed
        )
        if n_ext_tr > 0:
            tr = concat_batches(tr, take(ext, ext_tr_idx))
        ext_eval = take(ext, ext_ev_idx) if len(ext_ev_idx) else ext

        save_batch(tr, train_p)
        save_batch(va, val_p)
        save_batch(te, test_p)
        save_batch(ext_eval, ext_p)
        print(
            f"saved batches: train={len(tr.doy)} val={len(va.doy)} test={len(te.doy)} extreme={len(ext_eval.doy)}",
            flush=True,
        )
    else:
        print("=== loading existing batches ===", flush=True)
        tr = load_batch(train_p)
        va = load_batch(val_p)
        te = load_batch(test_p)
        ext_eval = load_batch(ext_p) if ext_p.exists() else te

    if not args.skip_train:
        print("=== training (short run) ===", flush=True)
        run_full_training_pipeline(tr, va, cfg, artifact_dir=art)

    if not args.skip_benchmark:
        print("=== benchmark ===", flush=True)
        run_benchmark(
            art, test_batch=te, extreme_batch=ext_eval, config=cfg, n_speed_samples=24
        )

    print("done.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
