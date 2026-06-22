"""Subprocess worker for IRI evaluations (single or batch)."""
from __future__ import annotations
import argparse
import json
import os
from datetime import datetime
from pathlib import Path

def _ensure_data_dir() -> None:
    if os.environ.get("IRI2020_DATA_DIR"):
        return
    here = Path(__file__).resolve()
    for c in [here.parents[2] / "data", here.parents[3] / "src" / "data",
              Path.cwd() / "src" / "data", Path.cwd() / "data"]:
        if c.is_dir() and (c / "ig_rz.dat").exists():
            os.environ["IRI2020_DATA_DIR"] = str(c.resolve())
            return

def _eval_one(req: dict) -> dict:
    from iri2020.base import IRI
    import numpy as np
    t = datetime(int(req["year"]), int(req["month"]), int(req["day"]),
                 int(req["hour"]), int(req.get("minute", 0)), 0)
    a = float(req["alt_km"])
    glat = float(req["glat"])
    glon = float(req["glon"])
    targets = list(req["targets"])
    try:
        ds = IRI(t, [a, a, 1.0], glat, glon)
        y = []
        for name in targets:
            if name in ds and "alt_km" in ds[name].dims:
                y.append(float(ds[name].values.ravel()[0]))
            else:
                y.append(float(np.asarray(ds[name].values).ravel()[0]))
        return {"ok": True, "f107": float(ds.attrs.get("f107", 100.0)),
                "ap": float(ds.attrs.get("ap", 10.0)), "y": y, "meta": req.get("meta")}
    except BaseException as e:
        return {"ok": False, "error": repr(e), "meta": req.get("meta")}

def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--req", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    _ensure_data_dir()
    with open(args.req) as f:
        req = json.load(f)
    if isinstance(req, dict) and "batch" in req:
        payload = {"ok": True, "results": [_eval_one(item) for item in req["batch"]]}
    else:
        payload = _eval_one(req)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(payload, f)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
