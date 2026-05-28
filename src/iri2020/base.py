from __future__ import annotations
from dateutil.parser import parse
from datetime import datetime
import xarray

from .iri2020 import run_iri_py  # type: ignore

SIMOUT = ["ne", "Tn", "Ti", "Te", "nO+", "nH+", "nHe+", "nO2+", "nNO+", "nCI", "nN+"]

__all__ = ["IRI"]


def IRI(time: str | datetime, altkmrange: list[float], glat: float, glon: float) -> xarray.Dataset:
    if isinstance(time, str):
        time = parse(time)

    assert len(altkmrange) == 3, "altitude (km) min, max, step"
    assert isinstance(glat, (int, float)) and isinstance(
        glon, (int, float)
    ), "glat, glon is scalar"

    dhour = time.hour + time.minute / 60.0 + time.second / 3600.0

    ret = run_iri_py(
        time.year,
        time.month,
        time.day,
        dhour,
        float(glat),
        float(glon),
        [float(a) for a in altkmrange],
    )

    altkm = ret["altkm"]
    outf = ret["outf"]
    oarr = ret["oarr"]

    dsf = {k: (("alt_km"), v) for (k, v) in zip(SIMOUT, outf[:11, :])}

    iono = xarray.Dataset(
        dsf,
        coords={"time": [time], "alt_km": altkm, "glat": glat, "glon": glon},
        attrs={"f107": oarr[40], "ap": oarr[51]},
    )

    for i, p in enumerate(["NmF2", "hmF2", "NmF1", "hmF1", "NmE", "hmE"]):
        iono[p] = (("time"), [oarr[i]])

    iono["TEC"] = (("time"), [oarr[36]])
    iono["EqVertIonDrift"] = (("time"), [oarr[43]])
    iono["foF2"] = (("time"), [oarr[99]])

    return iono
