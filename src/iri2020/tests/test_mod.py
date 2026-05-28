from pytest import approx
import pytest

import iri2020


def test_altitude_profile():
    time = "2015-12-13T10"
    altkmrange = (100, 1000, 10.0)
    glat = 65.1
    glon = -147.5

    iri = iri2020.IRI(time, altkmrange, glat, glon)

    # .item() necessary for stability across OS, pytest versions, etc.
    assert iri["ne"][10].item() == approx(24407842800.0, rel=0.001)
    assert iri.NmF2.item() == approx(77149454300.0, rel=0.001)
    assert iri.hmF2.item() == approx(265.249115, rel=0.001)
    assert iri.foF2.item() == approx(2.4943397, rel=0.001)


def test_input_validation_heistp():
    time = "2015-12-13T10"
    glat = 65.1
    glon = -147.5
    with pytest.raises(ValueError, match="heistp must be greater than 0.0"):
        iri2020.IRI(time, (100, 1000, 0.0), glat, glon)
    with pytest.raises(ValueError, match="heistp must be greater than 0.0"):
        iri2020.IRI(time, (100, 1000, -5.0), glat, glon)


def test_input_validation_num_alt():
    time = "2015-12-13T10"
    glat = 65.1
    glon = -147.5
    with pytest.raises(ValueError, match="num_alt must be less than or equal to 1000"):
        iri2020.IRI(time, (100, 1000, 0.5), glat, glon)


def test_output_shape():
    from iri2020.iri2020 import run_iri_py
    res = run_iri_py(2015, 12, 13, 10.0, 65.1, -147.5, [100.0, 200.0, 10.0])
    assert res["outf"].shape == (20, 11)

