"""Tests for the IOAPI reference reader (io/ioapi.py) on a synthetic Lambert file."""
import sys, os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import make_synthetic_ioapi as mki          # noqa: E402
from megan_sc.io import ioapi                # noqa: E402


@pytest.fixture()
def ioapi_file(tmp_path):
    return mki.write(str(tmp_path / "MEGANv31.SC_2019.CB6X.2019015.BDSNP.nc"))


def test_describe_reads_species_and_units(ioapi_file):
    d = ioapi.describe(ioapi_file)
    assert set(d["species"]) == {"ISOP", "TERP"}
    assert d["units"]["ISOP"] == "moles/s"
    assert d["ntimes"] == 2
    assert d["times"][0] == (2019015, 0)


def test_cell_centers_land_in_santa_catarina(ioapi_file):
    d = ioapi.describe(ioapi_file)
    lon, lat = ioapi.cell_centers_lonlat(d["grid"])
    assert lon.shape == (int(d["grid"]["NROWS"]), int(d["grid"]["NCOLS"]))
    # SC lies roughly within these bounds — verifies the LCC projection is set up right
    assert -54.5 < lon.min() and lon.max() < -47.5
    assert -30.0 < lat.min() and lat.max() < -25.5


def test_read_species_sum_and_instant(ioapi_file):
    both = ioapi.read_species(ioapi_file)                       # sum over 2 steps
    assert set(both["_available"]) == {"ISOP", "TERP"}
    assert np.allclose(both["ISOP"], 10.0)                      # 5.0 * 2 timesteps
    one = ioapi.read_species(ioapi_file, names=["isop"], tstep=0)
    assert np.allclose(one["ISOP"], 5.0)                        # single step
    assert one["_units"]["ISOP"] == "moles/s"


def test_unknown_species_skipped(ioapi_file):
    out = ioapi.read_species(ioapi_file, names=["ISOP", "NOTREAL"], tstep=0)
    assert "ISOP" in out and "NOTREAL" not in out
    assert "NOTREAL" not in out["_available"]
