"""Tests for the emission-factor processor.

Composition math is checked exactly with a tiny synthetic table; the real MEGAN
CSVs + the clipped SC ecotype grid are used for a shape/finiteness smoke test.
"""
import os
import numpy as np
import pandas as pd
import pytest

from megan_sc import efproc

DATA = os.path.join(os.path.dirname(__file__), "..", "data", "megan_inputs", "EFP240120")
EF_CSV = os.path.join(DATA, "EFv240120.csv")
SC_ECO = os.path.join(os.path.dirname(__file__), "..", "data", "sc_ecotype.nc")


def test_compose_ef_is_fraction_weighted_average():
    # Two species in one ecotype: A (0.7, EF1=10), B (0.3, EF1=20) -> EF1 = 13
    ef = pd.DataFrame({c: [0.0, 0.0] for c in efproc.EF_COLS}, index=["A", "B"])
    ef.loc["A", "EF1"] = 10.0
    ef.loc["B", "EF1"] = 20.0
    spec = pd.DataFrame({"EcoTypeID": [1, 1], "VegID": ["A", "B"], "frac": [0.7, 0.3]})
    out = efproc.compose_ef(ef, spec, "frac")
    assert out.loc[1, "EF1"] == pytest.approx(13.0)


@pytest.mark.skipif(not os.path.exists(EF_CSV), reason="EF data not present")
def test_load_ef_table_real():
    ef = efproc.load_ef_table(EF_CSV)
    assert all(c in ef.columns for c in efproc.EF_COLS)
    assert ef.index.is_unique
    assert np.all(np.isfinite(ef[efproc.EF_COLS].to_numpy()))


@pytest.mark.skipif(not (os.path.exists(EF_CSV) and os.path.exists(SC_ECO)),
                    reason="EF/ecotype data not present")
def test_sc_ecotype_to_ef_grid_smoke():
    from netCDF4 import Dataset
    ef = efproc.load_ef_table(EF_CSV)
    eco_tree = efproc.ef_by_ecotype(ef, DATA, "tree")
    assert list(eco_tree.columns) == efproc.EF_COLS
    assert np.all(eco_tree.to_numpy() >= 0.0)

    nc = Dataset(SC_ECO)
    grid = np.asarray(nc.variables["EVT"][:])
    nc.close()
    ef_grid = efproc.map_ecotype_ef(grid, eco_tree)
    assert ef_grid.shape == (19,) + grid.shape
    assert np.all(np.isfinite(ef_grid)) and np.all(ef_grid >= 0.0)


def test_assemble_ef_grid_weights_growth_forms():
    # synthetic: one ecotype (id 1), tree EF1=10, herb EF1=4 ; 60% tree / 40% herb
    tree = pd.DataFrame({c: [0.0] for c in efproc.EF_COLS}, index=pd.Index([1], name="EcoTypeID"))
    herb = tree.copy()
    tree.loc[1, "EF1"] = 10.0
    herb.loc[1, "EF1"] = 4.0
    grid = np.ones((2, 2), dtype=int)
    gf = {"tree": np.full((2, 2), 0.6), "herb": np.full((2, 2), 0.4)}
    out = efproc.assemble_ef_grid(grid, gf, {"tree": tree, "herb": herb})
    assert out.shape == (19, 2, 2)
    assert np.allclose(out[0], 0.6 * 10.0 + 0.4 * 4.0)   # EF1 == 8.4
