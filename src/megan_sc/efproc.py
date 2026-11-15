"""Emission-factor processor — builds gridded per-class EFs from MEGAN's vegetation data.

Reproduces the composition logic of MEGAN's EF processor using the datasets you
provided:
  * ``EFv240120.csv``            : per-species emission factors EF1..EF19 + LDF3-6
  * ``Speciation{Tree,Shrub,Herb,Crop}.csv`` : ecotype -> species fraction per growth form
  * ``EVT3d.nc``                 : global ecotype grid (uint16), clipped to the domain
  * ``GFv7nc``                   : growth-form fractions per cell (weights the growth forms)
These are the files of ``EFP240120.zip``, ``EVT3d.zip`` and ``GFv7nc.zip`` of the MEGAN 3.2 input files (Guenther 2024,
Zenodo, https://doi.org/10.5281/zenodo.10939297).

Per grid cell:
    EF_class = Σ_growthform  GF_frac · ( Σ_species  specfrac · EF_class(species) )

UNIT NOTE: the absolute EF unit in this table must be confirmed against MEGAN
documentation / the bundled EFMAP before trusting magnitudes; the composition
here is unit-preserving. Magnitudes will be validated against the 2019 totals.

Source cross-reference: MEGAN3 EF processor + mgn2mech.F EF handling.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

EF_COLS = [f"EF{i}" for i in range(1, 20)]      # 19 MEGAN classes
LDF_COLS = ["LDF3", "LDF4", "LDF5", "LDF6"]      # LDF for MT_PINE/ACYC/CAMP/SABI

# growth form -> (speciation filename, fraction column name)
GROWTH_FORMS = {
    "tree":  ("SpeciationTree230921.csv",  "treeSpecfrac"),
    "shrub": ("SpeciationShrub230822.csv", "ShrubSpecFrac"),
    "herb":  ("SpeciationHerb230822.csv",  "HerbSpecFrac"),
    "crop":  ("SpeciationCrop230822.csv",  "CropSpecFrac"),
}


def load_ef_table(path):
    """Load per-species EF table, indexed by VegID (duplicates dropped)."""
    df = pd.read_csv(path, encoding="latin-1")
    df = df.drop_duplicates(subset="VegID").set_index("VegID")
    for c in EF_COLS + [c for c in LDF_COLS if c in df.columns]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df


def load_speciation(dirpath, growth_form):
    """Load an ecotype->species fraction table for one growth form."""
    fname, fraccol = GROWTH_FORMS[growth_form]
    df = pd.read_csv(os.path.join(dirpath, fname), usecols=["EcoTypeID", "VegID", fraccol],
                     encoding="latin-1")
    df[fraccol] = pd.to_numeric(df[fraccol], errors="coerce").fillna(0.0)
    return df, fraccol


def compose_ef(ef_df, spec_df, fraccol):
    """Per-ecotype EF = Σ_species specfrac·EF_species (pure; testable in isolation).

    ``ef_df`` indexed by VegID with EF1..19; ``spec_df`` has EcoTypeID, VegID and
    the fraction column ``fraccol``. Returns DataFrame indexed by EcoTypeID.
    """
    merged = spec_df.merge(ef_df[EF_COLS], left_on="VegID", right_index=True, how="left")
    merged[EF_COLS] = merged[EF_COLS].fillna(0.0)
    weighted = merged[EF_COLS].to_numpy() * merged[[fraccol]].to_numpy()
    out = pd.DataFrame(weighted, columns=EF_COLS)
    out["EcoTypeID"] = merged["EcoTypeID"].to_numpy()
    return out.groupby("EcoTypeID")[EF_COLS].sum()


def ef_by_ecotype(ef_df, dirpath, growth_form):
    """Per-ecotype EF (EF1..19) for one growth form, from the data files."""
    spec, fraccol = load_speciation(dirpath, growth_form)
    return compose_ef(ef_df, spec, fraccol)


def map_ecotype_ef(ecotype_grid, ef_by_eco):
    """Map a 2-D ecotype-ID grid to a (19, ny, nx) EF array via table lookup.

    Ecotype IDs absent from the table (e.g. no-vegetation) get 0 EF.
    """
    ecotype_grid = np.asarray(ecotype_grid)
    ny, nx = ecotype_grid.shape
    lut = {int(eid): ef_by_eco.loc[eid].to_numpy() for eid in ef_by_eco.index}
    out = np.zeros((19, ny, nx), dtype=np.float64)
    flat = ecotype_grid.ravel()
    ef_flat = out.reshape(19, -1)
    for eid in np.unique(flat):
        vec = lut.get(int(eid))
        if vec is not None:
            ef_flat[:, flat == eid] = vec[:, None]
    return out


def assemble_ef_grid(ecotype_grid, gf_fractions, ef_by_eco_by_gf):
    """Growth-form-weighted gridded EFs (19, ny, nx).

    ``gf_fractions``     : {growth_form: 2-D fraction grid (0-1)}
    ``ef_by_eco_by_gf``  : {growth_form: DataFrame from :func:`ef_by_ecotype`}
    """
    ny, nx = np.asarray(ecotype_grid).shape
    total = np.zeros((19, ny, nx), dtype=np.float64)
    for gf, frac in gf_fractions.items():
        ef = map_ecotype_ef(ecotype_grid, ef_by_eco_by_gf[gf])
        total += np.asarray(frac, dtype=np.float64)[None, :, :] * ef
    return total
