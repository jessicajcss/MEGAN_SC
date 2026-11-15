"""Speciation, emission-factor application and unit conversion — port of mgn2mech.F.

Turns the 19-class activity factors (ER, from megvea) plus the per-class emission
factors (EF, nanomol m-2 s-1) into:
  * 201 speciated compounds (native flux, per-cell rate, tonnes)
  * roll-ups to the 19 MEGAN classes
  * roll-ups to reporting families (isoprene, monoterpenes, sesquiterpenes, ...)

Unit chain (mgn2mech.F): flux[nanomol m-2 s-1] -> mol s-1 per cell via
1e-9 * cell_area; tonnes via molar weight * 1e-6 * seconds. Soil NO (class 8)
uses the GAMNO activity factor instead of the megvea ER, as in the Fortran.

Source: ``reference/src/MGN2MECH/mgn2mech.F`` + ``INCLDIR/SPC_NOCONVER.EXT``.
The 201-species table is auto-generated in ``speciation_data.py``.

NOTE — faithful source quirk: SPC_NOCONVER.EXT contains a duplicated map line
(cadinene_d, species 48) and consequently leaves species 67 (selinene_b)
unmapped, so it emits zero. mgn2mech.F's assignment semantics do the same; we
reproduce it exactly rather than "correcting" it. Two classes' speciation
fractions therefore do not sum to 1 in the source: OXPROD=0.995, OTHER=1.005.
"""
from __future__ import annotations

import numpy as np

from . import species as sp
from . import speciation_data as sd

SPCA_NAMES = list(sd.SPCA_NAMES)
SPCA_MWT = np.asarray(sd.SPCA_MWT, dtype=np.float64)      # g/mol, per compound
MAP_MG20 = np.asarray(sd.MAP_MG20, dtype=int)            # parent class (0-based)
MAP_SPCA = np.asarray(sd.MAP_SPCA, dtype=int)            # target compound (0-based)
MAP_EFFS = np.asarray(sd.MAP_EFFS, dtype=np.float64)     # speciation fraction
N_SPCA = int(sd.N_SPCA)

# unit constants (mgn2mech.F)
NANO2MOL = 1e-9
G2TONNE = 1e-6
HR2SEC = 3600.0

# Per-compound parent class (last-write-wins, matching Fortran); -1 = unmapped.
SPEC_CLASS = np.full(N_SPCA, -1, dtype=int)
for _s in range(MAP_SPCA.size):
    SPEC_CLASS[MAP_SPCA[_s]] = MAP_MG20[_s]

# Reporting families as lists of 0-based MEGAN class indices.
FAMILIES = {
    "isoprene": [0],
    "MBO": [1],
    "monoterpenes": [2, 3, 4, 5, 6],       # MT_PINE, MT_ACYC, MT_CAMP, MT_SABI, MT_AROM
    "sesquiterpenes": [8, 9],              # SQT_HR, SQT_LR
    "OVOC": [10, 11, 12, 13, 15],          # MEOH, ACTO, ETOH, ACID, OXPROD
    "other_VOC": [14, 16, 17],             # LVOC, STRESS, OTHER
    "NO": [7],                             # soil nitric oxide (non-VOC)
    "CO": [18],                            # carbon monoxide (non-VOC)
}
# All VOC classes = everything except NO (7) and CO (18)
VOC_CLASSES = [c for c in range(sp.NCLASS) if c not in (7, 18)]


def speciate_to_compounds(er, ef, gamno=None):
    """Distribute the 19 class activity factors over 201 compounds.

    ``er``, ``ef`` have shape (19, ...) — ER is the megvea activity factor and EF
    the emission-factor map (nanomol m-2 s-1). ``gamno`` (soil-NO activity, per
    cell) replaces ER for the NO class if given. Returns compound flux array of
    shape (201, ...) in nanomol m-2 s-1. Assignment order matches mgn2mech.F,
    so unmapped compounds stay 0.
    """
    er = np.asarray(er, dtype=np.float64)
    ef = np.asarray(ef, dtype=np.float64)
    grid = er.shape[1:]
    flux = np.zeros((N_SPCA,) + grid, dtype=np.float64)
    for s in range(MAP_SPCA.size):
        c = MAP_MG20[s]
        if c == sp.IDX_NO and gamno is not None:
            src = np.asarray(gamno, dtype=np.float64)
        else:
            src = er[c]
        flux[MAP_SPCA[s]] = src * ef[c] * MAP_EFFS[s]
    return flux


def compounds_to_classes(compound_arr):
    """Roll a (201, ...) compound array up to (19, ...) by parent class."""
    compound_arr = np.asarray(compound_arr, dtype=np.float64)
    grid = compound_arr.shape[1:]
    out = np.zeros((sp.NCLASS,) + grid, dtype=np.float64)
    for c in range(sp.NCLASS):
        members = np.nonzero(SPEC_CLASS == c)[0]
        if members.size:
            out[c] = compound_arr[members].sum(axis=0)
    return out


def classes_to_families(class_arr):
    """Roll a (19, ...) class array up to named families (dict of arrays)."""
    class_arr = np.asarray(class_arr, dtype=np.float64)
    fam = {name: class_arr[idx].sum(axis=0) for name, idx in FAMILIES.items()}
    fam["total_BVOC"] = class_arr[VOC_CLASSES].sum(axis=0)
    return fam


def flux_to_mol_s(flux_nanomol_m2_s, cell_area_m2):
    """Native per-cell emission rate (mol s-1) from flux and geodesic cell area."""
    return np.asarray(flux_nanomol_m2_s, dtype=np.float64) * NANO2MOL * np.asarray(cell_area_m2, dtype=np.float64)


def compound_tonnes(mol_s, seconds):
    """Compound mass (tonnes) over ``seconds`` from per-compound mol s-1.

    ``mol_s`` shape (201, ...). Uses each compound's molar weight (SPCA_MWT).
    """
    mol_s = np.asarray(mol_s, dtype=np.float64)
    mwt = SPCA_MWT.reshape((N_SPCA,) + (1,) * (mol_s.ndim - 1))
    return mol_s * mwt * G2TONNE * float(seconds)
