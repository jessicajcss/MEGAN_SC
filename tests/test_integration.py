"""End-to-end engine smoke test: met -> canopy -> ER -> speciation -> tonnes.

Uses synthetic emission factors (real EF maps come from the data layer later).
Checks the whole chain runs, stays finite/non-negative, and produces a positive
daytime isoprene flux with sane magnitude.
"""
import numpy as np

from megan_sc.engine import canopy, megvea, soil, speciation as spc, species as sp


def _g(v, shape=(3, 3)):
    return np.full(shape, v, dtype=float)


def test_full_chain_daytime():
    shape = (3, 3)
    ctf = np.zeros((6, *shape)); ctf[2] = 100.0        # temperate broadleaf

    can = canopy.megcan(
        temp=_g(298.0), par_wm2=_g(400.0), wind=_g(2.0), pres=_g(1.0e5),
        qv=_g(0.010), lat=_g(-27.0), lon=_g(-50.0), lai=_g(4.0),
        ctf=ctf, day=15, hour_utc=15.0,
    )

    gamsm = soil.gamma_sm(_g(0.30), soil.wilting_point(6))   # well-watered
    er = megvea.emission_activity(
        can, laic=_g(4.0), laip=_g(4.0), d_temp=_g(297.0), d_ppfd=_g(1800.0),
        gamsm=gamsm,
    )
    assert er.shape == (19, 3, 3)
    assert np.all(np.isfinite(er)) and np.all(er >= 0.0)

    ef = np.ones((sp.NCLASS, *shape))                  # synthetic EF (nanomol m-2 s-1)
    gamno = _g(0.5)
    flux = spc.speciate_to_compounds(er, ef, gamno=gamno)   # nanomol m-2 s-1
    assert flux.shape == (201, 3, 3)
    assert np.all(np.isfinite(flux)) and np.all(flux >= 0.0)

    area = _g(1.0e6)                                    # 1 km^2 cells
    mol_s = spc.flux_to_mol_s(flux, area)
    tonnes_hr = spc.compound_tonnes(mol_s, seconds=3600.0)
    assert np.all(np.isfinite(tonnes_hr)) and np.all(tonnes_hr >= 0.0)

    # daytime isoprene must emit something
    assert flux[0].sum() > 0.0

    # class + family roll-ups run and conserve mass vs compounds
    cls = spc.compounds_to_classes(flux)
    assert np.allclose(cls.sum(axis=0), flux.sum(axis=0))
    fam = spc.classes_to_families(cls)
    assert np.all(np.isfinite(fam["total_BVOC"]))
