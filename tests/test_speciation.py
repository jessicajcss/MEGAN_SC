"""Tests for the speciation, EF-application and unit-conversion port (mgn2mech.F)."""
import numpy as np
import pytest

from megan_sc.engine import species as sp
from megan_sc.engine import speciation as spc


def test_table_sizes():
    assert spc.N_SPCA == 201
    assert len(spc.SPCA_NAMES) == 201
    assert spc.SPCA_MWT.shape == (201,)
    assert spc.MAP_MG20.shape == spc.MAP_SPCA.shape == spc.MAP_EFFS.shape == (201,)


def test_isoprene_is_species_zero():
    assert spc.SPCA_NAMES[0] == "isoprene"
    assert spc.SPCA_MWT[0] == pytest.approx(68.12)


def test_source_quirk_unmapped_species_emits_zero():
    # selinene_b (species index 66, 1-based 67) is never mapped -> always 0.
    er = np.ones((sp.NCLASS, 2, 2))
    ef = np.ones((sp.NCLASS, 2, 2))
    flux = spc.speciate_to_compounds(er, ef)
    assert np.all(flux[66] == 0.0)
    # cadinene_d (species 47) receives its 0.005 fraction (duplicated line, same value)
    assert np.allclose(flux[47], 0.005)


def test_speciation_fraction_applied():
    # With ER=EF=1 everywhere, each compound flux == its speciation fraction.
    er = np.ones((sp.NCLASS, 3))
    ef = np.ones((sp.NCLASS, 3))
    flux = spc.speciate_to_compounds(er, ef)
    assert np.allclose(flux[0], 1.0)                    # isoprene: fraction 1.0
    # pinene_a (species 2) is 100% of MT_PINE
    assert np.allclose(flux[2], 1.0)


def test_no_uses_gamno_not_er():
    er = np.ones((sp.NCLASS, 2, 2))
    ef = np.ones((sp.NCLASS, 2, 2))
    gamno = np.full((2, 2), 0.3)
    flux_with = spc.speciate_to_compounds(er, ef, gamno=gamno)
    flux_without = spc.speciate_to_compounds(er, ef)
    # find an NO compound (parent class == IDX_NO)
    no_species = np.nonzero(spc.SPEC_CLASS == sp.IDX_NO)[0]
    assert no_species.size >= 1
    j = no_species[0]
    assert np.allclose(flux_with[j], 0.3 * spc.MAP_EFFS[np.nonzero(spc.MAP_SPCA == j)[0][0]])
    assert np.allclose(flux_without[j], 1.0 * spc.MAP_EFFS[np.nonzero(spc.MAP_SPCA == j)[0][0]])


def test_rollup_conserves_and_matches_isoprene():
    er = np.ones((sp.NCLASS, 4))
    ef = np.ones((sp.NCLASS, 4))
    flux = spc.speciate_to_compounds(er, ef)
    cls = spc.compounds_to_classes(flux)
    assert cls.shape == (sp.NCLASS, 4)
    assert np.allclose(cls[sp.index_of("ISOP")], 1.0)     # isoprene class == 1
    assert np.allclose(cls[sp.index_of("MBO")], 1.0)
    # total across compounds == total across classes (rollup conserves mass)
    assert np.allclose(flux.sum(axis=0), cls.sum(axis=0))


def test_families_group_classes():
    cls = np.ones((sp.NCLASS, 2, 2))
    fam = spc.classes_to_families(cls)
    assert np.allclose(fam["monoterpenes"], 5.0)          # 5 MT classes
    assert np.allclose(fam["sesquiterpenes"], 2.0)        # 2 SQT classes
    assert np.allclose(fam["isoprene"], 1.0)
    assert np.allclose(fam["total_BVOC"], float(len(spc.VOC_CLASSES)))


def test_unit_conversions():
    # 1 nanomol m-2 s-1 over a 1 km^2 (1e6 m^2) cell -> 1e-3 mol s-1
    mol_s = spc.flux_to_mol_s(np.array([1.0]), np.array([1.0e6]))
    assert mol_s[0] == pytest.approx(1.0e-3)
    # tonnes of isoprene (MW 68.12) over 1 hour from that rate
    t = spc.compound_tonnes(np.array([1.0e-3]), seconds=3600.0)
    assert t[0] == pytest.approx(1.0e-3 * 68.12 * 1e-6 * 3600.0)
