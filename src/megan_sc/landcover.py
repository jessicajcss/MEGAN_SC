"""MapBiomas land-cover -> MEGAN vegetation structure (a documented, coherent crosswalk).

MapBiomas Brazil (30 m, annual) is far more SC-specific than the global 1 km
ecotype or 500 m MODIS land cover. We use it *coherently* — without displacing
MEGAN's calibrated science — as follows:

  * MapBiomas class  ->  growth-form fractions (tree/shrub/herb/crop) per cell
                         (replaces the GFv7nc growth-form file; better for Brazil)
  * growth-form      ->  MEGAN canopy-type fractions for megcan
  * vegetated?       ->  a 30 m vegetation mask (water/urban/bare -> 0 emission)

What MapBiomas does NOT replace, to protect data quality:
  * the ecotype (EVT) -> species -> emission-factor composition stays MEGAN's,
    now weighted by MapBiomas growth-form fractions instead of GFv7nc.
  * LAI stays the observed MODIS field (land cover is not a leaf-area measurement);
    MapBiomas is used to mask/sharpen it, not to invent LAI magnitudes.

Legend codes: MapBiomas Brazil Collection 9/10 (stable across recent collections).
Assumptions are called out in ``CROSSWALK_NOTES`` and are user-overridable.

Sources:
  Earth Engine asset ``projects/mapbiomas-public/assets/brazil/lulc/v1``;
  MapBiomas legend codes (brasil.mapbiomas.org).
"""
from __future__ import annotations

import numpy as np

# Growth-form order used throughout.
GROWTH_FORMS = ("tree", "shrub", "herb", "crop")
_GF_IDX = {g: i for i, g in enumerate(GROWTH_FORMS)}

# MEGAN canopy types (megcan.f / MEGCAN.EXT), 1-based:
#   1 needleleaf, 2 tropical broadleaf, 3 temperate broadleaf, 4 shrub, 5 herb, 6 crop
# SC forests are subtropical Atlantic broadleaf -> temperate broadleaf (3) by default.
TREE_CANOPY_TYPE_DEFAULT = 3

# MapBiomas class code -> growth-form weights (tree, shrub, herb, crop); rows sum to
# <=1, remainder is non-vegetated. Classes absent here are treated as non-veg (0).
CROSSWALK = {
    1:  (1.0, 0.0, 0.0, 0.0),   # Forest (generic)
    3:  (1.0, 0.0, 0.0, 0.0),   # Forest Formation
    4:  (0.5, 0.5, 0.0, 0.0),   # Savanna Formation (tree/shrub mix)  [ASSUMPTION]
    5:  (1.0, 0.0, 0.0, 0.0),   # Mangrove
    6:  (1.0, 0.0, 0.0, 0.0),   # Floodable Forest
    49: (1.0, 0.0, 0.0, 0.0),   # Wooded Sandbank / restinga arbórea
    9:  (1.0, 0.0, 0.0, 0.0),   # Forest Plantation (silviculture)     [broadleaf assumed]
    10: (0.0, 0.0, 1.0, 0.0),   # Non-Forest Natural Formation (parent)
    11: (0.0, 0.0, 1.0, 0.0),   # Wetland
    12: (0.0, 0.0, 1.0, 0.0),   # Grassland
    13: (0.0, 0.0, 1.0, 0.0),   # Other Non-Forest Natural
    50: (0.0, 0.0, 1.0, 0.0),   # Herbaceous Sandbank / restinga herbácea
    15: (0.0, 0.0, 1.0, 0.0),   # Pasture
    14: (0.0, 0.0, 0.5, 0.5),   # Farming (parent, mixed)              [ASSUMPTION]
    21: (0.0, 0.0, 0.5, 0.5),   # Mosaic of Uses (pasture+agri)        [ASSUMPTION]
    18: (0.0, 0.0, 0.0, 1.0),   # Agriculture (parent)
    19: (0.0, 0.0, 0.0, 1.0),   # Temporary Crop
    39: (0.0, 0.0, 0.0, 1.0),   # Soybean
    20: (0.0, 0.0, 0.0, 1.0),   # Sugar cane
    40: (0.0, 0.0, 0.0, 1.0),   # Rice
    62: (0.0, 0.0, 0.0, 1.0),   # Cotton
    41: (0.0, 0.0, 0.0, 1.0),   # Other Temporary Crops
    36: (0.0, 0.0, 0.0, 1.0),   # Perennial Crop
    46: (0.0, 0.0, 0.0, 1.0),   # Coffee
    47: (0.0, 0.0, 0.0, 1.0),   # Citrus
    35: (0.0, 0.0, 0.0, 1.0),   # Palm oil                             [treated as crop]
    48: (0.0, 0.0, 0.0, 1.0),   # Other Perennial Crops
    # --- non-vegetated / water / no-data (all zeros => masked) ---
    22: (0.0, 0.0, 0.0, 0.0),   # Non-vegetated area (parent)
    23: (0.0, 0.0, 0.0, 0.0),   # Beach/Dune/Sand
    24: (0.0, 0.0, 0.0, 0.0),   # Urban
    25: (0.0, 0.0, 0.0, 0.0),   # Other Non-vegetated
    30: (0.0, 0.0, 0.0, 0.0),   # Mining
    29: (0.0, 0.0, 0.0, 0.0),   # Rocky Outcrop
    32: (0.0, 0.0, 0.0, 0.0),   # Hypersaline Tidal Flat
    26: (0.0, 0.0, 0.0, 0.0),   # Water (parent)
    33: (0.0, 0.0, 0.0, 0.0),   # River/Lake/Ocean
    31: (0.0, 0.0, 0.0, 0.0),   # Aquaculture
    27: (0.0, 0.0, 0.0, 0.0),   # Not Observed
}

CROSSWALK_NOTES = (
    "Savanna (4) split 50/50 tree/shrub; Mosaic (21) and Farming (14) split 50/50 "
    "herb/crop; Silviculture (9) treated as broadleaf tree; Palm oil (35) as crop. "
    "Tree canopy type defaults to temperate broadleaf (3) for subtropical SC. "
    "All assumptions are overridable."
)

# --- MODIS MCD12Q1 'LC_Type5' (Plant Functional Type) -> growth-form weights ---
# Codes: 1 evergreen needleleaf, 2 evergreen broadleaf, 3 deciduous needleleaf,
# 4 deciduous broadleaf, 5 shrub, 6 grass, 7 cereal crop, 8 broadleaf crop trees;
# 0 water, 9 urban, 10 snow, 11 barren -> masked (all-zero). This band IS a PFT map,
# a direct fit for MEGAN's growth forms; at 500 m it aggregates cheaply to the grid.
MODIS_PFT_CROSSWALK = {
    1: (1.0, 0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0, 0.0),
    3: (1.0, 0.0, 0.0, 0.0), 4: (1.0, 0.0, 0.0, 0.0),
    5: (0.0, 1.0, 0.0, 0.0), 6: (0.0, 0.0, 1.0, 0.0),
    7: (0.0, 0.0, 0.0, 1.0), 8: (0.0, 0.0, 0.0, 1.0),
}

# --- WRF 'IVGTYP' (MODIFIED_IGBP_MODIS_NOAH, 21 classes) -> growth-form weights ---
# This is the land-use classification carried *inside* the WRF/MCIP files (attr MMINLU),
# so using it makes our run share the reference MEGAN run's exact land surface. Water(17),
# lake(21), urban(13), snow/ice(15), barren(16), tundra-barren(20) -> masked (all-zero).
IGBP_MODIS_CROSSWALK = {
    1:  (1.0, 0.0, 0.0, 0.0),   # Evergreen Needleleaf Forest -> tree
    2:  (1.0, 0.0, 0.0, 0.0),   # Evergreen Broadleaf Forest  -> tree
    3:  (1.0, 0.0, 0.0, 0.0),   # Deciduous Needleleaf Forest -> tree
    4:  (1.0, 0.0, 0.0, 0.0),   # Deciduous Broadleaf Forest  -> tree
    5:  (1.0, 0.0, 0.0, 0.0),   # Mixed Forests               -> tree
    6:  (0.0, 1.0, 0.0, 0.0),   # Closed Shrublands           -> shrub
    7:  (0.0, 1.0, 0.0, 0.0),   # Open Shrublands             -> shrub
    8:  (0.5, 0.5, 0.0, 0.0),   # Woody Savannas   [tree/shrub, ASSUMPTION]
    9:  (0.0, 0.5, 0.5, 0.0),   # Savannas         [shrub/herb, ASSUMPTION]
    10: (0.0, 0.0, 1.0, 0.0),   # Grasslands                  -> herb
    11: (0.0, 0.0, 1.0, 0.0),   # Permanent Wetlands          -> herb
    12: (0.0, 0.0, 0.0, 1.0),   # Croplands                   -> crop
    13: (0.0, 0.0, 0.0, 0.0),   # Urban and Built-Up          -> masked
    14: (0.0, 0.0, 0.5, 0.5),   # Cropland/Natural Mosaic [herb/crop, ASSUMPTION]
    15: (0.0, 0.0, 0.0, 0.0),   # Snow and Ice                -> masked
    16: (0.0, 0.0, 0.0, 0.0),   # Barren/Sparsely Vegetated   -> masked
    17: (0.0, 0.0, 0.0, 0.0),   # Water                       -> masked
    18: (0.5, 0.5, 0.0, 0.0),   # Wooded Tundra    [n/a in SC]
    19: (0.0, 0.5, 0.5, 0.0),   # Mixed Tundra     [n/a in SC]
    20: (0.0, 0.0, 0.0, 0.0),   # Barren Tundra               -> masked
    21: (0.0, 0.0, 0.0, 0.0),   # Lake                        -> masked
}

# --- SC-CORRECTED IGBP crosswalk (biogeographically valid for Santa Catarina) ---
# Santa Catarina is 100% Atlantic Forest biome (Mata Atlântica) — there is NO Cerrado/savanna
# in the state. The global IGBP-MODIS product nonetheless labels ~20% of SC (the north-central
# Planalto and western interior) as "Woody Savanna"/"Savanna". Evidence that this is a
# misclassification of FOREST, not real savanna: those pixels carry forest-level leaf area in
# the WRF file itself (mean LAI ≈ 3.8; 75% > 3.5), indistinguishable from the confirmed
# Evergreen Broadleaf Forest (LAI ≈ 4.9); true Cerrado/savanna is LAI ≈ 1-2. Geographically
# they coincide with the Araucária Forest (Floresta Ombrófila Mista) of the northern plateau
# and the seasonal deciduous forest of the west. So for SC we reclassify:
#   8 Woody Savanna -> tree (Araucária/seasonal forest)
#   9 Savanna       -> mostly tree, small campos (natural-grassland) component
# This is a documented, region-specific correction. The FULLY valid fix is to use MapBiomas
# (30 m, Brazil-specific), which carries ~0 Savanna Formation in SC — see growth_form_source:
# mapbiomas. Everything else is inherited from IGBP_MODIS_CROSSWALK.
IGBP_MODIS_CROSSWALK_SC = dict(IGBP_MODIS_CROSSWALK)
IGBP_MODIS_CROSSWALK_SC[8] = (1.0, 0.0, 0.0, 0.0)   # Woody Savanna -> tree
IGBP_MODIS_CROSSWALK_SC[9] = (0.75, 0.0, 0.25, 0.0)  # Savanna -> 75% tree / 25% herb (campos)


def _gf_lut(crosswalk):
    """Lookup table (max_code+1, 4) of growth-form weights indexed by class code."""
    max_code = max(crosswalk)
    lut = np.zeros((max_code + 1, 4), dtype=np.float64)
    for code, w in crosswalk.items():
        if 0 <= code <= max_code:
            lut[code] = w
    return lut, max_code


def growth_form_weights(class_grid, crosswalk=None):
    """Map a land-cover class grid to per-pixel growth-form weights.

    ``crosswalk`` selects the legend: ``CROSSWALK`` (MapBiomas, default) or
    ``MODIS_PFT_CROSSWALK``. Returns array of shape (..., 4) in GROWTH_FORMS order;
    unknown/no-data codes map to all-zero (non-vegetated).
    """
    crosswalk = CROSSWALK if crosswalk is None else crosswalk
    lut, max_code = _gf_lut(crosswalk)
    cg = np.asarray(class_grid)
    cg = np.where((cg >= 0) & (cg <= max_code), cg, 0).astype(int)
    return lut[cg]


def fractions_on_grid(class_grid, flon, flat, tlon, tlat, res, crosswalk=None):
    """Aggregate a fine land-cover class grid to growth-form fractions on the analysis grid.

    Returns dict with 4 arrays (tree/shrub/herb/crop), each (n_tlat, n_tlon) giving
    the fractional coverage of that growth form in the cell, plus 'veg' = their sum
    (the vegetated fraction / MEGAN TotalCT). ``crosswalk`` selects the legend.
    """
    from .domain import regrid_block_mean

    w = growth_form_weights(class_grid, crosswalk)           # (ny_f, nx_f, 4)
    out = {}
    for g, i in _GF_IDX.items():
        frac = regrid_block_mean(w[..., i], flon, flat, tlon, tlat, res)
        out[g] = np.nan_to_num(frac, nan=0.0)
    out["veg"] = sum(out[g] for g in GROWTH_FORMS)
    return out


def canopy_type_fractions(gf_fractions, tree_canopy_type=TREE_CANOPY_TYPE_DEFAULT):
    """MEGAN 6-type canopy fractions (PERCENT, 0-100) from growth-form fractions.

    tree -> ``tree_canopy_type`` (default 3, temperate broadleaf); shrub->4;
    herb->5; crop->6. Output matches megcan's CTF convention (percent; TotalCT =
    sum*0.01). Shape (6, ny, nx).
    """
    ny, nx = gf_fractions["tree"].shape
    ctf = np.zeros((6, ny, nx), dtype=np.float64)
    ctf[tree_canopy_type - 1] += gf_fractions["tree"] * 100.0
    ctf[3] += gf_fractions["shrub"] * 100.0   # type 4 shrub
    ctf[4] += gf_fractions["herb"] * 100.0    # type 5 herb
    ctf[5] += gf_fractions["crop"] * 100.0    # type 6 crop
    return ctf
