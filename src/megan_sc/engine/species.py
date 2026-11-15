"""MEGANv3.21 emission-class parameters — transcribed verbatim from MEGVEA.EXT.

Source: ``reference/src/MEGVEA/MEGVEA.EXT`` (Guenther & Huang, updated 2021).
NCLASS = 19. Arrays are 0-indexed here (Fortran class ``S`` == Python index
``S-1``). Do NOT edit values without re-checking against the Fortran header;
these are the ground truth the whole model depends on.

All arrays are float64 numpy arrays of length 19 (except NAMES/DESCR).
"""
from __future__ import annotations

import numpy as np

NCLASS: int = 19

# Class short names (MGN_SPC) and human descriptions, in Fortran order 1..19.
NAMES: tuple[str, ...] = (
    "ISOP", "MBO", "MT_PINE", "MT_ACYC", "MT_CAMP", "MT_SABI", "MT_AROM",
    "NO", "SQT_HR", "SQT_LR", "MEOH", "ACTO", "ETOH", "ACID", "LVOC",
    "OXPROD", "STRESS", "OTHER", "CO",
)
DESCR: tuple[str, ...] = (
    "isoprene",
    "2-methyl-3-buten-2-ol (MBO)",
    "monoterpenes: pines (alpha & beta pinene)",
    "monoterpenes: acyclic (e.g. myrcene, ocimenes)",
    "monoterpenes: carene, camphene, others",
    "monoterpenes: sabinene, limonene, terpinenes, others",
    "C10 aromatics (cymenes) & C8-C13 oxygenated (e.g. camphor)",
    "soil nitric oxide (NO)",
    "highly reactive sesquiterpenes (e.g. caryophyllene)",
    "less reactive sesquiterpenes (e.g. longifolene, copaene) + salates",
    "methanol",
    "acetone",
    "acetaldehyde and ethanol",
    "organic acids: formic, acetic, pyruvic",
    "C2-C4 hydrocarbons (e.g. ethene, ethane)",
    "oxidation products: aldehydes",
    "stress compounds (e.g. linalool)",
    "other VOC (e.g. indole, pentane, methyl bromide)",
    "carbon monoxide (CO)",
)

# Light-dependent fraction. NOTE: for classes 3,4,5,6 (Fortran S) the model reads
# LDF from a gridded LDFILE instead of this constant (megvea.f L651-666); these
# constants are the fallback/default only.
LDF = np.array(
    [1.0, 1.0, 0.6, 0.9, 0.2, 0.4, 0.1, 0.0, 0.5, 0.4,
     0.8, 0.2, 0.8, 0.8, 0.2, 0.2, 0.8, 0.1, 1.0], dtype=np.float64
)
LDF_FROM_MAP = (2, 3, 4, 5)  # 0-based indices (Fortran S=3,4,5,6) read LDF from map

# Light-dependent temperature coefficient CT1
CT1 = np.array(
    [95, 95, 80, 80, 80, 80, 80, 80, 130, 130,
     60, 80, 95, 95, 80, 80, 80, 80, 60], dtype=np.float64
)
# Cleo: coefficient in Eopt of the light-dependent temperature response
CLEO = np.array(
    [2.0, 2.0, 1.83, 1.83, 1.83, 1.83, 1.83, 1.83, 2.37, 2.37,
     1.6, 1.83, 2.0, 2.0, 1.83, 1.83, 1.83, 1.83, 1.6], dtype=np.float64
)
# beta: light-independent temperature coefficient
BETA = np.array(
    [0.13, 0.13, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.17, 0.17,
     0.08, 0.10, 0.13, 0.13, 0.10, 0.10, 0.10, 0.10, 0.08], dtype=np.float64
)

# Leaf-age relative activities
ANEW = np.array(
    [0.05, 0.05, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 0.4, 0.4,
     3.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float64
)
AGRO = np.array(
    [0.6, 0.6, 1.8, 1.8, 1.8, 1.8, 1.8, 1.0, 0.6, 0.6,
     3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float64
)
AMAT = np.ones(19, dtype=np.float64)  # all 1.0 in MEGVEA.EXT
AOLD = np.array(
    [0.9, 0.9, 1.05, 1.05, 1.05, 1.05, 1.05, 1.0, 0.95, 0.95,
     1.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float64
)

# Stress activity coefficients (applied only when the corresponding switch is ON;
# all default OFF in the run scripts).
CAQ = np.array([1, 1, 1, 5, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 1, 5, 1, 1], dtype=np.float64)
CHW = np.array([1, 1, 5, 5, 5, 5, 5, 1, 5, 5, 1, 1, 1, 1, 1, 1, 5, 1, 1], dtype=np.float64)
CHT = np.array([1, 1, 1, 5, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 1, 5, 1, 1], dtype=np.float64)
CLT = np.array([1, 1, 1, 5, 1, 1, 1, 1, 5, 5, 1, 1, 1, 1, 1, 1, 5, 1, 1], dtype=np.float64)
# Thresholds (uniform across classes in this version)
TAQ = np.full(19, 20.0);  THW = np.full(19, 12.0);  THT = np.full(19, 40.0);  TLT = np.full(19, 10.0)
DTAQ = np.full(19, 30.0); DTHW = np.full(19, 8.0);  DTHT = np.full(19, 8.0);  DTLT = np.full(19, 8.0)

# ---- scalar physical constants (MEGVEA.EXT) ----
CO2_STD = 400.0     # reference ambient CO2 (ppm); GAMCO2 == 1 at this value
ISMAX = 1.344       # CO2 inhibition parameter
CO2H = 1.4614       # CO2 inhibition exponent
CSTAR = 585.0       # CO2 inhibition half-saturation
PSTD = 200.0        # standard PPFD reference (unused in active GAMP; C1 fixed)
CCD1 = -0.2         # canopy-depth response slope
CCD2 = 1.3          # canopy-depth response intercept
CT2 = 230.0         # light-dependent temperature response denominator constant
TS_LI = 303.15      # reference temperature (K) for light-independent response
ALPHA_LIGHT = 0.004  # GAMP alpha
C1_LIGHT = 1.03      # GAMP C1 (fixed constant in this version)

# Leaf-age timestep length (days) between LAI periods, keyed by NLAI
TSTLEN = {46: 8, 37: 10, 12: 30}

# Canopy Gaussian-quadrature (5 layers). Nodes = fractional canopy depth; weights
# integrate over the canopy. For Layers != 5 the model uses uniform nodes/weights.
CDEPTH5 = np.array([0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899], dtype=np.float64)
VPGWT5 = np.array([0.1184635, 0.2393144, 0.284444444, 0.2393144, 0.1184635], dtype=np.float64)

# Indices of special-cased classes (0-based)
IDX_ISOP = 0     # isoprene: GAMCO2 applied
IDX_ETOH = 12    # acetaldehyde/ethanol: GAMBD applied
IDX_NO = 7       # soil NO
IDX_CO = 18      # carbon monoxide


def index_of(name: str) -> int:
    """Return the 0-based class index for a MGN_SPC name (case-insensitive)."""
    return NAMES.index(name.upper())


def canopy_nodes_weights(layers: int) -> tuple[np.ndarray, np.ndarray]:
    """Canopy depth nodes and quadrature weights (megvea.f GAMMA_CD & VPGWT).

    Layers == 5 uses the tabulated Gaussian nodes/weights; otherwise uniform.
    """
    if layers == 5:
        return CDEPTH5.copy(), VPGWT5.copy()
    k = np.arange(1, layers + 1, dtype=np.float64)
    cdepth = (k - 0.5) / layers
    vpgwt = np.full(layers, 1.0 / layers, dtype=np.float64)
    return cdepth, vpgwt
