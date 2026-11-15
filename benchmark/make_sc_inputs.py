"""Build MEGANv3.21 IOAPI inputs on the Santa Catarina analysis grid (0.025 deg, 144 x 224, south-first)
for one day, from the 2019 inputs of a run folder (MODIS LAI with Earth Engine default resampling, MapBiomas growth
forms, ERA5-Land monthly means).
The diurnal cycle is idealised (T: cosine, amplitude 5 K, max 15 h local; PAR: clear-sky cos(zenith)
shape scaled to the monthly-mean PAR).  These are TEST inputs for a code-equivalence benchmark, not weather.
usage: make_sc_inputs.py <engine src> <megan_sc_outputs dir> <MEGANv3.21 dir> <outdir> <month> <day>"""
import sys, os, datetime, numpy as np, rasterio
SRC, RAW, MGN, OUT, MON, DAY = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
sys.path.insert(0, SRC); sys.path.insert(0, os.path.dirname(__file__))
from megan_sc import landcover, fixes
from ioapi import write_ioapi
os.makedirs(OUT, exist_ok=True)
YEAR = 2019; RES = 0.025
doy = datetime.date(YEAR, MON, DAY).timetuple().tm_yday; SDATE = YEAR * 1000 + doy
def tif(path, band=None):
    with rasterio.open(path) as r:
        a = r.read().astype(np.float64); tr = r.transform
    return np.flip(a, axis=1), tr                                   # south-first rows
lai12 = np.stack([np.nan_to_num(tif(f'{RAW}/{YEAR}/inputs/lai_{YEAR}_{m:02d}.tif')[0][0], nan=0.0) for m in range(1, 13)])
met, tr = tif(f'{RAW}/{YEAR}/inputs/met_{YEAR}_{MON:02d}.tif')    # TEMP2, PAR, WINDSPD, PRES, QV, SOILM (monthly means)
ny, nx = lai12.shape[1:]
lon = tr.c + RES * (np.arange(nx) + 0.5); lat = (tr.f + tr.e * ny) + RES * (np.arange(ny) + 0.5)   # south-first centres
LON, LAT = np.meshgrid(lon, lat)
gf = dict(np.load(f'{RAW}/{YEAR}/inputs/growthform_{YEAR}.npz'))
ctf = landcover.canopy_type_fractions(gf, 3)                       # percent, (6, ny, nx); tree -> type 3 (production setting)
m = {'temp': met[0], 'par': met[1], 'wind': met[2], 'pres': met[3], 'qv': met[4], 'soilm': met[5]}
m, bad = fixes.fill_met_nearest(m)                                  # ocean / no-data cells -> nearest land (FIX-1 rule)
hours = np.arange(24)
loc = (hours[:, None, None] + LON[None] / 15.0) % 24.0
T = m['temp'][None] + 5.0 * np.cos(2 * np.pi * (loc - 15.0) / 24.0)
decl = np.deg2rad(23.45) * np.sin(2 * np.pi * (284 + doy) / 365.0)
cosz = np.sin(np.deg2rad(LAT))[None] * np.sin(decl) + np.cos(np.deg2rad(LAT))[None] * np.cos(decl) * np.cos(np.deg2rad(15.0 * (loc - 12.0)))
shape = np.maximum(cosz, 0.0)
PAR = shape * (m['par'] * 24.0 / np.maximum(shape.sum(0), 1e-6))[None]
const = lambda a: np.repeat(a[None], 24, axis=0)
common = dict(lon0=lon[0] - RES / 2, lat0=lat[0] - RES / 2, dlon=RES, dlat=RES)
TPL = f'{MGN}/Input/MGNMET/MET.MEGAN.SE_US.rad45.2018122.ncf'
write_ioapi(f'{OUT}/MET.ncf', TPL, {'TEMP2': T, 'PAR': PAR, 'WINDSPD': const(m['wind']), 'PRES': const(m['pres']), 'QV': const(m['qv']), 'PREC_ADJ': np.ones_like(T)},
            sdate=SDATE, stime=0, tstep=10000, units={'TEMP2': 'K', 'PAR': 'WATTS/M**2', 'WINDSPD': 'M/S', 'PRES': 'Pa', 'QV': 'KG/KG'}, **common)
laid = {f'LAI{k + 1:02d}': lai12[k][None] for k in range(12)}
# megvea.f builds the MGNERS header from leftover IOAPI description arrays and expects >= 20 variables in the
# LAI file (as in the 46-period official LAI3 file).  Pad with unused LAI13..LAI46 (never read with NLAI = 12).
laid.update({f'LAI{k + 1:02d}': np.zeros((1,) + lai12.shape[1:]) for k in range(12, 46)})
laid['LAT'] = LAT[None]; laid['LONG'] = LON[None]
write_ioapi(f'{OUT}/LAI.ncf', f'{MGN}/Input/MAP/LAI3_SE_US.ncf', laid, sdate=0, stime=0, tstep=10000, **common)
write_ioapi(f'{OUT}/CANTYP.ncf', f'{MGN}/Input/MAP/LDF_SE_US.2019b.ncf', {'CTS': ctf}, sdate=0, stime=0, tstep=10000, **common)
ldf = {'LDF03': 0.6, 'LDF04': 0.9, 'LDF05': 0.2, 'LDF06': 0.4}     # MEGVEA.EXT constants used by the production run
write_ioapi(f'{OUT}/LDF.ncf', f'{MGN}/Input/MAP/LDF_SE_US.2019b.ncf', {k: np.full((1, ny, nx), v) for k, v in ldf.items()}, sdate=0, stime=0, tstep=10000, **common)
with open(f'{OUT}/GRIDDESC', 'w') as f:
    f.write("' '\n'LATLON'\n  1  0.0  0.0  0.0  0.0  0.0\n' '\n'SC_BENCH'\n")
    f.write(f"'LATLON'  {lon[0] - RES / 2:.6f}  {lat[0] - RES / 2:.6f}  {RES:.6f}  {RES:.6f}  {nx}  {ny}  1\n' '\n")
np.savez_compressed(f'{OUT}/py_inputs.npz', T=T, PAR=PAR, wind=m['wind'], pres=m['pres'], qv=m['qv'], lai12=lai12, ctf=ctf, LAT=LAT, LON=LON, doy=doy, sdate=SDATE)
print('SC benchmark inputs:', OUT, 'SDATE', SDATE, 'grid', ny, nx, '| met cells filled', int(bad.sum()), '| LAI>0 cells', int((lai12[MON - 1] > 0).sum()))
