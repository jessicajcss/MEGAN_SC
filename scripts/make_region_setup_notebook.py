#!/usr/bin/env python3
"""Generator for run_your_region/region_setup.ipynb — the lay-person "Region Setup Wizard".

This produces a self-guided Colab notebook that lets someone with minimal computer
experience re-target the MEGAN pipeline to ANY region on Earth by (1) filling a few
plainly-labelled boxes and (2) uploading a boundary file. It reuses the SAME validated
`megan_sc` package functions the production notebook uses (build_domain, growth-form
pull, build_ef_grid_multiform, run_timestep, speciation), so none of the science is coded
twice, and if the built-in test runs, the region is ready for a full run.

Run:  python scripts/make_region_setup_notebook.py
Out:  run_your_region/region_setup.ipynb
"""
import json, os

CELLS = []
def md(text):   CELLS.append({"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)})
def code(text): CELLS.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": text.splitlines(keepends=True)})

# ─────────────────────────────────────────────────────────────────────────────
md(r"""# 🌎 Run MEGAN for **your** region — Setup Wizard

**You do not need to know how to code to use this.** You fill in a few boxes, press
one button, and this notebook does the technical work: it reads your region's outline,
pulls the weather, leaf-area and land-cover data for that exact area from Google Earth
Engine, builds the emission factors, and runs a short **test** to prove everything works
for your region. At the end it writes a ready-to-use configuration file and shows you a
sample emissions map of your own region.

This is an **independent helper**. It does not change or overwrite the original Santa
Catarina setup — it only *reads* the shared model code and *writes* new files named after
your region.

---
### Before you start you need three things (all free)
1. **A boundary file of your region** — a shapefile (`.shp` + its `.dbf/.shx/.prj`
   companions) **or** a GeoPackage (`.gpkg`). This is just the outline of the state,
   province, or country you want. Put it in your Google Drive inside the `megan-sc`
   folder, e.g. `MyDrive/megan-sc/data/boundaries/`. (Where to get one: see
   `START_HERE.md`, Step 1.)
2. **A Google account with Earth Engine turned on** — one-time, free sign-up at
   [earthengine.google.com](https://earthengine.google.com/). You will also need your
   *project id* (looks like `ee-yourname`). `START_HERE.md`, Step 2 shows exactly where
   to find it.
3. **The `megan-sc` folder already in your Google Drive** at `MyDrive/megan-sc` (the
   whole pipeline folder). If you can see it in Drive, you're set.

### How to run
In the menu: **Runtime ▸ Run all**. Then just watch — each step prints what it did. The
built-in test takes roughly **3–8 minutes** depending on your settings. If a step needs
something from you (like signing in to Earth Engine, or a missing file), it stops and
tells you in plain language what to do.
""")

# ── STEP 1 : settings ────────────────────────────────────────────────────────
md(r"""## ✏️ STEP 1 — Tell us about your region

Edit the values on the right of each `=` sign below. **Only this cell needs your input.**
Everything after it runs on its own. Read the comment (the text after `#`) on each line —
it explains what to type.
""")
code(r'''# ═══════════════ YOUR SETTINGS — edit these, then run every cell ═══════════════

REGION_NAME       = "my_region"        # a short nickname, NO SPACES (used in file names). e.g. "parana", "kenya"

SHAPEFILE_PATH    = "/content/drive/MyDrive/megan-sc/data/boundaries/my_region.shp"
                                       # full path to YOUR boundary file (.shp or .gpkg) inside Google Drive

COUNTRY_IS_BRAZIL = True               # True  -> uses MapBiomas 30 m land cover (best quality, Brazil only)
                                       # False -> uses global MODIS 500 m land cover (works anywhere on Earth)

# ── The period to simulate: first and last month, both included ──
START_YEAR        = 2020               # year of the first month
START_MONTH       = 1                  # first month: 1 = January ... 12 = December
END_YEAR          = 2020               # year of the last month
END_MONTH         = 12                 # last month. One month only: same year and month as the start

EE_PROJECT        = "ee-yourname"      # your Earth Engine project id  (START_HERE.md, Step 2)

RESOLUTION_DEG    = 0.05               # map cell size in degrees. 0.05 ≈ 5.5 km (fast).
                                       # 0.025 ≈ 2.8 km is sharper but ~4× slower. 0.1 for a quick look.

# ── Only needed if your region is OUTSIDE Santa Catarina (see STEP 5). Leave "" otherwise. ──
GLOBAL_ECOTYPE_PATH = ""               # path to the global MEGAN ecotype file (EVT), .nc or .tif, in your Drive

# ── Advanced (safe to leave as-is) ──
CANOPY_TYPE       = "auto"             # "auto" picks tropical vs temperate forest from your latitude; or set 2 or 3
TEST_MONTH        = "auto"             # "auto" picks a warm month of your period; or set 1–12
TEST_DAYS         = 2                  # how many days the quick test simulates (2 is plenty to prove it works)
# ═══════════════════════════════════════════════════════════════════════════════

# check the period and list its months (nothing to edit below this line)
import calendar
for _name, _v in (("START_YEAR", START_YEAR), ("END_YEAR", END_YEAR), ("START_MONTH", START_MONTH), ("END_MONTH", END_MONTH)):
    if not isinstance(_v, int) or isinstance(_v, bool):
        raise SystemExit(f"❌ {_name} must be a whole number without quotes (you typed {_v!r}).")
for _name, _v in (("START_MONTH", START_MONTH), ("END_MONTH", END_MONTH)):
    if not 1 <= _v <= 12:
        raise SystemExit(f"❌ {_name} must be a month number from 1 to 12 (you typed {_v}).")
if (END_YEAR, END_MONTH) < (START_YEAR, START_MONTH):
    raise SystemExit("❌ The last month comes before the first one. Check START_YEAR/START_MONTH and END_YEAR/END_MONTH.")
PERIOD_MONTHS = [(y, m) for y in range(START_YEAR, END_YEAR + 1) for m in range(1, 13)
                 if (START_YEAR, START_MONTH) <= (y, m) <= (END_YEAR, END_MONTH)]
PERIOD_DAYS = sum(calendar.monthrange(y, m)[1] for y, m in PERIOD_MONTHS)
PERIOD_START, PERIOD_END = f"{START_YEAR}-{START_MONTH:02d}", f"{END_YEAR}-{END_MONTH:02d}"
print("Settings recorded for region:", REGION_NAME)
print(f"Period: {PERIOD_START} to {PERIOD_END} = {len(PERIOD_MONTHS)} month(s), {PERIOD_DAYS} days "
      f"(months and days in UTC, the time of the ERA5-Land data)")
''')

# ── STEP 2 : environment ─────────────────────────────────────────────────────
md(r"""## ⚙️ STEP 2 — Set up the tools and connect your Google Drive

This installs the mapping libraries, connects your Google Drive (a pop-up will ask you to
allow it — click through), and finds the `megan-sc` folder. Nothing to edit here.
""")
code(r'''import sys, os, subprocess
print("Installing mapping libraries (about a minute)…")
subprocess.run([sys.executable, "-m", "pip", "-q", "install",
                "geopandas", "shapely", "rasterio", "netCDF4", "pyyaml", "earthengine-api"], check=False)

from google.colab import drive
drive.mount('/content/drive')

# locate the megan-sc folder in Drive
CANDIDATES = ["/content/drive/MyDrive/megan-sc", "/content/drive/MyDrive/MEGAN_072026/megan-sc",
              "/content/drive/My Drive/megan-sc"]
MEGAN_DIR = next((p for p in CANDIDATES if os.path.isdir(p)), None)
if MEGAN_DIR is None:
    raise SystemExit("❌ Could not find the 'megan-sc' folder in your Google Drive.\n"
                     "   Put the whole megan-sc folder at MyDrive/megan-sc and run this cell again.")
sys.path.insert(0, os.path.join(MEGAN_DIR, "src"))
print("✅ Found the pipeline at:", MEGAN_DIR)

import numpy as np, pandas as pd, calendar, datetime
import matplotlib.pyplot as plt
import geopandas as gpd
from megan_sc import domain as dom, run, landcover
from megan_sc.io import geotiff, netcdf as ncio
from megan_sc.engine import speciation as spc, species as sp
from megan_sc.gee import era5, lai as laimod, mapbiomas, modis_lc, region_from_bounds
print("✅ Model code loaded.")
''')

# ── STEP 3 : load + preview region ───────────────────────────────────────────
md(r"""## 🗺️ STEP 3 — Load your region and preview it

This reads your boundary file, lines it up to standard longitude/latitude (so any map
datum works — SIRGAS 2000, WGS84, UTM…), and lays a grid over it. The grey shape is your
region; the coloured cells are the ones the model will report on. Check that the shape
looks right before continuing.
""")
code(r'''if not os.path.exists(SHAPEFILE_PATH):
    raise SystemExit(f"❌ I can't find your boundary file:\n   {SHAPEFILE_PATH}\n"
                     "   Check the path in STEP 1. It must be inside your Google Drive.")

# Build a minimal config that points the grid at YOUR shapefile (shapefile mode = grid auto-fits the outline).
_cfg_p = os.path.join(MEGAN_DIR, "config", "example_config.yaml")
if not os.path.exists(_cfg_p):                       # the repository keeps it at its root
    _cfg_p = os.path.join(MEGAN_DIR, "example_config.yaml")
cfg = run.load_config(_cfg_p)
cfg["domain"]["mode"] = "shapefile"
cfg["domain"]["shapefile"] = SHAPEFILE_PATH
cfg["domain"]["resolution_deg"] = float(RESOLUTION_DEG)
# make the EF/speciation/ecotype paths absolute so they resolve from anywhere
for k in ("ef_table", "speciation_dir", "ecotype"):
    v = cfg["datasets"]["emission_factors"].get(k)
    if v and not os.path.isabs(v):
        cfg["datasets"]["emission_factors"][k] = os.path.join(MEGAN_DIR, v)

domain = dom.build_domain(cfg)                       # builds the grid AND a mask of cells inside your outline
lon, lat, area = domain["lon"], domain["lat"], domain["area"]
res = float(domain["resolution_deg"]); res_m = res * 111000.0
mask = domain["mask"]
if mask is None:
    raise SystemExit("❌ Your file loaded but no cells fell inside it — is it a valid polygon outline?")
inside = int(mask.sum())
mean_lat = float(np.mean(lat))

# auto choices
canopy_type = (2 if abs(mean_lat) < 23.5 else 3) if CANOPY_TYPE == "auto" else int(CANOPY_TYPE)
warm_month  = 1 if mean_lat < 0 else 7                 # January south of the equator, July north of it
want_month  = warm_month if TEST_MONTH == "auto" else int(TEST_MONTH)
# the test runs inside your period: the first month of the period with that number, else the first month of the period
test_year, test_month = next(((y, m) for y, m in PERIOD_MONTHS if m == want_month), PERIOD_MONTHS[0])
biome_txt   = "tropical broadleaf (type 2)" if canopy_type == 2 else "temperate broadleaf (type 3)"

print(f"Region: {REGION_NAME}")
print(f"Grid: {lat.size} rows × {lon.size} cols  @ {res}°  ({inside:,} cells inside your outline)")
print(f"Centre latitude {mean_lat:.1f}°  →  forest type: {biome_txt}")
print(f"Test month chosen: {test_year}-{test_month:02d}  "
      + ("(a warm month for your hemisphere)" if test_month == warm_month else "(inside your period)"))

# preview
gdf = gpd.read_file(SHAPEFILE_PATH).to_crs("EPSG:4326")
fig, ax = plt.subplots(figsize=(6, 6))
extent = [lon[0]-res/2, lon[-1]+res/2, lat[0]-res/2, lat[-1]+res/2]
ax.imshow(np.where(mask, 1.0, np.nan), origin="lower", extent=extent, cmap="YlGn", alpha=0.6)
gdf.boundary.plot(ax=ax, color="0.25", linewidth=0.8)
ax.set_title(f"{REGION_NAME} — grid cells the model will report on"); ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
plt.show()
''')

# ── STEP 4 : Earth Engine ────────────────────────────────────────────────────
md(r"""## 🔐 STEP 4 — Sign in to Google Earth Engine

A link/box appears — sign in with the Google account that has Earth Engine enabled, copy
the code back if asked. This is what lets the notebook pull weather, leaf-area and
land-cover data for your region. If you get a *project* error, fix `EE_PROJECT` in STEP 1.
""")
code(r'''import ee
try:
    ee.Authenticate()
    ee.Initialize(project=EE_PROJECT)
except Exception as e:
    raise SystemExit(f"❌ Earth Engine sign-in failed.\n   {e}\n"
                     "   Check EE_PROJECT in STEP 1, and that Earth Engine is enabled for your account.")
region = region_from_bounds(domain["bounds"])
print("✅ Earth Engine ready for project:", EE_PROJECT)

# are the weather and leaf-area data there for your whole period? (a few small requests)
def _n_images(collection, start, end):
    return int(ee.ImageCollection(collection).filterBounds(region).filterDate(start, end).size().getInfo())
_first, _last = PERIOD_MONTHS[0], PERIOD_MONTHS[-1]
_last_day = f"{_last[0]}-{_last[1]:02d}-{calendar.monthrange(*_last)[1]:02d}"
for _day in (f"{_first[0]}-{_first[1]:02d}-01", _last_day):
    if _n_images(era5.COLLECTION, ee.Date(_day), ee.Date(_day).advance(1, "day")) < 24:
        raise SystemExit(f"❌ Earth Engine has no complete ERA5-Land weather for {_day}.\n"
                         "   Choose a period inside the years the dataset covers (its page in the Earth Engine catalog "
                         "gives the last available date) and run again from STEP 1.")
_before = (_first[0] - (_first[1] == 1), (_first[1] - 2) % 12 + 1)   # leaf age needs the month before the period
for (_y, _m) in (_before, _last):
    _s = ee.Date.fromYMD(_y, _m, 1)
    if _n_images(laimod.COLLECTION, _s, _s.advance(1, "month")) == 0:
        raise SystemExit(f"❌ Earth Engine has no MODIS leaf-area images ({laimod.COLLECTION}) for {_y}-{_m:02d}.\n"
                         "   The model needs them for every month of the period and for the month before it. "
                         "Choose a later start or an earlier end and run again from STEP 1.")
print(f"✅ Weather and leaf-area data found for {PERIOD_START} to {PERIOD_END}.")
''')

# ── STEP 5 : land cover + ecotype ────────────────────────────────────────────
md(r"""## 🌱 STEP 5 — Land cover and vegetation type for your region

The model needs to know **what grows where** (trees, shrubs, grasses, crops) and MEGAN's
**ecotype** map (which sets the emission "recipe"). Two pieces:

**A) Growth-form fractions.** If your region is in **Brazil**, this uses MapBiomas (30 m,
very detailed). Anywhere else on Earth, it uses the global MODIS land-cover product
(500 m) — coarser, but it covers the whole planet. This is chosen automatically from your
`COUNTRY_IS_BRAZIL` setting.

**B) The ecotype map.** The pipeline ships an ecotype tile that covers **only Santa
Catarina**. If your region sits inside that box, it's reused automatically. If your region
is anywhere else, you need MEGAN's **global** ecotype file once: `EVT3d.nc`, in `EVT3d.zip` of
the MEGAN 3.2 input files on Zenodo (https://doi.org/10.5281/zenodo.10939297; `START_HERE.md`, Step 5). Put it in your Drive and set
`GLOBAL_ECOTYPE_PATH` in STEP 1; this cell will then **clip it to your region for you**.
""")
code(r'''VEG_YEAR = min(int(test_year), 2024)    # land cover of the test year (annual maps; MapBiomas coverage caps at 2024)

# ── A) growth-form fractions on your grid ──
if COUNTRY_IS_BRAZIL:
    print("Land cover: MapBiomas 30 m (Brazil). Pulling growth forms (1–2 min)…")
    gf_np = mapbiomas.growth_form_fractions_gridded(ee, region, domain, VEG_YEAR)
else:
    print("Land cover: global MODIS 500 m. Pulling growth forms…")
    img = modis_lc.growth_form_fractions(ee, region, VEG_YEAR, out_scale_m=res_m)
    bands = list(landcover.GROWTH_FORMS) + ["veg"]
    g = mapbiomas.sample_image_tiled(ee, img, domain, bands, ntiles=4)
    gf_np = {b: np.nan_to_num(g[b], nan=0.0) for b in bands}
print("   growth-form means:", {k: round(float(gf_np[k].mean()), 3) for k in landcover.GROWTH_FORMS})

# ── B) ecotype: reuse the SC tile if we're inside it, else clip a global file ──
sc_eco = os.path.join(MEGAN_DIR, "data", "sc_ecotype.nc")
import netCDF4 as _nc
def _coverage(path):
    d = _nc.Dataset(path); lo = d.variables["lon"][:]; la = d.variables["lat"][:]; d.close()
    return float(lo.min()), float(lo.max()), float(la.min()), float(la.max())
bx0, by0, bx1, by1 = domain["bounds"]                 # region extent (lon0, lat0, lon1, lat1)
elo0, elo1, ela0, ela1 = _coverage(sc_eco)
inside_sc = (bx0 >= elo0 and bx1 <= elo1 and by0 >= ela0 and by1 <= ela1)

def clip_global_ecotype(src, out_path, bbox, margin=0.15):
    """Clip a global MEGAN ecotype (EVT) file — NetCDF or GeoTIFF — to the region and save as NetCDF."""
    x0, y0, x1, y1 = bbox; x0-=margin; y0-=margin; x1+=margin; y1+=margin
    if src.lower().endswith((".nc", ".nc4", ".netcdf")):
        d = _nc.Dataset(src)
        vname = "EVT" if "EVT" in d.variables else [v for v in d.variables if v not in ("lon","lat","longitude","latitude")][0]
        lo = d.variables["lon"][:] if "lon" in d.variables else d.variables["longitude"][:]
        la = d.variables["lat"][:] if "lat" in d.variables else d.variables["latitude"][:]
        evt = np.asarray(d.variables[vname][:]); d.close()
        ix = np.where((lo >= x0) & (lo <= x1))[0]; iy = np.where((la >= y0) & (la <= y1))[0]
        lo2, la2, evt2 = lo[ix], la[iy], evt[np.ix_(iy, ix)]
    else:  # GeoTIFF
        import rasterio
        from rasterio.windows import from_bounds
        with rasterio.open(src) as r:
            win = from_bounds(x0, y0, x1, y1, r.transform)
            evt2 = r.read(1, window=win).astype("float64")
            tr = r.window_transform(win); h, w = evt2.shape
            lo2 = tr.c + tr.a*(np.arange(w)+0.5); la2 = tr.f + tr.e*(np.arange(h)+0.5)
            if la2[0] > la2[-1]:                          # make latitude ascending (south-first)
                la2 = la2[::-1]; evt2 = evt2[::-1, :]
    if evt2.size == 0:
        raise SystemExit("❌ The global ecotype file does not cover your region — is it truly global?")
    o = _nc.Dataset(out_path, "w"); o.createDimension("lat", la2.size); o.createDimension("lon", lo2.size)
    o.createVariable("lat", "f8", ("lat",))[:] = la2; o.createVariable("lon", "f8", ("lon",))[:] = lo2
    o.createVariable("EVT", "i4", ("lat", "lon"))[:] = np.rint(evt2).astype("int32"); o.close()
    return out_path

if inside_sc:
    print("Ecotype: your region is inside the Santa Catarina tile — reusing it. ✅")
    cfg["datasets"]["emission_factors"]["ecotype"] = sc_eco
elif GLOBAL_ECOTYPE_PATH:
    if not os.path.exists(GLOBAL_ECOTYPE_PATH):
        raise SystemExit(f"❌ GLOBAL_ECOTYPE_PATH not found: {GLOBAL_ECOTYPE_PATH}")
    out = os.path.join(MEGAN_DIR, "data", "boundaries", f"{REGION_NAME}_ecotype.nc")
    print("Ecotype: clipping the global MEGAN ecotype to your region…")
    clip_global_ecotype(GLOBAL_ECOTYPE_PATH, out, (bx0, by0, bx1, by1))
    cfg["datasets"]["emission_factors"]["ecotype"] = out
    print("   ✅ saved:", out)
else:
    raise SystemExit(
        "⚠️  Your region is OUTSIDE Santa Catarina, so you need MEGAN's GLOBAL ecotype file (EVT) once.\n"
        "    1) Download it (see START_HERE.md, Step 5) and put it in your Drive.\n"
        "    2) Set GLOBAL_ECOTYPE_PATH in STEP 1 to that file.\n"
        "    3) Run all cells again. This step will then clip it to your region automatically.")

# FIX-5: ecotype no-data pixels take the nearest valid ecotype, as in notebook 10
from megan_sc import fixes
_eco = cfg["datasets"]["emission_factors"]["ecotype"]
_eco_f = os.path.join(MEGAN_DIR, "data", "boundaries", f"{REGION_NAME}_ecotype_filled.nc")
print("Ecotype no-data pixels filled (FIX-5):", fixes.fill_ecotype_nearest(_eco, _eco_f))
cfg["datasets"]["emission_factors"]["ecotype"] = _eco_f
''')

# ── STEP 6 : EF + canopy ─────────────────────────────────────────────────────
md(r"""## 🧪 STEP 6 — Build the emission factors for your region

This combines MEGAN's ecotype→species→emission-factor recipe with your region's
growth-form mix — exactly the MEGAN structure, just fitted to your land cover.
""")
code(r'''ef  = run.build_ef_grid_multiform(cfg, domain, gf_np)
ctf = landcover.canopy_type_fractions(gf_np, canopy_type)
isop_i = sp.index_of("ISOP")
print(f"✅ Emission factors built. Mean isoprene EF over vegetated cells: "
      f"{float(np.nanmean(np.where(gf_np['veg']>0, ef[isop_i], np.nan))):.2f} (nanomol m⁻² s⁻¹)")
''')

# ── STEP 7 : the shared engine + quick test ──────────────────────────────────
md(r"""## ✅ STEP 7 — Quick test: run a couple of days and see your region's map

This runs the **real MEGAN engine** — the same code the full pipeline uses — over a few
days of your test month, then draws an isoprene emissions map of your region. If you see a
sensible map, **your region works** and you're ready for a full run (STEP 9). Nothing here
is a shortcut in the science; it's just a short time span so it finishes quickly.
""")
code(r'''from megan_sc import fixes
PATCHES = fixes.apply_engine_patches()              # FIX-6, FIX-7, FIX-8 (src/megan_sc/fixes.py)
assert fixes.patches_active()
MET_BANDS = ['TEMP2','PAR','WINDSPD','PRES','QV','SOILM']
CLS = list(sp.NAMES)
SAMPLE_HOURS = 6                                    # test: 4 steps/day (full runs use 3 -> 8 steps/day)
HOURS = list(range(0, 24, SAMPLE_HOURS))
SEC_PER_STEP = SAMPLE_HOURS * 3600.0

def sample_multi(img, bands):
    grid = {'dimensions': {'width': int(lon.size), 'height': int(lat.size)},
            'affineTransform': {'scaleX': res, 'shearX': 0.0, 'translateX': float(lon[0]-res/2),
                                'shearY': 0.0, 'scaleY': -res, 'translateY': float(lat[-1]+res/2)},
            'crsCode': 'EPSG:4326'}
    st_ = ee.data.computePixels({'expression': img.select(bands), 'grid': grid, 'fileFormat': 'NUMPY_NDARRAY'})
    return {b: np.flipud(np.asarray(st_[b], dtype='float64')) for b in bands}

_LAI_RAW = {}
def _shift(year, m, k):
    t = year * 12 + (m - 1) + k
    return t // 12, t % 12 + 1

def month_lai_raw(year, m, required=True):
    """Monthly MODIS LAI on the grid, non-vegetated cells set to 0 (cached)."""
    if (year, m) not in _LAI_RAW:
        try:
            L = sample_multi(laimod.monthly_lai(ee, region, year, m).rename('LAI'), ['LAI'])['LAI']
        except Exception as e:
            if required:
                raise
            print(f"   LAI {year}-{m:02d} not available ({str(e)[:60]}); not used for gap filling")
            L = np.zeros((lat.size, lon.size))
        veg = sum(gf_np[g] for g in landcover.GROWTH_FORMS)
        _LAI_RAW[(year, m)] = np.nan_to_num(L, nan=0.0, posinf=0.0, neginf=0.0) * (veg > 0)
    return _LAI_RAW[(year, m)]

def month_lai(year, m):
    """FIX-4: gaps of vegetated cells filled from the adjacent months (fixes.fill_lai_gaps, as notebook 10)."""
    veg = sum(gf_np[g] for g in landcover.GROWTH_FORMS)
    filled, _ = fixes.fill_lai_gaps(month_lai_raw(year, m), month_lai_raw(*_shift(year, m, -1), required=False),
                                    month_lai_raw(*_shift(year, m, +1), required=False), veg)
    return filled

def run_one_month(year, m, ndays):
    """Per-class tonnes for `ndays` days of month `m`, with the procedures FIX-1 to FIX-8 of notebook 10.

    FIX-1 (meteorology of cells without reanalysis data from the nearest valid cell), FIX-2 (daily means of the
    simulation day), FIX-3 (leaf age from the previous month's LAI) and FIX-4 (LAI gaps) act here; FIX-5 acted on the
    ecotype in STEP 5; FIX-6, FIX-7 and FIX-8 are the engine patches installed above. One difference from notebook
    10: the daily mean temperature is the mean of the sampled steps of the day, while notebook 10 averages the 24
    hourly ERA5-Land values (one extra Earth Engine request per day)."""
    ny, nx = lat.size, lon.size
    acc_t = np.zeros((len(CLS), ny, nx))
    lai = month_lai(year, m)                                # FIX-4
    laip = month_lai(*_shift(year, m, -1))                  # FIX-3
    for d in range(1, ndays + 1):
        doy = datetime.date(year, m, d).timetuple().tm_yday
        mets = {h: sample_multi(era5.megan_met_image(ee, region, ee.Date(f'{year}-{m:02d}-{d:02d}T{h:02d}:00')), MET_BANDS)
                for h in HOURS}
        tp = [fixes.fill_met_nearest({'temp': mets[h]['TEMP2'], 'par': mets[h]['PAR']}, keys=('temp', 'par'))[0] for h in HOURS]
        d_temp, d_ppfd = fixes.daily_means_from_samples([x['temp'] for x in tp], [x['par'] for x in tp])   # FIX-2
        for h in HOURS:
            met = mets[h]
            metd = {'temp': met['TEMP2'], 'par': met['PAR'], 'wind': met['WINDSPD'], 'pres': met['PRES'],
                    'qv': met['QV'], 'soilm': met['SOILM'], 'lai': lai}
            r = fixes.run_timestep_fixed(domain, ef, metd, day=doy, hour_utc=h, canopy_type_fractions=ctf,
                                         d_temp=d_temp, d_ppfd=d_ppfd, laip=laip, nlai=12, wilt=0.12,
                                         date_yyyyddd=year * 1000 + doy, strict=True)   # FIX-1 inside
            tph = spc.compound_tonnes(spc.flux_to_mol_s(r['compounds'], area[None, :, :]), SEC_PER_STEP)
            acc_t += spc.compounds_to_classes(tph)
        print(f"   day {d}/{ndays}")
    return acc_t

test_days = min(int(TEST_DAYS), calendar.monthrange(test_year, test_month)[1])
print(f"Running a {test_days}-day test for {REGION_NAME}, {test_year}-{test_month:02d}…")
tonnes = run_one_month(test_year, test_month, test_days)

# draw the isoprene map, clipped to your region
isop = np.where(mask, tonnes[isop_i], np.nan)
fig, ax = plt.subplots(figsize=(6.5, 6))
extent = [lon[0]-res/2, lon[-1]+res/2, lat[0]-res/2, lat[-1]+res/2]
im = ax.imshow(np.ma.masked_invalid(isop), origin='lower', extent=extent, cmap='Greens')
gdf.boundary.plot(ax=ax, color='0.3', linewidth=0.6)
ax.set_title(f"{REGION_NAME}: isoprene, {test_days}-day test ({test_year}-{test_month:02d})")
fig.colorbar(im, ax=ax, shrink=0.8, label='t / test-period')
plt.show()
print("✅ If that map looks sensible, your region is ready. Save the config in STEP 8.")
''')

# ── STEP 8 : write config ────────────────────────────────────────────────────
md(r"""## 💾 STEP 8 — Save your region's configuration file

This writes `config/<region>_run.yaml` — a complete settings file for your region. The
full-run notebook (or STEP 9 below) reads it. You can open it later in any text editor.
""")
code(r'''import yaml
cfg["run"]["name"] = f"{REGION_NAME}_{PERIOD_START}_{PERIOD_END}"
cfg["run"]["description"] = f"MEGAN biogenic inventory for {REGION_NAME}"
# write the region's box too, so BOTH shapefile-mode and bbox-mode runs target the right area
_x0, _y0, _x1, _y1 = domain["bounds"]
cfg["domain"]["bbox"] = {"lon_min": round(float(_x0), 4), "lon_max": round(float(_x1), 4),
                         "lat_min": round(float(_y0), 4), "lat_max": round(float(_y1), 4)}
cfg["period"]["start"] = PERIOD_START               # YYYY-MM, first month of the period
cfg["period"]["end"]   = PERIOD_END                 # YYYY-MM, last month (both included)
cfg["datasets"]["vegetation"]["growth_form_source"] = "mapbiomas" if COUNTRY_IS_BRAZIL else "modis"
cfg["datasets"]["vegetation"]["tree_canopy_type"] = int(canopy_type)
cfg["earth_engine"]["project"] = EE_PROJECT

os.makedirs(os.path.join(MEGAN_DIR, "config"), exist_ok=True)
cfg_path = os.path.join(MEGAN_DIR, "config", f"{REGION_NAME}_run.yaml")
with open(cfg_path, "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)
print("✅ Saved your configuration:", cfg_path)
print("   period      :", cfg["period"]["start"], "to", cfg["period"]["end"])
print("   land cover  :", cfg["datasets"]["vegetation"]["growth_form_source"])
print("   forest type :", cfg["datasets"]["vegetation"]["tree_canopy_type"])
print("   ecotype     :", cfg["datasets"]["emission_factors"]["ecotype"])
''')

# ── STEP 9 : full run ────────────────────────────────────────────────────────
md(r"""## 🚀 STEP 9 — Run the full inventory for your region (optional, longer)

The test in STEP 7 already used the real engine. To produce the **full monthly inventory**
for every month of the period you chose in STEP 1 (from `START_YEAR`/`START_MONTH` to
`END_YEAR`/`END_MONTH`), run the cell below. It reuses the exact same functions and writes,
per month, a NetCDF of tonnes and a statistics CSV into `megan_sc_outputs_<region>/<year>/`.
The time grows with the number of days and of grid cells: a full year at fine resolution
can take **hours**, so the first time choose one month and a coarse `RESOLUTION_DEG`.
Months are calendar months in UTC, the time of the ERA5-Land data.

> For the full-featured production run (checkpointing so it resumes after a disconnect,
> per-municipality tables, input GeoTIFFs, fix logs), use `notebooks/10_production_run.ipynb`
> with your region's configuration and boundary files; the edits are listed in
> `docs/tutorial/11_your_region.md`. This cell uses the same engine and procedures.
""")
code(r'''RUN_FULL = False    # ← change to True when you're ready for the (long) full run
SAMPLE_HOURS = 3                                    # full run: 8 steps per day, as the Santa Catarina production run
HOURS = list(range(0, 24, SAMPLE_HOURS)); SEC_PER_STEP = SAMPLE_HOURS * 3600.0

print(f"Period {PERIOD_START} to {PERIOD_END}: {len(PERIOD_MONTHS)} month(s), {PERIOD_DAYS} days, "
      f"{PERIOD_DAYS * len(HOURS)} model steps")
if RUN_FULL:
    OUT_ROOT = os.path.join(MEGAN_DIR, f"megan_sc_outputs_{REGION_NAME}")
    years = sorted({y for y, _ in PERIOD_MONTHS})
    for year in years:
        OUT = os.path.join(OUT_ROOT, str(year))
        for sub in ("emissions", "stats"): os.makedirs(os.path.join(OUT, sub), exist_ok=True)
        vy = min(year, 2024)
        # refresh land cover for this year (annual)
        if COUNTRY_IS_BRAZIL:
            gfy = mapbiomas.growth_form_fractions_gridded(ee, region, domain, vy)
        else:
            _img = modis_lc.growth_form_fractions(ee, region, vy, out_scale_m=res_m)
            _b = list(landcover.GROWTH_FORMS) + ["veg"]
            _g = mapbiomas.sample_image_tiled(ee, _img, domain, _b, ntiles=4)
            gfy = {b: np.nan_to_num(_g[b], nan=0.0) for b in _b}
        globals()['gf_np'] = gfy
        efy  = run.build_ef_grid_multiform(cfg, domain, gfy)
        cty  = landcover.canopy_type_fractions(gfy, canopy_type)
        globals()['ef'] = efy; globals()['ctf'] = cty
        for m in [mm for yy, mm in PERIOD_MONTHS if yy == year]:
            nd = calendar.monthrange(year, m)[1]
            print(f"{REGION_NAME} {year}-{m:02d}: full month ({nd} days)…")
            tm = run_one_month(year, m, nd)
            ncio.write_emissions(os.path.join(OUT, "emissions", f"emis_{year}_{m:02d}.nc"),
                                 tm, CLS, lon, lat, units="tonnes/month",
                                 global_attrs={"year": year, "month": m, "region": REGION_NAME})
            from megan_sc import stats as _st
            rows = _st.species_stats_rows(tm, CLS, year, m, mask=domain.get("mask"))
            pd.DataFrame(rows).to_csv(os.path.join(OUT, "stats", f"stats_{year}_{m:02d}.csv"), index=False)
        print(f"✅ {year} done → {OUT}")
    print("🎉 Full run complete for", REGION_NAME)
else:
    print("Full run is OFF. Set RUN_FULL = True in this cell when you want the complete inventory.")
''')

md(r"""---
### 🎉 That's it
When STEP 7 shows a sensible map, your region is set up. Your settings live in
`config/<region>_run.yaml`, your test map and (if you ran STEP 9) your inventory live in
`megan_sc_outputs_<region>/`. The original Santa Catarina files are untouched.

Questions or stuck on a step? See **`START_HERE.md`** in this folder — it walks through each
step with pictures-in-words, including where to get a boundary file and the ecotype file.
""")

# ─────────────────────────────────────────────────────────────────────────────
NB = {"cells": CELLS,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python"}, "colab": {"provenance": []}},
      "nbformat": 4, "nbformat_minor": 5}

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "run_your_region")
os.makedirs(OUT_DIR, exist_ok=True)
OUT = os.path.join(OUT_DIR, "region_setup.ipynb")
with open(OUT, "w") as f:
    json.dump(NB, f, indent=1)
print(f"wrote {OUT} with {len(CELLS)} cells")
