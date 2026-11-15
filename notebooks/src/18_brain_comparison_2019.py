# %% [markdown]
# # 18 · Comparison with the BRAIN biogenic emissions of 2019 (MEGAN, domain d02, 4 km)
#
# This notebook compares the 2019 biogenic emissions of this inventory (notebook 10) with the MEGAN emissions of the
# Brazilian Atmospheric Inventories (BRAIN; Hoinaski et al. 2024), domain d02 (southern Brazil, 4 km), published in the
# Science Data Bank (doi:10.57760/sciencedb.09886; Hoinaski et al. 2023). It serves three purposes:
#
# * compare domain, speciation, units and parametrizations of the two inventories before interpreting differences
#   between models;
# * give a technical basis to the discussion of the 7.6-fold (isoprene) and 18.9-fold (soil NO) differences from the
#   previous state inventory (Hoinaski et al. 2022), with BRAIN 2019 as a third estimate;
# * keep the wording conditional, because inventories built with different model versions and inputs cannot isolate
#   the effect of any one of them.
#
# 2019 is the base year of both BRAIN and the previous state inventory. It is not affected by the LAI correction of
# notebook 17 (2024-2025 only).
#
# ## 1. What the BRAIN files are
#
# Checked on a d02 file of the same series (`MEGANv31.SC_2019.CB6X.2020015.BDSNP.ncf`, 15 January 2020) and checked again
# by this notebook on every 2019 file (cells 6 and 7).
#
# | Item | Value | Source |
# |---|---|---|
# | format | I/O API netCDF (classic, 64-bit offset), one file per day, zipped; 25 hourly steps from 00:00 UTC of the day to 00:00 UTC of the next day | file header, `TFLAG` |
# | grid | `SC_2019_CROSS`: 233 columns × 243 rows; Mercator (`GDTYP` = 7) with true scale at the equator (`P_ALP` = 0), central meridian 54.857° W, cells of 4000 m in projected units; sphere of radius 6,370 km | header, `FILEDESC` |
# | extent | about 55.8° W to 47.4° W and 30.1° S to 22.2° S, so all of Santa Catarina | computed in cell 4 |
# | species | 34 CB6X species (`ISOP`, `TERP`, `SQT`, `MEOH`, `CO`, `NO`, ...), all in mol s⁻¹ per grid cell | `VAR-LIST`, `units` |
# | processing | MGN2MECH of MEGANv3.21 (path in `FILEDESC`), soil NO from BDSNP (`BDSNP` in the file name), meteorology from WRF ARW 4.0 through MCIP 5.1 | `FILEDESC`, `UPNAM`, file name |
#
# The article describes the BRAIN biogenic emissions as MEGAN 3.2 driven by WRF and MCIP, d01 over Brazil at 20 km and
# d02 over southern Brazil at 4 km, year 2019 (Hoinaski et al. 2024, Sect. 2.1.4 and 3). It does not report biogenic
# totals, and it does not name the LAI, land-cover or emission-factor inputs of MEGAN.
#
# ## 2. How the comparison is made
#
# 1. **Same region.** Santa Catarina is the IBGE polygon of the article (`sc_state_ibge.gpkg`). Every cell of both grids
#    is weighted by the fraction of its area inside the polygon, computed in an equal-area projection. For this
#    inventory, the cell-centre mask of the article is also reported, so the published totals can be traced.
# 2. **Same species.** This inventory stores the 19 MEGAN classes in mass. They are converted to the 34 CB6X species
#    with the two tables MGN2MECH uses, read from `MEGANv3.21.tar.gz` (Guenther 2024): `SPC_NOCONVER.EXT` (19 classes →
#    201 compounds, molar masses and molar fractions) and `MAP_CV2CB6X.EXT` (201 compounds → CB6X). MEGANv3.21 splits each
#    class into compounds with fixed molar fractions, so this conversion is exact. Where a class is the only source of a
#    CB6X species (isoprene, methanol, acetone, ethanol and acetaldehyde, organic acids, light hydrocarbons, oxidation
#    products, stress compounds, other VOC, CO and NO), the BRAIN class mass is recovered exactly too. Monoterpenes and
#    sesquiterpenes are compared as the CB6X species `TERP` and `SQT`, in moles and in mass with the CB6X molar masses
#    (136.23 and 204.35 g mol⁻¹), because several classes feed each of them.
# 3. **Same units.** Moles and tonnes. A BRAIN daily total is the sum of the 24 hourly rates of 00 to 23 UTC times
#    3600 s. The 25th step (00 UTC of the next day) belongs to the next file and is not counted. Soil NO is expressed as
#    mass of NO (30.01 g mol⁻¹), as in this inventory and in Hoinaski et al. (2022). The 46.00 g mol⁻¹ of `SPC_CB6X.EXT`
#    is used by MGN2MECH only for its tonnes-per-hour option.
# 4. **Same period.** 2019, calendar months in UTC (both models run in UTC). The download list holds 364 of the 365 daily
#    files (31 December is missing). A month with missing days is scaled by days in the month / days available and
#    enters the comparison only with at least `MIN_DAYS` days; the factors are written to the QA table.
# 5. **Same grid for maps and cell statistics.** This inventory (0.025°) is remapped conservatively onto the BRAIN cells
#    with exact area overlaps on the sphere (both grids are rectangles in longitude and latitude; Jones 1999). The two
#    datasets are then compared as fluxes per true area.
# 6. **Cell area used by MGN2MECH.** MGN2MECH turns fluxes (nmol m⁻² s⁻¹) into rates per cell with
#    `garea = xcell3d * ycell3d` (`mgn2mech.F` of MEGANv3.21; the line is checked in cell 2), the nominal area of 16 km².
#    On a Mercator grid true at the equator, the true area of a cell at latitude φ is about 16 cos²φ km² (Snyder 1987),
#    about 12.2 to 12.9 km² over Santa Catarina. The delivered rates are then MEGAN fluxes × 16 km², and the emission of
#    the true area of a cell is the delivered rate × A_true / A_nominal. The notebook reports the state totals both ways
#    ("BRAIN, true area" and "BRAIN, as delivered"). Maps and cell statistics use fluxes and do not depend on this choice.
# 7. **Co-registration.** The grid described by the BRAIN header is tested against this inventory (cell 9): correlation of
#    the isoprene fields for shifts of up to ±4 cells, and the position of the east coast row by row. On the
#    15 January 2020 file, the header grid gave the best correlation without shift and a median coast offset of 0 cells.
#    Two other readings of the header failed: a Mercator true at 27° S gave no correlation (r = −0.08), and cells of
#    4000 m true at 27° S also gave r = −0.08, with the coast displaced by a median of 8 cells.
#
# ## 3. Checks
#
# * one file per day; 25 hourly steps in order; identical grid, species list and units (mol s⁻¹) in all files;
# * CB6X species that `MAP_CV2CB6X.EXT` never produces (for example `XYL`, `TOL`, `BENZ`, `NR`, `NH3`, `KET`, `N2O`) must be
#   zero in BRAIN;
# * speciation fingerprints: CB6X species produced by a single MEGAN class keep fixed molar ratios (`ALD2`/`ETOH` = 1 from
#   class ETOH; `FACD`/`AACD` from class ACID; `CH4`, `TRS`, `CH3BR`, `CH3CL`, `CH3I` from class OTHER; four species from
#   class LVOC). The same ratios in BRAIN show that BRAIN used the same speciation tables, which the class-level
#   comparison requires;
# * daily totals from the 24 hourly rates against the trapezoid over the 25 steps (difference reported);
# * this inventory: the 2019 isoprene of the article (cell-centre mask, ρ = 1) is reproduced from the monthly files.
#
# ## 4. What the comparison can and cannot show
#
# Both inventories use the MEGANv3.21 code: BRAIN through the Fortran, this inventory through its Python port, which was
# checked against the Fortran (notebook 12). They differ in meteorology (WRF at 4 km against ERA5-Land at 0.1°), vegetation
# inputs (those of BRAIN are not named in the article; this inventory uses MapBiomas growth forms at 30 m and the MODIS
# MCD15A3H LAI), soil NO scheme (BDSNP against the MEGANv3.21 default), grid, time step (1 h against 3 h) and the canopy
# factor ρ of isoprene (not applied by the code; 0.96 in the published numbers of this inventory). The ratios measure the
# combined effect of all these differences. Assigning a difference to one input requires runs that change one input at
# a time (comment 38).
#
# ## 5. How to run
#
# 1. Copy `BRAIN_output_urls.txt` to Drive (`URL_FILE`, cell 1). Only its 364 lines of `/V1/MEGAN/SC_4X4/` are used; the
#    BRAVES, FINN and IND lines (vehicles, fires, industry) are ignored.
# 2. Keep `MEGANv3.21.tar.gz` in Drive (`MGN_TGZ`, the same file as notebook 12) and the production run of notebook 10 in
#    `ROOT`, with the boundary files of `data/boundaries` in `BND`.
# 3. Run all cells. The daily zips are downloaded from the Science Data Bank to the Colab disk, read and deleted. Each
#    month is saved in `OUT_DIR/brain_monthly/` when it is complete, so after a disconnection run all cells again. Each
#    zip holds one file of about 190 MB (34 species × 25 hours × 243 × 233 cells). If the downloads fail from Colab, save
#    the zips by hand in `BRAIN_DIR`; the notebook reads local files first.
# 4. No Earth Engine is needed.
#
# ## 6. Outputs (`OUT_DIR`)
#
# | File | Content |
# |---|---|
# | `brain_files_qa.csv` | one row per BRAIN day: steps, grid, units, species, trapezoid check, read time |
# | `brain_months.csv` | days available, scale factor and inclusion of each month |
# | `brain_speciation_checks.csv` | zero species and the speciation fingerprints |
# | `comparison_2019_cb6x_species.csv` | annual Santa Catarina totals of the 34 CB6X species (mol): BRAIN (true area, as delivered) and this inventory (area-weighted, cell-centre) |
# | `comparison_2019_classes_t.csv` | the same in tonnes by MEGAN class (or CB6X species for monoterpenes and sesquiterpenes), with the previous state inventory |
# | `comparison_2019_monthly_t.csv` | monthly totals of isoprene, monoterpenes, sesquiterpenes, methanol, CO and soil NO |
# | `comparison_2019_mesoregions_t.csv` | the same species by mesoregion |
# | `comparison_2019_cells.csv` | cell statistics on the 4-km grid, by dominant growth form |
# | `brain_daily_sc.csv` | daily BRAIN totals in Santa Catarina |
# | `parametrizations.csv` | inputs and settings of the two inventories, with the source of each entry |
# | `fig_brain_*.png` | co-registration, maps, monthly cycle, cell scatter, daily series, ratios |
# | `brain_monthly/brain_2019_MM.npz` | monthly BRAIN grids of the 34 species (mol per cell, as delivered) |
# | `provenance_brain_comparison.json` | settings, file counts, hashes of the MEGANv3.21 tables, grid header |
#
# ## References
#
# * Guenther A (2024) MEGANv3.21.tar.gz emission calculator code (3.21). Zenodo. https://doi.org/10.5281/zenodo.10526206
# * Guenther AB, Jiang X, Heald CL, Sakulyanontvittaya T, Duhl T, Emmons LK, Wang X (2012) The Model of Emissions of
#   Gases and Aerosols from Nature version 2.1 (MEGAN2.1). Geosci Model Dev 5:1471–1492. https://doi.org/10.5194/gmd-5-1471-2012
# * Hoinaski L et al (2022) Avaliação do impacto das emissões antropogênicas e naturais na qualidade do ar em Santa
#   Catarina: ano base 2021. LCQAr/UFSC, Florianópolis. ISBN 978-85-8328-156-6
# * Hoinaski L et al (2023) BRAIN emissions, domain d02 [dataset]. Science Data Bank. https://doi.org/10.57760/sciencedb.09886
# * Hoinaski L et al (2024) Brazilian Atmospheric Inventories – BRAIN: a comprehensive database of air quality in
#   Brazil. Earth Syst Sci Data 16(5):2385–2405. https://doi.org/10.5194/essd-16-2385-2024
# * Hudman RC, Moore NE, Mebust AK, Martin RV, Russell AR, Valin LC, Cohen RC (2012) Steps towards a mechanistic model
#   of global soil nitric oxide emissions: implementation and space based-constraints. Atmos Chem Phys 12:7779–7795.
#   https://doi.org/10.5194/acp-12-7779-2012
# * Jones PW (1999) First- and second-order conservative remapping schemes for grids in spherical coordinates. Mon
#   Weather Rev 127(9):2204–2210. https://doi.org/10.1175/1520-0493(1999)127<2204:FASOCR>2.0.CO;2
# * Otte TL, Pleim JE (2010) The Meteorology-Chemistry Interface Processor (MCIP) for the CMAQ modeling system: updates
#   through MCIPv3.4.1. Geosci Model Dev 3:243–256. https://doi.org/10.5194/gmd-3-243-2010
# * Snyder JP (1987) Map projections: a working manual. US Geological Survey Professional Paper 1395. https://doi.org/10.3133/pp1395
# * Yienger JJ, Levy II H (1995) Empirical model of global soil-biogenic NOx emissions. J Geophys Res 100(D6):11447–11464.
#   https://doi.org/10.1029/95JD00370

# %%
# 1 · Paths and settings ---------------------------------------------------------------------------------------------
import os, sys, glob, json, re, time, zipfile, tarfile, hashlib, shutil, warnings, calendar, urllib.parse, datetime as dt
try:
    from google.colab import drive; drive.mount('/content/drive'); IN_COLAB = True
except Exception:
    IN_COLAB = False
%pip -q install netCDF4 pyproj geopandas shapely rasterio requests
ROOT      = '/content/drive/MyDrive/megan_sc_outputs_v10'          # <-- production run of notebook 10
BND       = '/content/drive/MyDrive/megan-sc/data/boundaries'      # <-- sc_state_ibge.gpkg, sc_mesoregions_2024.gpkg
URL_FILE  = '/content/drive/MyDrive/BRAIN/BRAIN_output_urls.txt'   # <-- list of download links of the BRAIN d02 files
BRAIN_DIR = '/content/drive/MyDrive/BRAIN/MEGAN_SC_4X4'            # <-- zips saved by hand (optional; read before any download)
MGN_TGZ   = '/content/drive/MyDrive/MEGANv3.21.tar.gz'             # <-- same file as notebook 12 (Zenodo 10.5281/zenodo.10526206)
OUT_DIR   = os.path.join(ROOT, 'analysis_v10', 'brain_comparison')
WORK      = '/content/brain_tmp'                                   # Colab disk; zips and extracted files are deleted after reading
YEAR      = 2019
N_THREADS = 4                       # parallel downloads
KEEP_ZIPS = False                   # True: keep a copy of every downloaded zip in BRAIN_DIR
MIN_DAYS  = 25                      # a month enters the comparison only with at least this many BRAIN days
RHO_ISOP  = 0.96                    # isoprene canopy factor of the published numbers (Guenther et al. 2006); the code applies none
DOMINANT  = 0.6                     # growth-form fraction above which a cell is "dominated" (as in notebook 16)
FULL_CELL = 0.99                    # cells with at least this fraction inside Santa Catarina enter the cell statistics
PREV_2019_KG = dict(ISOP=33_280_584.0, CO=7_495_807.0, NO=571_106.94)   # Hoinaski et al. (2022), LCQAr report, Table 12
MAKE_FIGS = True

# %% [markdown]
# ## Cell 2 · MEGANv3.21 speciation tables
#
# The class → CB6X conversion is built from the tables MGN2MECH itself reads, with its semantics: each line of
# `SPC_NOCONVER.EXT` assigns a compound to a class with a molar fraction (a later line for the same compound replaces the
# earlier one, as the Fortran assignment does), and each line of `MAP_CV2CB6X.EXT` adds a compound, times a conversion
# factor, to a CB6X species. `K[c, j]` is the number of moles of CB6X species j per gram of class c. The cell also
# confirms the `garea` line of `mgn2mech.F` quoted in Section 2 of the header.

# %%
# 2 · MEGANv3.21 tables: 19 classes -> 201 compounds -> 34 CB6X species ------------------------------------------------
import numpy as np, pandas as pd
warnings.filterwarnings('ignore', category=RuntimeWarning)
assert os.path.exists(MGN_TGZ), f'MEGANv3.21.tar.gz not found at {MGN_TGZ} (Zenodo 10.5281/zenodo.10526206, the file of notebook 12)'
MEMBERS = {'SPC_NOCONVER.EXT': 'MEGANv3.21/src/MGN2MECH/INCLDIR/SPC_NOCONVER.EXT',
           'MAP_CV2CB6X.EXT': 'MEGANv3.21/src/MGN2MECH/INCLDIR/MAP_CV2CB6X.EXT',
           'SPC_CB6X.EXT': 'MEGANv3.21/src/MGN2MECH/INCLDIR/SPC_CB6X.EXT',
           'mgn2mech.F': 'MEGANv3.21/src/MGN2MECH/mgn2mech.F'}
TXT = {}
with tarfile.open(MGN_TGZ, 'r:gz') as tf:
    for k, m in MEMBERS.items():
        TXT[k] = tf.extractfile(m).read().decode('latin-1')
TABLE_SHA = {k: hashlib.sha256(v.encode('latin-1')).hexdigest()[:16] for k, v in TXT.items()}

def parse_noconver(txt):
    mw = {int(i): (n.strip(), float(m)) for i, n, m in
          re.findall(r"DATA\s+SPCA_SPC\(\s*(\d+)\)\s*,\s*SPCA_MWT\(\s*\d+\)\s*/\s*'([^']*)'\s*,\s*([0-9.Ee+-]+)\s*/", txt)}
    cls = {int(i): n.strip() for i, n in re.findall(r"DATA\s+MGN_SPC\(\s*(\d+)\)\s*/\s*'([^']*)'\s*/", txt)}
    maps = [(n.strip(), int(k), c.strip(), int(ci), float(e)) for n, k, c, ci, e in
            re.findall(r"/\s*'([^']*)'\s*,\s*(\d+)\s*,\s*'([^']*)'\s*,\s*(\d+)\s*,\s*([0-9.Ee+-]+)\s*/", txt)]
    return mw, cls, maps
def parse_mech_map(txt, mech='CB6X'):
    sp = re.findall(r"DATA\s+SPMH_NAM_%s\(\s*(\d+)\)\s*,\s*SPMH_MAP_%s\(\s*\d+\)[\s&]*/\s*'([^']*)'\s*,\s*(\d+)\s*/" % (mech, mech), txt)
    me = re.findall(r"DATA\s+MECH_NAM_%s\(\s*(\d+)\)\s*,\s*MECH_MAP_%s\(\s*\d+\)[\s&]*/\s*'([^']*)'\s*,\s*(\d+)\s*/" % (mech, mech), txt)
    cf = re.findall(r"DATA\s+CONV_FAC_%s\(\s*(\d+)\)[\s&]*/\s*([0-9.Ee+-]+)\s*/" % mech, txt)
    S = {int(i): (n.strip(), int(k)) for i, n, k in sp}; M = {int(i): (n.strip(), int(k)) for i, n, k in me}
    F = {int(i): float(v) for i, v in cf}
    assert set(S) == set(M) == set(F), (len(S), len(M), len(F))
    return [(S[i][0], S[i][1], M[i][0], M[i][1], F[i]) for i in sorted(S)]
def parse_mech_spc(txt, mech='CB6X'):
    return [(n.strip(), float(m)) for _, n, m in
            re.findall(r"DATA\s+MECH_SPC_%s\(\s*(\d+)\)\s*,\s*MECH_MWT_%s\(\s*\d+\)\s*/\s*'([^']*)'\s*,\s*([0-9.Ee+-]+)\s*/" % (mech, mech), txt)]

MWC, CLS_TAB, SMAP = parse_noconver(TXT['SPC_NOCONVER.EXT'])
MMAP = parse_mech_map(TXT['MAP_CV2CB6X.EXT']); MSPC = parse_mech_spc(TXT['SPC_CB6X.EXT'])
CLASSES = [CLS_TAB[i] for i in range(1, len(CLS_TAB) + 1)]
MECH = [n for n, _ in MSPC]; MECH_MW = dict(MSPC)
assert CLASSES == ['ISOP', 'MBO', 'MT_PINE', 'MT_ACYC', 'MT_CAMP', 'MT_SABI', 'MT_AROM', 'NO', 'SQT_HR', 'SQT_LR', 'MEOH',
                   'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO'], CLASSES
assert len(MWC) == 201 and len(SMAP) == 201 and len(MSPC) == 34, (len(MWC), len(SMAP), len(MSPC))
comp_cls, comp_eff = {}, {}
for _n, k, _c, ci, e in SMAP:                      # later lines replace earlier ones (Fortran assignment)
    comp_cls[k] = ci - 1; comp_eff[k] = e
W_CLASS = np.zeros(len(CLASSES))                   # g per mol of class "mixture": sum of fraction x molar mass
for k, ci in comp_cls.items():
    W_CLASS[ci] += comp_eff[k] * MWC[k][1]
K = np.zeros((len(CLASSES), len(MECH)))            # mol of CB6X species per g of class
for _s, k, _m, j, f in MMAP:
    if j == 999 or k not in comp_cls:
        continue
    ci = comp_cls[k]; K[ci, j - 1] += comp_eff[k] / W_CLASS[ci] * f
I_C = {c: i for i, c in enumerate(CLASSES)}; I_M = {m: j for j, m in enumerate(MECH)}
FEEDERS = {m: [CLASSES[i] for i in np.nonzero(K[:, j])[0]] for j, m in enumerate(MECH)}
SOLE = {}                                           # class -> CB6X species produced by that class alone
for m, fe in FEEDERS.items():
    if len(fe) == 1:
        SOLE.setdefault(fe[0], []).append(m)
NEVER = [m for m, fe in FEEDERS.items() if not fe]  # CB6X species no compound is mapped to
GARE = [l.strip() for l in TXT['mgn2mech.F'].splitlines() if re.search(r'^\s*garea\s*=\s*xcell3d\s*\*\s*ycell3d', l, re.I)]
assert GARE, 'mgn2mech.F: the line garea = xcell3d * ycell3d was not found'
print('MEGANv3.21 tables (SHA-256 prefixes):', TABLE_SHA)
print('mgn2mech.F:', GARE[0])
print(f'{len(CLASSES)} classes, {len(MWC)} compounds, {len(MMAP)} map lines, {len(MECH)} CB6X species')
print('CB6X species never produced by the map:', NEVER)
print('classes that are the only source of a CB6X species:', SOLE)
print(pd.DataFrame([{'class': c, 'g/mol': round(W_CLASS[i], 2),
                     'CB6X (mmol per g)': ', '.join(f'{MECH[j]} {1e3 * K[i, j]:.3f}' for j in np.nonzero(K[i])[0])}
                    for i, c in enumerate(CLASSES)]).to_string(index=False))

# %% [markdown]
# ## Cell 3 · The BRAIN file list
#
# The lines of `URL_FILE` whose `path` contains `/V1/MEGAN/` are kept. The day of the year comes from the file name
# (`MEGANv31.SC_2019.CB6X.<YYYYDDD>.BDSNP.zip`). Duplicates and missing days are reported.

# %%
# 3 · BRAIN file list -----------------------------------------------------------------------------------------------
assert os.path.exists(URL_FILE), f'URL list not found: {URL_FILE}'
rows = []
for u in (l.strip() for l in open(URL_FILE, encoding='utf-8', errors='replace')):
    if not u.startswith('http'):
        continue
    q = urllib.parse.parse_qs(urllib.parse.urlparse(u).query)
    path = q.get('path', [''])[0]; name = q.get('fileName', [os.path.basename(path)])[0]
    if '/MEGAN/' not in path:
        continue
    m = re.search(r'\.(\d{4})(\d{3})\.', name)
    if not m or int(m.group(1)) != YEAR:
        continue
    d = dt.date(YEAR, 1, 1) + dt.timedelta(days=int(m.group(2)) - 1)
    rows.append(dict(name=name, url=u, doy=int(m.group(2)), date=d, month=d.month, path=path))
FILES = pd.DataFrame(rows).sort_values('doy').reset_index(drop=True)
assert len(FILES), f'no MEGAN file of {YEAR} in {URL_FILE}'
DUP = FILES[FILES.doy.duplicated()].doy.tolist()
FILES = FILES.drop_duplicates('doy').reset_index(drop=True)
NDAYS_YEAR = 366 if calendar.isleap(YEAR) else 365
MISSING = sorted(set(range(1, NDAYS_YEAR + 1)) - set(FILES.doy))
_md = [str(dt.date(YEAR, 1, 1) + dt.timedelta(days=x - 1)) for x in MISSING]
print(f'{len(FILES)} MEGAN files of {YEAR}; duplicated days: {DUP}; missing days ({len(MISSING)}): '
      f'{_md[:10]}{" ..." if len(_md) > 10 else ""}')

# %% [markdown]
# ## Cell 4 · Reading the BRAIN files and their grid
#
# `fetch` returns a local copy (from `BRAIN_DIR`, from the Colab disk or by download, with retries and a check that the
# download is a zip). `read_day` opens the single netCDF file of a zip, checks the time flags, grid, species list and
# units against the first file, and returns the daily total of each species (mol per cell, as delivered).
#
# The grid comes from the I/O API header: Mercator (`GDTYP` = 7) with latitude of true scale `P_ALP`, central meridian
# `P_GAM`, origin (`XCENT`, `YCENT`), lower-left corner (`XORIG`, `YORIG`) and cell size (`XCELL`, `YCELL`) in projected
# metres, on the sphere named in `FILEDESC` (6,370 km for MCIP). Mercator cells are rectangles in longitude and latitude,
# so their edges and their true areas on the sphere follow exactly.

# %%
# 4 · Readers and grid ------------------------------------------------------------------------------------------------
import requests, netCDF4
from pyproj import CRS, Transformer
from concurrent.futures import ThreadPoolExecutor
for d_ in (OUT_DIR, WORK, os.path.join(WORK, 'zip'), os.path.join(WORK, 'nc')):
    os.makedirs(d_, exist_ok=True)
HDR_KEYS = ['GDTYP', 'P_ALP', 'P_BET', 'P_GAM', 'XCENT', 'YCENT', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 'NCOLS', 'NROWS']

def local_file(name):
    stem = re.sub(r'\.zip$', '', name)
    for d_ in (BRAIN_DIR, os.path.join(WORK, 'zip')):
        for c_ in (name, stem + '.ncf', stem + '.nc', stem):
            p_ = os.path.join(d_, c_)
            if os.path.isfile(p_) and os.path.getsize(p_) > 0:
                return p_
    return None

def fetch(row, tries=6):
    p_ = local_file(row['name'])
    if p_:
        return p_
    dst = os.path.join(WORK, 'zip', row['name']); tmp = dst + '.part'; err = None
    for k in range(tries):
        try:
            with requests.get(row['url'], stream=True, timeout=(30, 600), headers={'User-Agent': 'Mozilla/5.0 (MEGAN-SC notebook 18)'}) as r:
                r.raise_for_status()
                with open(tmp, 'wb') as f:
                    for chunk in r.iter_content(1 << 22):
                        f.write(chunk)
            if not zipfile.is_zipfile(tmp):
                head = open(tmp, 'rb').read(300)
                raise IOError(f'not a zip ({os.path.getsize(tmp)} bytes): {head[:200]!r}')
            os.replace(tmp, dst)
            if KEEP_ZIPS:
                os.makedirs(BRAIN_DIR, exist_ok=True); shutil.copy2(dst, os.path.join(BRAIN_DIR, row['name']))
            return dst
        except Exception as e:
            err = e; time.sleep(min(300, 15 * 2 ** k))
    raise RuntimeError(f"download failed after {tries} tries: {row['name']}: {err}")

def open_nc(path):
    """Path of the netCDF file: the file itself, or the single member of a zip extracted to WORK/nc (temporary)."""
    if not zipfile.is_zipfile(path):
        return path, False
    with zipfile.ZipFile(path) as z:
        mem = [i for i in z.infolist() if not i.is_dir()]
        assert len(mem) == 1, f'{path}: expected one file in the zip, found {[i.filename for i in mem]}'
        out = os.path.join(WORK, 'nc', os.path.basename(mem[0].filename))
        with z.open(mem[0]) as src, open(out, 'wb') as dst:
            shutil.copyfileobj(src, dst, 1 << 24)
    return out, True

def attrs_of(ds):
    a = {k: ds.getncattr(k) for k in ds.ncattrs()}
    a = {k: (v.item() if hasattr(v, 'item') and np.size(v) == 1 else v) for k, v in a.items()}
    a['SPECIES'] = [a['VAR-LIST'][i:i + 16].strip() for i in range(0, len(a['VAR-LIST']), 16) if a['VAR-LIST'][i:i + 16].strip()]
    return a

def earth_radius(filedesc):
    m = re.search(r'EARTH RADIUS ASSUMED IN MCIP:\s*([0-9.]+)\s*m', filedesc)
    return float(m.group(1)) if m else None

def mercator_grid(a):
    assert int(a['GDTYP']) == 7, f"GDTYP {a['GDTYP']}: this notebook handles the Mercator grid (7) of the BRAIN d02 files"
    R = earth_radius(a.get('FILEDESC', ''))
    if R is None:
        R = 6370000.0; print('FILEDESC gives no earth radius; 6,370 km (MCIP default) assumed')
    p4 = f"+proj=merc +lat_ts={a['P_ALP']} +lon_0={a['P_GAM']} +x_0=0 +y_0=0 +R={R} +units=m +no_defs"
    crs = CRS.from_proj4(p4)
    fwd = Transformer.from_crs('EPSG:4326', crs, always_xy=True); inv = Transformer.from_crs(crs, 'EPSG:4326', always_xy=True)
    xc, yc = fwd.transform(a['XCENT'], a['YCENT'])
    xe = a['XORIG'] + a['XCELL'] * np.arange(a['NCOLS'] + 1) + xc
    ye = a['YORIG'] + a['YCELL'] * np.arange(a['NROWS'] + 1) + yc
    lon_e = np.array(inv.transform(xe, np.full(xe.shape, yc))[0]); lat_e = np.array(inv.transform(np.full(ye.shape, xc), ye)[1])
    a_true = R ** 2 * np.deg2rad(lon_e[1] - lon_e[0]) * np.diff(np.sin(np.deg2rad(lat_e)))      # m2, per row (south first)
    return dict(R=R, crs=p4, lon_e=lon_e, lat_e=lat_e, lon_c=(lon_e[:-1] + lon_e[1:]) / 2,
                lat_c=(lat_e[:-1] + lat_e[1:]) / 2, a_true=a_true, a_nom=float(a['XCELL'] * a['YCELL']))

REF = {}
def read_day(path, row):
    """Daily totals (mol per cell, as delivered) of the 34 species, and the QA record of the file."""
    t0 = time.time(); ncp, tmp = open_nc(path)
    try:
        with netCDF4.Dataset(ncp) as ds:
            ds.set_auto_mask(False)
            a = attrs_of(ds)
            if not REF:
                REF.update({k: a[k] for k in HDR_KEYS}); REF['SPECIES'] = a['SPECIES']
            grid_ok = all(np.isclose(float(a[k]), float(REF[k]), rtol=0, atol=1e-6) for k in HDR_KEYS)
            spc_ok = a['SPECIES'] == REF['SPECIES'] and set(a['SPECIES']) == set(MECH)
            units = {v: str(getattr(ds.variables[v], 'units', '')).strip() for v in a['SPECIES']}
            units_ok = all(u == 'mol/s' for u in units.values())
            tf = np.asarray(ds.variables['TFLAG'][:])
            d0 = row['date']; d1 = d0 + dt.timedelta(days=1)
            exp = [(int(d0.strftime('%Y%j')), h * 10000) for h in range(24)] + [(int(d1.strftime('%Y%j')), 0)]
            tflag_ok = tf.shape[0] == 25 and bool((tf == tf[:, :1, :]).all()) and [tuple(map(int, x)) for x in tf[:, 0, :]] == exp
            out = np.zeros((len(MECH), int(a['NROWS']), int(a['NCOLS'])))
            trap = rect = np.nan
            for v in a['SPECIES']:
                x = np.asarray(ds.variables[v][:, 0, :, :], dtype='float64')
                out[I_M[v]] = x[:24].sum(0) * 3600.0
                if v == 'ISOP':
                    rect = x[:24].sum() * 3600.0; trap = (0.5 * x[0].sum() + x[1:24].sum() + 0.5 * x[24].sum()) * 3600.0
    finally:
        if tmp and os.path.exists(ncp):
            os.remove(ncp)
    qa = dict(date=str(row['date']), doy=row['doy'], file=row['name'], steps=int(tf.shape[0]), tflag_ok=tflag_ok, grid_ok=grid_ok,
              species_ok=spc_ok, units_ok=units_ok, isop_trapezoid_over_hourly=trap / rect if rect > 0 else np.nan,
              read_s=round(time.time() - t0, 1))
    return out, qa, a

# grid from the first file (saved, so a resumed run needs no download for it)
GRID_JSON = os.path.join(OUT_DIR, 'brain_grid_header.json')
if os.path.exists(GRID_JSON):
    HDR = json.load(open(GRID_JSON))
else:
    p0 = fetch(FILES.iloc[0]); ncp, tmp = open_nc(p0)
    with netCDF4.Dataset(ncp) as ds:
        HDR = {k: (v if not isinstance(v, (np.generic,)) else v.item()) for k, v in attrs_of(ds).items() if k in HDR_KEYS + ['FILEDESC', 'GDNAM', 'UPNAM', 'SPECIES', 'VAR-LIST']}
    if tmp:
        os.remove(ncp)
    json.dump({k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in HDR.items()}, open(GRID_JSON, 'w'), indent=1, default=str)
REF.update({k: HDR[k] for k in HDR_KEYS}); REF['SPECIES'] = HDR['SPECIES']
G = mercator_grid(HDR)
NR, NC = int(HDR['NROWS']), int(HDR['NCOLS'])
AREA_RATIO = G['a_true'] / G['a_nom']                            # per row: true area / nominal area used by MGN2MECH
print(f"BRAIN grid {HDR['GDNAM'].strip()}: {NC} x {NR} cells, Mercator lat_ts {HDR['P_ALP']} lon_0 {HDR['P_GAM']:.3f}, R {G['R']:.0f} m")
print(f"extent {G['lon_e'][0]:.3f} to {G['lon_e'][-1]:.3f} E, {G['lat_e'][0]:.3f} to {G['lat_e'][-1]:.3f} N")
print('FILEDESC:'); print('\n'.join('   ' + HDR['FILEDESC'][i:i + 80].rstrip() for i in range(0, len(HDR['FILEDESC']), 80) if HDR['FILEDESC'][i:i + 80].strip()))

# %% [markdown]
# ## Cell 5 · Santa Catarina and the mesoregions on both grids
#
# Fractions of each cell inside the polygons, from exact intersections in an Albers equal-area projection centred on the
# state (cell edges densified every 0.002°). The article's cell-centre mask of this inventory is built as in notebook 13.

# %%
# 5 · Polygon fractions -------------------------------------------------------------------------------------------------
import geopandas as gpd, shapely
from rasterio.features import rasterize
from rasterio.transform import from_origin
EA = CRS.from_proj4('+proj=aea +lat_1=-25.5 +lat_2=-29.5 +lat_0=-27.5 +lon_0=-51 +ellps=GRS80 +units=m +no_defs')
TO_EA = Transformer.from_crs('EPSG:4326', EA, always_xy=True)
def _ea(geoms):
    return shapely.transform(geoms, lambda xy: np.column_stack(TO_EA.transform(xy[:, 0], xy[:, 1])))
def _union(g):
    return g.geometry.union_all() if hasattr(g.geometry, 'union_all') else g.geometry.unary_union
sc_g = gpd.read_file(os.path.join(BND, 'sc_state_ibge.gpkg')).to_crs('EPSG:4326'); SC_POLY = _union(sc_g)
meso_g = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs('EPSG:4326')
MESO = list(meso_g['name_meso'])

def box_fractions(lon_e, lat_e, polys):
    """(len(polys), rows, cols) fraction of each lon/lat cell (rows south first) inside each polygon."""
    lo, la = np.meshgrid(lon_e, lat_e)
    boxes = shapely.box(lo[:-1, :-1].ravel(), la[:-1, :-1].ravel(), lo[1:, 1:].ravel(), la[1:, 1:].ravel())
    boxes = _ea(shapely.segmentize(boxes, 0.002)); area = shapely.area(boxes); tree = shapely.STRtree(boxes)
    out = np.zeros((len(polys), (lat_e.size - 1) * (lon_e.size - 1)))
    for k, p in enumerate(polys):
        pe = _ea(shapely.segmentize(p, 0.002)); idx = tree.query(pe, predicate='intersects')
        out[k, idx] = np.clip(shapely.area(shapely.intersection(boxes[idx], pe)) / area[idx], 0, 1)
    return out.reshape(len(polys), lat_e.size - 1, lon_e.size - 1)

FB = box_fractions(G['lon_e'], G['lat_e'], [SC_POLY] + list(meso_g.geometry))
F_SC_B, F_MESO_B = FB[0], FB[1:]
assert G['lon_e'][0] < SC_POLY.bounds[0] and G['lon_e'][-1] > SC_POLY.bounds[2] and G['lat_e'][0] < SC_POLY.bounds[1] \
    and G['lat_e'][-1] > SC_POLY.bounds[3], 'Santa Catarina is not inside the BRAIN grid'
A_SC_B = float((F_SC_B * G['a_true'][:, None]).sum())
print(f'BRAIN grid: {int((F_SC_B > 0).sum())} cells touch Santa Catarina, {int((F_SC_B >= FULL_CELL).sum())} are inside it; '
      f'area inside {A_SC_B / 1e6:,.0f} km2 (true areas); true/nominal cell area over the state '
      f'{AREA_RATIO[(F_SC_B > 0).any(1)].min():.3f} to {AREA_RATIO[(F_SC_B > 0).any(1)].max():.3f}')

# %% [markdown]
# ## Cell 6 · BRAIN, month by month
#
# For each month: the missing local files are downloaded in parallel, each day is read and added, and the month is
# saved as `brain_monthly/brain_<YEAR>_<MM>.npz` (mol per cell and month, as delivered, float32, with the days it holds).
# A month already saved with the same days is not read again. Daily state totals go to `brain_daily_sc.csv`.

# %%
# 6 · BRAIN monthly grids -------------------------------------------------------------------------------------------
MDIR = os.path.join(OUT_DIR, 'brain_monthly'); os.makedirs(MDIR, exist_ok=True)
W_TRUE = F_SC_B * AREA_RATIO[:, None]                 # weights giving the emission of the true area inside the state
W_DELV = F_SC_B                                       # weights giving the delivered rates inside the state
QA, DAILY = [], []
for mth in range(1, 13):
    sel = FILES[FILES.month == mth]
    fn = os.path.join(MDIR, f'brain_{YEAR}_{mth:02d}.npz')
    if os.path.exists(fn):
        z = np.load(fn, allow_pickle=True)
        if sorted(z['doys'].tolist()) == sorted(sel.doy.tolist()):
            QA += json.loads(str(z['qa'])); DAILY += json.loads(str(z['daily'])); print(f'{YEAR}-{mth:02d}: saved ({len(sel)} days)'); continue
    if not len(sel):
        print(f'{YEAR}-{mth:02d}: no file'); continue
    t0 = time.time(); acc = np.zeros((len(MECH), NR, NC)); qa_m, daily_m = [], []
    recs = sel.to_dict('records')
    with ThreadPoolExecutor(N_THREADS) as ex:
        futs = [ex.submit(fetch, r) for r in recs]
        for r, fu in zip(recs, futs):
            p_ = fu.result()
            day, qa, _ = read_day(p_, r)
            acc += day; qa_m.append(qa)
            daily_m.append(dict(date=str(r['date']), **{f'{s}_true_mol': float((day[I_M[s]] * W_TRUE).sum()) for s in MECH},
                                **{f'{s}_delivered_mol': float((day[I_M[s]] * W_DELV).sum()) for s in MECH}))
            if p_.startswith(os.path.join(WORK, 'zip')):
                os.remove(p_)
    np.savez_compressed(fn, mol=acc.astype('float32'), doys=np.array(sel.doy.tolist()), qa=json.dumps(qa_m, default=str),
                        daily=json.dumps(daily_m))
    QA += qa_m; DAILY += daily_m
    print(f'{YEAR}-{mth:02d}: {len(sel)} days in {time.time() - t0:.0f} s')
QA = pd.DataFrame(QA); QA.to_csv(os.path.join(OUT_DIR, 'brain_files_qa.csv'), index=False)
DAILY = pd.DataFrame(DAILY).sort_values('date'); DAILY.to_csv(os.path.join(OUT_DIR, 'brain_daily_sc.csv'), index=False)

# monthly grids, scale factors for missing days, months compared
BR = np.zeros((12, len(MECH), NR, NC)); MONTHS = []
for mth in range(1, 13):
    nd = calendar.monthrange(YEAR, mth)[1]; fn = os.path.join(MDIR, f'brain_{YEAR}_{mth:02d}.npz')
    nav = int((FILES.month == mth).sum()) if os.path.exists(fn) else 0
    use = nav >= MIN_DAYS
    if os.path.exists(fn) and nav:
        BR[mth - 1] = np.load(fn)['mol'].astype('float64') * (nd / nav)
    MONTHS.append(dict(month=mth, days_in_month=nd, days_available=nav, scale=nd / nav if nav else np.nan, compared=use))
MONTHS = pd.DataFrame(MONTHS); MONTHS.to_csv(os.path.join(OUT_DIR, 'brain_months.csv'), index=False)
MSEL = MONTHS.compared.values
print(MONTHS.to_string(index=False))
assert MSEL.any(), f'no month with at least MIN_DAYS = {MIN_DAYS} BRAIN days'

# %% [markdown]
# ## Cell 7 · Checks of the BRAIN files and of their speciation
#
# The file checks must all be true. Species that the CB6X map never produces must be zero. For each class that is the
# only source of two or more CB6X species, the molar ratios of those species in BRAIN must equal the ratios of `K`.

# %%
# 7 · BRAIN checks ----------------------------------------------------------------------------------------------------
bad = QA[~(QA.tflag_ok & QA.grid_ok & QA.species_ok & QA.units_ok)]
print(f'files read: {len(QA)}; files failing a check: {len(bad)}'); print(bad.to_string(index=False) if len(bad) else '')
print(f'daily isoprene, trapezoid over 25 steps / sum of 24 hourly rates: {QA.isop_trapezoid_over_hourly.min():.4f} to '
      f'{QA.isop_trapezoid_over_hourly.max():.4f} (median {QA.isop_trapezoid_over_hourly.median():.4f})')
ANN_B = BR[MSEL].sum(0)                                          # mol per cell, as delivered, months compared
TOT_B_DELV = {s: float((ANN_B[I_M[s]] * W_DELV).sum()) for s in MECH}
TOT_B_TRUE = {s: float((ANN_B[I_M[s]] * W_TRUE).sum()) for s in MECH}
chk = []
for s in NEVER:
    chk.append(dict(check='species never produced by MAP_CV2CB6X', item=s, expected=0.0, brain=float(ANN_B[I_M[s]].sum()),
                    ok=float(ANN_B[I_M[s]].sum()) == 0.0))
for c, sp in SOLE.items():
    for s in sp[1:]:
        exp = K[I_C[c], I_M[s]] / K[I_C[c], I_M[sp[0]]]
        got = TOT_B_DELV[s] / TOT_B_DELV[sp[0]] if TOT_B_DELV[sp[0]] > 0 else np.nan
        chk.append(dict(check=f'fingerprint of class {c}', item=f'{s}/{sp[0]}', expected=exp, brain=got,
                        ok=bool(np.isfinite(got) and abs(got / exp - 1) < 1e-3)))
CHK = pd.DataFrame(chk); CHK.to_csv(os.path.join(OUT_DIR, 'brain_speciation_checks.csv'), index=False)
print(CHK.to_string(index=False))
if not CHK.ok.all():
    print('WARNING: BRAIN does not follow the MEGANv3.21 CB6X speciation for the rows marked False; the class-level '
          'comparison (cell 10) is not exact for those classes.')

# %% [markdown]
# ## Cell 8 · This inventory, in CB6X species and on the BRAIN grid
#
# The monthly class grids of notebook 10 (tonnes per cell) are converted to CB6X moles with `K` (ρ = 1, as written by
# the code; the published isoprene is 0.96 of it). The conservative remapping uses the overlap of each 0.025° cell with
# each BRAIN cell: in longitude, the overlapping fraction of the cell width; in latitude, the overlapping fraction of
# sin(latitude), which is proportional to area on the sphere. Every 0.025° cell is fully inside the BRAIN grid, so the
# total of each species is preserved.

# %%
# 8 · This inventory ----------------------------------------------------------------------------------------------------
def read_ours(y):
    out = np.zeros((12, len(CLASSES), 0, 0)); lat = lon = None
    for mth in range(1, 13):
        fn = os.path.join(ROOT, str(y), 'emissions', f'emis_{y}_{mth:02d}.nc'); assert os.path.exists(fn), fn
        with netCDF4.Dataset(fn) as ds:
            ds.set_auto_mask(False)
            if lat is None:
                lat = np.asarray(ds.variables['lat'][:], 'float64'); lon = np.asarray(ds.variables['lon'][:], 'float64')
                out = np.zeros((12, len(CLASSES), lat.size, lon.size))
            for i, c in enumerate(CLASSES):
                a = np.asarray(ds.variables[c][:], dtype='float64'); out[mth - 1, i] = np.where(np.isfinite(a) & (a > -9990), a, 0.0)
    return out, lat, lon
OT, LAT, LON = read_ours(YEAR)                                   # tonnes per cell and month, 19 classes
assert np.all(np.diff(LAT) > 0), 'expected south-first latitudes'
RES = float(LON[1] - LON[0]); NY, NX = LAT.size, LON.size
OM = np.einsum('ci,mcyx->miyx', K, OT * 1e6)                     # mol per cell and month, 34 CB6X species (rho = 1)

# masks of this inventory: cell centre (article) and area fractions
TR = from_origin(LON[0] - RES / 2, LAT[-1] + RES / 2, RES, RES)
SC_CC = np.flipud(rasterize([(SC_POLY, 1)], out_shape=(NY, NX), transform=TR, fill=0, dtype='int32')) == 1
lon_eo = np.r_[LON - RES / 2, LON[-1] + RES / 2]; lat_eo = np.r_[LAT - RES / 2, LAT[-1] + RES / 2]
F_SC_O = box_fractions(lon_eo, lat_eo, [SC_POLY])[0]
stats_f = os.path.join(ROOT, str(YEAR), 'stats', f'annual_by_family_{YEAR}.csv')
if os.path.exists(stats_f):
    pub = pd.read_csv(stats_f).set_index('species')['annual_tonnes_SC']
    mine = float((OT[:, I_C['ISOP']].sum(0) * SC_CC).sum())
    print(f"isoprene {YEAR}, cell-centre mask, rho = 1: {mine:,.1f} t; notebook 10 stats: {pub['ISOP']:,.1f} t "
          f"(difference {100 * (mine / pub['ISOP'] - 1):+.4f} %)")

# conservative remapping onto the BRAIN cells
def overlap(a_lo, a_hi, b_lo, b_hi):
    lo = np.maximum(a_lo[:, None], b_lo[None, :]); hi = np.minimum(a_hi[:, None], b_hi[None, :])
    return np.clip(hi - lo, 0, None) / (a_hi - a_lo)[:, None]
s_ = lambda d: np.sin(np.deg2rad(d))
WY = overlap(s_(lat_eo[:-1]), s_(lat_eo[1:]), s_(G['lat_e'][:-1]), s_(G['lat_e'][1:]))      # (NY, NR)
WX = overlap(lon_eo[:-1], lon_eo[1:], G['lon_e'][:-1], G['lon_e'][1:])                      # (NX, NC)
assert np.allclose(WY.sum(1), 1) and np.allclose(WX.sum(1), 1), 'the 0.025 deg grid is not inside the BRAIN grid'
remap = lambda a: WY.T @ a @ WX
A_O = (G['R'] ** 2 * np.deg2rad(RES) * np.diff(s_(lat_eo)))[:, None] * np.ones((1, NX))
COVER = remap(A_O) / G['a_true'][:, None]                        # share of each BRAIN cell covered by the 0.025 deg grid
OB = np.stack([[remap(OM[mth, j]) for j in range(len(MECH))] for mth in range(12)])        # mol per BRAIN cell and month
GF = {}
gfz = sorted(glob.glob(os.path.join(ROOT, str(YEAR), 'inputs', 'growthform_*.npz')))
if gfz:
    z = np.load(gfz[0]); den = remap(A_O)
    GF = {k: np.where(den > 0, remap(z[k] * A_O) / den, np.nan) for k in ('tree', 'herb', 'crop')}
ANN_O = OM[MSEL].sum(0); ANN_OB = OB[MSEL].sum(0)
TOT_O_AW = {s: float((ANN_O[I_M[s]] * F_SC_O).sum()) for s in MECH}
TOT_O_CC = {s: float((ANN_O[I_M[s]] * SC_CC).sum()) for s in MECH}
TOT_OB = {s: float((ANN_OB[I_M[s]] * F_SC_B).sum()) for s in MECH}
print(f"isoprene, months compared: area-weighted {TOT_O_AW['ISOP'] * 68.12e-6:,.0f} t, cell-centre {TOT_O_CC['ISOP'] * 68.12e-6:,.0f} t, "
      f"on the BRAIN grid {TOT_OB['ISOP'] * 68.12e-6:,.0f} t (rho = 1)")
print(f'BRAIN cells inside the state with the 0.025 deg grid covering < 99.9 %: '
      f'{int(((F_SC_B > 0) & (COVER < 0.999)).sum())} (state area share {float((F_SC_B * G["a_true"][:, None])[(COVER < 0.999)].sum()) / A_SC_B:.2e})')

# %% [markdown]
# ## Cell 9 · Co-registration of the two grids
#
# Annual isoprene fluxes on the BRAIN cells covered by the 0.025° grid. The BRAIN field is shifted by up to ±4 cells in
# each direction; the correct georeferencing gives the highest correlation without shift. For each row, the eastern-most
# emitting column (the coast) is compared: a BRAIN cell counts as emitting when its isoprene is above zero, and for this
# inventory when more than half of the cell area is covered by emitting 0.025° cells (the remapping spreads coastal
# emission over partly covered cells). A failure stops the notebook, because every later result depends on it.

# %%
# 9 · Co-registration ---------------------------------------------------------------------------------------------------
fo = ANN_OB[I_M['ISOP']] / G['a_true'][:, None]; fb = ANN_B[I_M['ISOP']] / G['a_nom']
inside = COVER > 0.999
LANDF = remap((ANN_O[I_M['ISOP']] > 0) * A_O) / G['a_true'][:, None]     # share of each BRAIN cell covered by emitting 0.025 deg cells
res = []
for dy in range(-4, 5):
    for dx in range(-4, 5):
        sb = np.roll(np.roll(fb, dy, 0), dx, 1); m_ = inside & ((fo > 0) | (sb > 0))
        res.append(dict(dy=dy, dx=dx, r=float(np.corrcoef(fo[m_], sb[m_])[0, 1])))
SHIFT = pd.DataFrame(res).sort_values('r', ascending=False).reset_index(drop=True)
off = []
for r_ in range(NR):
    cols = np.nonzero(inside[r_])[0]
    if cols.size < 20:
        continue
    eo, eb = np.nonzero(LANDF[r_] > 0.5)[0], np.nonzero((fb[r_] > 0) & inside[r_])[0]
    if eo.size and eb.size and eo.max() < cols.max() - 2:
        off.append(int(eb.max() - eo.max()))
OFF = np.array(off)
print(SHIFT.head(5).to_string(index=False))
print(f'east coast, BRAIN minus this inventory (columns): median {np.median(OFF):.0f}, IQR {np.percentile(OFF, 25):.0f} to '
      f'{np.percentile(OFF, 75):.0f}, {OFF.size} rows')
assert SHIFT.loc[0, 'dy'] == 0 and SHIFT.loc[0, 'dx'] == 0 and abs(np.median(OFF)) <= 1, \
    'the BRAIN grid does not line up with this inventory: check the grid header (cell 4) before using any result'

# %% [markdown]
# ## Cell 10 · Annual totals in Santa Catarina
#
# Table A: the 34 CB6X species, in moles. Table B: tonnes by MEGAN class where a class is the only source of a CB6X
# species (the BRAIN class mass is then the moles of that species divided by `K`), and the CB6X species `TERP` and
# `SQT` for monoterpenes and sesquiterpenes (mass with the CB6X molar masses, for both datasets). Ratios are this
# inventory over BRAIN. The previous state inventory (Hoinaski et al. 2022) is added for isoprene, CO and soil NO; its
# values are annual, so they are compared only when all 12 months are compared.

# %%
# 10 · Annual totals ----------------------------------------------------------------------------------------------------
NCMP = int(MSEL.sum())
TA = pd.DataFrame([dict(species=s, sources=', '.join(FEEDERS[s]) or '(none)', brain_true_mol=TOT_B_TRUE[s],
                        brain_delivered_mol=TOT_B_DELV[s], this_inventory_mol=TOT_O_AW[s], this_inventory_cellcentre_mol=TOT_O_CC[s])
                   for s in MECH])
TA['ratio_vs_true'] = TA.this_inventory_mol / TA.brain_true_mol
TA['ratio_vs_delivered'] = TA.this_inventory_mol / TA.brain_delivered_mol
TA.to_csv(os.path.join(OUT_DIR, f'comparison_{YEAR}_cb6x_species.csv'), index=False)

rows = []
KEY_MASS = [('ISOP', 'ISOP'), ('MT (as TERP)', 'TERP'), ('SQT (as SQT)', 'SQT')]
for label, s in KEY_MASS:
    mw = MECH_MW[s]
    rows.append(dict(item=label, basis=f'CB6X {s} x {mw} g/mol', brain_true_t=TOT_B_TRUE[s] * mw * 1e-6,
                     brain_delivered_t=TOT_B_DELV[s] * mw * 1e-6, this_inventory_t=TOT_O_AW[s] * mw * 1e-6,
                     this_inventory_cellcentre_t=TOT_O_CC[s] * mw * 1e-6))
for c in ['MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER', 'CO', 'NO']:
    s = SOLE[c][0]; k = K[I_C[c], I_M[s]]
    rows.append(dict(item=c, basis=f'class mass from CB6X {s}', brain_true_t=TOT_B_TRUE[s] / k * 1e-6,
                     brain_delivered_t=TOT_B_DELV[s] / k * 1e-6,
                     this_inventory_t=float((OT[MSEL, I_C[c]].sum(0) * F_SC_O).sum()),
                     this_inventory_cellcentre_t=float((OT[MSEL, I_C[c]].sum(0) * SC_CC).sum())))
TB = pd.DataFrame(rows)
cmp_ = TB.item.isin(['ISOP', 'MT (as TERP)', 'SQT (as SQT)', 'MEOH', 'ACTO', 'ETOH', 'ACID', 'LVOC', 'OXPROD', 'STRESS', 'OTHER'])
TB = pd.concat([TB, pd.DataFrame([dict(item='VOC compared (sum of the rows above, MBO excluded)', basis='sum',
                                       **{k: TB.loc[cmp_, k].sum() for k in ['brain_true_t', 'brain_delivered_t', 'this_inventory_t',
                                                                             'this_inventory_cellcentre_t']})])], ignore_index=True)
TB['this_inventory_published_t'] = TB.this_inventory_t * np.where(TB.item == 'ISOP', RHO_ISOP, 1.0)
TB['ratio_vs_true'] = TB.this_inventory_t / TB.brain_true_t
TB['ratio_vs_delivered'] = TB.this_inventory_t / TB.brain_delivered_t
prev = {'ISOP': PREV_2019_KG['ISOP'] / 1e3, 'CO': PREV_2019_KG['CO'] / 1e3, 'NO': PREV_2019_KG['NO'] / 1e3}
TB['previous_inventory_t'] = TB.item.map(prev) if (NCMP == 12 and YEAR == 2019) else np.nan
TB['this_inventory_published_over_previous'] = TB.this_inventory_published_t / TB.previous_inventory_t
TB['brain_true_over_previous'] = TB.brain_true_t / TB.previous_inventory_t
TB['brain_delivered_over_previous'] = TB.brain_delivered_t / TB.previous_inventory_t
TB.to_csv(os.path.join(OUT_DIR, f'comparison_{YEAR}_classes_t.csv'), index=False)
pd.set_option('display.width', 250)
print(f'months compared: {NCMP} of 12 ({", ".join(str(m) for m in MONTHS.month[MSEL])})')
print(TB[['item', 'brain_true_t', 'brain_delivered_t', 'this_inventory_t', 'this_inventory_published_t', 'ratio_vs_true',
          'ratio_vs_delivered', 'previous_inventory_t']].round(3).to_string(index=False))

# %% [markdown]
# ## Cell 11 · Monthly totals and the daily BRAIN series

# %%
# 11 · Monthly -----------------------------------------------------------------------------------------------------------
KEY = [('ISOP', 'ISOP', MECH_MW['ISOP']), ('MT (as TERP)', 'TERP', MECH_MW['TERP']), ('SQT (as SQT)', 'SQT', MECH_MW['SQT']),
       ('MEOH', 'MEOH', MECH_MW['MEOH']), ('CO', 'CO', MECH_MW['CO']), ('NO', 'NO', W_CLASS[I_C['NO']])]
mrows = []
for mth in range(12):
    for label, s, mw in KEY:
        mrows.append(dict(month=mth + 1, item=label, compared=bool(MSEL[mth]),
                          brain_true_t=float((BR[mth, I_M[s]] * W_TRUE).sum()) * mw * 1e-6,
                          brain_delivered_t=float((BR[mth, I_M[s]] * W_DELV).sum()) * mw * 1e-6,
                          this_inventory_t=float((OM[mth, I_M[s]] * F_SC_O).sum()) * mw * 1e-6))
TM = pd.DataFrame(mrows); TM['ratio_vs_true'] = TM.this_inventory_t / TM.brain_true_t
TM.to_csv(os.path.join(OUT_DIR, f'comparison_{YEAR}_monthly_t.csv'), index=False)
for label, s, mw in KEY:
    q = TM[(TM.item == label) & TM.compared]
    r_ = np.corrcoef(q.this_inventory_t, q.brain_true_t)[0, 1] if len(q) > 2 else np.nan
    print(f'{label:14s} monthly ratio (this / BRAIN true area): {q.ratio_vs_true.min():.2f} to {q.ratio_vs_true.max():.2f}; '
          f'r of the monthly cycles {r_:.2f} (n = {len(q)})')

# %% [markdown]
# ## Cell 12 · Mesoregions
#
# Both datasets on the BRAIN grid, weighted by the fraction of each cell in each mesoregion (the 2024 municipal unions
# of the article). The share of the state total shows whether the spatial distribution agrees independently of the
# magnitude.

# %%
# 12 · Mesoregions ------------------------------------------------------------------------------------------------------
rr = []
for k, nm in enumerate(MESO):
    for label, s, mw in KEY:
        b = float((ANN_B[I_M[s]] * F_MESO_B[k] * AREA_RATIO[:, None]).sum()) * mw * 1e-6
        o = float((ANN_OB[I_M[s]] * F_MESO_B[k]).sum()) * mw * 1e-6
        rr.append(dict(mesoregion=nm, item=label, brain_true_t=b, this_inventory_t=o))
TR_ = pd.DataFrame(rr)
TR_['brain_share_pct'] = 100 * TR_.brain_true_t / TR_.groupby('item').brain_true_t.transform('sum')
TR_['this_share_pct'] = 100 * TR_.this_inventory_t / TR_.groupby('item').this_inventory_t.transform('sum')
TR_['ratio'] = TR_.this_inventory_t / TR_.brain_true_t
TR_.to_csv(os.path.join(OUT_DIR, f'comparison_{YEAR}_mesoregions_t.csv'), index=False)
print(TR_[TR_.item.isin(['ISOP', 'MT (as TERP)', 'NO'])].round(2).to_string(index=False))

# %% [markdown]
# ## Cell 13 · Cell statistics on the 4-km grid
#
# Cells at least `FULL_CELL` inside Santa Catarina and fully covered by the 0.025° grid. Fluxes per true area (this
# inventory: remapped tonnes / true cell area; BRAIN: delivered rate / nominal area, which is the MEGAN flux). Statistics:
# Pearson r of the fluxes and of their logarithms, Spearman ρ, the median and interquartile range of the ratio this /
# BRAIN, and the share of cells within a factor of 2. Strata: cells dominated (> `DOMINANT`) by trees, herbaceous
# vegetation or crops in this inventory's 2019 growth forms, and the remaining cells.

# %%
# 13 · Cell statistics -----------------------------------------------------------------------------------------------------
from scipy import stats as sst
CELL = (F_SC_B >= FULL_CELL) & (COVER > 0.999)
STRATA = {'all': CELL}
if GF:
    dom = {k: np.nan_to_num(GF[k]) > DOMINANT for k in GF}
    STRATA.update({f'{k}-dominated': CELL & dom[k] for k in dom})
    STRATA['other cells'] = CELL & ~np.any(np.stack(list(dom.values())), 0)
FLUX = {}
crow = []
for label, s, mw in KEY:
    fo_ = ANN_OB[I_M[s]] * mw * 1e-3 / (G['a_true'][:, None] * 1e-6)     # kg km-2 over the months compared
    fb_ = ANN_B[I_M[s]] * mw * 1e-3 / (G['a_nom'] * 1e-6)
    FLUX[label] = (fo_, fb_)
    for st, msk in STRATA.items():
        a, b = fo_[msk], fb_[msk]; pos = (a > 0) & (b > 0)
        if pos.sum() < 10:
            continue
        rt = a[pos] / b[pos]
        crow.append(dict(item=label, stratum=st, cells=int(msk.sum()), r=float(np.corrcoef(a, b)[0, 1]),
                         r_log=float(np.corrcoef(np.log10(a[pos]), np.log10(b[pos]))[0, 1]),
                         spearman=float(sst.spearmanr(a, b)[0]), ratio_median=float(np.median(rt)),
                         ratio_p25=float(np.percentile(rt, 25)), ratio_p75=float(np.percentile(rt, 75)),
                         within_factor2_pct=float(100 * np.mean((rt >= 0.5) & (rt <= 2))),
                         this_mean_kg_km2=float(a.mean()), brain_mean_kg_km2=float(b.mean())))
TC = pd.DataFrame(crow); TC.to_csv(os.path.join(OUT_DIR, f'comparison_{YEAR}_cells.csv'), index=False)
print(TC.round(3).to_string(index=False))

# %% [markdown]
# ## Cell 14 · Figures

# %%
# 14 · Figures ---------------------------------------------------------------------------------------------------------------
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
def outline(ax):
    for p in getattr(SC_POLY, 'geoms', [SC_POLY]):
        x, y = p.exterior.xy; ax.plot(x, y, color='k', lw=0.5)
    ax.set_xlim(SC_POLY.bounds[0] - 0.2, SC_POLY.bounds[2] + 0.2); ax.set_ylim(SC_POLY.bounds[1] - 0.2, SC_POLY.bounds[3] + 0.2)
    ax.set_aspect(1 / np.cos(np.deg2rad(27.5)))
if MAKE_FIGS:
    mask_sc = np.where(F_SC_B > 0.5, 1.0, np.nan)
    # (a) co-registration
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    piv = SHIFT.pivot(index='dy', columns='dx', values='r')
    im = ax[0].imshow(piv.values, origin='lower', extent=[-4.5, 4.5, -4.5, 4.5], cmap='viridis'); fig.colorbar(im, ax=ax[0], label='r')
    ax[0].set_xlabel('BRAIN shift, columns'); ax[0].set_ylabel('BRAIN shift, rows'); ax[0].set_title('Correlation of isoprene fluxes')
    ax[1].hist(OFF, bins=np.arange(OFF.min() - 0.5, OFF.max() + 1.5), color='0.4'); ax[1].set_xlabel('east coast, BRAIN minus this inventory (columns)')
    ax[1].set_ylabel('rows'); ax[1].set_title('Coastline position')
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, 'fig_brain_coregistration.png'), dpi=200); plt.show()
    # (b) maps
    items = ['ISOP', 'MT (as TERP)', 'NO']
    fig, axs = plt.subplots(len(items), 3, figsize=(13, 3.6 * len(items)))
    for i, it in enumerate(items):
        fo_, fb_ = FLUX[it]; a_, b_ = fo_ * mask_sc, fb_ * mask_sc
        v = np.concatenate([a_[np.isfinite(a_) & (a_ > 0)], b_[np.isfinite(b_) & (b_ > 0)]])
        nrm = LogNorm(vmin=np.percentile(v, 2), vmax=np.percentile(v, 98))
        for j, (arr, ttl) in enumerate([(a_, 'this inventory'), (b_, 'BRAIN')]):
            pc = axs[i, j].pcolormesh(G['lon_e'], G['lat_e'], np.where(arr > 0, arr, np.nan), norm=nrm, cmap='YlGn'); outline(axs[i, j])
            axs[i, j].set_title(f'{it}, {ttl} (kg km$^{{-2}}$)', fontsize=9); fig.colorbar(pc, ax=axs[i, j], shrink=0.8)
        lr = np.log2(np.where((a_ > 0) & (b_ > 0), a_ / b_, np.nan))
        pc = axs[i, 2].pcolormesh(G['lon_e'], G['lat_e'], lr, norm=TwoSlopeNorm(0, -3, 3), cmap='RdBu_r'); outline(axs[i, 2])
        axs[i, 2].set_title(f'{it}, log$_2$(this / BRAIN)', fontsize=9); fig.colorbar(pc, ax=axs[i, 2], shrink=0.8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, 'fig_brain_maps.png'), dpi=200); plt.show()
    # (c) monthly cycle
    fig, axs = plt.subplots(2, 3, figsize=(13, 6.5))
    for ax, (label, s, mw) in zip(axs.ravel(), KEY):
        q = TM[TM.item == label]
        ax.plot(q.month, q.this_inventory_t / 1e3, 'o-', color='tab:green', label='this inventory')
        ax.plot(q.month, q.brain_true_t / 1e3, 's-', color='tab:blue', label='BRAIN, true area')
        ax.plot(q.month, q.brain_delivered_t / 1e3, 's--', color='tab:blue', alpha=0.5, label='BRAIN, as delivered')
        ax.set_title(label); ax.set_xticks(range(1, 13)); ax.set_ylabel('kt per month')
    axs[0, 0].legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, 'fig_brain_monthly.png'), dpi=200); plt.show()
    # (d) cell scatter
    fig, axs = plt.subplots(2, 3, figsize=(13, 8))
    for ax, (label, s, mw) in zip(axs.ravel(), KEY):
        fo_, fb_ = FLUX[label]; a, b = fo_[CELL], fb_[CELL]; pos = (a > 0) & (b > 0)
        if pos.sum() > 10:
            ax.hexbin(np.log10(b[pos]), np.log10(a[pos]), gridsize=50, mincnt=1, cmap='Greys', bins='log')
            lo, hi = np.percentile(np.r_[np.log10(a[pos]), np.log10(b[pos])], [0.5, 99.5])
            for k_, ls in ((0, '-'), (np.log10(2), ':'), (-np.log10(2), ':')):
                ax.plot([lo, hi], [lo + k_, hi + k_], 'r' + ls, lw=0.8)
        ax.set_xlabel('log$_{10}$ BRAIN (kg km$^{-2}$)'); ax.set_ylabel('log$_{10}$ this inventory'); ax.set_title(label)
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, 'fig_brain_cells.png'), dpi=200); plt.show()
    # (e) daily BRAIN isoprene and the monthly mean daily rate of this inventory
    fig, ax = plt.subplots(figsize=(11, 3.6)); dd = pd.to_datetime(DAILY.date)
    ax.plot(dd, DAILY.ISOP_true_mol * 68.12e-9, lw=0.8, color='tab:blue', label='BRAIN, true area (daily)')
    q = TM[TM.item == 'ISOP']
    ax.step(pd.to_datetime([dt.date(YEAR, m_, 1) for m_ in q.month]), q.this_inventory_t.values / 1e3 / MONTHS.days_in_month.values,
            where='post', color='tab:green', label='this inventory (monthly mean of the day)')
    ax.set_ylabel('isoprene, kt per day'); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'fig_brain_daily_isoprene.png'), dpi=200); plt.show()
    # (f) ratios of the annual totals
    q = TB[TB.item.str.startswith('VOC compared') == False]
    fig, ax = plt.subplots(figsize=(9, 4.2)); y_ = np.arange(len(q))
    ax.barh(y_ - 0.2, q.ratio_vs_true, height=0.4, color='tab:blue', label='vs BRAIN, true area')
    ax.barh(y_ + 0.2, q.ratio_vs_delivered, height=0.4, color='tab:blue', alpha=0.45, label='vs BRAIN, as delivered')
    ax.set_xscale('log'); ax.axvline(1, color='k', lw=0.8); ax.set_yticks(y_); ax.set_yticklabels(q.item, fontsize=8)
    ax.set_xlabel('this inventory / BRAIN, Santa Catarina'); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'fig_brain_ratios.png'), dpi=200); plt.show()

# %% [markdown]
# ## Cell 15 · Parametrizations, provenance and summary
#
# The table records what differs between the two inventories and where each entry comes from: the BRAIN file header
# (read here), the BRAIN article (Hoinaski et al. 2024), or the methods of this inventory (article Sections 2.2-2.4 and
# notebook 10).

# %%
# 15 · Parametrizations, provenance, summary ---------------------------------------------------------------------------
fd = HDR['FILEDESC']; fdl = [fd[i:i + 80].strip() for i in range(0, len(fd), 80)]
def fd_get(key):
    hit = [l for l in fdl if l.upper().startswith(key)]
    return hit[0].split(':', 1)[1].strip() if hit else 'not in FILEDESC'
P = [
 ('Model code', 'MEGANv3.21, Python port checked against the Fortran (notebook 12)',
  f"MEGAN 3.2 (article); MGN2MECH of MEGANv3.21 (FILEDESC: {fdl[1][:60]}...)", 'article Sect. 2.2; BRAIN header and Hoinaski et al. 2024'),
 ('Grid', f'0.025 deg lon-lat ({NX} x {NY} cells)',
  f"{HDR['GDNAM'].strip()}: Mercator, true at {HDR['P_ALP']} deg, {HDR['XCELL']:.0f} m projected cells ({NC} x {NR}); sphere R = {G['R']:.0f} m",
  'notebook 10; BRAIN header'),
 ('Cell area for per-cell rates', 'true area on the sphere (megan_sc.domain)',
  f"nominal XCELL x YCELL = {G['a_nom'] / 1e6:.0f} km2 (mgn2mech.F: {GARE[0]})", 'notebook 10; mgn2mech.F of MEGANv3.21'),
 ('Meteorology', 'ERA5-Land, 0.1 deg, hourly fields sampled every 3 h (00-21 UTC)',
  f"WRF ARW 4.0 via MCIP 5.1 (FILEDESC), 4 km; GFS 0.25 deg, 36-h runs, first 12 h discarded (article); "
  f"PBL {fd_get('PBL SCHEME')}; land surface {fd_get('LAND-SURFACE SCHEME')}; microphysics {fd_get('MICROPHYSICS')}",
  'article Sect. 2.3.4; BRAIN header and Hoinaski et al. 2024'),
 ('Land use of the meteorology', 'not applicable (ERA5-Land)', fd_get('LAND USE CLASSIFICATION'), 'BRAIN header'),
 ('Vegetation for MEGAN', 'MapBiomas Collection 10 (30 m) growth forms of the year; MEGAN v3 ecotype map and emission factors of the MEGAN 3.2 input files',
  'not named in the article (MEGAN input data from bai.ess.uci.edu/megan/data-and-code)', 'article Sect. 2.3.3; Hoinaski et al. 2024'),
 ('Leaf area index', 'MODIS MCD15A3H C6.1, monthly means of main-algorithm retrievals on the 0.025 deg grid', 'not named in the article',
  'article Sect. 2.3.2; Hoinaski et al. 2024'),
 ('Soil NO', 'MEGANv3.21 default (activity factor GAMNO, Yienger and Levy 1995 type), growing season set for the Southern Hemisphere',
  'BDSNP (Hudman et al. 2012), from the file name', 'notebook 10; BRAIN file name'),
 ('Chemical speciation', '19 MEGAN classes in mass (converted here to CB6X with the MEGANv3.21 tables)', '34 CB6X species in mol/s',
  'notebook 10; BRAIN header'),
 ('Isoprene canopy factor rho', f'1 in the code; {RHO_ISOP} in the published numbers', '1 (not applied by MEGANv3.21)', 'article Sect. 2.2'),
 ('Time', 'UTC; monthly totals of 3-hourly steps', 'UTC; 25 hourly steps per daily file', 'notebook 10; BRAIN header'),
 ('Year', str(YEAR), str(YEAR), ''),
]
PT = pd.DataFrame(P, columns=['item', 'this inventory', 'BRAIN d02', 'source'])
PT.to_csv(os.path.join(OUT_DIR, 'parametrizations.csv'), index=False)
prov = dict(notebook='18_brain_comparison_2019', run_utc=dt.datetime.utcnow().isoformat(timespec='seconds'), year=YEAR,
            url_file=URL_FILE, files_listed=int(len(FILES)), missing_days=MISSING, duplicated_days=DUP, files_read=int(len(QA)),
            files_failing=int(len(bad)), months_compared=[int(m) for m in MONTHS.month[MSEL]], min_days=MIN_DAYS,
            megan_tables_sha256=TABLE_SHA, megan_tgz=MGN_TGZ, grid={k: HDR[k] for k in HDR_KEYS}, earth_radius_m=G['R'],
            rho_isop=RHO_ISOP, dominant=DOMINANT, full_cell=FULL_CELL, root=ROOT,
            versions=dict(numpy=np.__version__, pandas=pd.__version__, netCDF4=netCDF4.__version__, shapely=shapely.__version__))
json.dump(prov, open(os.path.join(OUT_DIR, 'provenance_brain_comparison.json'), 'w'), indent=1, default=str)
print(PT.to_string(index=False))
print('\nSummary (Santa Catarina, months compared:', NCMP, ')')
for _, r_ in TB.iterrows():
    print(f"  {r_['item']:52s} this {r_['this_inventory_t']:>12,.0f} t | BRAIN true {r_['brain_true_t']:>12,.0f} t "
          f"(x{r_['ratio_vs_true']:.2f}) | BRAIN delivered {r_['brain_delivered_t']:>12,.0f} t (x{r_['ratio_vs_delivered']:.2f})")
print('Outputs in', OUT_DIR)
