# %% [markdown]
# # 21 · ERA5-Land temperature and radiation against the INMET stations: bias, its dependence on elevation, lapse rate
#
# **Why.** In MEGANv3.21 the leaf-level isoprene response (`gamma.gamtld`) changes by d ln γ / dT ≈ 0.11-0.13 K⁻¹ at
# 20-30 °C when only the hourly temperature changes, and by ≈ 0.155-0.16 K⁻¹ (about +17 % per kelvin) when the daily
# mean changes with it, so a temperature bias of the driving meteorology moves the inventory directly. ERA5-Land (≈9 km) represents the mean elevation of its cells; in the Serra
# Geral, the coastal escarpment and the valleys of the Itajaí and Tubarão rivers the terrain within one cell spans
# hundreds of metres. This notebook measures, with INMET observations:
#
# * **A. Normals 1991-2020.** ERA5-Land monthly means 1991-2020 against the INMET climatological normals of mean
#   temperature (`Normal-Climatologica-TMEDSECA.xlsx`, downloaded from the INMET portal; the INMET booklet of the
#   1991-2020 normals lists the compensated mean temperature among its variables). The INMET mean is computed with the
#   formula of the conventional network, not as the mean of 24 hourly values, so part of a station bias here is a
#   difference of definition; it is common to all stations to first order and goes into the intercept of the fit.
# * **B. 2019, hourly, automatic stations** (INMET historical data, `dadoshistoricos/2019.zip`): the same definition on
#   both sides (instantaneous 2 m temperature at the full hour, UTC). This is the main estimate. Bias by station, month
#   and hour; daytime (12-21 UTC, 9-18 h local) separately, because isoprene is emitted by day. Global radiation is
#   compared as daily totals (the timestamp convention of the hourly station values is checked by the lag of the
#   diurnal cycle).
#
# **Lapse rate.** The bias b = T(ERA5-Land) − T(station) is regressed on Δz = z(station) − z(ERA5-Land cell): if the
# whole difference came from elevation, b = Γ Δz. The slope is an empirical lapse rate (K km⁻¹) with an ordinary least
# squares fit (heteroscedasticity-robust HC3 standard error), the Theil-Sen slope, and a bootstrap interval over
# stations (5,000 resamples). The rule that picks the value used by notebook 20 is fixed here, before any result: the
# 2019 daytime estimate if it has at least 15 stations and a bootstrap interval that excludes zero; else the 2019
# all-hours estimate under the same rule; else the normals estimate under the same rule; else 6.5 K km⁻¹.
#
# **Elevations.** Station elevation from the INMET station lists you have (1981-2010 list; 1961-1990 list for stations
# absent from it). Cell elevation: mean of SRTM 90 m v4 over the 0.1° ERA5-Land pixel (Earth Engine), as notebook 20.
# Each station is compared with the ERA5-Land pixel that contains it; a station whose pixel is sea in ERA5-Land takes the
# mean of the land pixels around it (flagged `era5_filled`, the analogue of FIX-1) and is excluded from the fits.
#
# **Area.** Stations inside the box 30.5-25.0° S, 54.5-47.5° W (Santa Catarina and the neighbouring parts of Paraná and
# Rio Grande do Sul): the slope is a property of ERA5-Land in this terrain, so neighbouring stations add information.
# Tables give the Santa Catarina stations separately, with their mesoregion.

# %%
# 1 · Settings -----------------------------------------------------------------------------------------------------------------
%pip -q install earthengine-api pandas numpy scipy statsmodels openpyxl xlrd geopandas matplotlib requests
import os, sys, io, re, json, time, zipfile, unicodedata, datetime, warnings
from concurrent.futures import ThreadPoolExecutor
try:
    from google.colab import drive; drive.mount('/content/drive')
except Exception:
    pass
ROOT_V10    = '/content/drive/MyDrive/megan_sc_outputs_v10'
OUT_DIR     = ROOT_V10 + '/analysis_v10/era5land_vs_inmet'                 # lapse_rate.json is read by notebook 20
BND         = '/content/drive/MyDrive/MEGAN_SC/data/boundaries'            # <-- sc_state_ibge.gpkg, sc_mesoregions_2024.gpkg
LIST_8110   = '/content/drive/MyDrive/INMET/Estações-Normal-Climatoógica-1981-2010.xls'   # <-- INMET station list 1981-2010
LIST_6190   = '/content/drive/MyDrive/INMET/Relac_Est_Meteo_NC.xls'                       # <-- INMET station list 1961-1990
NORMALS_URL = 'https://portal.inmet.gov.br/uploads/normais/Normal-Climatologica-TMEDSECA.xlsx'
NORMALS_XLSX = None          # local copy of that file if the download fails (portal.inmet.gov.br/normais)
HOURLY_URL  = 'https://portal.inmet.gov.br/uploads/dadoshistoricos/2019.zip'
HOURLY_ZIP  = None           # local copy of that file if the download fails (portal.inmet.gov.br/dadoshistoricos)
YEAR        = 2019
BOX         = dict(lat_min=-30.5, lat_max=-25.0, lon_min=-54.5, lon_max=-47.5)
DAY_HOURS   = list(range(12, 22))   # UTC hours 12-21 = 9-18 h local (UTC-3)
MIN_VALID   = 0.7                   # share of valid hourly values a station needs in the year
N_BOOT      = 5000
WORK        = '/content/inmet_tmp'
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, requests
warnings.filterwarnings('ignore', category=RuntimeWarning)
def ascii_(s):
    return unicodedata.normalize('NFKD', str(s)).encode('ascii', 'ignore').decode().strip()

# %%
# 2 · Earth Engine ---------------------------------------------------------------------------------------------------------------
EE_PROJECT = 'your-ee-project'         # <-- your Earth Engine (Google Cloud) project id
import ee
ee.Authenticate(); ee.Initialize(project=EE_PROJECT)
HOURLY = 'ECMWF/ERA5_LAND/HOURLY'
PJ = ee.ImageCollection(HOURLY).first().select('temperature_2m').projection().getInfo()
print('ERA5-Land projection:', PJ['crs'], PJ['transform'])
assert abs(PJ['transform'][0] - 0.1) < 1e-9, PJ

# %% [markdown]
# ## Cell 3 · Station lists (coordinates and elevations)

# %%
# 3 · Station lists ----------------------------------------------------------------------------------------------------------------
def dms(s):
    """'27°35'S' -> -27.583 (degrees and minutes of the 1961-1990 list)."""
    m = re.match(r"\s*(\d+)\D+(\d+)\D*([NSEWO])", str(s))
    if not m:
        return np.nan
    v = int(m[1]) + int(m[2]) / 60.0
    return -v if m[3] in 'SWO' else v
def read_lists(p8110, p6190):
    a = pd.read_excel(p8110, header=None)
    h = [i for i in range(min(10, len(a))) if any(ascii_(x).lower() == 'codigo' for x in a.iloc[i].tolist())][0]
    a.columns = [ascii_(x) for x in a.iloc[h]]; a = a.iloc[h + 1:]
    a = a.rename(columns={'Codigo': 'code', 'Nome da Estacao': 'name', 'UF': 'uf', 'Latitude': 'lat', 'Longitude': 'lon',
                          'Atitude': 'z', 'Altitude': 'z'})
    a = a[pd.to_numeric(a.code, errors='coerce').notna()].assign(source='INMET list 1981-2010')
    b = pd.read_excel(p6190)
    b.columns = [ascii_(x) for x in b.columns]
    b = b.rename(columns={'Codigo': 'code', 'Nome': 'name', 'UF': 'uf', 'Altitude (m)': 'z'})
    b['lat'] = b['Latitude'].map(dms); b['lon'] = b['Longitude'].map(dms)
    b = b[pd.to_numeric(b.code, errors='coerce').notna()].assign(source='INMET list 1961-1990')
    keep = ['code', 'name', 'uf', 'lat', 'lon', 'z', 'source']
    s = pd.concat([a[keep], b[keep]], ignore_index=True)
    s['code'] = s.code.astype(int); s[['lat', 'lon', 'z']] = s[['lat', 'lon', 'z']].apply(pd.to_numeric, errors='coerce')
    return s.drop_duplicates('code', keep='first').set_index('code')
LISTS = read_lists(LIST_8110, LIST_6190)
inbox = lambda d: d[(d.lat >= BOX['lat_min']) & (d.lat <= BOX['lat_max']) & (d.lon >= BOX['lon_min']) & (d.lon <= BOX['lon_max'])]
print(f'{len(LISTS)} conventional stations in the two lists; {len(inbox(LISTS))} inside the box')

# %% [markdown]
# ## Cell 4 · Elevation of the ERA5-Land pixel and ERA5-Land values at points
#
# `era5_points` returns, for each point, the value of an image in the ERA5-Land pixel that contains it (sampling on the
# ERA5-Land grid itself), and, where that pixel has no value (sea), the mean of the valid pixels among its eight
# neighbours, with a flag.

# %%
# 4 · ERA5-Land sampling helpers ------------------------------------------------------------------------------------------------------
DEM = ee.Image('CGIAR/SRTM90_V4').select('elevation')
def fc_of(df):
    return ee.FeatureCollection([ee.Feature(ee.Geometry.Point([float(r.lon), float(r.lat)]), {'sid': str(i)}) for i, r in df.iterrows()])
def reduce_at(img, df, bands):
    out = img.reduceRegions(collection=fc_of(df), reducer=ee.Reducer.first().forEachBand(img) if len(bands) > 1 else ee.Reducer.first(),
                            crs=PJ['crs'], crsTransform=PJ['transform']).getInfo()['features']
    rows = {f['properties']['sid']: {b: f['properties'].get(b if len(bands) > 1 else 'first') for b in bands} for f in out}
    return pd.DataFrame.from_dict(rows, orient='index').reindex([str(i) for i in df.index]).astype(float).set_axis(df.index)
def era5_points(img, df, bands):
    """Values in the containing ERA5-Land pixel; NaN rows filled with the mean of the eight neighbouring pixel centres."""
    v = reduce_at(img, df, bands); filled = v.isna().any(axis=1)
    if filled.any():
        nb = []
        for i in df.index[filled]:
            c0, r0 = round(df.loc[i, 'lon'], 1), round(df.loc[i, 'lat'], 1)
            nb += [dict(sid=f'{i}_{k}', owner=i, lon=c0 + dx, lat=r0 + dy) for k, (dx, dy) in
                   enumerate([(a, b) for a in (-0.1, 0, 0.1) for b in (-0.1, 0, 0.1) if (a, b) != (0, 0)])]
        nb = pd.DataFrame(nb).set_index('sid'); w = reduce_at(img, nb, bands)
        w['owner'] = nb.owner
        v.loc[filled] = w.groupby('owner')[bands].mean().reindex(df.index[filled]).values
    return v, filled
Z_ERA_IMG = DEM.reduceResolution(ee.Reducer.mean(), maxPixels=65536).rename('z_era5')

# %% [markdown]
# ## Cell 5 · A. Normals 1991-2020

# %%
# 5 · Normals 1991-2020 ------------------------------------------------------------------------------------------------------------
MONTHS_PT = ['janeiro', 'fevereiro', 'marco', 'abril', 'maio', 'junho', 'julho', 'agosto', 'setembro', 'outubro', 'novembro', 'dezembro']
def read_normals(path):
    raw = pd.read_excel(path, header=None)
    title = ' | '.join(ascii_(x) for x in raw.iloc[:3].values.ravel() if str(x) != 'nan')
    h = [i for i in range(min(15, len(raw))) if any(ascii_(x).lower() == 'janeiro' for x in raw.iloc[i].tolist())][0]
    cols = [ascii_(x).lower() for x in raw.iloc[h]]
    d = raw.iloc[h + 1:].copy(); d.columns = cols
    code_c = [c for c in cols if 'codigo' in c][0]; uf_c = [c for c in cols if c == 'uf'][0]
    d = d[pd.to_numeric(d[code_c], errors='coerce').notna()]
    out = pd.DataFrame({'code': d[code_c].astype(int).values, 'uf': d[uf_c].values})
    for k, m in enumerate(MONTHS_PT, 1):
        out[f'm{k:02d}'] = pd.to_numeric(d[m].astype(str).str.replace(',', '.'), errors='coerce').values
    return out.set_index('code'), title
if NORMALS_XLSX is None:
    NORMALS_XLSX = os.path.join(WORK, 'Normal-Climatologica-TMEDSECA.xlsx')
    r = requests.get(NORMALS_URL, timeout=120, headers={'User-Agent': 'MEGAN_SC notebook 21'}); r.raise_for_status()
    open(NORMALS_XLSX, 'wb').write(r.content)
NOR, TITLE = read_normals(NORMALS_XLSX)
print('normals file title rows:', TITLE)
NS = NOR.join(LISTS[['name', 'lat', 'lon', 'z', 'source']], how='inner')
print(f'{len(NOR)} stations in the normals; {len(NS)} with coordinates from the two lists; '
      f'{len(NOR) - len(NS)} without (listed in normals_stations_without_coordinates.csv)')
NOR[~NOR.index.isin(NS.index)].to_csv(os.path.join(OUT_DIR, 'normals_stations_without_coordinates.csv'))
NS = inbox(NS); NS = NS[NS[[f'm{k:02d}' for k in range(1, 13)]].notna().sum(1) == 12]
print(f'{len(NS)} stations inside the box with all twelve months')
MON = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR').filterDate('1991-01-01', '2021-01-01')
assert 'temperature_2m' in MON.first().bandNames().getInfo(), 'ERA5_LAND/MONTHLY_AGGR has no band temperature_2m'
assert MON.size().getInfo() == 360, MON.size().getInfo()
MON = MON.select('temperature_2m')
CLIM = ee.Image.cat([MON.filter(ee.Filter.calendarRange(k, k, 'month')).mean().subtract(273.15).rename(f'e{k:02d}') for k in range(1, 13)])
EV, FILLED_N = era5_points(CLIM, NS, [f'e{k:02d}' for k in range(1, 13)])
ZE, _ = era5_points(Z_ERA_IMG, NS, ['z_era5'])
ZS = reduce_at(DEM.rename('z_srtm_point'), NS, ['z_srtm_point'])
NS = NS.join(EV).join(ZE).join(ZS); NS['era5_filled'] = FILLED_N
for k in range(1, 13):
    NS[f'b{k:02d}'] = NS[f'e{k:02d}'] - NS[f'm{k:02d}']
NS['bias_annual'] = NS[[f'b{k:02d}' for k in range(1, 13)]].mean(1); NS['dz'] = NS.z - NS.z_era5
NS['z_list_minus_srtm'] = NS.z - NS.z_srtm_point
NS.to_csv(os.path.join(OUT_DIR, 'normals_1991_2020_stations.csv'))
print(NS[['name', 'uf', 'z', 'z_era5', 'dz', 'bias_annual', 'era5_filled', 'z_list_minus_srtm']].round(2).to_string())

# %% [markdown]
# ## Cell 6 · B. 2019 hourly, automatic stations
#
# The annual archive holds one CSV per station: eight metadata lines (region, state, name, WMO code, latitude,
# longitude, altitude, foundation date), a header line and hourly rows; `;` separated, decimal comma, -9999 for missing
# values, latin-1. The parser finds the columns by name (date, `HORA`, `BULBO SECO`, `RADIACAO GLOBAL`) and stops if
# the hour is not in UTC.

# %%
# 6 · 2019 hourly stations --------------------------------------------------------------------------------------------------------------
def parse_inmet_csv(raw_bytes):
    txt = raw_bytes.decode('latin-1').replace('\r\n', '\n').replace('\r', '\n')
    lines = txt.split('\n')
    meta = {}
    for l in lines[:8]:
        k, _, v = l.partition(';'); meta[ascii_(k).rstrip(':').upper()] = v.strip().strip(';')
    num = lambda s: float(str(s).replace(',', '.')) if str(s).strip() not in ('', '-9999') else np.nan
    info = dict(uf=meta.get('UF'), name=meta.get('ESTACAO'), code=meta.get('CODIGO (WMO)'), lat=num(meta.get('LATITUDE')),
                lon=num(meta.get('LONGITUDE')), z=num(meta.get('ALTITUDE')))
    df = pd.read_csv(io.StringIO('\n'.join(lines[8:])), sep=';', decimal=',', na_values=['-9999', '-9999,0', ''], dtype=str)
    cols = {c: ascii_(c).upper() for c in df.columns}
    c_date = [c for c, a in cols.items() if a.startswith('DATA')][0]
    c_hour = [c for c, a in cols.items() if a.startswith('HORA')][0]
    c_t = [c for c, a in cols.items() if 'BULBO SECO' in a][0]
    c_r = [c for c, a in cols.items() if 'RADIACAO' in a]
    sample = str(df[c_hour].dropna().iloc[0])
    assert 'UTC' in cols[c_hour] or 'UTC' in sample.upper(), f'hour column {c_hour!r} ({sample!r}) is not marked UTC'
    hh = df[c_hour].astype(str).str.extract(r'(\d{1,2})[:]?(\d{2})?')[0].astype(int)
    day = pd.to_datetime(df[c_date].astype(str).str.replace('/', '-'), errors='coerce')
    out = pd.DataFrame({'time': day + pd.to_timedelta(hh, unit='h'),
                        't_obs': pd.to_numeric(df[c_t].astype(str).str.replace(',', '.'), errors='coerce')})
    out['sw_obs_kj'] = pd.to_numeric(df[c_r[0]].astype(str).str.replace(',', '.'), errors='coerce') if c_r else np.nan
    out.loc[(out.t_obs < -10) | (out.t_obs > 45) | (out.t_obs == -9999), 't_obs'] = np.nan
    out.loc[out.sw_obs_kj <= -9999, 'sw_obs_kj'] = np.nan
    return info, out.dropna(subset=['time']).set_index('time')
if HOURLY_ZIP is None:
    HOURLY_ZIP = os.path.join(WORK, f'inmet_{YEAR}.zip')
    if not os.path.exists(HOURLY_ZIP):
        with requests.get(HOURLY_URL, stream=True, timeout=(30, 900), headers={'User-Agent': 'MEGAN_SC notebook 21'}) as r:
            r.raise_for_status()
            with open(HOURLY_ZIP, 'wb') as f:
                for ch in r.iter_content(1 << 22):
                    f.write(ch)
ST, OBS = {}, {}
with zipfile.ZipFile(HOURLY_ZIP) as z:
    names = [n for n in z.namelist() if n.lower().endswith('.csv')]
    for n in names:
        info, o = parse_inmet_csv(z.read(n))
        if not (np.isfinite(info['lat']) and BOX['lat_min'] <= info['lat'] <= BOX['lat_max'] and BOX['lon_min'] <= info['lon'] <= BOX['lon_max']):
            continue
        o = o[o.index.year == YEAR]
        if o.t_obs.notna().mean() * len(o) < MIN_VALID * 8760:
            continue
        ST[info['code']] = {**info, 'file': n, 'valid_t_share': round(o.t_obs.notna().sum() / 8760, 3)}; OBS[info['code']] = o
SA = pd.DataFrame.from_dict(ST, orient='index'); SA.index.name = 'code'
print(f'{len(names)} station files in the archive; {len(SA)} automatic stations inside the box with at least {MIN_VALID:.0%} valid hours')

# %% [markdown]
# ## Cell 7 · ERA5-Land hourly series at the automatic stations
#
# One `getRegion` request per station and month on the ERA5-Land grid (2 m temperature and the accumulated
# downward shortwave radiation), with the hour before the month for the de-accumulation (FIX-8 rule of `fixes.py`:
# the 01 UTC field is already a one-hour accumulation; every other hour, 00 UTC included, subtracts the previous hour).

# %%
# 7 · ERA5-Land at the stations ------------------------------------------------------------------------------------------------------------
def deaccumulate(acc_cur, acc_prev, hour_utc):                 # numpy reference of FIX-8 (src/megan_sc/fixes.py)
    prev = np.where(np.asarray(hour_utc) == 1, 0.0, np.asarray(acc_prev, dtype='float64'))
    return np.maximum(np.asarray(acc_cur, dtype='float64') - prev, 0.0) / 3600.0
def series_at(lon, lat, t0, t1):
    col = ee.ImageCollection(HOURLY).filterDate(t0, t1).select(['temperature_2m', 'surface_solar_radiation_downwards'])
    for k in range(4):
        try:
            rows = col.getRegion(ee.Geometry.Point([lon, lat]), crs=PJ['crs'], crsTransform=PJ['transform']).getInfo(); break
        except Exception as e:
            if k == 3: raise
            time.sleep(10 * (k + 1))
    d = pd.DataFrame(rows[1:], columns=rows[0])
    d['time'] = pd.to_datetime(d['time'], unit='ms'); d = d.set_index('time').sort_index()
    return d[['temperature_2m', 'surface_solar_radiation_downwards']].astype(float)
def station_era5(code):
    r = SA.loc[code]; parts = []
    for m in range(1, 13):
        t0 = (pd.Timestamp(YEAR, m, 1) - pd.Timedelta(hours=1)).strftime('%Y-%m-%dT%H:%M'); t1 = pd.Timestamp(YEAR, m, 1) + pd.offsets.MonthBegin(1)
        parts.append(series_at(float(r.lon), float(r.lat), t0, t1.strftime('%Y-%m-%dT%H:%M')))
    d = pd.concat(parts); d = d[~d.index.duplicated()]
    filled = d.temperature_2m.isna().mean() > 0.5
    if filled:                                                  # sea pixel: mean of the land pixels around it
        c0, r0 = round(float(r.lon), 1), round(float(r.lat), 1); nbs = []
        for dx in (-0.1, 0, 0.1):
            for dy in (-0.1, 0, 0.1):
                if (dx, dy) != (0, 0):
                    nbs.append(pd.concat([series_at(c0 + dx, r0 + dy, (pd.Timestamp(YEAR, m, 1) - pd.Timedelta(hours=1)).strftime('%Y-%m-%dT%H:%M'),
                                                    (pd.Timestamp(YEAR, m, 1) + pd.offsets.MonthBegin(1)).strftime('%Y-%m-%dT%H:%M')) for m in range(1, 13)]))
        d = pd.concat(nbs).groupby(level=0).mean()
    sw = deaccumulate(d.surface_solar_radiation_downwards.values, d.surface_solar_radiation_downwards.shift(1).values, d.index.hour)
    out = pd.DataFrame({'t_era5': d.temperature_2m.values - 273.15, 'sw_era5_wm2': sw}, index=d.index)
    return out[out.index.year == YEAR], filled
t0 = time.time(); E5, FL = {}, {}
with ThreadPoolExecutor(6) as ex:
    for code, (e, f) in zip(SA.index, ex.map(station_era5, SA.index)):
        E5[code], FL[code] = e, f
print(f'ERA5-Land series for {len(E5)} stations in {time.time() - t0:.0f} s; sea pixels filled from neighbours: {sum(FL.values())}')
ZA, _ = era5_points(Z_ERA_IMG, SA[['lat', 'lon']], ['z_era5'])
SA = SA.join(ZA); SA['dz'] = SA.z - SA.z_era5; SA['era5_filled'] = pd.Series(FL)

# %% [markdown]
# ## Cell 8 · Bias statistics and lapse-rate fits

# %%
# 8 · Bias and lapse rate --------------------------------------------------------------------------------------------------------------------
import geopandas as gpd, statsmodels.api as sm
from scipy import stats
meso = gpd.read_file(os.path.join(BND, 'sc_mesoregions_2024.gpkg')).to_crs('EPSG:4326')
def meso_of(lon, lat):
    p = gpd.points_from_xy([lon], [lat])[0]
    hit = meso[meso.contains(p)]
    return hit.name_meso.iloc[0] if len(hit) else None
rows, monthly, hourly = [], [], []
for code in SA.index:
    j = OBS[code].join(E5[code], how='inner'); j = j[j.t_obs.notna() & j.t_era5.notna()]
    j['bias'] = j.t_era5 - j.t_obs; day = j.index.hour.isin(DAY_HOURS)
    rows.append(dict(code=code, n=len(j), bias_all=j.bias.mean(), bias_day=j.bias[day].mean(), bias_night=j.bias[~day].mean(),
                     rmse_all=np.sqrt((j.bias ** 2).mean()), r_hourly=j[['t_era5', 't_obs']].corr().iloc[0, 1],
                     t_obs_day=j.t_obs[day].mean(), t_era5_day=j.t_era5[day].mean()))
    for m, g in j.groupby(j.index.month):
        monthly.append(dict(code=code, month=m, bias_all=g.bias.mean(), bias_day=g.bias[g.index.hour.isin(DAY_HOURS)].mean()))
    for h, g in j.groupby(j.index.hour):
        hourly.append(dict(code=code, hour_utc=h, bias=g.bias.mean()))
B = SA.join(pd.DataFrame(rows).set_index('code'))
B['uf_meso'] = [meso_of(r.lon, r.lat) or r.uf for _, r in B.iterrows()]
B.to_csv(os.path.join(OUT_DIR, f'hourly_{YEAR}_stations_bias.csv'))
pd.DataFrame(monthly).to_csv(os.path.join(OUT_DIR, f'hourly_{YEAR}_bias_by_month.csv'), index=False)
HB = pd.DataFrame(hourly); HB.to_csv(os.path.join(OUT_DIR, f'hourly_{YEAR}_bias_by_hour.csv'), index=False)
def fit(df, y, x='dz', label=''):
    d = df[[x, y]].dropna(); n = len(d)
    if n < 5:
        return dict(label=label, n=n)
    X = sm.add_constant(d[x].values / 1000.0); o = sm.OLS(d[y].values, X).fit(cov_type='HC3')
    ts = stats.theilslopes(d[y].values, d[x].values / 1000.0)
    rng = np.random.default_rng(21); bs = []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        if np.ptp(d[x].values[i]) > 0:
            bs.append(np.polyfit(d[x].values[i] / 1000.0, d[y].values[i], 1)[0])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(label=label, n=n, gamma_ols_K_per_km=o.params[1], se_hc3=o.bse[1], p_hc3=o.pvalues[1], intercept_K=o.params[0],
                r2=o.rsquared, gamma_theil_sen=ts[0], ts_lo=ts[2], ts_hi=ts[3], boot_lo=lo, boot_hi=hi)
use = lambda df: df[~df.era5_filled.astype(bool)]
FITS = pd.DataFrame([fit(use(B), 'bias_day', label=f'{YEAR} daytime (12-21 UTC)'), fit(use(B), 'bias_all', label=f'{YEAR} all hours'),
                     fit(use(NS), 'bias_annual', label='normals 1991-2020, annual'),
                     fit(use(B[B.uf == 'SC']), 'bias_day', label=f'{YEAR} daytime, SC stations only')])
FITS.to_csv(os.path.join(OUT_DIR, 'lapse_rate_fits.csv'), index=False)
pd.set_option('display.width', 220); print(FITS.round(3).to_string(index=False))
def ok_(r):
    return r.get('n', 0) >= 15 and np.isfinite(r.get('boot_lo', np.nan)) and (r['boot_lo'] > 0 or r['boot_hi'] < 0)
choice = next((r for _, r in FITS.iloc[:3].iterrows() if ok_(r)), None)
GAMMA_N20 = float(choice['gamma_ols_K_per_km']) if choice is not None else 6.5
json.dump(dict(rule='first of: 2019 daytime, 2019 all hours, normals, with n >= 15 and a bootstrap 95 % interval excluding 0; else 6.5',
               chosen=(choice['label'] if choice is not None else 'standard atmosphere 6.5 K/km'), gamma_for_notebook20_K_per_km=GAMMA_N20,
               fits=FITS.to_dict(orient='records'), year=YEAR, day_hours_utc=DAY_HOURS, box=BOX,
               normals_file=NORMALS_URL, normals_title=TITLE, hourly_file=HOURLY_URL, stations_hourly=int(len(B)), stations_normals=int(len(NS)),
               run_utc=datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z'),
          open(os.path.join(OUT_DIR, 'lapse_rate.json'), 'w'), indent=1, default=float)
print(f'lapse rate for notebook 20: {GAMMA_N20:.2f} K/km ({choice["label"] if choice is not None else "6.5 default"})')
# diagnostics only (the rule above and the value passed to notebook 20 are not changed by them):
# weight of the single most influential station on the two 2019 fits, and the mean bias over the stations
d_ = use(B)[['dz', 'bias_day', 'bias_all']].dropna(); X_ = sm.add_constant(d_.dz.values / 1000.0); INF = []
for ycol in ('bias_day', 'bias_all'):
    yv = d_[ycol].values; o_ = sm.OLS(yv, X_).fit(); cook = o_.get_influence().cooks_distance[0]; k = int(np.argmax(cook))
    o2 = sm.OLS(np.delete(yv, k), np.delete(X_, k, 0)).fit(cov_type='HC3'); ci = o2.conf_int()
    ts2 = stats.theilslopes(np.delete(yv, k), np.delete(X_[:, 1], k))
    loo = [sm.OLS(np.delete(yv, i), np.delete(X_, i, 0)).fit().params[1] for i in range(len(yv))]
    INF.append(dict(response=ycol, n=len(yv), most_influential=d_.index[k], cook_d=cook[k], gamma_all=o_.params[1],
                    gamma_without=o2.params[1], ols_hc3_lo_without=ci[1, 0], ols_hc3_hi_without=ci[1, 1], theil_sen_without=ts2[0],
                    ts_lo_without=ts2[2], ts_hi_without=ts2[3], leave_one_out_min=min(loo), leave_one_out_max=max(loo)))
INF = pd.DataFrame(INF); INF.to_csv(os.path.join(OUT_DIR, 'lapse_rate_influence.csv'), index=False)
print('\nmost influential station of each 2019 fit (Cook distance) and the fit without it:'); print(INF.round(3).to_string(index=False))
BS = []
for lab, s_ in (('box', use(B)), ('Santa Catarina', use(B[B.uf == 'SC']))):
    for ycol in ('bias_day', 'bias_all'):
        b_ = s_[ycol].dropna(); h = stats.t.ppf(0.975, len(b_) - 1) * b_.std(ddof=1) / np.sqrt(len(b_))
        BS.append(dict(stations=lab, response=ycol, n=len(b_), mean_K=b_.mean(), lo=b_.mean() - h, hi=b_.mean() + h, median_K=b_.median(),
                       wilcoxon_p=stats.wilcoxon(b_).pvalue))
BS = pd.DataFrame(BS); BS.to_csv(os.path.join(OUT_DIR, f'bias_summary_{YEAR}.csv'), index=False)
print(f'\nmean bias ERA5-Land minus station, {YEAR} (K; t interval, Wilcoxon signed-rank p):'); print(BS.round(3).to_string(index=False))
print('\nSanta Catarina stations, 2019 (bias = ERA5-Land minus station, K):')
print(B[B.uf == 'SC'][['name', 'uf_meso', 'z', 'z_era5', 'dz', 'bias_all', 'bias_day', 'bias_night', 'r_hourly', 'era5_filled']].sort_values('uf_meso').round(2).to_string())
print(B[B.uf == 'SC'].groupby('uf_meso')[['bias_all', 'bias_day']].agg(['mean', 'count']).round(2).to_string())

# %% [markdown]
# ## Cell 9 · Radiation (daily totals) and figures
#
# Daily global radiation: ERA5-Land de-accumulated hourly flux summed over the day (MJ m⁻²) against the sum of the
# hourly station values (kJ m⁻² → MJ m⁻²), on days where the station has a value for every hour with ERA5-Land
# radiation above 10 W m⁻². The lag of the station diurnal cycle with respect to ERA5-Land (−2 to +2 h) shows the
# timestamp convention of the station values; daily totals do not depend on it.

# %%
# 9 · Radiation and figures --------------------------------------------------------------------------------------------------------------------
import matplotlib; import matplotlib.pyplot as plt
rad, lags = [], []
for code in B.index:
    j = OBS[code].join(E5[code], how='inner')
    if j.sw_obs_kj.notna().sum() < 1000:
        continue
    j['sw_obs_wm2'] = j.sw_obs_kj.clip(lower=0) * 1000.0 / 3600.0
    dc_o = j.groupby(j.index.hour).sw_obs_wm2.mean(); dc_e = j.groupby(j.index.hour).sw_era5_wm2.mean()
    lags.append(dict(code=code, **{f'r_lag{k:+d}h': np.corrcoef(np.roll(dc_o.values, -k), dc_e.values)[0, 1] for k in range(-2, 3)}))
    sun = j.sw_era5_wm2 > 10
    g = j[sun].groupby(j[sun].index.date)
    dd = pd.DataFrame({'obs_MJ': g.sw_obs_wm2.sum() * 3600 / 1e6, 'era5_MJ': g.sw_era5_wm2.sum() * 3600 / 1e6, 'n_obs': g.sw_obs_wm2.count(), 'n_sun': g.size()})
    dd = dd[dd.n_obs == dd.n_sun]
    if len(dd) >= 200:
        rad.append(dict(code=code, name=B.loc[code, 'name'], uf_meso=B.loc[code, 'uf_meso'], days=len(dd), obs_MJ=dd.obs_MJ.mean(),
                        era5_MJ=dd.era5_MJ.mean(), ratio=dd.era5_MJ.mean() / dd.obs_MJ.mean(), r_daily=dd[['obs_MJ', 'era5_MJ']].corr().iloc[0, 1]))
RAD = pd.DataFrame(rad); RAD.to_csv(os.path.join(OUT_DIR, f'radiation_{YEAR}_daily_totals.csv'), index=False)
LAG = pd.DataFrame(lags); LAG.to_csv(os.path.join(OUT_DIR, f'radiation_{YEAR}_diurnal_lag.csv'), index=False)
if len(RAD):
    print(f'daily global radiation, {len(RAD)} stations: ERA5-Land / station median {RAD.ratio.median():.3f} '
          f'(IQR {RAD.ratio.quantile(0.25):.3f} to {RAD.ratio.quantile(0.75):.3f}); median r of daily totals {RAD.r_daily.median():.2f}')
    print('diurnal-cycle correlation by lag (median over stations):', LAG.drop(columns='code').median().round(3).to_dict())
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for a, (df, y, t) in zip(ax[:2], ((use(B), 'bias_day', f'{YEAR}, 12-21 UTC'), (use(NS), 'bias_annual', 'normals 1991-2020'))):
    a.scatter(df.dz, df[y], s=18, c=np.where(df.uf == 'SC', 'C0', 'C7')); a.axhline(0, color='k', lw=0.5)
    r = FITS[FITS.label.str.startswith(t.split(',')[0] if 'normals' not in t else 'normals')].iloc[0]
    xx = np.linspace(df.dz.min(), df.dz.max(), 2)
    if np.isfinite(r.get('gamma_ols_K_per_km', np.nan)):
        a.plot(xx, r['intercept_K'] + r['gamma_ols_K_per_km'] * xx / 1000, 'C3', label=f"{r['gamma_ols_K_per_km']:.1f} K/km")
    a.set_xlabel('station minus ERA5-Land cell elevation (m)'); a.set_ylabel('ERA5-Land minus station (K)'); a.set_title(t); a.legend()
hb = HB.merge(B[['uf']], left_on='code', right_index=True)
hb[hb.uf == 'SC'].groupby('hour_utc').bias.median().plot(ax=ax[2], marker='o'); ax[2].axhline(0, color='k', lw=0.5)
ax[2].set_xlabel('hour (UTC)'); ax[2].set_ylabel('median bias, SC stations (K)'); ax[2].set_title(f'diurnal bias {YEAR}')
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, 'era5land_vs_inmet.png'), dpi=200); plt.show()
print('outputs in', OUT_DIR)
