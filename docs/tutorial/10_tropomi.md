# 10. The TROPOMI evaluation

This is the only comparison of the inventory with observations. This chapter explains why formaldehyde columns can test
isoprene emissions, how the test was built and fixed in advance, what came out, and how to compute the magnitude
check by hand.

## 10.1 Why formaldehyde

Isoprene reacts with OH within about an hour in daylight (Atkinson; Arey, 2003, give 1.4 h against OH), and its
oxidation produces formaldehyde (HCHO), which lives a few hours. Without transport, the column of HCHO over a region is
proportional to the local emission (Palmer et al., 2003, their equation 1):

```
Ω = (1 / k_HCHO) · Σ_i  Y_i · E_i          Ω: HCHO column (molecules cm⁻²)
                                           E_i: emission of compound i (atoms C cm⁻² s⁻¹)
                                           Y_i: HCHO yield per carbon atom; 1/k_HCHO = τ, the HCHO lifetime
```

For isoprene, Palmer et al. (2003) give ultimate yields of 0.43 per carbon atom at 1 ppb of NOx and 0.29 at 0.1 ppb.
The test is therefore weighted towards isoprene, and cannot check monoterpenes, sesquiterpenes or soil NO separately.

## 10.2 Data and processing (notebook 14)

**Satellite side.**

* TROPOMI tropospheric HCHO, L3 OFFL product of Earth Engine (De Smedt et al., 2018), 2019-2025, in which pixels with
  qa_value below 0.5 are already removed;
* screening: cloud fraction < 0.3, solar zenith angle < 70°, and removal of pixels within 20 km of a MODIS active-fire
  detection (FIRMS) of the same or the previous UTC day;
* monthly means on the 0.025° model grid (one Earth Engine batch export per month), with the count of kept
  observations;
* background: the monthly median over an Atlantic box east of the state (45° W to 41° W, 30° S to 25° S), subtracted
  to give ΔHCHO.

**Model side.** Isoprene at 16:30 UTC, near the overpass (13:30 local solar time; 12:55 to 13:17 local mean solar
time across the state): the monthly mean of each cell multiplied by the diurnal profile of the domain, interpolated
between the 15 and 18 UTC steps, and expressed in carbon atoms (atoms C cm⁻² s⁻¹).

## 10.3 The four tests fixed in advance

The tests were written into the notebook before it was run, so that their outcome could not shape their definition:

| Tier | Question | Estimator | Result |
|---|---|---|---|
| 1 spatial | does ΔHCHO sit where the model puts isoprene? | Pearson r of summer (December-February) means on 0.1° and 0.2° blocks, model smoothed by Gaussian kernels of 0 to 40 km, block-bootstrap interval | r = 0.12 and 0.18 (95 % CI −0.07 to 0.29 and −0.03 to 0.37); with smoothing, −0.11 to 0.16 |
| 2 seasonal | same annual cycle? | r of the 12 monthly climatologies | r = −0.26 (p = 0.41); model peak in January, ΔHCHO peak in September |
| 3 anomalies | same month-to-month departures? | r of the deseasonalised monthly series (n = 84), with the effective sample size of Bretherton et al. (1999) | r = 0.45, N_eff = 72, p < 0.001; ΔHCHO with temperature r = 0.21 |
| 4 slope | is the column per unit emission plausible? | slope S of ΔHCHO on E over 0.2° blocks; S/Y is an implied lifetime | r = 0.18, too weak to interpret |

Files: `results/v10/notebook14/tier1_*` to `tier4_*` and `monthly_series_*`.

## 10.4 Diagnostics after the tests

Labelled post hoc in the documents (`results/v10/diagnostics/nb14/`, `results/v10/final_review/`):

* **The spatial pattern follows the terrain.** The summer ΔHCHO pattern is reproducible between halves of the period
  (split-half r = 0.76 and 0.86), a stable feature of the observations, and it follows ERA5-Land surface pressure (r = 0.82 and
  0.88): larger columns on the coast and in the lowlands, smaller ones on the plateau. Modelled isoprene is nearly
  unrelated to pressure (r ≤ 0.10); with pressure held fixed, the partial correlation of ΔHCHO with isoprene is 0.15
  to 0.21. The documents list the possible explanations, which these data cannot separate: a shorter air column over
  the plateau, for which a single oceanic background is too large; anthropogenic VOC sources concentrated on the coast;
  or a wrong spatial pattern of the modelled emissions.
* **The background shares the seasonal cycle.** The Atlantic background itself follows the seasonal cycle of modelled
  isoprene (r = 0.67 between monthly means), so subtracting it removes most of the summer rise (the raw state column has
  r = 0.37 with isoprene). ΔHCHO peaks in August and September with the fraction of observations removed by the
  screening (r = 0.80).
* **The anomaly test survives its controls.** Without August-October the anomaly correlation is 0.41 (N_eff = 58, p = 0.001). The background anomalies also follow those of isoprene (r = 0.43), a regional signal common to land and
  ocean boxes; with the background anomaly held fixed the partial correlation is 0.65, and with the temperature held fixed 0.48. By season: 0.54 (December-February), 0.45 (March-May), −0.21 (June-August), 0.67 (September-November).

## 10.5 The magnitude, by hand

The summer state means are ΔHCHO = 2.19 × 10¹⁵ molecules cm⁻², or 3.43 × 10¹⁵ molecules cm⁻² with the bias correction,
and E = 1.78 × 10¹² atoms C cm⁻² s⁻¹ (`diagnostics/nb14/nb14_diagnostics.csv`, D2). From Ω = Y · E · τ:

```python
E = 1.78e12                                     # modelled isoprene at the overpass, atoms C cm-2 s-1
for name, dHCHO in (('raw', 2.19e15), ('bias-corrected', 3.43e15)):
    S = dHCHO / E                               # column per unit emission, s: 1230 and 1927
    for Y in (0.29, 0.43):                      # HCHO yield per carbon atom (Palmer et al., 2003)
        tau = S / Y / 3600                      # implied lifetime, h: 1.18 and 0.79 (raw), 1.85 and 1.24 (corrected)
        print(name, Y, round(tau, 2), [round(tau / L, 2) for L in (1.48, 1.76)])
```

The bias correction inverts the regression of Vigouroux et al. (2020), TROPOMI = 1.10 × 10¹⁵ + 0.64 × FTIR. Applied to
the column and to the background, the offset cancels in ΔHCHO, which is divided by 0.64 (2.19/0.64 = 3.42; D2 gives
3.43 from the unrounded mean), and so is the implied lifetime.

The last list compares the implied lifetime with the HCHO column lifetimes of 1.48 to 1.76 h that Palmer et al. (2003,
Table 3) simulated for summer over North America. The ratio of observed to expected ΔHCHO is the ratio of the two
lifetimes: 0.45 to 0.80 without the bias correction and 0.71 to 1.25 with it (`final_review/final_review_numbers.csv`,
key `ratio_obs_expected_palmer_tau`; `final_review/tropomi_magnitude_grid.csv` for a grid of lifetimes and yields).
From the unrounded means, `diagnostics/nb14/nb14_diagnostics.log` (D2) gives S = 1235 s and lifetimes of 1.18 and
0.80 h. The range spans 1, so the magnitude check leaves open whether the modelled isoprene is too high.

## 10.6 What the evaluation supports

* The month-to-month departures of modelled isoprene are seen in the satellite columns, beyond what a regional
  background and the temperature explain.
* The spatial pattern and the seasonal cycle are not reproduced; for the seasonal cycle, the background subtraction
  limits the test itself.
* The magnitude is open. Local flux or concentration measurements in Santa Catarina are what would settle it.

## 10.7 Rerunning and varying the test

Notebook 14 runs from the shipped cache (`reproduce.py nb14`) or from Earth Engine (chapter 6, section 6.5). Its
settings cell exposes every choice: `CF_MAX`, `SZA_MAX`, `FIRE_MASK`, `FIRE_BUF_KM`, `BG_BOX`, `OVERPASS_UTC` and
`BIAS_CORRECT`. A change of screening builds a new cache under a new `TAG` and needs Earth Engine.
