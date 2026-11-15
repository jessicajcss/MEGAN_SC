# Time conventions of the inputs and of the model

Every statement names its source: the documentation of the dataset, the MEGANv3.21 code
(Guenther 2024; file and line), the code of this repository, or a published result file.

## 1. Time reference of each input

| Input | Time reference | Source |
|---|---|---|
| ERA5-Land hourly 2 m temperature, dew point, 10 m wind, surface pressure, soil water | instantaneous, valid at the hour of the timestamp (UTC) | ECMWF, *ERA5-Land: data documentation*: parameters described as instantaneous "are valid at the specified time" |
| ERA5-Land surface solar radiation downwards | accumulated from 00 UTC to the validity time; the field valid at 00 UTC holds the 24 h of the previous day. After de-accumulation, the value of hour *h* is the mean flux over the hour ending at *h* | ECMWF, same page ("For the CDS time, or validity time, of 00 UTC, the accumulations are over the 24 hours ending at 00 UTC"); Earth Engine catalog page of `ECMWF/ERA5_LAND/HOURLY` (accumulations reset at midnight); `fixes.deaccumulate` and `fixes.hourly_flux_fixed` (FIX-8) |
| MODIS MCD15A3H LAI | 4-day composites, averaged over the images dated in each calendar month | Myneni et al. (2021); `src/megan_sc/gee/lai.py` |
| VIIRS LAI (VNP15A2H, VJ115A2H, VJ215A2H; diagnostics only) | 8-day composites | Myneni (2023a, 2023b, 2026) |
| MapBiomas land cover | annual maps | Souza Junior et al. (2020) |
| TROPOMI HCHO (evaluation only) | overpass near 13:30 local solar time; images averaged by month; fire mask on the same or the previous UTC day | De Smedt et al. (2018); `notebooks/src/14_tropomi_hcho_evaluation.py` |
| INMET hourly station data (notebook 21) | UTC | `notebooks/src/21_era5land_vs_inmet.py` |

LAI and land cover have no time of day. A shift of the day boundary by three hours is far below their 4-day, 8-day or
annual steps.

## 2. How the model uses time

The production run keeps the conventions of the MEGANv3.21 code, which the Python port reproduces.

* **Solar position.** `megcan.f` (lines 500-508) computes the hour as `ITIME/10000 + LONG/15`, the UTC hour of the
  step plus the longitude of the cell divided by 15, moves negative hours to the previous day, and passes it to
  `Calcbeta`, whose formula has no equation-of-time term. `engine/canopy.py` (`megcan`) does the same. Light at the
  leaves is therefore placed in local mean solar time in every cell, whatever the time zone. Across Santa Catarina
  (cell centres from 53.89° W to 48.31° W), local mean solar time is UTC − 3.59 h to UTC − 3.22 h.
* **Daily mean temperature.** `daymet.f` averages the 24 hourly records of the meteorological file of the day starting
  at 00 (TEMP2), and `megvea.f` reads that daily value for the same date. In `GAMTLD`, `T240 = T24`, so the daily mean
  of the UTC day stands for both the 24-h and the 240-h temperatures. Guenther et al. (2012, Sect. 2.2, Eqs. 9-10)
  define T24 and T240 as the average leaf temperature of the past 24 and 240 h. Notebook 10 averages the 24 hourly
  ERA5-Land 2 m temperatures from 00 to 23 UTC of the day of the step (FIX-2). The region wizard averages its sampled
  steps of the day instead.
* **Daily mean PPFD.** In `GAMP` of `megvea.f`, the coefficient C1 is the constant 1.03 (the form that depends on
  PPFD24 is commented out), so the daily PPFD only switches light-dependent emission off where it is below 0.01. The port
  keeps this (`engine/gamma.py`, `gamp`).
* **Pairing within a step.** At step *h* the engine uses the temperature at *h*, the radiation averaged over the hour
  ending at *h*, and the solar elevation at *h*. The canopy model also splits radiation into direct and diffuse parts
  from its ratio to the clear-sky value at that solar elevation (`canopy.solar_fractions`).
* **Integration.** Each 3-h step *h* (00, 03, ..., 21 UTC) stands for [*h*, *h* + 3 h) in the monthly sums, so days,
  months and years of the output are UTC periods.

## 3. Effect on the Santa Catarina results

* **Month and year boundaries.** A UTC month starts at 00 UTC, about 20:24-20:47 local mean solar time on the last day
  of the previous month. Counting months from 03 UTC (local midnight at UTC−3) instead would remove the 00 UTC step of
  the first day and add that of the first day of the next month. Computed from the diurnal profile of each month of the
  run (`analysis/v10/review/month_boundary.py`, `results/v10/final_review/month_boundary_effect.csv`), this changes
  monthly totals by at most 0.25% for monoterpenes, 0.28% for soil NO and 0.17% for all VOC classes in moles (95th
  percentiles 0.14%, 0.21% and 0.10%), leaves isoprene unchanged (the 00 UTC step is at night), and changes annual totals
  by at most 0.01%. The first-day values are approximated by the monthly mean of the 00 UTC step.
* **Half-hour offset between radiation and solar position.** The radiation of step *h* is centred 30 min before the
  solar position and the temperature it is paired with. This was not tested. Test S1 (3-hourly against hourly steps,
  −0.06%) does not measure it, because both runs pair the inputs in the same way.
* **Window of the daily temperature.** The UTC day runs from about 20:30 local solar time of the previous day to the same
  time of the day, and it includes hours after the step, as in MEGANv3.21. A local-day window or a past-24-h window was
  not tested.
* **TROPOMI comparison.** The model was sampled at 16:30 UTC (12:54-13:17 local mean solar time), interpolated between
  the 15 and 18 UTC steps of the diurnal profile (notebook 14).

## 4. A test that would measure the two untested items

Both need Earth Engine. In the 28 days of test S1 (the first seven days of January, April, July and October 2019),
one variant would compute the solar elevation at *h* − 0.5 h, the centre of the radiation hour, and another would
average the daily temperature over the 24 h before each step. Each variant compared with the base run gives the size of
one item. These runs were not made.

## Sources

* ECMWF (2025) ERA5-Land: data documentation. https://confluence.ecmwf.int/spaces/CKB/pages/140385202/ERA5-Land+data+documentation
  (read on 15 November 2026)
* Google Earth Engine, catalog page of `ECMWF/ERA5_LAND/HOURLY` (read on 15 November 2026)
* Guenther A (2024) MEGANv3.21.tar.gz emission calculator code (3.21). Zenodo. https://doi.org/10.5281/zenodo.10526206
  (`src/DAYMET/daymet.f`, `src/MEGCAN/megcan.f`, `src/MEGVEA/megvea.f`)
* Guenther AB et al. (2012) Geoscientific Model Development 5:1471-1492. https://doi.org/10.5194/gmd-5-1471-2012
* Myneni R, Knyazikhin Y, Park T (2021) MCD15A3H V061. https://doi.org/10.5067/MODIS/MCD15A3H.061
