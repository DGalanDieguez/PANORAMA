# `panorama.obs` — Observation Planning and Timing

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](../LICENSE)

> Coordinates, barycentric corrections, exposure times and SNR, FITS observing logs, night and annual visibility, periodic events and orbital phases, telescope pointing and the BinarAlt planning interface.

[← Back to the main README](../README.md)

---
<br>

## Table of Contents
1. [Overview](#1-overview)
2. [Function Summary](#2-function-summary)
3. [Function Reference](#3-function-reference)
4. [Dependencies](#4-dependencies)
5. [Internal Structure](#5-internal-structure)
6. [Demo](#6-demo)
7. [References](#7-references)

---
<br>

## 1. Overview

`obs` covers the preparation, execution and bookkeeping of the spectroscopic monitoring of (massive) stars, with special attention to periodic systems (eclipsing and spectroscopic binaries, pulsators):

* **Coordinates:** angular separations and position angles, and a position in decimal, sexagesimal and Galactic coordinates.
* **Time:** barycentric and heliocentric radial-velocity corrections, and JD, MJD, HJD and BJD_TDB time stamps of the observations.
* **Signal-to-noise:** exposure time needed for a target SNR, SNR for a new exposure time, and SNR of co-added spectra.
* **Data:** CSV catalog of the FITS files of an observing campaign, with their header metadata and barycentric corrections.
* **Visibility:** Sun and Moon ephemerides of a night, altitude of several targets along a night (StarAlt-like plot), and altitude and observable hours of a target on every night of a period.
* **Periodic objects:** upcoming eclipses (or any other phase) and whether they are observable, phases at given times, and phase coverage of an observing run.
* **Observations:** azimuth, elevation, airmass, hour and parallactic angles of a target, and the **BinarAlt** browser-based planning interface.

<br>

**Conventions shared by the module:**
* **Coordinates** (`coords`, `targets`) are accepted as:
  * a string: `"05:35:17.3 -05:23:28"`, `"05h35m17.3s -05d23m28s"`, `"05 35 17.3 -05 23 28"` (as copied from SIMBAD) or `"83.82208 -5.39111"`;
  * a tuple or list `(RA, Dec)` whose elements use the same notations (`("05 35 17.3", "-05 23 28")`, `("83.82208", "-5.39111")`), plain numbers (`[83.82208, -5.39111]`) or `astropy` quantities (`(5.5881389 * u.hourangle, -5.39111 * u.deg)`);
  * a dictionary with `'ra'` (or `'right_ascension'`) and `'dec'` (or `'declination'`), and optional `'unit'` and `'frame'` keys;
  * a `SkyCoord` object (ICRS unless another frame is given).

  RA values with `:` or `h` are read in hours; plain numbers, decimal strings and `d`-style strings are read in degrees unless explicit units are given. Dec is always read in degrees.
* **Observatories** (`observatory`) are given by a pre-defined code (case-insensitive) — `"OT"`, `"ORM"`, `"CAHA"`, `"OHP"`, `"LASILLA"`, `"PARANAL"`, `"SAAO"`, `"MAUNAKEA"` — or by a dictionary with `'longitude'` [deg, East +], `'latitude'` [deg, North +], `'elevation'` [m], `'timezone'` (IANA name, e.g. `"Atlantic/Canary"`) and an optional `'name'`.
* Times are in UTC unless stated otherwise. **Nights** are identified by the date on which they start; their local (true solar) midnight, with the Sun at lower culmination, is found from the mean local midnight (00:00 UTC of the next day minus longitude / 15 h), which selects the right night at every longitude.
* **Night definitions** are the same in every function: sunset and sunrise with the Sun at −0.833° (upper limb plus standard refraction); civil, nautical and astronomical twilights at −6°, −12° and −18° (geometric center of the Sun); **dark time** with the Sun below −18°; and a target is **observable** during dark time with altitude ≥ `min_alt`. Crossings are interpolated linearly between grid points (2-min grid in `visibility` and `obs_phase_night`, 5-min grid in `annual_visibility`), so durations and windows do not depend on the grid and agree between functions to a few seconds.
* Periods are in days, and reference epochs (`hjd0`, `jd0`) are (Heliocentric) Julian Dates in UTC; `hjd0` is used if both are given. HJD values follow the classical (UTC) convention; use BJD_TDB from `bary_corr` for timing analyses.
* Angles (altitudes, azimuths, separations) are in degrees, velocities in km/s, exposure times in seconds and the Moon illumination in percent (0–100). Altitudes are geometric (no refraction), and airmass is sec z.
* File arguments (inputs and outputs) accept a `str` or a `pathlib.Path`. Output directories are never created: they must exist.
* Plotting functions take a `show` argument and return the figure in `'fig'`, so that it can be saved afterward (e.g. with [`savefig_pkl`](tex_README.md#savefig_pkl)).
* All functions validate their inputs before computing (`TypeError` for wrong types, `FileNotFoundError` / `IsADirectoryError` / `NotADirectoryError` / `PermissionError` for paths, `ValueError` for invalid values, including empty or non-finite ones and coordinates, observatories or times that cannot be parsed). Error messages name the physical quantity, the argument and its unit (e.g. `Observatory longitude 'observatory['longitude']' must lie in [-180, 180] deg (East positive).`).

---
<br>

## 2. Function Summary

### Coordinates (`COORDINATES`)

| Function | Description |
|---|---|
| [`sky_dist`](#sky_dist) | Angular separation and position angle between two positions. |
| [`coord_info`](#coord_info) | Decimal, sexagesimal and Galactic coordinates of a position. |

<br>

### Time (`TIME`)

| Function | Description |
|---|---|
| [`bary_corr`](#bary_corr) | Barycentric and heliocentric RV corrections, and JD, MJD, HJD and BJD_TDB of the observations. |

<br>

### Signal-to-noise (`SIGNAL-TO-NOISE`)

| Function | Description |
|---|---|
| [`estimate_exposure`](#estimate_exposure) | Exposure time needed to reach a target SNR. |
| [`estimate_SNR`](#estimate_snr) | SNR expected for a new exposure time. |
| [`combine_SNR`](#combine_snr) | SNR of the optimal combination of several spectra. |

<br>

### Data (`DATA`)

| Function | Description |
|---|---|
| [`fits_catalog`](#fits_catalog) | CSV catalog of the FITS files of an instrument, with header metadata and barycentric corrections. |

<br>

### Visibility (`VISIBILITY`)

| Function | Description |
|---|---|
| [`night_info`](#night_info) | Sunset, twilights, midnight, sunrise, night and dark time, LST and Moon of a night. |
| [`visibility`](#visibility) | Altitude of one or several targets along a night, with the night events (plot). |
| [`annual_visibility`](#annual_visibility) | Altitude and observable dark hours of a target on every night of a period (plot). |

<br>

### Periodic objects (`PERIODIC OBJECTS`)

| Function | Description |
|---|---|
| [`periodic_event`](#periodic_event) | Next periodic events (eclipses, quadratures...) and their observability. |
| [`obs_phase`](#obs_phase) | Phases of a periodic object at given times. |
| [`obs_phase_night`](#obs_phase_night) | Phases covered on each night of an observing run (plot). |

<br>

### Observations (`OBSERVATIONS`)

| Function | Description |
|---|---|
| [`telescope_pointing`](#telescope_pointing) | Azimuth, elevation, airmass, hour and parallactic angles of a target (plot). |
| [`binaralt`](#binaralt) | Launch the BinarAlt observation-planning interface in the web browser. |

---
<br>

## 3. Function Reference

Each function lists its arguments in a **Parameters** table and its outputs in a **Returns** table. For functions that return a dictionary, the first row is the dictionary itself and the following rows are its keys (keys of a nested dictionary are written as `'outer'[name]['key']`).

### `sky_dist`
```python
sky_dist(
         coords1 : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
         coords2 : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
    ) -> dict[str, float]
```
Angular separation between two positions (Vincenty formula, accurate at all separations) and position angle of the second one with respect to the first, measured from North through East. If the positions are given in different frames, `coords2` is transformed into the frame of `coords1`. Typical use: checking for visual companions or contaminating sources within the fiber or slit aperture, or matching targets between catalogs.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords1` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Reference position. |
| `coords2` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Second position. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float]` | Dictionary with the keys below. |
| `'arcsec'` | `float` | Angular separation [arcsec]. |
| `'deg'` | `float` | Angular separation [deg]. |
| `'pa'` | `float` | Position angle of `coords2` with respect to `coords1`, from North through East [deg, in [0, 360)]. |

**Raises:** `TypeError` for unsupported types and `ValueError` for empty or unparsable coordinates.

<br>

### `coord_info`
```python
coord_info(
           coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
    ) -> dict[str, float | str]
```
Expresses a position in ICRS decimal degrees, sexagesimal notation and Galactic coordinates (coordinates in another frame are first transformed into ICRS). OB stars are strongly concentrated toward the Galactic plane, so a large |b| for an OB star of known distance can point to a runaway origin (Blaauw 1961).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float \| str]` | Dictionary with the keys below. |
| `'ra'` | `float` | Right ascension (ICRS) [deg]. |
| `'dec'` | `float` | Declination (ICRS) [deg]. |
| `'ra_hms'` | `str` | Right ascension as `hh:mm:ss.ss`. |
| `'dec_dms'` | `str` | Declination as `±dd:mm:ss.s`. |
| `'l'` | `float` | Galactic longitude [deg, in [0, 360)]. |
| `'b'` | `float` | Galactic latitude [deg, in [−90, 90]]. |

**Raises:** `TypeError` for unsupported types and `ValueError` for empty or unparsable coordinates.

<br>

### `bary_corr`
```python
bary_corr(
          coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
          observatory : str | dict,
          obstime     : str | datetime | Time | list[str | datetime | Time],
          exptime     : int | float | list[int | float] | None = None,
          rv          : int | float | list[int | float] | None = None
    ) -> dict[str, SkyCoord | Observer | np.ndarray | list[str] | float]
```
Barycentric (and heliocentric) radial-velocity correction of one or several observations, computed with `astropy` at mid-exposure if `exptime` is given, together with their JD, MJD, HJD (UTC) and BJD_TDB, the recommended time stamp for time-series analysis (Eastman et al. 2010). Measured velocities are corrected as rv_bary = rv + v_bary + rv v_bary / c; to first order, a spectrum is moved to the barycentric frame with λ_bary = λ_obs (1 + v_bary / c). The correction reaches ~30 km/s (Earth's orbit) plus ~0.5 km/s (Earth's rotation) and changes by up to ~0.1 km/s per hour, so the mid-exposure time matters for long exposures. Many reduction pipelines already apply the correction or store it in the header: check before correcting, to avoid applying it twice.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `obstime` | `str \| datetime \| Time \| list[str \| datetime \| Time]` | — | Start of the exposure(s) [UTC], e.g. `DATE-OBS`. |
| `exptime` | `int \| float \| list[int \| float] \| None` | `None` | Exposure time(s) [s]; one value for all observations or one per observation; non-negative. |
| `rv` | `int \| float \| list[int \| float] \| None` | `None` | Measured (topocentric) radial velocity(ies) [km/s]; one value for all or one per observation. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, SkyCoord \| Observer \| np.ndarray \| list[str] \| float]` | Dictionary with the keys below; the arrays have one value per observation. |
| `'coords'` | `SkyCoord` | Input coordinates converted to a `SkyCoord` object. |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'obstime'` | `list[str]` | Times at which the correction is evaluated (mid-exposure if `exptime` is given), as `YYYY-MM-DD HH:MM:SS.sss` [UTC]. |
| `'jd'` | `np.ndarray[float]` | Julian Dates [UTC]. |
| `'mjd'` | `np.ndarray[float]` | Modified Julian Dates, MJD = JD − 2400000.5 [UTC]. |
| `'hjd'` | `np.ndarray[float]` | Heliocentric Julian Dates [UTC]. |
| `'bjd_tdb'` | `np.ndarray[float]` | Barycentric Julian Dates [TDB]. |
| `'v_bary'` | `np.ndarray[float]` | Barycentric radial-velocity correction [km/s]. |
| `'v_helio'` | `np.ndarray[float]` | Heliocentric radial-velocity correction [km/s]. |
| `'rv_bary'` | `np.ndarray[float] \| float` | Barycentric radial velocities (`rv` corrected) [km/s]; `NaN` if `rv` is not given. |

**Raises:** `TypeError` for wrong types and `ValueError` for unparsable coordinates, observatories or times, non-finite values, negative exposure times or lists whose length differs from `obstime`.

<br>

### `estimate_exposure`
```python
estimate_exposure(
                  snr_init   : int | float,
                  t_init     : int | float,
                  snr_target : int | float
    ) -> float
```
Exposure time for a target SNR in the photon-noise (source-dominated) regime, SNR ∝ √t: t_final = t_init (SNR_target / SNR_init)². For faint targets, bright sky or short exposures (read-out and sky-noise regimes) the time is underestimated. The reference exposure should be taken under similar conditions (airmass, seeing, Moon, setup). Inverse of [`estimate_SNR`](#estimate_snr).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snr_init` | `int \| float` | — | SNR of a reference exposure; positive. |
| `t_init` | `int \| float` | — | Exposure time of the reference exposure [any time unit, e.g. s]; positive. |
| `snr_target` | `int \| float` | — | Desired SNR; positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `t_final` | `float` | Exposure time needed to reach `snr_target` [same units as `t_init`]. |

**Raises:** `TypeError` if any input is not a real number and `ValueError` if any input is not finite or not positive.

<br>

### `estimate_SNR`
```python
estimate_SNR(
             snr_init : int | float,
             t_init   : int | float,
             t_final  : int | float
    ) -> float
```
SNR for a new exposure time in the photon-noise regime: SNR_final = SNR_init √(t_final / t_init). In read-out or sky-noise dominated regimes the SNR is overestimated. Inverse of [`estimate_exposure`](#estimate_exposure).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snr_init` | `int \| float` | — | SNR of a reference exposure; positive. |
| `t_init` | `int \| float` | — | Exposure time of the reference exposure [any time unit]; positive. |
| `t_final` | `int \| float` | — | New exposure time [same units as `t_init`]; positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `snr_final` | `float` | Expected SNR for `t_final` [dimensionless]. |

**Raises:** `TypeError` if any input is not a real number and `ValueError` if any input is not finite or not positive.

<br>

### `combine_SNR`
```python
combine_SNR(
            snr_list : list[int | float] | tuple[int | float, ...] | np.ndarray
    ) -> float
```
SNR of the inverse-variance (1/σ²) weighted co-addition of spectra with the same continuum level and independent noise, SNR_comb = √(Σ SNR_i²); N spectra of equal SNR give SNR √N. A plain (unweighted) average of spectra with different SNR gives a lower value, so the result is an upper limit in that case.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snr_list` | `list[int \| float] \| tuple[int \| float, ...] \| np.ndarray` | — | SNR of each spectrum; finite and positive values. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `snr_comb` | `float` | SNR of the combined spectrum [dimensionless]. |

**Raises:** `TypeError` if `snr_list` is not a list, tuple or array or contains non-real values, and `ValueError` if it is empty or contains non-finite or non-positive values.

<br>

### `fits_catalog`
```python
fits_catalog(
             folder      : str | Path,
             instrument  : str,
             output_file : str | Path = "fits_catalog.csv"
    ) -> pd.DataFrame
```
Scans a folder recursively for the reduced spectra of an instrument (HERMES: `*_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fit*`), reads their headers and writes one row per file, sorted by observing night and file name. Exposures taken before the instrument split hour (09 UTC for HERMES) are assigned to the night that started on the previous date. The airmass is read from `AIRMASS` or computed from the telescope altitude (`TELALT`) as sec z. The barycentric correction is computed for every file with [`bary_corr`](#bary_corr), at mid-exposure (`DATE-OBS` + `EXPTIME` / 2), from the header `RA` and `DEC` [deg] and the instrument site; if they are missing, the barycentric values are `NaN` and a warning is issued (a missing `EXPTIME` is taken as 0 s). Unreadable files issue a warning and are skipped, and an existing `output_file` is overwritten with a warning.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `folder` | `str \| Path` | — | Root directory scanned recursively. |
| `instrument` | `str` | — | Instrument (case-insensitive); currently `"HERMES"` (Mercator, ORM). |
| `output_file` | `str \| Path` | `"fits_catalog.csv"` | Output `.csv` file (its directory must exist). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `df` | `pd.DataFrame` | Catalog with one row per FITS file (also written to `output_file`), with the columns `file`, `telescope`, `instrument`, `object`, `ra`, `dec`, `airmass`, `exptime_s`, `filter`, `date_obs`, `night_obs`, `time_obs`, `obstime` (mid-exposure, UTC), `jd`, `mjd`, `hjd`, `bjd_tdb`, `v_bary_kms` and `v_helio_kms`. Missing header values are `N/A` (the object is `Not specified`) and missing barycentric values `NaN`. |

**Raises:** `TypeError` for wrong types; `FileNotFoundError`, `NotADirectoryError`, `IsADirectoryError` or `PermissionError` for invalid paths; `ValueError` if `instrument` is empty or not supported, or `output_file` is not a `.csv` file.

<br>

### `night_info`
```python
night_info(
           observatory : str | dict,
           night       : str | datetime | Time | None = None
    ) -> dict[str, Observer | str | float | None]
```
Sun and Moon ephemerides of an observing night: sunset and sunrise (−0.833°), end and start of the civil (−6°), nautical (−12°) and astronomical (−18°) twilights, local midnight, night length (sunset to sunrise), dark time (Sun below −18°), local sidereal time at midnight (objects with RA ≈ LST transit at midnight), and Moon illumination, altitude, rise and set. Moonrise and moonset are only reported if they happen during the night, as in [`visibility`](#visibility). Bright time (Moon illumination ≳ 50–70 % with the Moon up) is usually enough for high-SNR spectroscopy of bright OB stars, while faint targets need dark or gray time.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `night` | `str \| datetime \| Time \| None` | `None` | Night, identified by its start date (`"YYYY-MM-DD"`); `None` = night starting on the current UTC date. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, Observer \| str \| float \| None]` | Dictionary with the keys below; event times are `YYYY-MM-DD HH:MM:SS` [UTC], or `None` if the event does not happen. |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'night'` | `str` | Night of observation (`YYYY-MM-DD`). |
| `'utc_offset'` | `float` | Local time − UTC at midnight, with daylight saving time [h] (`NaN` if it cannot be computed). |
| `'sunset'` | `str \| None` | Sunset (Sun at −0.833°). |
| `'civil_evening'` | `str \| None` | End of the evening civil twilight (Sun at −6°). |
| `'nautical_evening'` | `str \| None` | End of the evening nautical twilight (Sun at −12°). |
| `'astro_evening'` | `str \| None` | End of the evening astronomical twilight (Sun at −18°). |
| `'midnight'` | `str \| None` | Local (true solar) midnight. |
| `'astro_morning'` | `str \| None` | Start of the morning astronomical twilight (Sun at −18°). |
| `'nautical_morning'` | `str \| None` | Start of the morning nautical twilight (Sun at −12°). |
| `'civil_morning'` | `str \| None` | Start of the morning civil twilight (Sun at −6°). |
| `'sunrise'` | `str \| None` | Sunrise (Sun at −0.833°). |
| `'night_hours'` | `float` | Night length, sunset to sunrise [h] (`NaN` in polar day or night). |
| `'dark_hours'` | `float` | Dark time, between the astronomical twilights [h] (0.0 if the Sun never reaches −18°; 24.0 if it stays below −18° all day). |
| `'lst_midnight'` | `str` | Local sidereal time at midnight (`hh:mm:ss`). |
| `'moon_illumination'` | `float` | Illuminated fraction of the Moon at midnight [%]. |
| `'moon_alt_midnight'` | `float` | Altitude of the Moon at midnight [deg]. |
| `'moonrise'` | `str \| None` | Moonrise within the night (`None` if the Moon does not rise between sunset and sunrise). |
| `'moonset'` | `str \| None` | Moonset within the night (`None` if the Moon does not set between sunset and sunrise). |

**Raises:** `TypeError` for wrong types and `ValueError` if the observatory or `night` cannot be parsed.

<br>

### `visibility`
```python
visibility(
           targets     : dict | str | tuple[str | float, str | float] | list[str | float] | SkyCoord,
           observatory : str | dict,
           night       : str | datetime | Time | None = None,
           min_alt     : int | float | u.Quantity = 30,
           show        : bool = True
    ) -> dict[str, Observer | str | float | np.ndarray | dict | Figure]
```
Altitude of one or several targets from sunset to sunrise (2-min grid; ±12 h around midnight in polar day or night), plotted against UTC with an airmass axis, the Moon and `min_alt`. The day plus civil, nautical and astronomical twilights are shaded (without gaps between areas), the Sun events are marked on top of the plot, and a table below the legend gives the times of the night events, the night length and the dark time. For each target it gives the time of maximum altitude, the observable window (dark time with altitude ≥ `min_alt`, interpolated between grid points) and its duration, and the separation from the Moon at midnight. A dictionary is read as a single target (named `'Target'`) when all its keys are coordinate keys (`'ra'`, `'dec'`, `'right_ascension'`, `'declination'`, `'unit'`, `'frame'`), and as `{name: coordinates}` otherwise. Airmass sec z is accurate to ~1 % for z < 70°.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `targets` | `dict \| str \| tuple[str \| float, str \| float] \| list[str \| float] \| SkyCoord` | — | Coordinates of one target, or `{name: coordinates}` for several. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `night` | `str \| datetime \| Time \| None` | `None` | Night (start date); `None` = night starting on the current UTC date. |
| `min_alt` | `int \| float \| u.Quantity` | `30` | Minimum altitude [deg, in [0, 90]] (30° ≈ airmass 2). |
| `show` | `bool` | `True` | Display the figure with `plt.show()`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, Observer \| str \| float \| np.ndarray \| dict \| Figure]` | Dictionary with the keys below. |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'night'` | `str` | Night of observation (`YYYY-MM-DD`). |
| `'times'` | `np.ndarray[str]` | UTC times of the 2-min grid (`YYYY-MM-DD HH:MM`), from sunset to sunrise. |
| `'moon_alt'` | `np.ndarray[float]` | Altitude of the Moon along the grid [deg]. |
| `'moon_illumination'` | `float` | Illuminated fraction of the Moon at midnight [%]. |
| `'events'` | `dict[str, str \| None]` | Same event keys as `night_info` (`'sunset'` ... `'moonset'`); each value is its UTC time (`YYYY-MM-DD HH:MM`, rounded to the minute), or `None` if it does not happen within the plotted night. |
| `'night_hours'` | `float` | Night length, sunset to sunrise [h] (as in `night_info`). |
| `'dark_hours'` | `float` | Dark time, between the astronomical twilights [h] (as in `night_info`). |
| `'targets'` | `dict[str, dict]` | Dictionary keyed by target name (`name`; `'Target'` for a single unnamed target); each value is a dictionary with the keys below. |
| `'targets'[name]['altitude']` | `np.ndarray[float]` | Altitude along the grid [deg]. |
| `'targets'[name]['airmass']` | `np.ndarray[float]` | Airmass sec z along the grid (`NaN` below the horizon). |
| `'targets'[name]['transit']` | `str` | UTC time of maximum altitude within the night (`YYYY-MM-DD HH:MM`). |
| `'targets'[name]['max_alt']` | `float` | Maximum altitude within the night [deg]. |
| `'targets'[name]['start']` | `str \| None` | UTC time at which the target becomes observable (`None` if never observable). |
| `'targets'[name]['end']` | `str \| None` | UTC time at which the target stops being observable (`None` if never observable). |
| `'targets'[name]['hours']` | `float` | Dark time with altitude ≥ `min_alt` [h]. |
| `'targets'[name]['moon_sep']` | `float` | Separation from the Moon at midnight [deg]. |
| `'fig'` | `Figure` | The `matplotlib` visibility figure. |

**Raises:** `TypeError` for wrong types and `ValueError` for unparsable coordinates, observatories or nights, empty target names, or `min_alt` not finite or outside [0, 90]°.

<br>

### `annual_visibility`
```python
annual_visibility(
                  coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                  observatory : str | dict,
                  start_date  : str | datetime | Time | None = None,
                  n_nights    : int = 365,
                  min_alt     : int | float | u.Quantity = 30,
                  show        : bool = True
    ) -> dict[str, SkyCoord | Observer | list[str] | np.ndarray | str | int | Figure]
```
Altitude and observable dark hours of a target on each night of a period, sampled every 5 min from mean local noon to the next mean local noon (crossings interpolated, so the hours agree with [`visibility`](#visibility) and [`night_info`](#night_info) to a few seconds). The upper panel shows the maximum altitude in dark time and the altitude at local midnight; the lower panel, the observable hours and the total dark time; both show the Moon illumination at local midnight. The date axis runs from `start_date` to `start_date` + `n_nights` days (the morning after the last night), with both limits labeled. Typical use: choosing the semester of a proposal and estimating how many epochs of a binary can be secured. With 365 nights it takes about 20 s.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `start_date` | `str \| datetime \| Time \| None` | `None` | First night (`"YYYY-MM-DD"`); `None` = night starting on the current UTC date. |
| `n_nights` | `int` | `365` | Number of consecutive nights; positive. |
| `min_alt` | `int \| float \| u.Quantity` | `30` | Minimum altitude [deg, in [0, 90]]. |
| `show` | `bool` | `True` | Display the figure with `plt.show()`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, SkyCoord \| Observer \| list[str] \| np.ndarray \| str \| int \| Figure]` | Dictionary with the keys below; the arrays have one value per night. |
| `'coords'` | `SkyCoord` | Input coordinates converted to a `SkyCoord` object. |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'nights'` | `list[str]` | Nights (`YYYY-MM-DD`). |
| `'hours'` | `np.ndarray[float]` | Dark hours with the target above `min_alt` [h]. |
| `'dark_hours'` | `np.ndarray[float]` | Total dark time [h] (0.0 without dark time). |
| `'alt_midnight'` | `np.ndarray[float]` | Altitude of the target at local (true solar) midnight [deg]. |
| `'max_alt'` | `np.ndarray[float]` | Maximum altitude of the target during dark time [deg] (`NaN` without dark time). |
| `'moon_illumination'` | `np.ndarray[float]` | Illuminated fraction of the Moon at local midnight [%]. |
| `'best_night'` | `str` | Night with the most observable hours (the first one in case of ties). |
| `'n_observable'` | `int` | Number of nights with at least 1 observable hour. |
| `'fig'` | `Figure` | The `matplotlib` figure (altitude and observable-hours panels). |

**Raises:** `TypeError` for wrong types and `ValueError` for unparsable coordinates, observatories or dates, `n_nights` not positive, or `min_alt` not finite or outside [0, 90]°.

<br>

### `periodic_event`
```python
periodic_event(
               coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
               period      : int | float,
               observatory : str | dict,
               hjd0        : int | float | None = None,
               jd0         : int | float | None = None,
               start_time  : str | datetime | Time | None = None,
               n_events    : int = 10,
               min_alt     : int | float | u.Quantity = 30,
               phase       : int | float = 0.
    ) -> dict[str, SkyCoord | float | Observer | str | list[str | float | bool | int]]
```
Upcoming events of the linear ephemeris T(E) = T0 + (phase + E) P and whether each one is observable (dark time with the target above `min_alt`). With `hjd0` the ephemeris is linear in HJD and each event is converted to UTC with its own heliocentric light-travel time (up to ±8.3 min); with `jd0` it is applied directly in JD. The uncertainty of the predictions grows as σ_T(E) ≈ √(σ_T0² + E² σ_P²).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |
| `period` | `int \| float` | — | Period [d]; positive. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `hjd0` | `int \| float \| None` | `None` | Reference epoch of phase 0 [HJD, UTC]; used if both epochs are given. |
| `jd0` | `int \| float \| None` | `None` | Reference epoch of phase 0 [JD, UTC]; one of the two epochs is required. |
| `start_time` | `str \| datetime \| Time \| None` | `None` | Start of the search [UTC]; `None` = now. |
| `n_events` | `int` | `10` | Number of events; positive. |
| `min_alt` | `int \| float \| u.Quantity` | `30` | Minimum altitude [deg, in [0, 90]]. |
| `phase` | `int \| float` | `0.` | Phase of the events, in [0, 1): 0 primary eclipse, 0.5 secondary eclipse (circular orbits), 0.25/0.75 quadratures. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, SkyCoord \| float \| Observer \| str \| list[str \| float \| bool \| int]]` | Dictionary with the keys below; the lists have one value per event. |
| `'coords'` | `SkyCoord` | Input coordinates converted to a `SkyCoord` object. |
| `'period'` | `float` | Period [d] (input value). |
| `'hjd0'` | `float` | Reference epoch [HJD, UTC]. |
| `'jd0'` | `float` | Reference epoch [JD, UTC]. |
| `'phase'` | `float` | Phase of the events (input value). |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'start_time'` | `str` | Start of the search (`YYYY-MM-DD HH:MM:SS` UTC). |
| `'cycles'` | `list[int]` | Cycle number E of each event, counted from the reference epoch. |
| `'events'` | `list[str]` | Times of the events (`YYYY-MM-DD HH:MM:SS` UTC). |
| `'events_hjd'` | `list[float]` | Times of the events [HJD, UTC]. |
| `'observable'` | `list[bool]` | `True` if the event happens in dark time (Sun below −18°) with the target above `min_alt`. |
| `'altitude'` | `list[float]` | Altitude of the target at each event [deg]. |

**Raises:** `TypeError` for wrong types and `ValueError` for non-finite or out-of-range values, a missing reference epoch, or unparsable coordinates, observatories or times.

<br>

### `obs_phase`
```python
obs_phase(
          coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
          period      : int | float,
          observatory : str | dict,
          times       : str | datetime | Time | list[str | datetime | Time],
          hjd0        : int | float | None = None,
          jd0         : int | float | None = None
    ) -> dict[str, SkyCoord | float | Observer | list]
```
Phase φ = ((t − T0) mod P) / P of a periodic object at one or several times, using HJD for both t and T0 if `hjd0` is given and JD otherwise.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |
| `period` | `int \| float` | — | Period [d]; positive. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `times` | `str \| datetime \| Time \| list[str \| datetime \| Time]` | — | Observation time(s) [UTC]. |
| `hjd0` | `int \| float \| None` | `None` | Reference epoch [HJD, UTC]; used if both epochs are given. |
| `jd0` | `int \| float \| None` | `None` | Reference epoch [JD, UTC]; one of the two epochs is required. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, SkyCoord \| float \| Observer \| list]` | Dictionary with the keys below; the lists have one value per observation time. |
| `'coords'` | `SkyCoord` | Input coordinates converted to a `SkyCoord` object. |
| `'period'` | `float` | Period [d] (input value). |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'hjd0'` | `float` | Reference epoch [HJD, UTC]. |
| `'jd0'` | `float` | Reference epoch [JD, UTC]. |
| `'times'` | `list[str]` | Observation times (`YYYY-MM-DD HH:MM:SS` UTC). |
| `'times_jd'` | `list[float]` | Observation times [JD, UTC]. |
| `'times_hjd'` | `list[float]` | Observation times [HJD, UTC]. |
| `'phases'` | `list[float]` | Phases, in [0, 1). |

**Raises:** `TypeError` for wrong types and `ValueError` for non-finite or out-of-range values, a missing reference epoch, or unparsable coordinates, observatories or times.

<br>

### `obs_phase_night`
```python
obs_phase_night(
                coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                period      : int | float,
                observatory : str | dict,
                nights      : str | datetime | Time | list[str | datetime | Time],
                hjd0        : int | float | None = None,
                jd0         : int | float | None = None,
                min_alt     : int | float | u.Quantity = 30,
                show        : bool = True
    ) -> dict[str, SkyCoord | float | Observer | list | Figure]
```
Range of phases covered on each night while the target is observable (between the astronomical twilights, 2-min grid, with the `min_alt` crossings interpolated), plotted per night with the UTC start and end of each window, and as the cumulative coverage of the run. Nights whose window is longer than the period are flagged as full coverage. Where the Sun never reaches −18°, a warning is issued and the night is reported as not observable. Useful to check whether a campaign samples the whole orbit (e.g. for SB2 disentangling or light-curve coverage).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |
| `period` | `int \| float` | — | Period [d]; positive. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `nights` | `str \| datetime \| Time \| list[str \| datetime \| Time]` | — | Night(s), identified by their start date. |
| `hjd0` | `int \| float \| None` | `None` | Reference epoch [HJD, UTC]; used if both epochs are given. |
| `jd0` | `int \| float \| None` | `None` | Reference epoch [JD, UTC]; one of the two epochs is required. |
| `min_alt` | `int \| float \| u.Quantity` | `30` | Minimum altitude [deg, in [0, 90]]. |
| `show` | `bool` | `True` | Display the figure with `plt.show()`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, SkyCoord \| float \| Observer \| list \| Figure]` | Dictionary with the keys below; the lists have one value per night. |
| `'coords'` | `SkyCoord` | Input coordinates converted to a `SkyCoord` object. |
| `'period'` | `float` | Period [d] (input value). |
| `'hjd0'` | `float` | Reference epoch [HJD, UTC]. |
| `'jd0'` | `float` | Reference epoch [JD, UTC]. |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'nights'` | `list[str]` | Nights of observation (`YYYY-MM-DD`). |
| `'phases'` | `list[list[float \| None]]` | Phases `[start, end]` covered on each night while the target is observable; wraps through 0 if start > end; `[None, None]` if not observable. |
| `'full_coverage'` | `list[bool]` | `True` if the observable window is longer than the period (every phase covered). |
| `'times'` | `list[list[str \| float]]` | UTC times `[start, end]` of the observable window (`YYYY-MM-DD HH:MM:SS`); `[NaN, NaN]` if not observable. |
| `'fig'` | `Figure` | The `matplotlib` phase-coverage figure. |

**Raises:** `TypeError` for wrong types and `ValueError` for non-finite or out-of-range values, a missing reference epoch, or unparsable coordinates, observatories or nights.

<br>

### `telescope_pointing`
```python
telescope_pointing(
                   coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                   observatory : str | dict,
                   obstime     : str | datetime | Time | None = None,
                   show        : bool = True
    ) -> dict[str, SkyCoord | Time | Observer | float | str | Figure]
```
Azimuth, elevation, cardinal direction, airmass, hour angle, local sidereal time and parallactic angle of a target, plus its separation from the Moon, shown in a compass diagram and an elevation diagram. The hour angle uses the apparent right ascension of date (TETE frame), so precession and nutation since J2000 are taken into account. With a slit spectrograph, align the slit with the parallactic angle at large airmass (≳ 1.3) to avoid wavelength-dependent losses from atmospheric differential refraction (Filippenko 1982).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `coords` | `str \| tuple[str \| float, str \| float] \| list[str \| float] \| dict \| SkyCoord` | — | Target coordinates. |
| `observatory` | `str \| dict` | — | Observatory code or dictionary. |
| `obstime` | `str \| datetime \| Time \| None` | `None` | Observation time [UTC]; `None` = now. |
| `show` | `bool` | `True` | Display the figure with `plt.show()`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, SkyCoord \| Time \| Observer \| float \| str \| Figure]` | Dictionary with the keys below. |
| `'coords'` | `SkyCoord` | Input coordinates converted to a `SkyCoord` object. |
| `'obstime'` | `Time` | Observation time. |
| `'observatory'` | `Observer` | Observatory object (location and timezone). |
| `'elevation'` | `float` | Elevation (altitude) above the horizon [deg]. |
| `'azimuth'` | `float` | Azimuth [deg, 0 = N, 90 = E]. |
| `'direction'` | `str` | Cardinal direction of the pointing (`'N'`, `'NE'`, ...). |
| `'airmass'` | `float` | Airmass sec z (`NaN` below the horizon). |
| `'hour_angle'` | `float` | Hour angle [h, in [−12, 12), negative = East of the meridian]. |
| `'lst'` | `str` | Local apparent sidereal time (`hh:mm:ss`). |
| `'parallactic_angle'` | `float` | Parallactic angle [deg]. |
| `'moon_sep'` | `float` | Angular separation from the Moon [deg]. |
| `'moon_illumination'` | `float` | Illuminated fraction of the Moon [%]. |
| `'fig'` | `Figure` | The `matplotlib` pointing figure (compass and elevation diagrams). |

**Raises:** `TypeError` for wrong types and `ValueError` for unparsable coordinates, observatories or times, or more than one time.

<br>

### `binaralt`
```python
binaralt(
         block : bool = False
    ) -> None
```
Launches **BinarAlt**, a browser-based planning tool with target visibility, ephemerides, lunar information, live observatory weather, exposure scheduling, an interactive exposure timeline and an observation-plan table and logbook. The interface (`binaralt.html`) is served on 127.0.0.1, on a free port chosen by the operating system, and opened in a new tab of the default web browser (not inline in the console or the notebook). It can be launched from a script, the Python console or a local Jupyter notebook, but not from a remote Jupyter server (JupyterHub, Colab, SSH) unless the port is forwarded. It requires a modern web browser with JavaScript enabled. The astronomical calculations run in the browser (Astronomy Engine, bundled in the page, so they work offline); the local server only relays the live weather requests to the observatory websites, which do not allow direct requests from the page. Targets can be imported from and exported to CSV files with the columns `Name`, `RA`, `Dec`, `Period (d)`, `HJD0`, `JD0`, `Nº exposures`, `Exp. (s)`, `Own overhead (s)`, `Own min. alt. (°)` and `Priority (1-5)` (example: [`Demos/data_and_results/obs/data/northcircumpolar_targets.csv`](../Demos/data_and_results/obs/data/northcircumpolar_targets.csv)). A detailed user guide is available in [`BinarAlt_UserGuide.pdf`](BinarAlt_UserGuide.pdf); a copy is installed with the package next to the interface (`src/panorama/support/obs/`), from where it is opened with the question-mark icon in the corner of the interface.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `block` | `bool` | `False` | Keep the server in the foreground until `Ctrl+C` (needed in scripts), or run it in a background thread that stops when the Python session ends. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `None` | `None` | The interface is opened in the default web browser. |

**Raises:** `TypeError` if `block` is not a boolean and `FileNotFoundError` if `binaralt.html` is missing from the installation.

---
<br>

## 4. Dependencies

| Package | Used for |
|---|---|
| `astroplan` | Observer, Sun and Moon events, twilights and target altitudes. |
| `astropy` | Coordinates, times, barycentric corrections and FITS headers. |
| `matplotlib` | Visibility, phase-coverage and pointing figures. |
| `numpy` | Numerical operations and grids. |
| `pandas` | FITS catalogs (`fits_catalog`). |

`binaralt` uses only the Python standard library on the server side; its live weather panel needs internet access. All the other functions work offline (with `astropy`'s bundled ephemerides and Earth-orientation tables).

---
<br>

## 5. Internal Structure

| Path | Content |
|---|---|
| `src/panorama/obs.py` | Public functions, grouped in seven sections: `COORDINATES`, `TIME`, `SIGNAL-TO-NOISE`, `DATA`, `VISIBILITY`, `PERIODIC OBJECTS`, `OBSERVATIONS`. |
| `src/panorama/support/obs/obs_support.py` | Pre-defined observatories (`_OBSERVATORIES`) and instruments (`_INSTRUMENTS`), coordinate parsing (`_coord`), observer creation (`_get_observer`), validation of the arguments shared by several functions (`_check_bool`, `_check_period`, `_check_epochs`, `_check_min_alt`, `_check_night`, `_check_time`, `_check_times`, `_check_per_observation`), conversion between the HJD and JD reference epochs shared by `periodic_event`, `obs_phase` and `obs_phase_night` (`_reference_epochs`), check of astroplan event times (`_is_event`) and their formatting and durations (`_event_iso` and `_event_hours`, used by `night_info`; `_event_within`, used by `night_info` and `_event_minute`; `_event_minute` and `_format_hhmm`, used by `visibility`), phases of a periodic object at given times shared by `obs_phase` and `obs_phase_night` (`_phases`), drawing of phase ranges wrapped through phase 0 for `obs_phase_night` (`_draw_phase_range`), observing night of an exposure (`_compute_obs_night`), FITS header reading (`_get_header_value`), the interpolated observable intervals of a grid shared by `visibility` and `annual_visibility` (`_above`, `_observable`), the mean and true (apparent) solar midnight of a night shared by `night_info`, `visibility`, `annual_visibility` and `obs_phase_night` (`_mean_midnight`, `_solar_midnight`), and the weather service (`_WeatherService`, queried with `_WeatherService._get`) and HTTP request handler (`_BinarAltRequestHandler`, whose `do_GET` and `log_message` keep the names imposed by `http.server`) of BinarAlt. |
| `src/panorama/support/obs/binaralt.html` | BinarAlt web interface. |
| `READMEs/BinarAlt_UserGuide.pdf` | BinarAlt user guide (documentation copy, linked from this README). |
| `src/panorama/support/units/units_support.py` | Validation of real numbers (`_check_real`), reused for exposure times, signal-to-noise ratios, periods and phases. |
| `src/panorama/support/tex/tex_support.py` | Validation of output files (`_check_output_file`), reused by `fits_catalog`. |

---
<br>

## 6. Demo

The Jupyter notebook [`Demos/obs_demo.ipynb`](../Demos/obs_demo.ipynb) runs every function of the module in 13 examples, one section per example:

* **Examples 1-2:** every accepted coordinate format (all giving the same position), the Galactic coordinates of the runaway star ζ Oph, and angular separations and position angles, also between positions given in different frames (ICRS and Galactic).
* **Examples 3-4:** barycentric correction and BJD_TDB of three spectra of ρ Leo (one exposure time per spectrum, mock topocentric velocities corrected back to a constant barycentric value, and a single observation without velocities), and exposure-time and SNR scaling, with co-added spectra of equal and different SNR.
* **Example 5:** catalog of *mock* HERMES FITS files (headers only, distributed in `Demos/data_and_results/obs/data/observation_fits/`) with barycentric corrections: an early-morning exposure assigned to the previous night, an airmass computed from `TELALT`, a header without `DATE-OBS` (no correction, with its warning), an unreadable file (skipped, with its warning) and a calibration file ignored by the pattern; the warnings issued by `fits_catalog` are printed.
* **Examples 6-8:** the ephemerides of a night at the ORM (all twilights, moonrise only within the night) and at a high-latitude site without astronomical night (dark time 0 h), the visibility of three OB stars over that night (with a target never observable) and of a single target given as a coordinate dictionary with `min_alt` as a `Quantity`, and the visibility of ρ Leo over a year (365 nights).
* **Examples 9-11:** for a short-period eclipsing binary, the upcoming primary and secondary eclipses (with the epoch given as HJD and as JD), its phases at given times (HJD and JD epochs compared), and the phase coverage of a 6-night run, plus a night with full coverage (short period) and a night without astronomical night (with its warning).
* **Examples 12-13:** a telescope pointing, and the BinarAlt input check; the interface is launched only if `RUN_BINARALT = True` in the setup cell (it opens a browser tab; the server runs in the background, `block=False`, and stops with the kernel). The file `Demos/data_and_results/obs/data/northcircumpolar_targets.csv` is an example target list for the *Import CSV* button of BinarAlt.

Input files (the mock HERMES FITS files and the BinarAlt example target list) are read from `Demos/data_and_results/obs/data/` and outputs (CSV catalog, and figures saved with `savefig_pkl` as PDF + pickle) to `Demos/data_and_results/obs/results/`. No internet access is needed (only the live weather panel of BinarAlt uses it). The visibility, annual-visibility, phase-coverage and pointing figures are shown inline (`%matplotlib inline`). The notebook takes about 40 s (mostly example 8).

```bash
pip install -e ".[demos]"     # PANORAMA + JupyterLab and PyQt6 (demo extras)
cd Demos
jupyter lab obs_demo.ipynb    # Or: jupyter notebook obs_demo.ipynb
```

---
<br>

## 7. References

* Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022, *ApJ*, 935, 167 — `astropy` (coordinates, times, barycentric corrections, FITS).
* Blaauw, A. 1961, *BAN*, 15, 265 — runaway OB stars.
* Cross, D. — *Astronomy Engine* (astronomical calculations of BinarAlt), [github.com/cosinekitty/astronomy](https://github.com/cosinekitty/astronomy).
* Eastman, J., Siverd, R. & Gaudi, B. S. 2010, *PASP*, 122, 935 — BJD_TDB as the standard time stamp.
* Filippenko, A. V. 1982, *PASP*, 94, 715 — atmospheric differential refraction and the parallactic angle.
* Hunter, J. D. 2007, *CSE*, 9, 90 — `matplotlib`.
* Morris, B. M., Tollerud, E., Sipőcz, B. et al. 2018, *AJ*, 155, 128 — `astroplan`.
* Raskin, G., van Winckel, H., Hensberge, H. et al. 2011, *A&A*, 526, A69 — the HERMES spectrograph.
* Wright, J. T. & Eastman, J. D. 2014, *PASP*, 126, 838 — barycentric corrections at the cm/s level.
