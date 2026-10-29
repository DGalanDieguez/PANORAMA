# `panorama.units` — Physical Conversions and Derived Stellar Parameters

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](../LICENSE)

> Doppler shifts, vacuum/air wavelengths, resolution limits, characteristic stellar velocities and derived stellar parameters (luminosities, radii, masses, centrifugal and radiative corrections of log g, and Eddington/ΩΓ limits).

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

`units` gathers the small physical conversions that appear constantly in the spectroscopic analysis of (massive) stars:

* **Wavelength:** wavelength shifts ↔ radial velocities (non-relativistic or relativistic Doppler formula) and vacuum ↔ air wavelengths (Edlén 1966).
* **Velocity:** resolving power ↔ velocity resolution limit, thermal and microturbulent Doppler velocities, and Newtonian and radiation-corrected escape velocities.
* **Stellar parameters:** Stefan-Boltzmann luminosities and radii, spectroscopic masses, the centrifugal and radiative corrections of an effective log g, the spectroscopic luminosity ℒ = Teff⁴/g, bolometric magnitudes, and the Eddington factor with the radiative and centrifugal accelerations (Eddington and ΩΓ limits).

<br>

**Conventions shared by the module:**
* Wavelengths in Å; velocities in km/s, with positive radial velocities for a receding source (redshift).
* Radii, masses and luminosities in solar units (R☉, M☉, L☉), using the IAU 2015 nominal solar values provided by `astropy`.
* Surface gravities as log10(g / [cm s⁻²]); accelerations in cm s⁻².
* **Gravity notation:** `g_eff` is the effective gravity (reduced by rotation and/or radiation), `g_grav` the gravitational (Newtonian) acceleration, and `g_cen` and `g_rad` the centrifugal and radiative accelerations, so that g_grav = g_eff + g_cen + g_rad. Each has its logarithm (`logg_eff`, `logg_grav`, `logg_cen`, `logg_rad`).
* **Newtonian log g:** `loggradius2mass`, `spec_luminosity` and `eddington` expect the gravitational log g, not reduced by rotation or radiation. A fitted (effective) gravity must be corrected first with `logg_centrifugal`, `logg_radiative` or `logg_rad_centrifugal`, passing their `'logg_grav'`.
* All functions work on scalar inputs. Functions with a single result return a `float`; functions with several results return a `dict` with the inputs and the results.
* All functions validate their inputs before computing: `TypeError` for non-numeric inputs (including `bool`), `ValueError` for non-finite (`nan`, `inf`) or out-of-range values. Error messages name the physical quantity, the argument and its unit (e.g. `Rest wavelength 'wl0' must be positive.`).
* Warnings (`UserWarning`) flag results that are computed but may be inaccurate: non-relativistic Doppler shifts with |v| ≥ 1 % of c, Edlén conversions outside 2000–10000 Å, and stars beyond the ΩΓ limit in `eddington`.

The treatment of rotation is first order and spherical (statistical v sin i correction of log g; ΩΓ limit at the equator). The full Roche geometry with gravity darkening will be provided by the `rot` module (under development).

---
<br>

## 2. Function Summary

### Wavelength (`WAVELENGTH`)

| Function | Description |
|---|---|
| [`wl2kmps`](#wl2kmps) | Wavelength shift → radial velocity (non-relativistic Doppler). |
| [`kmps2wl`](#kmps2wl) | Radial velocity → wavelength shift (non-relativistic Doppler). |
| [`wl2kmps_rel`](#wl2kmps_rel) | Wavelength shift → radial velocity (relativistic Doppler). |
| [`kmps2wl_rel`](#kmps2wl_rel) | Radial velocity → wavelength shift (relativistic Doppler). |
| [`vac2air`](#vac2air) | Vacuum → air wavelength (Edlén 1966). |
| [`air2vac`](#air2vac) | Air → vacuum wavelength (Edlén 1966, exact inverse of `vac2air`). |

<br>

### Velocity (`VELOCITY`)

| Function | Description |
|---|---|
| [`R2kmps`](#r2kmps) | Resolving power → velocity resolution limit. |
| [`kmps2R`](#kmps2r) | Velocity resolution limit → required resolving power. |
| [`vthermal`](#vthermal) | Thermal velocity of an ion and total Doppler velocity with microturbulence. |
| [`vescape`](#vescape) | Newtonian and effective (radiation-corrected) escape velocities. |

<br>

### Stellar parameters (`STELLAR PARAMETERS`)

| Function | Description |
|---|---|
| [`teffradius2logL`](#teffradius2logl) | Teff and radius → log(L/L☉) (Stefan-Boltzmann law). |
| [`logLteff2radius`](#loglteff2radius) | log(L/L☉) and Teff → radius (Stefan-Boltzmann law). |
| [`loggradius2mass`](#loggradius2mass) | log g and radius → spectroscopic mass. |
| [`logg_centrifugal`](#logg_centrifugal) | Centrifugal correction of an effective log g (rotation). |
| [`logg_radiative`](#logg_radiative) | Radiative correction of an effective log g (electron scattering). |
| [`logg_rad_centrifugal`](#logg_rad_centrifugal) | Combined radiative and centrifugal correction of an effective log g. |
| [`spec_luminosity`](#spec_luminosity) | Spectroscopic luminosity ℒ = Teff⁴/g in solar units. |
| [`logL2Mbol`](#logl2mbol) | log(L/L☉) → absolute bolometric magnitude. |
| [`Mbol2logL`](#mbol2logl) | Absolute bolometric magnitude → log(L/L☉). |
| [`eddington`](#eddington) | Eddington factor, radiative and centrifugal accelerations, effective gravity and Eddington/ΩΓ limits. |

---
<br>

## 3. Function Reference

Each function lists its arguments in an **Inputs** table and its outputs in a **Returns** table. For functions that return a dictionary, the first row is the dictionary itself and the following rows are its keys.

### `wl2kmps`
```python
wl2kmps(
        dwl : int | float,
        wl0 : int | float
    ) -> float
```
Converts a wavelength shift into a radial velocity with the first-order Doppler formula, v = (Δλ / λ0) c. Positive shifts are redshifts (receding source). Δλ and λ0 must refer to the same medium (both in air or both in vacuum). A warning is issued when |v| ≥ 1 % of c; in that regime use `wl2kmps_rel`. Inverse of `kmps2wl`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dwl` | `int \| float` | — | Wavelength shift Δλ = λ_obs − λ0 [Å]; λ0 + Δλ must be positive. |
| `wl0` | `int \| float` | — | Rest wavelength λ0 [Å]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `dkmps` | `float` | Radial velocity [km/s]; positive for a receding source. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `wl0` is not positive or the observed wavelength λ0 + Δλ is not positive.

<br>

### `kmps2wl`
```python
kmps2wl(
        dkmps : int | float,
        wl0   : int | float
    ) -> float
```
Converts a radial velocity into a wavelength shift with the first-order Doppler formula, Δλ = (v / c) λ0. Also useful to express velocity widths (v sin i, macroturbulence, instrumental FWHM) in wavelength units at a given line. A warning is issued when |v| ≥ 1 % of c; in that regime use `kmps2wl_rel`. Inverse of `wl2kmps`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dkmps` | `int \| float` | — | Radial velocity [km/s]. |
| `wl0` | `int \| float` | — | Rest wavelength λ0 [Å]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `dwl` | `float` | Wavelength shift Δλ = λ_obs − λ0 [Å]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite or `wl0` is not positive.

<br>

### `wl2kmps_rel`
```python
wl2kmps_rel(
            dwl : int | float,
            wl0 : int | float
    ) -> float
```
Converts a wavelength shift into a radial velocity with the special-relativistic Doppler formula for radial motion, β = ((1 + z)² − 1) / ((1 + z)² + 1), with z = Δλ / λ0 and v = β c. Any positive observed wavelength gives −1 < β < 1. Exact inverse of `kmps2wl_rel`; reduces to `wl2kmps` for |v| ≪ c.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dwl` | `int \| float` | — | Wavelength shift Δλ = λ_obs − λ0 [Å]; λ0 + Δλ must be positive. |
| `wl0` | `int \| float` | — | Rest wavelength λ0 [Å]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `dkmps` | `float` | Radial velocity [km/s]; positive for a receding source. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `wl0` is not positive, the observed wavelength is not positive, or |β| ≥ 1 in floating-point arithmetic.

<br>

### `kmps2wl_rel`
```python
kmps2wl_rel(
            dkmps : int | float,
            wl0   : int | float
    ) -> float
```
Converts a radial velocity into a wavelength shift with the special-relativistic Doppler formula, z = √((1 + β) / (1 − β)) − 1 and Δλ = z λ0. Exact inverse of `wl2kmps_rel`; reduces to `kmps2wl` for |v| ≪ c.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dkmps` | `int \| float` | — | Radial velocity [km/s]; must satisfy \|v\| < c. |
| `wl0` | `int \| float` | — | Rest wavelength λ0 [Å]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `dwl` | `float` | Wavelength shift Δλ = λ_obs − λ0 [Å]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `wl0` is not positive or \|v\| ≥ c.

<br>

### `vac2air`
```python
vac2air(
        wl_vac : int | float
    ) -> float
```
Converts a vacuum wavelength into its value in standard air (15 °C, 101325 Pa, dry, 0.03 % CO₂), λ_air = λ_vac / n, with the Edlén (1966) refractive index evaluated at the vacuum wavenumber. Accurate between 2000 and 10000 Å; a warning is issued outside this range (the formula diverges near 1603 Å). Inverse of `air2vac`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `wl_vac` | `int \| float` | — | Vacuum wavelength [Å]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `wl_air` | `float` | Wavelength in standard air [Å]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if `wl_vac` is not finite or not positive.

<br>

### `air2vac`
```python
air2vac(
        wl_air : int | float
    ) -> float
```
Converts an air wavelength into its vacuum value, λ_vac = λ_air n(λ_vac), with the Edlén (1966) refractive index. Since n depends on the unknown λ_vac, the equation is solved by fixed-point iteration (five iterations), so that `air2vac` is the inverse of `vac2air` to floating-point precision in 2000–10000 Å. A warning is issued outside this range.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `wl_air` | `int \| float` | — | Air wavelength [Å]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `wl_vac` | `float` | Wavelength in vacuum [Å]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if `wl_air` is not finite or not positive.

<br>

### `R2kmps`
```python
R2kmps(
       R : int | float
    ) -> float
```
Converts a spectral resolving power into a velocity resolution limit, 1.32 c / R. Here c / R is the velocity FWHM of the instrumental profile and 1.32 is an empirical factor giving the smallest broadening that can be considered resolved (e.g. a lower limit for a measurable v sin i). For R = 46000 it gives ≈ 8.6 km/s.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `R` | `int \| float` | — | Resolving power λ/Δλ [dimensionless]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `limit_kmps` | `float` | Velocity resolution limit, 1.32 c / R [km/s]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if `R` is not finite or not positive.

<br>

### `kmps2R`
```python
kmps2R(
       dkmps : int | float
    ) -> float
```
Inverse of `R2kmps` with the same empirical factor, R = 1.32 c / v: the resolving power required to resolve a given velocity (e.g. ≈ 39600 for v sin i = 10 km/s).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `dkmps` | `int \| float` | — | Velocity resolution limit [km/s]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `R` | `float` | Required resolving power λ/Δλ [dimensionless]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if `dkmps` is not finite or not positive.

<br>

### `vthermal`
```python
vthermal(
         teff        : int | float,
         atomic_mass : int | float,
         vmic        : int | float = 0.
    ) -> dict[str, float]
```
Computes the thermal (most probable) velocity of an ion, v_th = √(2 k_B T / (A m_u)), and the total Doppler velocity with microturbulence added in quadrature, v_D = √(v_th² + ξ²). The Doppler width is Δλ_D = (v_D / c) λ0 (see `kmps2wl`) and the Gaussian FWHM is 2 √(ln 2) Δλ_D. Macroturbulence is deliberately not included: like rotation, it is applied as a convolution of the emergent flux profile and its radial-tangential profile is not Gaussian. In OB stars thermal broadening is only relevant for H and He lines (v_th(H, 30 kK) ≈ 22 km/s vs. v_th(Si, 30 kK) ≈ 4 km/s).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `teff` | `int \| float` | — | Temperature [K], usually the effective temperature; must be positive. |
| `atomic_mass` | `int \| float` | — | Atomic mass of the absorbing species [u] (e.g. 1.008 H, 4.003 He, 28.086 Si); must be positive. |
| `vmic` | `int \| float` | `0.` | Microturbulent velocity ξ [km/s]; must be non-negative. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float]` | Dictionary with the keys below. |
| `'teff'` | `float` | Temperature [K] (input value). |
| `'atomic_mass'` | `float` | Atomic mass [u] (input value). |
| `'vmic'` | `float` | Microturbulent velocity [km/s] (input value). |
| `'v_th'` | `float` | Thermal (most probable) velocity [km/s]. |
| `'v_doppler'` | `float` | Total Doppler velocity, thermal and microturbulent motions added in quadrature [km/s]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `teff` or `atomic_mass` is not positive, or `vmic` is negative.

<br>

### `vescape`
```python
vescape(
        mass   : int | float,
        radius : int | float,
        gamma  : int | float = 0.
    ) -> dict[str, float]
```
Computes the Newtonian photospheric escape velocity, v_esc = √(2 G M / R), and the effective escape velocity reduced by the radiative acceleration on free electrons, v_esc,eff = v_esc √(1 − Γ_e). For line-driven winds, v_∞ / v_esc,eff ≈ 2.6 for Teff ≳ 21 kK and ≈ 1.3 for cooler stars (bi-stability jump; Lamers et al. 1995).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mass` | `int \| float` | — | Stellar mass [M☉]; must be positive. |
| `radius` | `int \| float` | — | Stellar radius [R☉]; must be positive. |
| `gamma` | `int \| float` | `0.` | Electron-scattering Eddington factor Γ_e [dimensionless], in [0, 1) (see `eddington`). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float]` | Dictionary with the keys below. |
| `'mass'` | `float` | Stellar mass [M☉] (input value). |
| `'radius'` | `float` | Stellar radius [R☉] (input value). |
| `'gamma'` | `float` | Eddington factor Γ_e [dimensionless] (input value). |
| `'v_esc'` | `float` | Newtonian escape velocity [km/s]. |
| `'v_esc_eff'` | `float` | Effective escape velocity, reduced by the radiative acceleration on free electrons [km/s]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `mass` or `radius` is not positive, or `gamma` is outside [0, 1).

<br>

### `teffradius2logL`
```python
teffradius2logL(
                teff   : int | float,
                radius : int | float
    ) -> float
```
Bolometric luminosity from the Stefan-Boltzmann law, L = 4π R² σ Teff⁴, with the IAU 2015 nominal L☉ = 3.828 × 10²⁶ W (which implies Teff☉ ≈ 5772 K). Inverse of `logLteff2radius`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `teff` | `int \| float` | — | Effective temperature [K]; must be positive. |
| `radius` | `int \| float` | — | Stellar radius [R☉]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `logL` | `float` | Bolometric luminosity, log(L/L☉). |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite or not positive.

<br>

### `logLteff2radius`
```python
logLteff2radius(
                logL : int | float,
                teff : int | float
    ) -> float
```
Stellar radius from the Stefan-Boltzmann law, R = √(L / (4π σ Teff⁴)). Typical use: a photometric luminosity (Gaia distance + bolometric correction) combined with a spectroscopic Teff. Inverse of `teffradius2logL`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logL` | `int \| float` | — | Bolometric luminosity, log(L/L☉). |
| `teff` | `int \| float` | — | Effective temperature [K]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `radius` | `float` | Stellar radius [R☉]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite or `teff` is not positive.

<br>

### `loggradius2mass`
```python
loggradius2mass(
                logg   : int | float,
                radius : int | float
    ) -> float
```
Spectroscopic mass, M = g R² / G. `logg` must be the gravitational (Newtonian) gravity: a gravity reduced by the centrifugal acceleration (rotating stars, see `logg_centrifugal`) or by radiation pressure (if not included in the model atmosphere, see `logg_radiative`) underestimates the mass. Spectroscopic masses of OB stars are often systematically lower than evolutionary ones (mass discrepancy; Herrero et al. 1992).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logg` | `int \| float` | — | Gravitational (Newtonian) surface gravity, log10(g / [cm s⁻²]); a fitted gravity must be corrected first (`logg_centrifugal`, `logg_radiative`, `logg_rad_centrifugal`). |
| `radius` | `int \| float` | — | Stellar radius [R☉]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `mass` | `float` | Spectroscopic mass [M☉]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite or `radius` is not positive.

<br>

### `logg_centrifugal`
```python
logg_centrifugal(
                 logg   : int | float,
                 vsini  : int | float,
                 radius : int | float
    ) -> dict[str, float]
```
Corrects a spectroscopic (effective) gravity for the centrifugal acceleration with the statistical first-order correction g_grav = g_eff + (v sin i)² / R (Herrero et al. 1992; Repolust et al. 2004). Since v sin i ≤ v_eq, the correction is a lower limit; it is ~0.01 dex for slow rotators and can reach ~0.1 dex for v sin i ≳ 250 km/s in giants and supergiants. Use `'logg_grav'` for spectroscopic masses and as input of `eddington`. Only the centrifugal term is corrected; see `logg_radiative` and `logg_rad_centrifugal` when the gravity is also reduced by radiation pressure.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logg` | `int \| float` | — | Spectroscopic (effective) surface gravity, log10(g / [cm s⁻²]). |
| `vsini` | `int \| float` | — | Projected rotational velocity [km/s]; must be non-negative. |
| `radius` | `int \| float` | — | Stellar radius [R☉]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float]` | Dictionary with the keys below. |
| `'logg'` | `float` | Input (effective) surface gravity, log10(g / [cm s⁻²]). |
| `'vsini'` | `float` | Projected rotational velocity [km/s] (input value). |
| `'radius'` | `float` | Stellar radius [R☉] (input value). |
| `'g_eff'` | `float` | Effective surface gravity, 10^logg [cm s⁻²]. |
| `'g_cen'` | `float` | Centrifugal acceleration, (v sin i)² / R [cm s⁻²]. |
| `'logg_cen'` | `float` | log10(g_cen / [cm s⁻²]) (−inf if v sin i = 0). |
| `'g_grav'` | `float` | Gravitational acceleration, g_eff + g_cen [cm s⁻²]. |
| `'logg_grav'` | `float` | Centrifugally corrected surface gravity, log10(g_grav / [cm s⁻²]). |
| `'dlogg'` | `float` | Correction applied, logg_grav − logg [dex, non-negative]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `vsini` is negative, `radius` is not positive, or the computed accelerations are not finite.

<br>

### `logg_radiative`
```python
logg_radiative(
               logg   : int | float,
               teff   : int | float,
               H_frac : int | float = 0.7
    ) -> dict[str, float]
```
Corrects an effective gravity for the radiative acceleration on free electrons, g_grav = g_eff + g_rad with g_rad = κ_e σ Teff⁴ / c and κ_e = 0.20 (1 + X) cm² g⁻¹ (the same expression as in `eddington`). g_rad depends only on Teff and X, not on the mass or the radius, and is a lower limit to the total radiative acceleration. **Apply it only when `logg` does not already include radiation pressure:** model atmospheres that include it in the hydrostatic equilibrium (e.g. FASTWIND, CMFGEN, TLUSTY, ATLAS) already return the gravitational gravity with respect to radiation, and correcting them again counts g_rad twice. The correction is large for hot, luminous stars (g_rad ≈ 520 cm s⁻² at 30 kK, ≈ +0.18 dex at log g = 3.0).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logg` | `int \| float` | — | Effective surface gravity (reduced by radiation pressure), log10(g / [cm s⁻²]). |
| `teff` | `int \| float` | — | Effective temperature [K]; must be positive. |
| `H_frac` | `int \| float` | `0.7` | Hydrogen mass fraction X [dimensionless], in [0, 1]. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float]` | Dictionary with the keys below. |
| `'logg'` | `float` | Input (effective) surface gravity, log10(g / [cm s⁻²]). |
| `'teff'` | `float` | Effective temperature [K] (input value). |
| `'H_frac'` | `float` | Hydrogen mass fraction X [dimensionless] (input value). |
| `'g_eff'` | `float` | Effective surface gravity, 10^logg [cm s⁻²]. |
| `'g_rad'` | `float` | Radiative acceleration on free electrons, κ_e σ Teff⁴ / c [cm s⁻²]. |
| `'logg_rad'` | `float` | log10(g_rad / [cm s⁻²]). |
| `'g_grav'` | `float` | Gravitational acceleration, g_eff + g_rad [cm s⁻²]. |
| `'logg_grav'` | `float` | Radiatively corrected surface gravity, log10(g_grav / [cm s⁻²]). |
| `'dlogg'` | `float` | Correction applied, logg_grav − logg [dex, non-negative]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `teff` is not positive, `H_frac` is outside [0, 1], or the computed accelerations are not finite.

<br>

### `logg_rad_centrifugal`
```python
logg_rad_centrifugal(
                     logg   : int | float,
                     teff   : int | float,
                     vsini  : int | float,
                     radius : int | float,
                     H_frac : int | float = 0.7
    ) -> dict[str, float]
```
Corrects an effective gravity for both the radiative and the centrifugal accelerations, g_grav = g_eff + g_rad + (v sin i)² / R, reusing `logg_radiative` and `logg_centrifugal` so that the three corrections are always consistent. The individual corrections in dex are not additive; compare `'g_rad'` and `'g_cen'` instead. Both terms are lower limits. The same applicability caveat as `logg_radiative` holds: for gravities fitted with FASTWIND, CMFGEN, TLUSTY or ATLAS, use `logg_centrifugal` alone. Passing `'logg_grav'` to `eddington` with the same `vsini` recovers g_eff = 10^logg.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logg` | `int \| float` | — | Effective surface gravity (reduced by radiation and rotation), log10(g / [cm s⁻²]). |
| `teff` | `int \| float` | — | Effective temperature [K]; must be positive. |
| `vsini` | `int \| float` | — | Projected rotational velocity [km/s]; must be non-negative. |
| `radius` | `int \| float` | — | Stellar radius [R☉]; must be positive. |
| `H_frac` | `int \| float` | `0.7` | Hydrogen mass fraction X [dimensionless], in [0, 1]. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float]` | Dictionary with the keys below. |
| `'logg'` | `float` | Input (effective) surface gravity, log10(g / [cm s⁻²]). |
| `'teff'` | `float` | Effective temperature [K] (input value). |
| `'vsini'` | `float` | Projected rotational velocity [km/s] (input value). |
| `'radius'` | `float` | Stellar radius [R☉] (input value). |
| `'H_frac'` | `float` | Hydrogen mass fraction X [dimensionless] (input value). |
| `'g_eff'` | `float` | Effective surface gravity, 10^logg [cm s⁻²]. |
| `'g_rad'` | `float` | Radiative acceleration on free electrons, κ_e σ Teff⁴ / c [cm s⁻²]. |
| `'logg_rad'` | `float` | log10(g_rad / [cm s⁻²]). |
| `'g_cen'` | `float` | Centrifugal acceleration, (v sin i)² / R [cm s⁻²]. |
| `'logg_cen'` | `float` | log10(g_cen / [cm s⁻²]) (−inf if v sin i = 0). |
| `'g_grav'` | `float` | Gravitational acceleration, g_eff + g_rad + g_cen [cm s⁻²]. |
| `'logg_grav'` | `float` | Radiatively and centrifugally corrected surface gravity, log10(g_grav / [cm s⁻²]). |
| `'dlogg'` | `float` | Total correction applied, logg_grav − logg [dex, non-negative]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `teff` or `radius` is not positive, `vsini` is negative, `H_frac` is outside [0, 1], or the computed accelerations are not finite.

<br>

### `spec_luminosity`
```python
spec_luminosity(
                teff : int | float,
                logg : int | float
    ) -> float
```
Spectroscopic luminosity ℒ = Teff⁴ / g ∝ L / M (Langer & Kudritzki 2014), in solar units: log(ℒ/ℒ☉) = 4 log Teff − log g − log ℒ☉, with log ℒ☉ ≈ 10.61 derived from the IAU 2015 nominal solar values. It is distance independent, proportional to the electron-scattering Eddington factor, and the vertical axis of the spectroscopic HR diagram (sHRD).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `teff` | `int \| float` | — | Effective temperature [K]; must be positive. |
| `logg` | `int \| float` | — | Gravitational (Newtonian) surface gravity, log10(g / [cm s⁻²]); a fitted gravity must be corrected first (`logg_centrifugal`, `logg_radiative`, `logg_rad_centrifugal`). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `logLspec` | `float` | Spectroscopic luminosity, log(ℒ/ℒ☉). |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite or `teff` is not positive.

<br>

### `logL2Mbol`
```python
logL2Mbol(
          logL : int | float
    ) -> float
```
Absolute bolometric magnitude, Mbol = 4.74 − 2.5 log(L/L☉), with the IAU 2015 Resolution B2 zero point (Mbol☉ = 4.74 mag). Combined with a bolometric correction, Mbol = M_V + BC_V. Inverse of `Mbol2logL`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `logL` | `int \| float` | — | Bolometric luminosity, log(L/L☉). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `Mbol` | `float` | Absolute bolometric magnitude [mag]. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if `logL` is not finite.

<br>

### `Mbol2logL`
```python
Mbol2logL(
          Mbol : int | float
    ) -> float
```
Bolometric luminosity, log(L/L☉) = (4.74 − Mbol) / 2.5, with the IAU 2015 Resolution B2 zero point. Inverse of `logL2Mbol`.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `Mbol` | `int \| float` | — | Absolute bolometric magnitude [mag]. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `logL` | `float` | Bolometric luminosity, log(L/L☉). |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if `Mbol` is not finite.

<br>

### `eddington`
```python
eddington(
          teff   : int | float,
          logg   : int | float,
          radius : int | float,
          mass   : int | float | None = None,
          H_frac : int | float = 0.7,
          vsini  : int | float = 0.
    ) -> dict[str, float | bool]
```
Computes the Eddington luminosity L_edd = 4π G M c / κ_e, the bolometric luminosity, the electron-scattering Eddington factor Γ_e = L_bol / L_edd, the radiative acceleration g_rad = κ_e σ Teff⁴ / c, the centrifugal acceleration g_cen = (v sin i)² / R (same expression as `logg_centrifugal`) and the effective gravity g_eff = g_grav − g_rad − g_cen. The opacity is that of electron scattering in a fully ionized H/He plasma, κ_e = 0.20 (1 + X) cm² g⁻¹, so Γ_e is a lower limit to the true Eddington factor. If no mass is given, the spectroscopic mass M = g_grav R² / G is used; in that case g_eff ≤ 0 is equivalent to the ΩΓ limit, Γ_e + (v sin i / v_kep)² ≥ 1 (Langer 1997, 1998). A warning is issued beyond the ΩΓ limit, where `'logg_eff'` is NaN.

`logg` must be the **gravitational (Newtonian)** gravity. A fitted gravity is already reduced by the centrifugal acceleration (and, if the model atmosphere does not include radiation pressure, by g_rad): correct it first with `logg_centrifugal` (or `logg_rad_centrifugal`) and pass its `'logg_grav'`, otherwise g_cen (or g_rad) is subtracted twice. With that workflow `'g_eff'` recovers the fitted gravity. The treatment is spherical (`radius` taken as the equatorial radius, uniform flux) and v sin i ≤ v_eq, so g_cen is a lower limit to the equatorial centrifugal acceleration.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `teff` | `int \| float` | — | Effective temperature [K]; must be positive. |
| `logg` | `int \| float` | — | Gravitational (Newtonian) surface gravity, log10(g / [cm s⁻²]); a fitted gravity must be corrected first (`logg_centrifugal`, `logg_radiative`, `logg_rad_centrifugal`). |
| `radius` | `int \| float` | — | Stellar radius [R☉]; must be positive. |
| `mass` | `int \| float \| None` | `None` | Stellar mass [M☉]; if `None`, the spectroscopic mass M = g R² / G. |
| `H_frac` | `int \| float` | `0.7` | Hydrogen mass fraction X [dimensionless], in [0, 1]. |
| `vsini` | `int \| float` | `0.` | Projected rotational velocity [km/s]; must be non-negative. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, float \| bool]` | Dictionary with the keys below. |
| `'teff'` | `float` | Effective temperature [K] (input value). |
| `'logg'` | `float` | Newtonian surface gravity, log10(g / [cm s⁻²]) (input value). |
| `'radius'` | `float` | Stellar radius [R☉] (input value). |
| `'mass'` | `float` | Stellar mass [M☉] (input value, or the spectroscopic mass if `None` was given). |
| `'H_frac'` | `float` | Hydrogen mass fraction X [dimensionless] (input value). |
| `'vsini'` | `float` | Projected rotational velocity [km/s] (input value). |
| `'L_edd'` | `float` | Eddington luminosity [L☉]. |
| `'L_bol'` | `float` | Bolometric luminosity [L☉]. |
| `'gamma_e'` | `float` | Electron-scattering Eddington factor Γ_e = L_bol / L_edd [dimensionless]. |
| `'g_grav'` | `float` | Gravitational (Newtonian) acceleration, 10^logg [cm s⁻²]. |
| `'g_cen'` | `float` | Centrifugal acceleration, (v sin i)² / R [cm s⁻²] (0 if v sin i = 0). |
| `'g_rad'` | `float` | Radiative acceleration on free electrons [cm s⁻²]. |
| `'g_eff'` | `float` | Effective gravity, g_grav − g_rad − g_cen [cm s⁻²] (negative beyond the ΩΓ limit). |
| `'logg_eff'` | `float` | log10(g_eff / [cm s⁻²]) (NaN if g_eff ≤ 0). |
| `'eddington_limit'` | `bool` | `True` if Γ_e ≥ 1, `False` otherwise. |
| `'omega_gamma_limit'` | `bool` | `True` if g_eff ≤ 0 (ΩΓ limit reached), `False` otherwise. |

**Raises:** `TypeError` if any input is not a real number (booleans excluded); `ValueError` if any input is not finite, `teff`, `radius` or `mass` is not positive, `H_frac` is outside [0, 1] or `vsini` is negative.

---
<br>

## 4. Dependencies

| Package | Used for |
|---|---|
| `astropy` | Physical constants (IAU 2015 nominal solar values) and unit conversions. |
| `numpy` | Numerical operations and finiteness checks. |

---
<br>

## 5. Internal Structure

| Path | Content |
|---|---|
| `src/panorama/units.py` | Public functions, grouped in three sections: `WAVELENGTH`, `VELOCITY`, `STELLAR PARAMETERS`. |
| `src/panorama/support/units/units_support.py` | Validation of real-valued inputs (`_check_real`: type, finiteness and range, with the physical quantity and unit in the error message), shared by every function of the module and reused by `panorama.obs` (exposure times, signal-to-noise ratios, periods and phases) and `panorama.spec` (wavelength limits and sample regions; `allow_inf=True` accepts open limits). `eddington` and `logg_rad_centrifugal` reuse `logg_centrifugal` and `logg_radiative` directly. |

---
<br>

## 6. Demo

The Jupyter notebook [`Demos/units_demo.ipynb`](../Demos/units_demo.ipynb) runs every function of the module in 13 examples, one section per example:

* **Examples 1-4:** non-relativistic Doppler conversions and velocity widths in Å; relativistic vs. non-relativistic shifts (with the ≥ 1 % c warning and the `ValueError` for unphysical shifts); vacuum/air conversion of optical lines (round trip and out-of-range warning); and resolution limits.
* **Examples 5-6:** thermal and Doppler velocities of several ions and the corresponding Doppler widths, and Newtonian vs. radiation-corrected escape velocities.
* **Examples 7-11:** Stefan-Boltzmann and bolometric-magnitude conversions, a solar sanity check, the centrifugal correction of log g and its effect on the spectroscopic mass, the radiative and combined corrections (with their applicability caveat), and the spectroscopic luminosity.
* **Examples 12-13:** the Eddington factor of a non-rotating star (spectroscopic and independent mass), and the full workflow for a rotating star (fitted log g → `logg_centrifugal` → `eddington`, recovering the fitted gravity as `'g_eff'`), including a star beyond the ΩΓ limit (warning, `'logg_eff'` = NaN, `'omega_gamma_limit'` = `True`).

The notebook needs no input data and no internet access, and writes no files. Warnings are captured and printed, so the output is self-explanatory.

```bash
pip install -e ".[demos]"     # PANORAMA + JupyterLab and PyQt6 (demo extras)
cd Demos
jupyter lab units_demo.ipynb  # Or: jupyter notebook units_demo.ipynb
```

---
<br>

## 7. References

* Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022, *ApJ*, 935, 167 — `astropy` (physical constants and units).
* Edlén, B. 1966, *Metrologia*, 2, 71 — refractive index of standard air (`vac2air`, `air2vac`).
* Herrero, A., Kudritzki, R. P., Vílchez, J. M. et al. 1992, *A&A*, 261, 209 — centrifugal correction of log g; mass discrepancy.
* Lamers, H. J. G. L. M., Snow, T. P. & Lindholm, D. M. 1995, *ApJ*, 455, 269 — v_∞ / v_esc ratio and bi-stability jump (`vescape`).
* Langer, N. 1997, *ASP Conf. Ser.*, 120, 83 — the ΩΓ limit (`eddington`).
* Langer, N. 1998, *A&A*, 329, 551 — critical rotation of luminous stars (`eddington`).
* Langer, N. & Kudritzki, R. P. 2014, *A&A*, 564, A52 — spectroscopic HR diagram and ℒ = Teff⁴/g (`spec_luminosity`).
* Mamajek, E. E., Torres, G., Prša, A. et al. 2015, arXiv:1510.06262 — IAU 2015 Resolution B2, bolometric magnitude scale (`logL2Mbol`, `Mbol2logL`).
* Repolust, T., Puls, J. & Herrero, A. 2004, *A&A*, 415, 349 — centrifugal correction of log g (`logg_centrifugal`).
