"""
Physical conversions and derived stellar quantities for the spectroscopic
analysis of (massive) stars.

Conventions shared by every function in this module:
    - Wavelengths are in Angstroms (Å).
    - Velocities are in km/s; positive radial velocities correspond to a
      receding source (redshift).
    - Radii, masses and luminosities are in solar units (R☉, M☉, L☉),
      using the IAU 2015 nominal solar values provided by astropy.
    - Surface gravities are given as log10(g / [cm s^-2]). Functions that
      use 'logg' to derive masses, luminosities or Eddington quantities
      ('loggradius2mass', 'spec_luminosity', 'eddington') expect the
      Newtonian gravity, i.e. not reduced by rotation or radiation; use
      'logg_centrifugal', 'logg_radiative' or 'logg_rad_centrifugal' to
      correct an effective gravity first.
    - All functions work on scalar inputs and return Python floats or
      dictionaries of Python floats.
    - All functions validate their inputs before computing (with
      '_check_real', see 'support/units/units_support.py'): wrong types
      raise TypeError, and non-finite or out-of-range values raise
      ValueError.
"""

from __future__ import annotations

# ========================================================================
# ============================== WAVELENGTH ==============================
# ========================================================================

# ==============================
def wl2kmps(
            dwl : int | float,
            wl0 : int | float
    ) -> float:
    """
    Convert a wavelength shift into a radial velocity (non-relativistic Doppler formula).

    ----------
    Parameters
    ----------
        dwl : int | float
            Wavelength shift Δλ = λ_obs - λ0, in Angstroms. Positive values
            correspond to redshifts (receding source). Must satisfy
            λ_obs = λ0 + Δλ > 0.

        wl0 : int | float
            Rest (reference) wavelength λ0, in Angstroms. Must be positive.

    -------
    Returns
    -------
        dkmps : float
            Radial velocity in km/s. Positive values correspond to a
            receding source.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'wl0' is not positive, or the
            observed wavelength λ0 + Δλ is not positive.

    -----
    Notes
    -----
        - Uses the first-order (non-relativistic) Doppler formula:
                v = (Δλ / λ0) * c
          where c = 299792.458 km/s is the speed of light.

        - The relative error of the approximation grows as ~β/2 (β = v/c).
          A warning is issued when |v| >= 1% of c; in that regime use
          'wl2kmps_rel'.

        - Δλ and λ0 must refer to the same medium (both in air or both in
          vacuum; see 'air2vac' and 'vac2air').

        - Inverse: 'kmps2wl'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    import numpy as np
    import warnings


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    dwl = _check_real(dwl, "dwl", "Wavelength shift", "in Å")
    #
    wl0 = _check_real(wl0, "wl0", "Rest wavelength", "in Å", "positive")
    #
    if wl0 + dwl <= 0:  # Observed wavelength must be physical
        raise ValueError(f"Wavelength shift 'dwl' must be greater than -wl0 (the observed wavelength wl0 + dwl must be positive). Got dwl = {dwl} Å and wl0 = {wl0} Å.")


    # ####################  COMPUTATION  ####################
    # ##########  FORMULA  ##########
    dkmps = (dwl / wl0) * cte.c.value / 1000  # v = (Δλ / λ0) * c, with c converted from m/s to km/s


    # ##########  CHECKS  ##########
    # Warn when the non-relativistic approximation is no longer accurate (|v| >= 1% of c)
    if np.abs(dkmps) / (cte.c.value / 1000) >= 0.01:
        warnings.warn(
                      f"Computed radial velocity {dkmps} km/s "
                      "is >= 1% of the speed of light. Consider using "
                      "the relativistic conversion ('wl2kmps_rel').",
                      stacklevel=2  # Point the warning to the caller
                     )


    # ####################  OUTPUT  ####################
    return float(dkmps)




# ==============================
def kmps2wl(
            dkmps : int | float,
            wl0   : int | float
    ) -> float:
    """
    Convert a radial velocity into a wavelength shift (non-relativistic Doppler formula).

    ----------
    Parameters
    ----------
        dkmps : int | float
            Radial velocity in km/s. Positive values correspond to a
            receding source (redshift).

        wl0 : int | float
            Rest (reference) wavelength λ0, in Angstroms. Must be positive.

    -------
    Returns
    -------
        dwl : float
            Wavelength shift Δλ = λ_obs - λ0, in Angstroms.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, or 'wl0' is not positive.

    -----
    Notes
    -----
        - Uses the first-order (non-relativistic) Doppler formula:
                Δλ = (v / c) * λ0
          where c = 299792.458 km/s is the speed of light.

        - A warning is issued when |v| >= 1% of c; in that regime use
          'kmps2wl_rel'.

        - Useful to translate velocity widths (e.g. v sin i, macroturbulence,
          instrumental FWHM) into wavelength units at a given line.

        - Inverse: 'wl2kmps'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    import numpy as np
    import warnings


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    dkmps = _check_real(dkmps, "dkmps", "Radial velocity", "in km/s")
    #
    wl0 = _check_real(wl0, "wl0", "Rest wavelength", "in Å", "positive")


    # ##########  CHECKS  ##########
    if np.abs(dkmps) / (cte.c.value / 1000) >= 0.01:  # Non-relativistic approximation no longer accurate
        warnings.warn(
                      f"Radial velocity 'dkmps' = {dkmps} km/s "
                      "is >= 1% of the speed of light. Consider "
                      "using the relativistic conversion ('kmps2wl_rel').",
                      stacklevel=2  # Point the warning to the caller
                     )


    # ####################  COMPUTATION  ####################
    dwl = (dkmps * 1000 * wl0) / cte.c.value  # Δλ = (v / c) * λ0, with v converted from km/s to m/s


    # ####################  OUTPUT  ####################
    return float(dwl)




# ==============================
def wl2kmps_rel(
                dwl : int | float,
                wl0 : int | float
    ) -> float:
    """
    Convert a wavelength shift into a radial velocity (relativistic Doppler formula).

    ----------
    Parameters
    ----------
        dwl : int | float
            Wavelength shift Δλ = λ_obs - λ0, in Angstroms. Positive values
            correspond to redshifts (receding source). Must satisfy
            λ_obs = λ0 + Δλ > 0.

        wl0 : int | float
            Rest (reference) wavelength λ0, in Angstroms. Must be positive.

    -------
    Returns
    -------
        dkmps : float
            Radial velocity in km/s. Positive values correspond to a
            receding source.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'wl0' is not positive, the observed
            wavelength λ0 + Δλ is not positive, or the shift implies
            |β| >= 1 in floating-point arithmetic.

    -----
    Notes
    -----
        - Uses the special-relativistic Doppler formula for purely radial
          motion. With z = Δλ / λ0:
                β = ((1 + z)^2 - 1) / ((1 + z)^2 + 1)
          and the velocity is v = β * c.

        - Any positive observed wavelength gives -1 < β < 1; |β| = 1 can
          only be reached numerically for extreme shifts.

        - It is the exact inverse of 'kmps2wl_rel'. For |v| << c it reduces
          to 'wl2kmps'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    dwl = _check_real(dwl, "dwl", "Wavelength shift", "in Å")
    #
    wl0 = _check_real(wl0, "wl0", "Rest wavelength", "in Å", "positive")
    #
    if wl0 + dwl <= 0:  # Otherwise (1 + z)^2 would map unphysical shifts onto valid velocities
        raise ValueError(f"Wavelength shift 'dwl' must be greater than -wl0 (the observed wavelength wl0 + dwl must be positive). Got dwl = {dwl} Å and wl0 = {wl0} Å.")


    # ####################  COMPUTATION  ####################
    z    = dwl / wl0  # Redshift
    beta = ((1 + z)**2 - 1) / ((1 + z)**2 + 1)  # β = v / c from the relativistic Doppler relation

    if abs(beta) >= 1:  # Only reachable through floating-point rounding for extreme shifts
        raise ValueError(f"Wavelength shift 'dwl' = {dwl} Å at rest wavelength 'wl0' = {wl0} Å implies |β| = {abs(beta)} >= 1 (faster than or equal to the speed of light); the relativistic Doppler formula is undefined.")
    dkmps = beta * cte.c.value / 1000  # Convert from m/s to km/s


    # ####################  OUTPUT  ####################
    return float(dkmps)




# ==============================
def kmps2wl_rel(
                dkmps : int | float,
                wl0   : int | float
    ) -> float:
    """
    Convert a radial velocity into a wavelength shift (relativistic Doppler formula).

    ----------
    Parameters
    ----------
        dkmps : int | float
            Radial velocity in km/s. Positive values correspond to a
            receding source (redshift). Must satisfy |v| < c.

        wl0 : int | float
            Rest (reference) wavelength λ0, in Angstroms. Must be positive.

    -------
    Returns
    -------
        dwl : float
            Wavelength shift Δλ = λ_obs - λ0, in Angstroms.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'wl0' is not positive, or |v| >= c.

    -----
    Notes
    -----
        - Uses the special-relativistic Doppler formula for purely radial
          motion:
                z = sqrt((1 + β) / (1 - β)) - 1
          where β = v / c.

        - The wavelength shift is then Δλ = z * λ0.

        - It is the exact inverse of 'wl2kmps_rel'. For |v| << c it reduces
          to 'kmps2wl'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    dkmps = _check_real(dkmps, "dkmps", "Radial velocity", "in km/s")
    #
    wl0 = _check_real(wl0, "wl0", "Rest wavelength", "in Å", "positive")


    # ##########  CHECKS  ##########
    beta = (dkmps * 1000) / cte.c.value  # β = v / c, with v converted from km/s to m/s
    if abs(beta) >= 1:  # Unphysical: at or above the speed of light
        raise ValueError(f"Radial velocity 'dkmps' = {dkmps} km/s implies |β| = {abs(beta)} >= 1 (faster than or equal to the speed of light); the relativistic Doppler formula is undefined.")


    # ####################  COMPUTATION  ####################
    z   = np.sqrt((1 + beta) / (1 - beta)) - 1  # Relativistic redshift
    dwl = z * wl0  # Δλ = z * λ0


    # ####################  OUTPUT  ####################
    return float(dwl)




# ==============================
def vac2air(
            wl_vac : int | float
    ) -> float:
    """
    Convert a vacuum wavelength into an air wavelength (Edlén 1966 refractive index).

    ----------
    Parameters
    ----------
        wl_vac : int | float
            Wavelength in vacuum, in Angstroms. Must be positive.

    -------
    Returns
    -------
        wl_air : float
            Wavelength in standard air, in Angstroms.

    ------
    Raises
    ------
        TypeError
            If 'wl_vac' is not a real number.

        ValueError
            If 'wl_vac' is not finite or not positive.

    -----
    Notes
    -----
        - Uses the Edlén (1966) dispersion formula for the refractive index
          of standard air (15 °C, 101325 Pa, dry, 0.03% CO2):
                n = 1 + 1e-8 * (8342.13 + 2406030 / (130 - s^2) + 15997 / (38.9 - s^2))
          where s = 1e4 / λ_vac is the vacuum wavenumber in μm^-1.

        - The air wavelength is λ_air = λ_vac / n.

        - Accurate in the range 2000-10000 Å; a warning is issued outside it.
          Below 2000 Å wavelengths are conventionally quoted in vacuum, and
          the formula diverges near 1603 Å (s^2 = 38.9).

        - Typical use: comparing line lists or synthetic spectra (often in
          vacuum) with ground-based spectra (usually in air), or vice versa.

        - Inverse: 'air2vac'.

    ----------
    References
    ----------
        - Edlén, B. 1966, Metrologia, 2, 71.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import warnings


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    wl_vac = _check_real(wl_vac, "wl_vac", "Vacuum wavelength", "in Å", "positive")


    # ##########  CHECKS  ##########
    if wl_vac < 2000 or wl_vac > 10000:  # Validity range of the Edlén (1966) formula
        warnings.warn(
                      f"Vacuum wavelength 'wl_vac' = {wl_vac} Å is outside "
                      "the range 2000-10000 Å; the conversion may be inaccurate.",
                      stacklevel=2  # Point the warning to the caller
                     )


    # ####################  COMPUTATION  ####################
    s      = 1e4 / wl_vac  # Vacuum wavenumber in μm^-1
    n      = 1 + ((8342.13 + 2406030 / (130 - s**2) + 15997 / (38.9 - s**2)) / 1e8)  # Refractive index of air, Edlén (1966)
    wl_air = wl_vac / n


    # ####################  OUTPUT  ####################
    return float(wl_air)




# ==============================
def air2vac(
            wl_air : int | float
    ) -> float:
    """
    Convert an air wavelength into a vacuum wavelength (Edlén 1966 refractive index).

    ----------
    Parameters
    ----------
        wl_air : int | float
            Wavelength in standard air, in Angstroms. Must be positive.

    -------
    Returns
    -------
        wl_vac : float
            Wavelength in vacuum, in Angstroms.

    ------
    Raises
    ------
        TypeError
            If 'wl_air' is not a real number.

        ValueError
            If 'wl_air' is not finite or not positive.

    -----
    Notes
    -----
        - Uses the Edlén (1966) dispersion formula for the refractive index
          of standard air (15 °C, 101325 Pa, dry, 0.03% CO2):
                n = 1 + 1e-8 * (8342.13 + 2406030 / (130 - s^2) + 15997 / (38.9 - s^2))
          where s = 1e4 / λ_vac is the vacuum wavenumber in μm^-1.

        - The vacuum wavelength is λ_vac = λ_air * n(λ_vac). Since n depends
          on the unknown λ_vac, the equation is solved by fixed-point
          iteration starting from λ_vac = λ_air. Each iteration reduces the
          error by a factor ~1e-5, so five iterations make 'air2vac' the
          inverse of 'vac2air' to floating-point precision in 2000-10000 Å.

        - Accurate in the range 2000-10000 Å; a warning is issued outside it.

        - Inverse: 'vac2air'.

    ----------
    References
    ----------
        - Edlén, B. 1966, Metrologia, 2, 71.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import warnings


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    wl_air = _check_real(wl_air, "wl_air", "Air wavelength", "in Å", "positive")


    # ##########  CHECKS  ##########
    if wl_air < 2000 or wl_air > 10000:  # Validity range of the Edlén (1966) formula
        warnings.warn(
                      f"Air wavelength 'wl_air' = {wl_air} Å is outside "
                      "the range 2000-10000 Å; the conversion may be inaccurate.",
                      stacklevel=2  # Point the warning to the caller
                     )


    # ####################  COMPUTATION  ####################
    wl_vac = wl_air  # Initial guess (n - 1 ≈ 3e-4)
    for _ in range(5):  # Fixed-point iteration: s must be the vacuum wavenumber
        s      = 1e4 / wl_vac  # Vacuum wavenumber in μm^-1
        n      = 1 + ((8342.13 + 2406030 / (130 - s**2) + 15997 / (38.9 - s**2)) / 1e8)  # Refractive index of air, Edlén (1966)
        wl_vac = wl_air * n


    # ####################  OUTPUT  ####################
    return float(wl_vac)




# ======================================================================
# ============================== VELOCITY ==============================
# ======================================================================

# ==============================
def R2kmps(
           R : int | float
    ) -> float:
    """
    Convert a spectral resolving power into a velocity resolution limit.

    ----------
    Parameters
    ----------
        R : int | float
            Spectral resolving power R = λ / Δλ (dimensionless). Must be
            positive.

    -------
    Returns
    -------
        limit_kmps : float
            Velocity resolution limit in km/s, 1.32 * c / R.

    ------
    Raises
    ------
        TypeError
            If 'R' is not a real number.

        ValueError
            If 'R' is not finite or not positive.

    -----
    Notes
    -----
        - c / R is the velocity FWHM of the instrumental profile.

        - The factor 1.32 is an empirical scaling applied on top of c / R to
          give the smallest velocity that can be considered resolved
          (e.g. a lower limit for measurable line broadening such as v sin i).
          Use 'kmps2R' for the inverse relation with the same convention.

        - Example: R = 46000 (FIES/HERMES-like) gives ≈ 8.6 km/s.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    R = _check_real(R, "R", "Spectral resolving power", "dimensionless", "positive")


    # ####################  COMPUTATION  ####################
    limit_kmps = 1.32 * (cte.c.value / 1000) / R  # Empirical factor 1.32 times the velocity FWHM c / R


    # ####################  OUTPUT  ####################
    return float(limit_kmps)




# ==============================
def kmps2R(
           dkmps : int | float
    ) -> float:
    """
    Convert a velocity resolution limit into the spectral resolving power required to reach it.

    ----------
    Parameters
    ----------
        dkmps : int | float
            Velocity resolution limit in km/s. Must be positive.

    -------
    Returns
    -------
        R : float
            Spectral resolving power R = λ / Δλ (dimensionless).

    ------
    Raises
    ------
        TypeError
            If 'dkmps' is not a real number.

        ValueError
            If 'dkmps' is not finite or not positive.

    -----
    Notes
    -----
        - Exact inverse of 'R2kmps', with the same empirical factor:
                R = 1.32 * c / v

        - Useful to check which instrument can resolve a given broadening,
          e.g. kmps2R(10) ≈ 39600 for v sin i = 10 km/s.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    dkmps = _check_real(dkmps, "dkmps", "Velocity resolution limit", "in km/s", "positive")


    # ####################  COMPUTATION  ####################
    R = 1.32 * (cte.c.value / 1000) / dkmps  # Inverse of 'R2kmps'


    # ####################  OUTPUT  ####################
    return float(R)




# ==============================
def vthermal(
             teff        : int | float,
             atomic_mass : int | float,
             vmic        : int | float = 0.
    ) -> dict[str, float]:
    """
    Compute the thermal velocity of an ion and its combination with microturbulence.

    ----------
    Parameters
    ----------
        teff : int | float
            Temperature in Kelvin (K). The effective temperature is usually
            adopted as representative of the line-forming region. Must be
            positive.

        atomic_mass : int | float
            Atomic mass of the absorbing species in atomic mass units (u),
            e.g. 1.008 (H), 4.003 (He), 12.011 (C), 14.007 (N), 15.999 (O),
            24.305 (Mg), 28.086 (Si), 55.845 (Fe). Must be positive.

        vmic : int | float, optional
            Microturbulent velocity ξ in km/s. Must be non-negative.
            Default is 0.

    -------
    Returns
    -------
        results_dict : dict[str, float]
            A dictionary containing the following key-value pairs:
                - 'teff' : float
                    Temperature in K (input value).

                - 'atomic_mass' : float
                    Atomic mass in u (input value).

                - 'vmic' : float
                    Microturbulent velocity in km/s (input value).

                - 'v_th' : float
                    Thermal (most probable) velocity in km/s.

                - 'v_doppler' : float
                    Total Doppler velocity in km/s, combining thermal and
                    microturbulent motions in quadrature.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'teff' or 'atomic_mass' is not
            positive, or 'vmic' is negative.

    -----
    Notes
    -----
        - Thermal velocity (most probable speed of a Maxwellian distribution):
                v_th = sqrt(2 * k_B * T / (A * m_u))
          where A is the atomic mass in u, m_u is the atomic mass unit, and
          k_B is Boltzmann's constant.

        - Total Doppler velocity, with microturbulence added in quadrature:
                v_D = sqrt(v_th^2 + ξ^2)

        - The Doppler width of a line is Δλ_D = (v_D / c) * λ0
          (see 'kmps2wl'), and the Gaussian FWHM is 2 * sqrt(ln 2) * Δλ_D.

        - Macroturbulence is deliberately not included. Thermal and
          microturbulent motions broaden the line absorption coefficient
          (they enter the radiative transfer and change the equivalent
          width), whereas macroturbulence, like rotation, is applied as a
          convolution of the emergent flux profile and leaves the equivalent
          width unchanged. Its usual radial-tangential profile is not
          Gaussian, so it must not be added in quadrature to v_D.

        - In OB stars thermal broadening is only relevant for light
          elements (H, He): e.g. v_th(H, 30 kK) ≈ 22 km/s, while
          v_th(Si, 30 kK) ≈ 4 km/s, so metal lines are dominated by
          rotation and macroturbulence.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")
    #
    atomic_mass = _check_real(atomic_mass, "atomic_mass", "Atomic mass", "in u", "positive")
    #
    vmic = _check_real(vmic, "vmic", "Microturbulent velocity", "in km/s", "non-negative")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'teff'        : teff,
                    'atomic_mass' : atomic_mass,
                    'vmic'        : vmic,
                    'v_th'        : np.nan,
                    'v_doppler'   : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  THERMAL VELOCITY  ##########
    v_th = np.sqrt(2 * cte.k_B.value * teff / (atomic_mass * cte.u.value)) / 1000  # SI (m/s) → km/s


    # ##########  TOTAL DOPPLER VELOCITY  ##########
    v_doppler = np.sqrt(v_th**2 + vmic**2)  # Thermal and microturbulent motions added in quadrature


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'v_th'      : float(v_th),
                         'v_doppler' : float(v_doppler),
                        })
    return results_dict




# ==============================
def vescape(
            mass   : int | float,
            radius : int | float,
            gamma  : int | float = 0.
    ) -> dict[str, float]:
    """
    Compute the Newtonian and effective (radiation-corrected) photospheric escape velocities of a star.

    ----------
    Parameters
    ----------
        mass : int | float
            Stellar mass in solar masses. Must be positive.

        radius : int | float
            Stellar radius in solar radii. Must be positive.

        gamma : int | float, optional
            Electron-scattering Eddington factor Γ_e (dimensionless; see
            'eddington'). Must lie in [0, 1). Default is 0 (pure Newtonian
            escape velocity).

    -------
    Returns
    -------
        results_dict : dict[str, float]
            A dictionary containing the following key-value pairs:
                - 'mass' : float
                    Stellar mass in solar masses (input value).

                - 'radius' : float
                    Stellar radius in solar radii (input value).

                - 'gamma' : float
                    Eddington factor Γ_e (input value).

                - 'v_esc' : float
                    Newtonian escape velocity in km/s.

                - 'v_esc_eff' : float
                    Effective escape velocity in km/s, reduced by the
                    radiative acceleration on free electrons.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'mass' or 'radius' is not positive,
            or 'gamma' is not in [0, 1).

    -----
    Notes
    -----
        - Newtonian escape velocity:
                v_esc = sqrt(2 * G * M / R)

        - Effective escape velocity, accounting for the radiative
          acceleration on free electrons (relevant for hot, luminous stars):
                v_esc,eff = sqrt(2 * G * M * (1 - Γ_e) / R)

        - The terminal velocity of line-driven winds scales with v_esc,eff:
          v_∞ / v_esc,eff ≈ 2.6 for Teff ≳ 21 kK and ≈ 1.3 for cooler
          stars (bi-stability jump; Lamers et al. 1995).

    ----------
    References
    ----------
        - Lamers, H. J. G. L. M., Snow, T. P. & Lindholm, D. M. 1995,
          ApJ, 455, 269.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    mass = _check_real(mass, "mass", "Stellar mass", "in M☉", "positive")
    #
    radius = _check_real(radius, "radius", "Stellar radius", "in R☉", "positive")
    #
    gamma = _check_real(gamma, "gamma", "Eddington factor", "dimensionless", (0., 1., "[)"))


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'mass'      : mass,
                    'radius'    : radius,
                    'gamma'     : gamma,
                    'v_esc'     : np.nan,
                    'v_esc_eff' : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  NEWTONIAN  ##########
    v_esc = np.sqrt(2 * cte.G * mass * u.Msun / (radius * u.Rsun)).to(u.km / u.s)  # v_esc = sqrt(2GM / R)


    # ##########  EFFECTIVE (RADIATION-CORRECTED)  ##########
    v_esc_eff = v_esc * np.sqrt(1 - gamma)  # Gravity reduced by a factor (1 - Γ_e)


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'v_esc'     : float(v_esc.value),
                         'v_esc_eff' : float(v_esc_eff.value),
                        })
    return results_dict




# ================================================================================
# ============================== STELLAR PARAMETERS ==============================
# ================================================================================

# ==============================
def teffradius2logL(
                    teff   : int | float,
                    radius : int | float
    ) -> float:
    """
    Compute the bolometric luminosity from the effective temperature and radius (Stefan-Boltzmann law).

    ----------
    Parameters
    ----------
        teff : int | float
            Effective temperature in Kelvin (K). Must be positive.

        radius : int | float
            Stellar radius in solar radii. Must be positive.

    -------
    Returns
    -------
        logL : float
            Bolometric luminosity, log10(L / L☉).

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite or not positive.

    -----
    Notes
    -----
        - Stefan-Boltzmann law:
                L = 4 * π * R^2 * σ * Teff^4

        - L☉ is the IAU 2015 nominal solar luminosity (3.828e26 W), which
          together with the nominal R☉ implies Teff☉ ≈ 5772 K.

        - Inverse: 'logLteff2radius'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")
    #
    radius = _check_real(radius, "radius", "Stellar radius", "in R☉", "positive")


    # ####################  COMPUTATION  ####################
    lum  = (4 * np.pi * (radius * u.Rsun)**2 * cte.sigma_sb * (teff * u.K)**4).to(u.Lsun)  # Stefan-Boltzmann law
    logL = np.log10(lum.value)


    # ####################  OUTPUT  ####################
    return float(logL)




# ==============================
def logLteff2radius(
                    logL : int | float,
                    teff : int | float
    ) -> float:
    """
    Compute the stellar radius from the bolometric luminosity and effective temperature (Stefan-Boltzmann law).

    ----------
    Parameters
    ----------
        logL : int | float
            Bolometric luminosity, log10(L / L☉).

        teff : int | float
            Effective temperature in Kelvin (K). Must be positive.

    -------
    Returns
    -------
        radius : float
            Stellar radius in solar radii.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, or 'teff' is not positive.

    -----
    Notes
    -----
        - Stefan-Boltzmann law solved for the radius:
                R = sqrt(L / (4 * π * σ * Teff^4))

        - Typical use: deriving R from a photometric luminosity (e.g. from
          Gaia distances + bolometric corrections) and a spectroscopic Teff.

        - Inverse: 'teffradius2logL'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    logL = _check_real(logL, "logL", "Bolometric luminosity", "log10(L / L☉)")
    #
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")


    # ####################  COMPUTATION  ####################
    lum    = (10**logL) * u.Lsun
    radius = np.sqrt(lum / (4 * np.pi * cte.sigma_sb * (teff * u.K)**4)).to(u.Rsun)  # Stefan-Boltzmann law solved for R


    # ####################  OUTPUT  ####################
    return float(radius.value)




# ==============================
def loggradius2mass(
                    logg   : int | float,
                    radius : int | float
    ) -> float:
    """
    Compute the spectroscopic mass from the surface gravity and the radius.

    ----------
    Parameters
    ----------
        logg : int | float
            Newtonian surface gravity, log10(g / [cm s^-2]): the gravity
            due to the stellar mass alone, not reduced by rotation or
            radiation. Correct a fitted gravity first with
            'logg_centrifugal' (rotation) and, only if it is an effective
            gravity that does not include radiation pressure, with
            'logg_radiative' or 'logg_rad_centrifugal' (radiation).

        radius : int | float
            Stellar radius in solar radii. Must be positive.

    -------
    Returns
    -------
        mass : float
            Spectroscopic mass in solar masses.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, or 'radius' is not positive.

    -----
    Notes
    -----
        - Newtonian surface gravity solved for the mass:
                M = g * R^2 / G

        - 'logg' must be the Newtonian gravity. A gravity reduced by the
          centrifugal acceleration (rotating stars, see 'logg_centrifugal')
          or by radiation pressure (if not included in the model atmosphere,
          see 'logg_radiative') underestimates the mass.

        - Spectroscopic masses are known to be systematically lower than
          evolutionary masses for part of the OB-star population (the
          "mass discrepancy"; Herrero et al. 1992).

    ----------
    References
    ----------
        - Herrero, A., Kudritzki, R. P., Vílchez, J. M. et al. 1992,
          A&A, 261, 209.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    logg = _check_real(logg, "logg", "Surface gravity", "log10(g / [cm s^-2])")
    #
    radius = _check_real(radius, "radius", "Stellar radius", "in R☉", "positive")


    # ####################  COMPUTATION  ####################
    g    = (10**logg) * (u.cm / u.s**2)
    mass = (g * (radius * u.Rsun)**2 / cte.G).to(u.Msun)  # M = g * R^2 / G


    # ####################  OUTPUT  ####################
    return float(mass.value)




# ==============================
def logg_centrifugal(
                     logg   : int | float,
                     vsini  : int | float,
                     radius : int | float
    ) -> dict[str, float]:
    """
    Correct a spectroscopic surface gravity for the centrifugal acceleration due to rotation.

    ----------
    Parameters
    ----------
        logg : int | float
            Spectroscopic (effective) surface gravity, log10(g / [cm s^-2]).

        vsini : int | float
            Projected rotational velocity in km/s. Must be non-negative.

        radius : int | float
            Stellar radius in solar radii. Must be positive.

    -------
    Returns
    -------
        results_dict : dict[str, float]
            A dictionary containing the following key-value pairs:
                - 'logg' : float
                    Input (effective) surface gravity, log10(g / [cm s^-2]).

                - 'vsini' : float
                    Projected rotational velocity in km/s (input value).

                - 'radius' : float
                    Stellar radius in solar radii (input value).

                - 'g_eff' : float
                    Effective surface gravity, 10**logg, in cm/s^2.

                - 'g_cen' : float
                    Centrifugal acceleration (v sin i)^2 / R in cm/s^2.

                - 'logg_cen' : float
                    Logarithm of the centrifugal acceleration,
                    log10(g_cen / [cm s^-2]) (-inf if vsini = 0).

                - 'g_grav' : float
                    Gravitational acceleration, g_eff + g_cen, in cm/s^2.

                - 'logg_grav' : float
                    Centrifugally corrected surface gravity,
                    log10(g_grav / [cm s^-2]).

                - 'dlogg' : float
                    Correction applied, logg_grav - logg, in dex
                    (non-negative).

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'vsini' is negative, 'radius' is
            not positive, or the computed accelerations are not finite.

    -----
    Notes
    -----
        - Plane-parallel model atmospheres fitted to a rotating star yield an
          effective gravity reduced by the centrifugal acceleration. A
          statistical first-order correction is:
                g_grav = g_eff + (v sin i)^2 / R
          (Herrero et al. 1992; Repolust et al. 2004).

        - Since v sin i <= v_eq, the correction is a lower limit to the
          equatorial centrifugal acceleration. It is ~0.01 dex for slow
          rotators, but can reach ~0.1 dex for v sin i ≳ 250 km/s in
          giants/supergiants.

        - Use 'logg_grav' to compute spectroscopic masses
          (see 'loggradius2mass') and as input of 'eddington'.

        - Only the centrifugal term is corrected. See 'logg_radiative' and
          'logg_rad_centrifugal' when the gravity is also reduced by
          radiation pressure.

    ----------
    References
    ----------
        - Herrero, A., Kudritzki, R. P., Vílchez, J. M. et al. 1992,
          A&A, 261, 209.
        - Repolust, T., Puls, J. & Herrero, A. 2004, A&A, 415, 349.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    logg = _check_real(logg, "logg", "Surface gravity", "log10(g / [cm s^-2])")
    #
    vsini = _check_real(vsini, "vsini", "Projected rotational velocity", "in km/s", "non-negative")
    #
    radius = _check_real(radius, "radius", "Stellar radius", "in R☉", "positive")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'logg'      : logg,
                    'vsini'     : vsini,
                    'radius'    : radius,
                    'g_eff'     : np.nan,
                    'g_cen'     : np.nan,
                    'logg_cen'  : np.nan,
                    'g_grav'    : np.nan,
                    'logg_grav' : np.nan,
                    'dlogg'     : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  ACCELERATIONS  ##########
    # ##### COMPUTATION #####
    with np.errstate(over="ignore", invalid="ignore"):  # Overflow is reported below as a ValueError
        g_eff = float(np.power(10., logg))  # Effective (fitted) gravity in cm/s^2
        g_cen = float(((vsini * u.km / u.s)**2 / (radius * u.Rsun)).to(u.cm / u.s**2).value)  # (v sin i)^2 / R


    # ##### CHECKS #####
    if not np.isfinite(g_eff):
        raise ValueError(f"Effective surface gravity 'logg' = {logg} gives a non-finite effective gravity 10**logg; it must be finite.")
    elif not np.isfinite(g_cen):
        raise ValueError(f"Projected rotational velocity 'vsini' = {vsini} km/s and stellar radius 'radius' = {radius} R☉ give a non-finite centrifugal acceleration; it must be finite.")


    # ##########  CORRECTED GRAVITY  ##########
    logg_cen  = np.log10(g_cen) if g_cen > 0 else -np.inf  # log10(0) = -inf without a numpy warning
    g_grav    = g_eff + g_cen  # Accelerations add linearly, not in log space
    logg_grav = np.log10(g_grav)  # Newtonian gravity corrected for the centrifugal support


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'g_eff'     : g_eff,
                         'g_cen'     : g_cen,
                         'logg_cen'  : float(logg_cen),
                         'g_grav'    : float(g_grav),
                         'logg_grav' : float(logg_grav),
                         'dlogg'     : float(logg_grav - logg),
                        })
    return results_dict




# ==============================
def logg_radiative(
                   logg   : int | float,
                   teff   : int | float,
                   H_frac : int | float = 0.7
    ) -> dict[str, float]:
    """
    Correct an effective surface gravity for the radiative acceleration on free electrons.

    ----------
    Parameters
    ----------
        logg : int | float
            Effective surface gravity (reduced by radiation pressure),
            log10(g / [cm s^-2]).

        teff : int | float
            Effective temperature in Kelvin (K). Must be positive.

        H_frac : int | float, optional
            Hydrogen mass fraction X, in [0, 1]. Default is 0.7
            (approximately solar).

    -------
    Returns
    -------
        results_dict : dict[str, float]
            A dictionary containing the following key-value pairs:
                - 'logg' : float
                    Input (effective) surface gravity, log10(g / [cm s^-2]).

                - 'teff' : float
                    Effective temperature in K (input value).

                - 'H_frac' : float
                    Hydrogen mass fraction X (input value).

                - 'g_eff' : float
                    Effective surface gravity, 10**logg, in cm/s^2.

                - 'g_rad' : float
                    Radiative acceleration on free electrons
                    κ_e * σ * Teff^4 / c in cm/s^2.

                - 'logg_rad' : float
                    Logarithm of the radiative acceleration,
                    log10(g_rad / [cm s^-2]).

                - 'g_grav' : float
                    Gravitational acceleration, g_eff + g_rad, in cm/s^2.

                - 'logg_grav' : float
                    Radiatively corrected surface gravity,
                    log10(g_grav / [cm s^-2]).

                - 'dlogg' : float
                    Correction applied, logg_grav - logg, in dex
                    (non-negative).

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'teff' is not positive, 'H_frac' is
            outside [0, 1], or the computed accelerations are not finite.

    -----
    Notes
    -----
        1. CORRECTION
            - The radiative acceleration on free electrons (Thomson
              scattering of a fully ionized H/He plasma) is added to the
              effective gravity:
                    g_grav = g_eff + g_rad,   g_rad = κ_e * σ * Teff^4 / c
              with κ_e = 0.20 * (1 + X) cm^2/g. It depends only on Teff and
              X, not on the mass or the radius.

            - g_rad is a lower limit to the total radiative acceleration,
              which also includes line and bound-free opacity.

        2. APPLICABILITY
            - Apply it only when 'logg' is an effective gravity that does not
              include the radiative support (e.g. g_eff = g_grav * (1 - Γ_e)
              from a model or relation parameterized by the effective
              gravity).

            - Model atmospheres that include radiation pressure in the
              hydrostatic equilibrium (e.g. FASTWIND, CMFGEN, TLUSTY, ATLAS)
              already return the gravitational acceleration with respect to
              radiation: correcting their 'logg' again would count g_rad
              twice. For those, only the centrifugal correction
              ('logg_centrifugal') applies.

            - The correction is large for hot, luminous stars:
              g_rad ≈ 520 cm/s^2 at Teff = 30 kK (≈ +0.18 dex at
              log g = 3.0).

            - See 'logg_rad_centrifugal' for the combined radiative and
              centrifugal correction, and 'eddington' for Γ_e.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    logg = _check_real(logg, "logg", "Surface gravity", "log10(g / [cm s^-2])")
    #
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")
    #
    H_frac = _check_real(H_frac, "H_frac", "Hydrogen mass fraction", "dimensionless", (0., 1., "[]"))


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'logg'      : logg,
                    'teff'      : teff,
                    'H_frac'    : H_frac,
                    'g_eff'     : np.nan,
                    'g_rad'     : np.nan,
                    'logg_rad'  : np.nan,
                    'g_grav'    : np.nan,
                    'logg_grav' : np.nan,
                    'dlogg'     : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  ACCELERATIONS  ##########
    # ##### COMPUTATION #####
    kappa = 0.20 * (1 + H_frac) * (u.cm**2 / u.g)  # Electron-scattering opacity: κ_e = 0.20 * (1 + X) cm^2/g
    with np.errstate(over="ignore", invalid="ignore"):  # Overflow is reported below as a ValueError
        g_eff = float(np.power(10., logg))  # Effective gravity in cm/s^2
        g_rad = float((kappa * cte.sigma_sb * (teff * u.K)**4 / cte.c).to(u.cm / u.s**2).value)  # Same expression as in 'eddington'


    # ##### CHECKS #####
    if not np.isfinite(g_eff):
        raise ValueError(f"Effective surface gravity 'logg' = {logg} gives a non-finite effective gravity 10**logg; it must be finite.")
    elif not np.isfinite(g_rad):
        raise ValueError(f"Effective temperature 'teff' = {teff} K gives a non-finite radiative acceleration; it must be finite.")


    # ##########  CORRECTED GRAVITY  ##########
    logg_rad  = np.log10(g_rad)  # g_rad > 0 since teff > 0
    g_grav    = g_eff + g_rad  # Accelerations add linearly, not in log space
    logg_grav = np.log10(g_grav)  # Newtonian gravity corrected for the radiative support


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'g_eff'     : g_eff,
                         'g_rad'     : g_rad,
                         'logg_rad'  : float(logg_rad),
                         'g_grav'    : float(g_grav),
                         'logg_grav' : float(logg_grav),
                         'dlogg'     : float(logg_grav - logg),
                        })
    return results_dict




# ==============================
def logg_rad_centrifugal(
                         logg   : int | float,
                         teff   : int | float,
                         vsini  : int | float,
                         radius : int | float,
                         H_frac : int | float = 0.7
    ) -> dict[str, float]:
    """
    Correct an effective surface gravity for both the radiative and the centrifugal accelerations.

    ----------
    Parameters
    ----------
        logg : int | float
            Effective surface gravity (reduced by radiation pressure and
            rotation), log10(g / [cm s^-2]).

        teff : int | float
            Effective temperature in Kelvin (K). Must be positive.

        vsini : int | float
            Projected rotational velocity in km/s. Must be non-negative.

        radius : int | float
            Stellar radius in solar radii. Must be positive.

        H_frac : int | float, optional
            Hydrogen mass fraction X, in [0, 1]. Default is 0.7
            (approximately solar).

    -------
    Returns
    -------
        results_dict : dict[str, float]
            A dictionary containing the following key-value pairs:
                - 'logg' : float
                    Input (effective) surface gravity, log10(g / [cm s^-2]).

                - 'teff' : float
                    Effective temperature in K (input value).

                - 'vsini' : float
                    Projected rotational velocity in km/s (input value).

                - 'radius' : float
                    Stellar radius in solar radii (input value).

                - 'H_frac' : float
                    Hydrogen mass fraction X (input value).

                - 'g_eff' : float
                    Effective surface gravity, 10**logg, in cm/s^2.

                - 'g_rad' : float
                    Radiative acceleration on free electrons
                    κ_e * σ * Teff^4 / c in cm/s^2.

                - 'logg_rad' : float
                    Logarithm of the radiative acceleration,
                    log10(g_rad / [cm s^-2]).

                - 'g_cen' : float
                    Centrifugal acceleration (v sin i)^2 / R in cm/s^2.

                - 'logg_cen' : float
                    Logarithm of the centrifugal acceleration,
                    log10(g_cen / [cm s^-2]) (-inf if vsini = 0).

                - 'g_grav' : float
                    Gravitational acceleration, g_eff + g_rad + g_cen, in
                    cm/s^2.

                - 'logg_grav' : float
                    Radiatively and centrifugally corrected surface gravity,
                    log10(g_grav / [cm s^-2]).

                - 'dlogg' : float
                    Total correction applied, logg_grav - logg, in dex
                    (non-negative).

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, 'teff' or 'radius' is not positive,
            'vsini' is negative, 'H_frac' is outside [0, 1], or the computed
            accelerations are not finite.

    -----
    Notes
    -----
        - Both accelerations add linearly to the effective gravity:
                g_grav = g_eff + g_rad + (v sin i)^2 / R
          with g_rad from 'logg_radiative' and g_cen from 'logg_centrifugal'
          (both functions are reused, so the three corrections are always
          consistent).

        - The individual corrections in dex are not additive; use 'g_rad'
          and 'g_cen' (or 'logg_rad' and 'logg_cen') to compare them.

        - Since v sin i <= v_eq and g_rad only includes electron scattering,
          both terms are lower limits.

        - Apply the radiative term only when 'logg' does not already include
          radiation pressure (see 'logg_radiative', Notes 2). For gravities
          fitted with FASTWIND, CMFGEN, TLUSTY or ATLAS, use
          'logg_centrifugal' alone.

        - Passing 'logg_grav' to 'eddington' with the same 'vsini' gives
          g_eff = 10**logg in 'eddington', i.e. the input effective gravity
          is recovered.

    ----------
    References
    ----------
        - Herrero, A., Kudritzki, R. P., Vílchez, J. M. et al. 1992,
          A&A, 261, 209.
        - Repolust, T., Puls, J. & Herrero, A. 2004, A&A, 415, 349.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    logg = _check_real(logg, "logg", "Surface gravity", "log10(g / [cm s^-2])")
    #
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")
    #
    vsini = _check_real(vsini, "vsini", "Projected rotational velocity", "in km/s", "non-negative")
    #
    radius = _check_real(radius, "radius", "Stellar radius", "in R☉", "positive")
    #
    H_frac = _check_real(H_frac, "H_frac", "Hydrogen mass fraction", "dimensionless", (0., 1., "[]"))


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'logg'      : logg,
                    'teff'      : teff,
                    'vsini'     : vsini,
                    'radius'    : radius,
                    'H_frac'    : H_frac,
                    'g_eff'     : np.nan,
                    'g_rad'     : np.nan,
                    'logg_rad'  : np.nan,
                    'g_cen'     : np.nan,
                    'logg_cen'  : np.nan,
                    'g_grav'    : np.nan,
                    'logg_grav' : np.nan,
                    'dlogg'     : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  ACCELERATIONS  ##########
    radiative   = logg_radiative(logg, teff, H_frac)  # g_eff and g_rad (also checks that they are finite)
    centrifugal = logg_centrifugal(logg, vsini, radius)  # g_cen (also checks that it is finite)


    # ##########  CORRECTED GRAVITY  ##########
    g_grav    = radiative['g_eff'] + radiative['g_rad'] + centrifugal['g_cen']  # Accelerations add linearly, not in log space
    logg_grav = np.log10(g_grav)  # Newtonian gravity corrected for both supports


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'g_eff'     : radiative['g_eff'],
                         'g_rad'     : radiative['g_rad'],
                         'logg_rad'  : radiative['logg_rad'],
                         'g_cen'     : centrifugal['g_cen'],
                         'logg_cen'  : centrifugal['logg_cen'],
                         'g_grav'    : float(g_grav),
                         'logg_grav' : float(logg_grav),
                         'dlogg'     : float(logg_grav - logg),
                        })
    return results_dict




# ==============================
def spec_luminosity(
                    teff : int | float,
                    logg : int | float
    ) -> float:
    """
    Compute the spectroscopic luminosity ℒ = Teff^4 / g in solar units.

    ----------
    Parameters
    ----------
        teff : int | float
            Effective temperature in Kelvin (K). Must be positive.

        logg : int | float
            Newtonian surface gravity, log10(g / [cm s^-2]): the gravity
            due to the stellar mass alone, not reduced by rotation or
            radiation. Correct a fitted gravity first with
            'logg_centrifugal' (rotation) and, only if it is an effective
            gravity that does not include radiation pressure, with
            'logg_radiative' or 'logg_rad_centrifugal' (radiation).

    -------
    Returns
    -------
        logLspec : float
            Spectroscopic luminosity, log10(ℒ / ℒ☉).

    ------
    Raises
    ------
        TypeError
            If any input is not a real number.

        ValueError
            If any input is not finite, or 'teff' is not positive.

    -----
    Notes
    -----
        - Definition (Langer & Kudritzki 2014):
                ℒ = Teff^4 / g  ∝  L / M
          so it depends only on spectroscopic quantities (distance-independent).

        - In solar units:
                log(ℒ / ℒ☉) = 4 * log(Teff) - log(g) - [4 * log(Teff☉) - log(g☉)]
          with Teff☉ and g☉ derived from the IAU 2015 nominal solar values
          (Teff☉ ≈ 5772 K, log g☉ ≈ 4.438), i.e. log ℒ☉ ≈ 10.61
          (ℒ in K^4 cm^-1 s^2).

        - It is directly related to the electron-scattering Eddington factor:
                Γ_e = κ_e * σ * Teff^4 / (c * g)  ∝  ℒ
          and is the vertical axis of the spectroscopic HR diagram (sHRD).

    ----------
    References
    ----------
        - Langer, N. & Kudritzki, R. P. 2014, A&A, 564, A52.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")
    #
    logg = _check_real(logg, "logg", "Surface gravity", "log10(g / [cm s^-2])")


    # ####################  COMPUTATION  ####################
    # ##########  SOLAR REFERENCE  ##########
    teff_sun     = ((cte.L_sun / (4 * np.pi * cte.R_sun**2 * cte.sigma_sb))**0.25).to(u.K).value  # ≈ 5772 K
    logg_sun     = np.log10((cte.G * cte.M_sun / cte.R_sun**2).to(u.cm / u.s**2).value)  # ≈ 4.438
    logLspec_sun = 4 * np.log10(teff_sun) - logg_sun  # ≈ 10.61


    # ##########  SPECTROSCOPIC LUMINOSITY  ##########
    logLspec = 4 * np.log10(teff) - logg - logLspec_sun  # log(L/L☉) with L = Teff^4 / g, in solar units


    # ####################  OUTPUT  ####################
    return float(logLspec)




# ==============================
def logL2Mbol(
              logL : int | float
    ) -> float:
    """
    Convert a bolometric luminosity into an absolute bolometric magnitude.

    ----------
    Parameters
    ----------
        logL : int | float
            Bolometric luminosity, log10(L / L☉).

    -------
    Returns
    -------
        Mbol : float
            Absolute bolometric magnitude (mag).

    ------
    Raises
    ------
        TypeError
            If 'logL' is not a real number.

        ValueError
            If 'logL' is not finite.

    -----
    Notes
    -----
        - Uses the IAU 2015 Resolution B2 zero point, Mbol☉ = 4.74 mag:
                Mbol = 4.74 - 2.5 * log(L / L☉)

        - Combined with a bolometric correction, Mbol = M_V + BC_V.

        - Inverse: 'Mbol2logL'.

    ----------
    References
    ----------
        - Mamajek, E. E., Torres, G., Prša, A. et al. 2015, arXiv:1510.06262.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    logL = _check_real(logL, "logL", "Bolometric luminosity", "log10(L / L☉)")


    # ####################  COMPUTATION  ####################
    Mbol = 4.74 - 2.5 * logL  # IAU 2015 B2 zero point (Mbol☉ = 4.74 mag)


    # ####################  OUTPUT  ####################
    return float(Mbol)




# ==============================
def Mbol2logL(
              Mbol : int | float
    ) -> float:
    """
    Convert an absolute bolometric magnitude into a bolometric luminosity.

    ----------
    Parameters
    ----------
        Mbol : int | float
            Absolute bolometric magnitude (mag).

    -------
    Returns
    -------
        logL : float
            Bolometric luminosity, log10(L / L☉).

    ------
    Raises
    ------
        TypeError
            If 'Mbol' is not a real number.

        ValueError
            If 'Mbol' is not finite.

    -----
    Notes
    -----
        - Uses the IAU 2015 Resolution B2 zero point, Mbol☉ = 4.74 mag:
                log(L / L☉) = (4.74 - Mbol) / 2.5

        - Inverse: 'logL2Mbol'.

    ----------
    References
    ----------
        - Mamajek, E. E., Torres, G., Prša, A. et al. 2015, arXiv:1510.06262.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    Mbol = _check_real(Mbol, "Mbol", "Absolute bolometric magnitude", "in mag")


    # ####################  COMPUTATION  ####################
    logL = (4.74 - Mbol) / 2.5  # IAU 2015 B2 zero point (Mbol☉ = 4.74 mag)


    # ####################  OUTPUT  ####################
    return float(logL)




# ==============================
def eddington(
              teff   : int | float,
              logg   : int | float,
              radius : int | float,
              mass   : int | float | None = None,
              H_frac : int | float = 0.7,
              vsini  : int | float = 0.
    ) -> dict[str, float | bool]:
    """
    Compute the Eddington luminosity, Eddington factor, radiative and centrifugal accelerations and effective gravity of a (rotating) star.

    ----------
    Parameters
    ----------
        teff : int | float
            Effective temperature in Kelvin (K). Must be positive.

        logg : int | float
            Newtonian surface gravity, log10(g / [cm s^-2]): the gravity
            due to the stellar mass alone, not reduced by rotation or
            radiation. Correct a fitted gravity first with
            'logg_centrifugal' (rotation) and, only if it is an effective
            gravity that does not include radiation pressure, with
            'logg_radiative' or 'logg_rad_centrifugal' (radiation).

        radius : int | float
            Stellar radius in solar radii. Must be positive.

        mass : int | float | None, optional
            Stellar mass in solar masses. Must be positive if given. If None,
            the spectroscopic mass is derived from 'logg' and 'radius'.
            Default is None.

        H_frac : int | float, optional
            Hydrogen mass fraction X, in [0, 1]. Default is 0.7
            (approximately solar).

        vsini : int | float, optional
            Projected rotational velocity in km/s. Must be non-negative.
            Default is 0 (no centrifugal acceleration).

    -------
    Returns
    -------
        results_dict : dict[str, float | bool]
            A dictionary containing the following key-value pairs:
                - 'teff' : float
                    Effective temperature in K (input value).

                - 'logg' : float
                    Newtonian surface gravity, log10(g / [cm s^-2])
                    (input value).

                - 'radius' : float
                    Stellar radius in solar radii (input value).

                - 'mass' : float
                    Stellar mass in solar masses (input value, or the
                    spectroscopic mass from 'logg' if 'mass' was None).

                - 'H_frac' : float
                    Hydrogen mass fraction X (input value).

                - 'vsini' : float
                    Projected rotational velocity in km/s (input value).

                - 'L_edd' : float
                    Eddington luminosity in solar luminosities.

                - 'L_bol' : float
                    Bolometric luminosity in solar luminosities.

                - 'gamma_e' : float
                    Electron-scattering Eddington factor Γ_e = L_bol / L_edd
                    (dimensionless).

                - 'g_grav' : float
                    Gravitational (Newtonian) acceleration, 10**logg, in
                    cm/s^2.

                - 'g_cen' : float
                    Centrifugal acceleration (v sin i)^2 / R in cm/s^2
                    (0 if vsini = 0).

                - 'g_rad' : float
                    Radiative acceleration on free electrons in cm/s^2.

                - 'g_eff' : float
                    Effective gravity, including radiation and rotation,
                    g_grav - g_rad - g_cen, in cm/s^2 (negative above the
                    ΩΓ limit; equal to g_grav - g_rad if vsini = 0).

                - 'logg_eff' : float
                    Effective surface gravity, log10(g_eff / [cm s^-2])
                    (NaN if g_eff <= 0).

                - 'eddington_limit' : bool
                    True if the Eddington limit is reached (Γ_e >= 1),
                    False otherwise.

                - 'omega_gamma_limit' : bool
                    True if the ΩΓ limit is reached (g_eff <= 0), False
                    otherwise.

    ------
    Raises
    ------
        TypeError
            If any input is not a real number ('mass' may also be None).

        ValueError
            If any input is not finite, 'teff', 'radius' or 'mass' is not
            positive, 'H_frac' is outside [0, 1], or 'vsini' is negative.

    -----
    Notes
    -----
        1. OPACITY
            - The opacity is assumed to be dominated by electron (Thomson)
              scattering of a fully ionized H/He plasma:
                    κ_e = 0.20 * (1 + X) cm^2/g
              (≈ 0.34 cm^2/g for X = 0.7).

            - This is a good approximation in the atmospheres of hot, massive
              stars (OB, WR, LBV). Γ_e is therefore a lower limit to the true
              Eddington factor, which also includes line and bound-free
              opacity.

        2. MASS
            - If no mass is provided, the spectroscopic mass is used:
                    M = g_grav * R^2 / G
              (see 'loggradius2mass').

        3. LUMINOSITIES
            - Eddington luminosity:
                    L_edd = 4 * π * G * M * c / κ_e

            - Bolometric luminosity (Stefan-Boltzmann law):
                    L_bol = 4 * π * R^2 * σ * Teff^4

        4. ACCELERATIONS
            - Newtonian gravity:
                    g_grav = 10^logg  (in cm/s^2)

            - Radiative acceleration on free electrons:
                    g_rad = κ_e * σ * Teff^4 / c

            - Centrifugal acceleration, with the same expression as
              'logg_centrifugal':
                    g_cen = (v sin i)^2 / R

            - Effective gravity, including radiation and rotation:
                    g_eff = g_grav - g_rad - g_cen

            - 'logg' must be the Newtonian gravity. The gravity fitted with
              plane-parallel model atmospheres to a rotating star is already
              reduced by the centrifugal acceleration; passing it here
              would subtract g_cen twice. Correct it first with
              'logg_centrifugal' and pass its 'logg_grav'. Likewise, an
              effective gravity that does not include radiation pressure
              must be corrected with 'logg_rad_centrifugal' (or
              'logg_radiative', passing their 'logg_grav'); otherwise g_rad
              is also subtracted twice.

            - Since v sin i <= veq, g_cen is a lower limit to the equatorial
              centrifugal acceleration, and g_eff an upper limit to the
              equatorial effective gravity.

        5. LIMITS
            - The Eddington limit is reached when Γ_e >= 1 ('eddington_limit'
              is True).

            - The ΩΓ limit is reached when g_eff <= 0 ('omega_gamma_limit'
              is True). A warning is then issued and 'logg_eff' is NaN.

            - When the mass is derived from 'logg', g_eff <= 0 is
              equivalent to the ΩΓ limit (Langer 1997, 1998):
                    Γ_e + (v sin i)^2 / v_kep^2 >= 1
              with v_kep = sqrt(G * M / R).

            - If 'mass' is provided and is inconsistent with 'logg' and
              'radius', 'gamma_e' (based on the mass) and 'g_eff' (based on
              'logg') will not satisfy g_rad / g_grav = Γ_e exactly.

            - This is a first-order, spherical approximation: 'radius' is
              used as the equatorial radius and the radiative flux is assumed
              uniform. Roche distortion and gravity darkening are not
              included (they are left for the planned 'rot' module).

    ----------
    References
    ----------
        - Herrero, A., Kudritzki, R. P., Vílchez, J. M. et al. 1992,
          A&A, 261, 209.
        - Langer, N. 1997, ASP Conf. Ser., 120, 83.
        - Langer, N. 1998, A&A, 329, 551.
        - Repolust, T., Puls, J. & Herrero, A. 2004, A&A, 415, 349.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    from astropy import units as u
    import numpy as np
    import warnings


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    teff = _check_real(teff, "teff", "Effective temperature", "in K", "positive")
    #
    logg = _check_real(logg, "logg", "Surface gravity", "log10(g / [cm s^-2])")
    #
    radius = _check_real(radius, "radius", "Stellar radius", "in R☉", "positive")
    #
    mass = _check_real(mass, "mass", "Stellar mass", "in M☉", "positive", allow_none=True)
    #
    H_frac = _check_real(H_frac, "H_frac", "Hydrogen mass fraction", "dimensionless", (0., 1., "[]"))
    #
    vsini = _check_real(vsini, "vsini", "Projected rotational velocity", "in km/s", "non-negative")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'teff'              : teff,
                    'logg'              : logg,
                    'radius'            : radius,
                    'mass'              : np.nan,
                    'H_frac'            : H_frac,
                    'vsini'             : vsini,
                    'L_edd'             : np.nan,
                    'L_bol'             : np.nan,
                    'gamma_e'           : np.nan,
                    'g_grav'            : np.nan,
                    'g_cen'             : np.nan,
                    'g_rad'             : np.nan,
                    'g_eff'             : np.nan,
                    'logg_eff'          : np.nan,
                    'eddington_limit'   : None,
                    'omega_gamma_limit' : None,
                   }


    # ####################  PARAMETERS  ####################
    # ##########  CENTRIFUGAL ACCELERATION  ##########
    g_cen = logg_centrifugal(logg, vsini, radius)['g_cen'] * (u.cm / u.s**2)  # (v sin i)^2 / R, same expression as 'logg_centrifugal'


    # ##########  UNITS  ##########
    teff   = teff * u.K  # Quantities from here on, so that astropy handles the unit conversions
    radius = radius * u.Rsun
    g_grav = (10**logg) * (u.cm / u.s**2)  # Newtonian surface gravity
    kappa  = 0.20 * (1 + H_frac) * (u.cm**2 / u.g)  # Electron-scattering opacity: κ_e = 0.20 * (1 + X) cm^2/g


    # ##########  MASS  ##########
    if mass is None:
        mass = (g_grav * radius**2 / cte.G).to(u.Msun)  # Spectroscopic mass: M = g * R^2 / G
    else:
        mass = mass * u.Msun  # Given mass (e.g. evolutionary)


    # ####################  COMPUTATION  ####################
    # ##########  EDDINGTON LUMINOSITY  ##########
    l_edd = ((4 * np.pi * cte.G * mass * cte.c) / kappa).to(u.Lsun)  # L_edd = 4πGMc / κ_e


    # ##########  BOLOMETRIC LUMINOSITY  ##########
    l_bol = (4 * np.pi * radius**2 * cte.sigma_sb * teff**4).to(u.Lsun)  # Stefan-Boltzmann law


    # ##########  EDDINGTON FACTOR (GAMMA)  ##########
    gamma_e = (l_bol / l_edd).decompose()  # Γ_e = L_bol / L_edd (dimensionless)


    # ##########  GRAVITIES  ##########
    g_rad = (kappa * cte.sigma_sb * teff**4 / cte.c).to(u.cm / u.s**2)  # Radiative acceleration on free electrons
    g_eff = g_grav - g_rad - g_cen  # Gravity minus radiative and centrifugal support (negative above the ΩΓ limit)

    if g_eff.value > 0:
        logg_eff = np.log10(g_eff.value)  # g_eff in cm/s^2
    else:  # Logarithm undefined: the star is beyond the ΩΓ limit
        warnings.warn(
                      f"Effective gravity g_eff = {g_eff.value:.3e} cm/s^2 is non-positive. "
                      "The star is beyond the ΩΓ limit; 'logg_eff' is undefined and set to NaN.",
                      stacklevel=2  # Point the warning to the caller
                     )
        logg_eff = np.nan


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'mass'              : float(mass.value),
                         'L_edd'             : float(l_edd.value),
                         'L_bol'             : float(l_bol.value),
                         'gamma_e'           : float(gamma_e.value),
                         'g_grav'            : float(g_grav.value),
                         'g_cen'             : float(g_cen.value),
                         'g_rad'             : float(g_rad.value),
                         'g_eff'             : float(g_eff.value),
                         'logg_eff'          : float(logg_eff),
                         'eddington_limit'   : bool(gamma_e.value >= 1),  # Python bool, not numpy.bool_
                         'omega_gamma_limit' : bool(g_eff.value <= 0),
                        })
    return results_dict
