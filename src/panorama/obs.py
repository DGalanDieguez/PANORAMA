"""
Observation planning, timing and bookkeeping utilities for the spectroscopic
monitoring of (massive) stars.

Conventions shared by every function in this module:
    - Coordinates ('coords', 'targets') are accepted as strings, (RA, Dec)
      sequences, dictionaries or astropy SkyCoord objects (ICRS unless
      another frame is given); plain numbers are read in degrees.
    - Observatories ('observatory') are given by a pre-defined code ('OT',
      'ORM', 'CAHA', 'OHP', 'LASILLA', 'PARANAL', 'SAAO', 'MAUNAKEA') or by a
      dictionary with longitude, latitude, elevation and timezone.
    - Times are in UTC unless stated otherwise, and observing nights are
      identified by the date on which they start ('YYYY-MM-DD').
    - Periods are in days, and reference epochs ('hjd0', 'jd0') are
      (Heliocentric) Julian Dates in UTC.
    - Angles (altitudes, azimuths, separations) are in degrees, velocities
      in km/s and exposure times in seconds.
    - Altitudes are geometric (no atmospheric refraction) and airmasses are
      computed as sec(z).
    - Plotting functions return the matplotlib figure; they do not save it
      (use 'panorama.tex.savefig_pkl' to save it).
    - File arguments are paths (path + filename), given either as strings
      or as pathlib.Path objects.
    - All functions validate their inputs before computing; the arguments
      shared by several functions are validated by the same helpers (see
      'support/obs/obs_support.py'), so they give the same error messages.
"""

from __future__ import annotations

# =========================================================================
# ============================== COORDINATES ==============================
# =========================================================================

# ==============================
def sky_dist(
             coords1 : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
             coords2 : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
    ) -> dict[str, float]:
    """
    Compute the angular separation and position angle between two sky positions.

    ----------
    Parameters
    ----------
        coords1 : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the first (reference) position. See Notes.

        coords2 : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the second position. See Notes.

    -------
    Returns
    -------
        results_dict : dict[str, float]
            A dictionary containing the following key-value pairs:
                - 'arcsec' : float
                    Angular separation in arcseconds.

                - 'deg' : float
                    Angular separation in degrees.

                - 'pa' : float
                    Position angle of 'coords2' with respect to 'coords1', in
                    degrees, measured from North through East, in [0, 360).

    ------
    Raises
    ------
        TypeError
            If any coordinate has an unsupported type.

        ValueError
            If any coordinate is empty or cannot be parsed.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. COMPUTATION
            - The separation is computed on the sphere with the Vincenty
              formula (astropy), which is accurate at all separations.

            - If the two positions are given in different frames, 'coords2'
              is transformed into the frame of 'coords1'.

            - Typical use: checking for visual companions or contaminating
              sources within the fiber or slit aperture, or matching targets
              between catalogs.

    ----------
    References
    ----------
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import astropy.units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _coord


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    c1 = _coord(coords1, "coords1")  # Reference position (its frame is used)
    #
    c2 = _coord(coords2, "coords2")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'arcsec' : np.nan,
                    'deg'    : np.nan,
                    'pa'     : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    separation     = c1.separation(c2)  # Angle object, convertible to any angular unit
    position_angle = c1.position_angle(c2).wrap_at(360 * u.deg)  # East of North, wrapped to [0, 360)


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'arcsec' : float(separation.arcsec),
                         'deg'    : float(separation.degree),
                         'pa'     : float(position_angle.degree),
                        })
    return results_dict




# ==============================
def coord_info(
               coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
    ) -> dict[str, float | str]:
    """
    Express a sky position in decimal, sexagesimal and Galactic coordinates.

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the target. See Notes.

    -------
    Returns
    -------
        results_dict : dict[str, float | str]
            A dictionary containing the following key-value pairs:
                - 'ra' : float
                    Right ascension (ICRS) in degrees.

                - 'dec' : float
                    Declination (ICRS) in degrees.

                - 'ra_hms' : str
                    Right ascension as 'hh:mm:ss.ss'.

                - 'dec_dms' : str
                    Declination as '±dd:mm:ss.s'.

                - 'l' : float
                    Galactic longitude in degrees, in [0, 360).

                - 'b' : float
                    Galactic latitude in degrees, in [-90, 90].

    ------
    Raises
    ------
        TypeError
            If the coordinates have an unsupported type.

        ValueError
            If the coordinates are empty or cannot be parsed.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. GALACTIC COORDINATES
            - Coordinates given in another frame are first transformed into
              ICRS, so every output refers to the same position.

            - OB stars are strongly concentrated toward the Galactic plane
              (scale height of a few tens of pc). A large |b| for an OB star
              of known distance can point to a runaway origin (Blaauw 1961).

    ----------
    References
    ----------
        - Blaauw, A. 1961, BAN, 15, 265.
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import astropy.units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _coord


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    c = _coord(coords, "coords").icrs  # Always report ICRS values


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'ra'      : np.nan,
                    'dec'     : np.nan,
                    'ra_hms'  : '',
                    'dec_dms' : '',
                    'l'       : np.nan,
                    'b'       : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  SEXAGESIMAL  ##########
    ra_hms  = str(c.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True))  # 'hh:mm:ss.ss'
    dec_dms = str(c.dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True))  # '±dd:mm:ss.s'


    # ##########  GALACTIC  ##########
    gal = c.galactic  # Transformation ICRS -> Galactic


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'ra'      : float(c.ra.deg),
                         'dec'     : float(c.dec.deg),
                         'ra_hms'  : ra_hms,
                         'dec_dms' : dec_dms,
                         'l'       : float(gal.l.deg),
                         'b'       : float(gal.b.deg),
                        })
    return results_dict




# ==================================================================
# ============================== TIME ==============================
# ==================================================================

# ==============================
def bary_corr(
              coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
              observatory : str | dict,
              obstime     : str | datetime | Time | list[str | datetime | Time],
              exptime     : int | float | list[int | float] | None = None,
              rv          : int | float | list[int | float] | None = None
    ) -> dict[str, SkyCoord | Observer | np.ndarray | list[str] | float]:
    """
    Compute the barycentric radial-velocity correction and the barycentric/heliocentric Julian dates of one or several observations.

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the target. See Notes.

        observatory : str | dict
            Observatory parameters. See Notes.

        obstime : str | datetime | Time | list[str | datetime | Time]
            Observation time(s) in UTC, in any format accepted by
            astropy.time.Time (e.g. '2026-11-01T01:30:00'). It is
            taken as the start of the exposure (see 'exptime').

        exptime : int | float | list[int | float] | None, optional
            Exposure time(s) in seconds. If given, the correction is computed
            at mid-exposure, obstime + exptime / 2. A single value is applied
            to every observation; a list must have the same length as
            'obstime'. Must be non-negative. Default is None ('obstime' is
            used as given).

        rv : int | float | list[int | float] | None, optional
            Measured (topocentric) radial velocity(ies) in km/s to be
            corrected to the barycentric frame. A single value is applied to
            every observation; a list must have the same length as 'obstime'.
            Default is None (no velocity is corrected).

    -------
    Returns
    -------
        results_dict : dict[str, SkyCoord | Observer | np.ndarray | list[str] | float]
            A dictionary containing the following key-value pairs:
                - 'coords' : SkyCoord
                    The input coordinates converted to a SkyCoord object.

                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'obstime' : list[str]
                    Times at which the correction is evaluated (mid-exposure
                    if 'exptime' is given), as 'YYYY-MM-DD HH:MM:SS.sss' UTC.

                - 'jd' : np.ndarray[float]
                    Julian Dates (UTC).

                - 'mjd' : np.ndarray[float]
                    Modified Julian Dates (UTC), MJD = JD - 2400000.5.

                - 'hjd' : np.ndarray[float]
                    Heliocentric Julian Dates (UTC).

                - 'bjd_tdb' : np.ndarray[float]
                    Barycentric Julian Dates in the TDB time scale.

                - 'v_bary' : np.ndarray[float]
                    Barycentric radial-velocity correction in km/s.

                - 'v_helio' : np.ndarray[float]
                    Heliocentric radial-velocity correction in km/s.

                - 'rv_bary' : np.ndarray[float] | float
                    Barycentric radial velocities in km/s ('rv' corrected),
                    or NaN if 'rv' was not given.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the coordinates, the observatory or any time cannot be
            parsed, a value is not finite, an exposure time is negative, or
            a list has a length different from 'obstime'.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. VELOCITY CORRECTION
            - 'v_bary' is the projection of the velocity of the observatory
              with respect to the Solar System barycenter onto the line of
              sight (astropy SkyCoord.radial_velocity_correction). Its
              amplitude is up to ~30 km/s (Earth's orbital motion) plus
              ~0.5 km/s (Earth's rotation).

            - The measured radial velocity is corrected as:
                    rv_bary = rv + v_bary + rv * v_bary / c
              where the last (relativistic) cross term is ~10 m/s for
              rv = 100 km/s and v_bary = 30 km/s.

            - To first order, a spectrum is moved to the barycentric frame
              by scaling its wavelengths:
                    λ_bary = λ_obs * (1 + v_bary / c)

            - The precision of the astropy correction (a few m/s or better)
              is far below what is needed for massive stars, whose lines are
              broadened by tens to hundreds of km/s.

            - Many reduction pipelines already apply this correction or store
              it in the FITS header. Check this before applying it, to avoid
              correcting the spectra twice.

        4. TIME STAMPS
            - BJD_TDB is the recommended time stamp for time-series analysis
              (Eastman et al. 2010). HJD (UTC) differs from it by up to
              ~4 s (barycenter vs. heliocenter) plus the TDB - UTC offset
              (~69 s), which is irrelevant for orbital phases of days but
              matters for timing analyses.

            - The barycentric correction changes by up to ~0.1 km/s per
              hour (Earth's rotation), so for long exposures using the
              mid-exposure time matters at the level of tens of m/s.

    ----------
    References
    ----------
        - Eastman, J., Siverd, R. & Gaudi, B. S. 2010, PASP, 122, 935.
        - Wright, J. T. & Eastman, J. D. 2014, PASP, 126, 838.
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy import constants as cte
    import astropy.units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_per_observation, _check_times, _coord, _get_observer


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coords_obs = _coord(coords, "coords")
    #
    observatory_obs = _get_observer(observatory)
    #
    times   = _check_times(obstime, "obstime", "Observation time")  # Always an array-valued Time (UTC)
    n_times = len(times)
    #
    exptime_arr = np.zeros(n_times) if exptime is None else _check_per_observation(exptime, n_times, "exptime", "Exposure time", "in s", "non-negative")
    #
    rv_arr = None if rv is None else _check_per_observation(rv, n_times, "rv", "Radial velocity", "in km/s")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'coords'      : coords_obs,
                    'observatory' : observatory_obs,
                    'obstime'     : [],
                    'jd'          : np.nan,
                    'mjd'         : np.nan,
                    'hjd'         : np.nan,
                    'bjd_tdb'     : np.nan,
                    'v_bary'      : np.nan,
                    'v_helio'     : np.nan,
                    'rv_bary'     : np.nan,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  MID-EXPOSURE TIMES  ##########
    t_mid    = times + exptime_arr * u.s / 2.  # Equal to 'times' if no exposure time is given
    location = observatory_obs.location  # EarthLocation of the observatory (topocentric corrections)


    # ##########  VELOCITY CORRECTIONS  ##########
    v_bary  = coords_obs.radial_velocity_correction(kind="barycentric", obstime=t_mid, location=location).to(u.km / u.s).value
    v_helio = coords_obs.radial_velocity_correction(kind="heliocentric", obstime=t_mid, location=location).to(u.km / u.s).value


    # ##########  JULIAN DATES  ##########
    ltt_bary  = t_mid.light_travel_time(coords_obs, "barycentric", location=location)  # Light-travel time to the barycenter
    ltt_helio = t_mid.light_travel_time(coords_obs, "heliocentric", location=location)  # Light-travel time to the heliocenter

    jd      = t_mid.utc.jd  # JD (UTC)
    mjd     = t_mid.utc.mjd  # MJD = JD - 2400000.5
    hjd     = jd + ltt_helio.to_value(u.day)  # HJD in UTC (classical convention)
    bjd_tdb = (t_mid.tdb + ltt_bary).jd  # BJD in TDB (Eastman et al. 2010)


    # ##########  CORRECTED RADIAL VELOCITIES  ##########
    if rv_arr is not None:
        c_kms   = cte.c.to(u.km / u.s).value  # Speed of light in km/s
        rv_bary = rv_arr + v_bary + rv_arr * v_bary / c_kms  # Includes the (small) relativistic cross term
    else:
        rv_bary = None  # Reported as NaN


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'obstime' : [t.iso for t in t_mid.utc],  # 'YYYY-MM-DD HH:MM:SS.sss'
                         'jd'      : np.atleast_1d(jd).astype(float),
                         'mjd'     : np.atleast_1d(mjd).astype(float),
                         'hjd'     : np.atleast_1d(hjd).astype(float),
                         'bjd_tdb' : np.atleast_1d(bjd_tdb).astype(float),
                         'v_bary'  : np.atleast_1d(v_bary).astype(float),
                         'v_helio' : np.atleast_1d(v_helio).astype(float),
                         'rv_bary' : np.nan if rv_bary is None else np.atleast_1d(rv_bary).astype(float),
                        })
    return results_dict




# =============================================================================
# ============================== SIGNAL-TO-NOISE ==============================
# =============================================================================

# ==============================
def estimate_exposure(
                      snr_init   : int | float,
                      t_init     : int | float,
                      snr_target : int | float
    ) -> float:
    """
    Estimate the exposure time required to reach a target signal-to-noise ratio.

    ----------
    Parameters
    ----------
        snr_init : int | float
            Signal-to-noise ratio of a reference exposure (dimensionless).
            Must be positive.

        t_init : int | float
            Exposure time of the reference exposure (any time unit, e.g.
            seconds). Must be positive.

        snr_target : int | float
            Desired signal-to-noise ratio (dimensionless). Must be positive.

    -------
    Returns
    -------
        t_final : float
            Estimated exposure time needed to reach 'snr_target', in the same
            units as 't_init'.

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
        - Assumes the photon-noise (source-dominated) regime, where
          SNR ∝ sqrt(t):
                t_final = t_init * (snr_target / snr_init)^2

        - For faint targets, bright sky or short exposures (read-out and
          sky-noise dominated regimes) SNR grows more slowly than sqrt(t),
          and the exposure time is underestimated.

        - The reference exposure should be taken under similar conditions
          (airmass, seeing, Moon, instrument setup); otherwise scale it
          accordingly.

        - Inverse: 'estimate_SNR'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    snr_init = _check_real(snr_init, "snr_init", "Reference signal-to-noise ratio", "dimensionless", "positive")
    #
    t_init = _check_real(t_init, "t_init", "Reference exposure time", "any time unit", "positive")
    #
    snr_target = _check_real(snr_target, "snr_target", "Target signal-to-noise ratio", "dimensionless", "positive")


    # ####################  COMPUTATION  ####################
    t_final = t_init * (snr_target / snr_init)**2  # SNR ∝ sqrt(t)


    # ####################  OUTPUT  ####################
    return float(t_final)




# ==============================
def estimate_SNR(
                 snr_init : int | float,
                 t_init   : int | float,
                 t_final  : int | float
    ) -> float:
    """
    Estimate the signal-to-noise ratio expected for a new exposure time.

    ----------
    Parameters
    ----------
        snr_init : int | float
            Signal-to-noise ratio of a reference exposure (dimensionless).
            Must be positive.

        t_init : int | float
            Exposure time of the reference exposure (any time unit, e.g.
            seconds). Must be positive.

        t_final : int | float
            New exposure time, in the same units as 't_init'. Must be
            positive.

    -------
    Returns
    -------
        snr_final : float
            Estimated signal-to-noise ratio for 't_final'.

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
        - Assumes the photon-noise (source-dominated) regime, where
          SNR ∝ sqrt(t):
                snr_final = snr_init * sqrt(t_final / t_init)

        - In read-out or sky-noise dominated regimes the SNR is
          overestimated.

        - Inverse: 'estimate_exposure'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    snr_init = _check_real(snr_init, "snr_init", "Reference signal-to-noise ratio", "dimensionless", "positive")
    #
    t_init = _check_real(t_init, "t_init", "Reference exposure time", "any time unit", "positive")
    #
    t_final = _check_real(t_final, "t_final", "New exposure time", "same unit as 't_init'", "positive")


    # ####################  COMPUTATION  ####################
    snr_final = snr_init * (t_final / t_init)**0.5  # SNR ∝ sqrt(t)


    # ####################  OUTPUT  ####################
    return float(snr_final)




# ==============================
def combine_SNR(
                snr_list : list[int | float] | tuple[int | float, ...] | np.ndarray
    ) -> float:
    """
    Compute the signal-to-noise ratio of the optimal combination of several spectra.

    ----------
    Parameters
    ----------
        snr_list : list[int | float] | tuple[int | float, ...] | np.ndarray
            Signal-to-noise ratios of the individual spectra (dimensionless).
            Every value must be finite and positive. A NumPy array is also
            accepted and flattened.

    -------
    Returns
    -------
        snr_comb : float
            Signal-to-noise ratio of the combined spectrum.

    ------
    Raises
    ------
        TypeError
            If 'snr_list' is not a list, tuple or array, or any element is
            not a real number.

        ValueError
            If 'snr_list' is empty or contains non-finite or non-positive
            values.

    -----
    Notes
    -----
        - For inverse-variance (1 / σ²) weighted co-addition of spectra with
          the same continuum level and independent noise:
                SNR_comb = sqrt(Σ SNR_i^2)
          so N spectra with equal SNR give SNR * sqrt(N).

        - A plain (unweighted) average of spectra with different SNR gives a
          lower value; the formula is an upper limit in that case.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numpy as np


    # ##### PANORAMA #####
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if isinstance(snr_list, np.ndarray):  # Arrays are accepted and flattened
        snr_list = snr_list.ravel().tolist()
    if not isinstance(snr_list, (list, tuple)):
        raise TypeError(f"Signal-to-noise ratios 'snr_list' must be a list, tuple or numpy array of real numbers. Got {type(snr_list).__name__}.")
    elif len(snr_list) == 0:
        raise ValueError("Signal-to-noise ratios 'snr_list' must contain at least one value. Got an empty sequence.")
    snr_values = [_check_real(snr, "snr_list", "Every element of the signal-to-noise ratios", "dimensionless", "positive") for snr in snr_list]


    # ####################  COMPUTATION  ####################
    snr_comb = np.sqrt(np.sum(np.array(snr_values, dtype=float)**2))  # Quadrature sum (inverse-variance weighting)


    # ####################  OUTPUT  ####################
    return float(snr_comb)




# ==================================================================
# ============================== DATA ==============================
# ==================================================================

# ==============================
def fits_catalog(
                 folder      : str | Path,
                 instrument  : str,
                 output_file : str | Path = "fits_catalog.csv"
    ) -> pd.DataFrame:
    """
    Build a catalog of the FITS files of an instrument with their header metadata and barycentric corrections, and save it to a CSV file.

    ----------
    Parameters
    ----------
        folder : str | Path
            Root directory scanned recursively for FITS files.

        instrument : str
            Instrument the data come from (case-insensitive). Available
            options:
                >>> "HERMES" : HERMES @ Mercator (ORM, La Palma)

        output_file : str | Path, optional
            Path and filename of the output CSV file. Default is
            "fits_catalog.csv".

    -------
    Returns
    -------
        df : pd.DataFrame
            Catalog with one row per FITS file, also written to
            'output_file'. Columns include the original header metadata and
            the following barycentric correction outputs:
            'obstime', 'jd', 'mjd', 'hjd', 'bjd_tdb', 'v_bary_kms' and
            'v_helio_kms'. Missing header values are 'N/A' (the object name is
            'Not specified'), and missing barycentric values are NaN.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        FileNotFoundError
            If 'folder' or the directory of 'output_file' does not exist.

        NotADirectoryError
            If 'folder' or the parent of 'output_file' exists but is not a
            directory.

        IsADirectoryError
            If 'output_file' is an existing directory.

        PermissionError
            If 'folder' cannot be read.

        ValueError
            If 'instrument' is empty or not supported, or 'output_file' is
            not a .csv file.

    -----
    Notes
    -----
        1. FILES
            - Only the files matching the reduced-spectrum pattern of the
              instrument are read (HERMES:
              '*_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fit*').

            - Files that cannot be read issue a warning and are skipped.

            - An existing 'output_file' is overwritten (with a warning).

        2. OBSERVING NIGHT
            - 'night_obs' is the date on which the night started: exposures
              taken before the instrument split hour (e.g., 09 UTC for HERMES)
              are assigned to the previous calendar date. 'time_obs' is the
              UTC time of DATE-OBS.

            - Rows are sorted by observing night and file name.

        3. AIRMASS
            - Taken from the AIRMASS keyword if present; otherwise computed
              from the telescope altitude (TELALT) as sec(z) = 1 / sin(alt).

        4. BARYCENTRIC CORRECTION
            - Computed for every FITS file using 'bary_corr', at mid-exposure
              (DATE-OBS + EXPTIME / 2), from the RA and DEC in the header
              (in degrees) and the instrument observatory.

            - The catalog includes the UTC observation time at mid-exposure,
              JD and MJD in UTC, HJD in UTC, BJD in TDB, and the barycentric
              and heliocentric radial-velocity corrections in km/s.

            - If the coordinates or DATE-OBS are missing or invalid, the
              barycentric outputs are set to NaN and a warning is issued.
              If EXPTIME is missing, zero seconds is assumed.

            - Check whether the reduction pipeline has already applied the
              barycentric correction before correcting spectra, to avoid
              applying it twice.

    ----------
    References
    ----------
        - Raskin, G., van Winckel, H., Hensberge, H. et al. 2011,
          A&A, 526, A69.
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy.io import fits
    import numpy as np
    import os
    import pandas as pd
    from pathlib import Path
    import warnings


    # ##### PANORAMA #####
    from .support.obs.obs_support import _compute_obs_night, _get_header_value, _INSTRUMENTS
    from .support.tex.tex_support import _check_output_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(folder, (str, os.PathLike)):
        raise TypeError(f"Data directory 'folder' must be a path (str or pathlib.Path) to a directory. Got {type(folder).__name__}.")
    folder = Path(folder)  # Converted to Path right after the type check
    if not folder.exists():
        raise FileNotFoundError(f"Data directory not found: '{folder}'.")
    elif not folder.is_dir():
        raise NotADirectoryError(f"Expected a data directory but found a file: '{folder}'.")
    elif not os.access(folder, os.R_OK):
        raise PermissionError(f"Permission denied to read data directory: '{folder}'.")
    #
    if not isinstance(instrument, str):
        raise TypeError(f"Instrument name 'instrument' must be a string. Got {type(instrument).__name__}.")
    elif not instrument.strip():
        raise ValueError(f"Instrument name 'instrument' must not be empty. Available options: {list(_INSTRUMENTS.keys())}. Got '{instrument}'.")
    instrument = instrument.strip().upper()  # Case-insensitive (e.g. 'hermes' -> 'HERMES')
    if instrument not in _INSTRUMENTS:
        raise ValueError(f"Instrument name 'instrument' is not supported. Available options: {list(_INSTRUMENTS.keys())}. Got '{instrument}'.")
    #
    output_file = _check_output_file(output_file)
    if output_file.suffix.lower() != ".csv":
        raise ValueError(f"Output catalog 'output_file' must have a .csv extension. Got '{output_file.name}'.")
    elif output_file.exists():
        warnings.warn(f"Output catalog '{output_file}' already exists and will be overwritten.", stacklevel=2)


    # ####################  PARAMETERS  ####################
    instrument_config = _INSTRUMENTS[instrument]
    pattern           = instrument_config['pattern']  # Glob pattern of the reduced spectra
    night_split_hour  = instrument_config['night_split_hour']  # UTC hour at which the observing night changes date
    header_keywords   = instrument_config['header_keywords']  # Catalog column -> FITS keywords to try
    columns = [
               'file', 'telescope', 'instrument', 'object', 'ra', 'dec',
               'airmass', 'exptime_s', 'filter', 'date_obs', 'night_obs',
               'time_obs', 'obstime', 'jd', 'mjd', 'hjd', 'bjd_tdb',
               'v_bary_kms', 'v_helio_kms'
              ]  # Keys of every record, in the column order of the catalog
    missing = ("N/A", None, '')  # Header values treated as missing


    # ####################  COMPUTATION  ####################
    # ##########  FITS FILES  ##########
    fits_files = sorted(set(folder.rglob(pattern)))  # Recursive search; set() removes duplicates, sorted() gives a reproducible order

    if not fits_files:
        warnings.warn(f"No FITS files matching the pattern '{pattern}' were found in '{folder}'.", stacklevel=2)


    # ##########  HEADERS  ##########
    records = []  # One dictionary per readable file

    for f in fits_files:
        try:
            with fits.open(f) as hdul:
                header = hdul[0].header  # Metadata are in the primary header; the file is closed right after


            # ##### DATES #####
            date_obs_val                = _get_header_value(header, header_keywords['date_obs'], "N/A")  # ISO date and time (UTC)
            night_obs_val, time_obs_val = _compute_obs_night(date_obs_val if date_obs_val != "N/A" else None, night_split_hour)  # Night on which the exposure was taken, and UTC time


            # ##### AIRMASS #####
            airmass_val = _get_header_value(header, header_keywords['airmass'], None)
            if airmass_val is None:  # No AIRMASS keyword: computed from the telescope altitude
                alt_val = _get_header_value(header, header_keywords['altitude'], None)
                try:
                    alt_val     = float(alt_val)  # None or non-numeric values raise and give 'N/A'
                    airmass_val = f"{1. / np.sin(np.radians(alt_val)):.4f}" if alt_val > 0 else "N/A"  # sec(z) = 1 / sin(alt)
                except (TypeError, ValueError):
                    airmass_val = "N/A"


            # ##### RECORD #####
            record = {
                      'file'        : f.name,
                      'telescope'   : _get_header_value(header, header_keywords['telescope'], "N/A"),
                      'instrument'  : _get_header_value(header, header_keywords['instrument'], "N/A"),
                      'object'      : _get_header_value(header, header_keywords['object'], "Not specified"),
                      'ra'          : _get_header_value(header, header_keywords['ra'], "N/A"),
                      'dec'         : _get_header_value(header, header_keywords['dec'], "N/A"),
                      'airmass'     : airmass_val,
                      'exptime_s'   : _get_header_value(header, header_keywords['exptime_s'], "N/A"),
                      'filter'      : _get_header_value(header, header_keywords['filter'], "N/A"),
                      'date_obs'    : date_obs_val,
                      'night_obs'   : night_obs_val,
                      'time_obs'    : time_obs_val,
                      'obstime'     : np.nan,
                      'jd'          : np.nan,
                      'mjd'         : np.nan,
                      'hjd'         : np.nan,
                      'bjd_tdb'     : np.nan,
                      'v_bary_kms'  : np.nan,
                      'v_helio_kms' : np.nan,
                     }


            # ##### BARYCENTRIC CORRECTION #####
            try:
                if record['ra'] in missing or record['dec'] in missing or date_obs_val in missing:
                    raise ValueError("Missing RA, DEC or DATE-OBS in the FITS header.")

                coords  = (float(record['ra']), float(record['dec']))  # Header RA and DEC, in degrees
                exptime = float(record['exptime_s']) if record['exptime_s'] not in missing else 0.  # Missing EXPTIME: start of the exposure used

                bary_results = bary_corr(  # At mid-exposure (DATE-OBS + EXPTIME / 2)
                                         coords=coords,
                                         observatory=instrument_config['observatory'],
                                         obstime=date_obs_val,
                                         exptime=exptime
                                        )

                record['obstime']     = bary_results['obstime'][0]  # Single observation: first element of each array
                record['jd']          = float(bary_results['jd'][0])
                record['mjd']         = float(bary_results['mjd'][0])
                record['hjd']         = float(bary_results['hjd'][0])
                record['bjd_tdb']     = float(bary_results['bjd_tdb'][0])
                record['v_bary_kms']  = float(bary_results['v_bary'][0])
                record['v_helio_kms'] = float(bary_results['v_helio'][0])

            except Exception as e:  # Non-fatal: the record is kept with NaN barycentric values
                warnings.warn(
                              f"Could not compute barycentric correction for '{f.name}' (barycentric values set to NaN): {e}",
                              stacklevel=2
                             )
            records.append(record)

        except Exception as e:  # Unreadable or corrupt file: skipped, the loop goes on
            warnings.warn(f"Error reading '{f.name}' (file skipped): {e}", stacklevel=2)


    # ##########  CATALOG  ##########
    df = pd.DataFrame(records, columns=columns)  # Fixed column order (also for an empty catalog)

    if not df.empty:  # Sorted by observing night and file name for readability
        df = df.sort_values(by=['night_obs', 'file'], na_position="last").reset_index(drop=True)


    # ####################  OUTPUT  ####################
    df.to_csv(output_file, index=False)  # Without the DataFrame index column
    return df




# ========================================================================
# ============================== VISIBILITY ==============================
# ========================================================================

# ==============================
def night_info(
               observatory : str | dict,
               night       : str | datetime | Time | None = None
    ) -> dict[str, Observer | str | float | None]:
    """
    Compute the Sun and Moon ephemerides of an observing night: sunset, twilights, sunrise, dark time and Moon phase.

    ----------
    Parameters
    ----------
        observatory : str | dict
            Observatory parameters. See Notes.

        night : str | datetime | Time | None, optional
            Night of observation, identified by the date on which it starts
            ('YYYY-MM-DD'). Default is None (the night starting on the
            current UTC date).

    -------
    Returns
    -------
        results_dict : dict[str, Observer | str | float | None]
            A dictionary containing the following key-value pairs:
                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'night' : str
                    Night of observation ('YYYY-MM-DD').

                - 'utc_offset' : float
                    Local time minus UTC at midnight, in hours (NaN if it
                    cannot be computed).

                - 'sunset' : str | None
                    Sunset (upper limb, with refraction), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'civil_evening' : str | None
                    End of the evening civil twilight (Sun at -6 deg), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'nautical_evening' : str | None
                    End of the evening nautical twilight (Sun at -12 deg), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'astro_evening' : str | None
                    End of the evening astronomical twilight (Sun at -18 deg),
                    as 'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'midnight' : str | None
                    Local solar midnight, as 'YYYY-MM-DD HH:MM:SS' UTC.

                - 'astro_morning' : str | None
                    Start of the morning astronomical twilight (Sun at -18
                    deg), as 'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not
                    occur.

                - 'nautical_morning' : str | None
                    Start of the morning nautical twilight (Sun at -12 deg),
                    as 'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'civil_morning' : str | None
                    Start of the morning civil twilight (Sun at -6 deg), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'sunrise' : str | None
                    Sunrise (upper limb, with refraction), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if it does not occur.

                - 'night_hours' : float
                    Time between sunset and sunrise, in hours (NaN if any of
                    them does not occur).

                - 'dark_hours' : float
                    Time between the astronomical twilights, in hours (0.0 if
                    the Sun never reaches -18 deg, 24.0 if it stays below -18
                    deg during the whole day).

                - 'lst_midnight' : str
                    Local sidereal time at midnight ('hh:mm:ss'); objects
                    with RA ≈ LST transit at midnight.

                - 'moon_illumination' : float
                    Illuminated fraction of the Moon at midnight, in percent
                    (0-100).

                - 'moon_alt_midnight' : float
                    Altitude of the Moon at midnight, in degrees.

                - 'moonrise' : str | None
                    Moonrise during the night (between sunset and sunrise), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if the Moon does not
                    rise within the night.

                - 'moonset' : str | None
                    Moonset during the night (between sunset and sunrise), as
                    'YYYY-MM-DD HH:MM:SS' UTC, or None if the Moon does not
                    set within the night.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the observatory or 'night' cannot be parsed.

    -----
    Notes
    -----
        1. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        2. DEFINITIONS
            - Sunset and sunrise use a horizon of -0.833 deg (solar radius
              plus standard refraction). Twilights use the geometric
              position of the Sun's center.

            - The night is the one that starts in the evening of 'night'
              (local date) and ends the next morning: its local midnight is
              searched around the mean local midnight (see Notes of
              'visibility').

            - The night runs from sunset to sunrise; where the Sun does not
              set or rise (polar day or night), it covers ±12 h around
              midnight. Moonrise and moonset are only reported if they
              happen within these limits, as in 'visibility'.

            - Dark time is the time with the Sun below -18 deg (between the
              astronomical twilights), the same definition used by
              'visibility', 'annual_visibility', 'periodic_event' and
              'obs_phase_night'.

            - Bright time is commonly defined by a Moon illumination
              ≳ 50-70 % with the Moon above the horizon. Bright time is
              usually enough for high-SNR spectroscopy of bright OB stars,
              while faint targets need dark or gray time.

    ----------
    References
    ----------
        - Morris, B. M., Tollerud, E., Sipőcz, B. et al. 2018, AJ, 155, 128.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import astropy.units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_night, _event_hours, _event_iso, _event_within, _get_observer, _is_event, _mean_midnight, _solar_midnight


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    observatory_obs = _get_observer(observatory)
    #
    night_str = _check_night(night, "night", "Observing night")  # 'YYYY-MM-DD'; None -> current UTC date


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'observatory'       : observatory_obs,
                    'night'             : night_str,
                    'utc_offset'        : np.nan,
                    'sunset'            : None,
                    'civil_evening'     : None,
                    'nautical_evening'  : None,
                    'astro_evening'     : None,
                    'midnight'          : None,
                    'astro_morning'     : None,
                    'nautical_morning'  : None,
                    'civil_morning'     : None,
                    'sunrise'           : None,
                    'night_hours'       : np.nan,
                    'dark_hours'        : np.nan,
                    'lst_midnight'      : '',
                    'moon_illumination' : np.nan,
                    'moon_alt_midnight' : np.nan,
                    'moonrise'          : None,
                    'moonset'           : None,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  MIDNIGHT  ##########
    midnight = _solar_midnight(observatory_obs, _mean_midnight(observatory_obs, night_str))  # True solar midnight (Sun at lower culmination) of the night starting on 'night_str'


    # ##########  SUN  ##########
    sunset  = observatory_obs.sun_set_time(midnight, which="previous", horizon=-0.8333 * u.deg)  # Upper limb + refraction
    sunrise = observatory_obs.sun_rise_time(midnight, which="next", horizon=-0.8333 * u.deg)
    civil_e = observatory_obs.twilight_evening_civil(midnight, which="previous")  # Sun at -6 deg; masked if it does not occur
    civil_m = observatory_obs.twilight_morning_civil(midnight, which="next")
    naut_e  = observatory_obs.twilight_evening_nautical(midnight, which="previous")  # Sun at -12 deg
    naut_m  = observatory_obs.twilight_morning_nautical(midnight, which="next")
    astro_e = observatory_obs.twilight_evening_astronomical(midnight, which="previous")  # Sun at -18 deg
    astro_m = observatory_obs.twilight_morning_astronomical(midnight, which="next")


    # ##########  NIGHT LIMITS  ##########
    # Same limits as in 'visibility': sunset to sunrise, or ±12 h around midnight in polar day or night
    if _is_event(sunset) and _is_event(sunrise):
        night_start, night_end = sunset, sunrise
    else:
        night_start, night_end = midnight - 12 * u.hour, midnight + 12 * u.hour


    # ##########  NIGHT LENGTH  ##########
    night_hours = _event_hours(sunset, sunrise)  # NaN if the Sun does not set or rise
    if _is_event(astro_e) and _is_event(astro_m):
        dark_hours = _event_hours(astro_e, astro_m)  # Between the astronomical twilights
    else:  # No -18 deg crossing: no dark time, or dark during the whole ±12 h window (polar night)
        dark_hours = 24. if observatory_obs.sun_altaz(midnight).alt.deg < -18. else 0.


    # ##########  SIDEREAL TIME  ##########
    lst = observatory_obs.local_sidereal_time(midnight)  # Apparent local sidereal time


    # ##########  MOON  ##########
    moon_illum = observatory_obs.moon_illumination(midnight)  # Illuminated fraction (0-1)
    moon_alt   = observatory_obs.moon_altaz(midnight).alt.deg  # Geometric altitude, in degrees
    moonrise   = _event_within(observatory_obs.moon_rise_time(night_start, which="next"), night_start, night_end)  # None if not within the night
    moonset    = _event_within(observatory_obs.moon_set_time(night_start, which="next"), night_start, night_end)


    # ##########  LOCAL TIME  ##########
    try:
        utc_offset = midnight.to_datetime(timezone=observatory_obs.timezone).utcoffset().total_seconds() / 3600.  # Includes daylight saving time
    except Exception:  # Unknown timezone offset: reported as NaN
        utc_offset = np.nan


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'utc_offset'        : float(utc_offset),
                         'sunset'            : _event_iso(sunset),
                         'civil_evening'     : _event_iso(civil_e),
                         'nautical_evening'  : _event_iso(naut_e),
                         'astro_evening'     : _event_iso(astro_e),
                         'midnight'          : _event_iso(midnight),
                         'astro_morning'     : _event_iso(astro_m),
                         'nautical_morning'  : _event_iso(naut_m),
                         'civil_morning'     : _event_iso(civil_m),
                         'sunrise'           : _event_iso(sunrise),
                         'night_hours'       : night_hours,
                         'dark_hours'        : float(dark_hours),
                         'lst_midnight'      : str(lst.to_string(unit=u.hourangle, sep=":", precision=0, pad=True)),
                         'moon_illumination' : float(moon_illum) * 100,
                         'moon_alt_midnight' : float(moon_alt),
                         'moonrise'          : _event_iso(moonrise),
                         'moonset'           : _event_iso(moonset),
                        })
    return results_dict




# ==============================
def visibility(
               targets     : dict | str | tuple[str | float, str | float] | list[str | float] | SkyCoord,
               observatory : str | dict,
               night       : str | datetime | Time | None = None,
               min_alt     : int | float | u.Quantity = 30,
               show        : bool = True
    ) -> dict[str, Observer | str | float | np.ndarray | dict | Figure]:
    """
    Compute and plot the altitude of one or several targets along an observing night (StarAlt-like visibility plot).

    ----------
    Parameters
    ----------
        targets : dict | str | tuple[str | float, str | float] | list[str | float] | SkyCoord
            Either the coordinates of a single target, or a dictionary
            {name: coordinates} with several targets, e.g.
                >>> {'HD 37042': "05:35:17.3 -05:23:28", 'HD 91316': (158.0, 9.3)}
            See Notes for the accepted coordinate formats.

        observatory : str | dict
            Observatory parameters. See Notes.

        night : str | datetime | Time | None, optional
            Night of observation, identified by the date on which it starts
            ('YYYY-MM-DD'). Default is None (the night starting on the
            current UTC date).

        min_alt : int | float | u.Quantity, optional
            Minimum altitude for a target to be considered observable, in
            [0, 90] deg. If int or float, it is read in degrees. Default is
            30 (airmass 2).

        show : bool, optional
            If True, display the figure with plt.show(). Default is True.

    -------
    Returns
    -------
        results_dict : dict[str, Observer | str | float | np.ndarray | dict | Figure]
            A dictionary containing the following key-value pairs:
                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'night' : str
                    Night of observation ('YYYY-MM-DD').

                - 'times' : np.ndarray[str]
                    UTC times of the grid ('YYYY-MM-DD HH:MM'), from sunset
                    to sunrise in steps of 2 minutes.

                - 'moon_alt' : np.ndarray[float]
                    Altitude of the Moon along the grid, in degrees.

                - 'moon_illumination' : float
                    Illuminated fraction of the Moon at midnight, in percent
                    (0-100).

                - 'events' : dict[str, str | None]
                    A dictionary keyed by night event ('sunset',
                    'civil_evening', 'nautical_evening', 'astro_evening',
                    'midnight', 'astro_morning', 'nautical_morning',
                    'civil_morning', 'sunrise', 'moonrise', 'moonset'); each
                    value is its UTC time ('YYYY-MM-DD HH:MM'), or None if
                    the event does not happen within the plotted night.

                - 'night_hours' : float
                    Length of the night (sunset to sunrise), in hours (NaN
                    if the Sun does not set or rise).

                - 'dark_hours' : float
                    Dark time (between the evening and morning astronomical
                    twilights), in hours (0.0 if the Sun never reaches
                    -18 deg, 24.0 if it stays below -18 deg during the whole
                    day).

                - 'targets' : dict[str, dict]
                    A dictionary keyed by target name ('Target' for a single
                    target given without a name); each value is a dictionary
                    containing the following key-value pairs:
                        - 'altitude' : np.ndarray[float]
                            Altitude along the grid, in degrees.

                        - 'airmass' : np.ndarray[float]
                            Airmass sec(z) along the grid (NaN below the
                            horizon).

                        - 'transit' : str
                            UTC time of maximum altitude within the night
                            ('YYYY-MM-DD HH:MM').

                        - 'max_alt' : float
                            Maximum altitude within the night, in degrees.

                        - 'start' : str | None
                            UTC time ('YYYY-MM-DD HH:MM') at which the target
                            becomes observable (Sun below -18 deg and
                            altitude >= 'min_alt'), or None if never
                            observable.

                        - 'end' : str | None
                            UTC time ('YYYY-MM-DD HH:MM') at which the target
                            stops being observable, or None if never
                            observable.

                        - 'hours' : float
                            Hours of dark time (Sun below -18 deg) with
                            altitude >= 'min_alt'.

                        - 'moon_sep' : float
                            Separation from the Moon at midnight, in degrees.

                - 'fig' : Figure
                    The matplotlib visibility figure.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If any coordinate, the observatory or the night cannot be
            parsed, a target name is empty, or 'min_alt' is not finite or
            lies outside [0, 90] deg.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

            - A dictionary is interpreted as a single target when all its
              keys are coordinate keys ('ra', 'dec', 'right_ascension',
              'declination', 'unit', 'frame'), and as {name: coordinates}
              otherwise.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. NIGHT
            - The grid runs from sunset to sunrise (Sun at -0.833 deg)
              around the local midnight of the night starting on 'night'.
              Where the Sun does not set or rise (polar day or night), it
              covers ±12 h around midnight.

            - The local midnight is searched around the mean local midnight
              (00:00 UTC of the next day minus longitude / 15 h), which
              gives the right night at every longitude.

            - Dark time is the time with the Sun below -18 deg. The
              observable window of a target (Sun below -18 deg and altitude
              >= 'min_alt') is interpolated linearly between grid points, so
              that 'start', 'end' and 'hours' do not depend on the 2-min
              grid (the same step as in 'obs_phase_night'). The same definition
              is used by 'night_info', 'annual_visibility', 'periodic_event'
              and 'obs_phase_night'.

            - Times are rounded to the nearest minute.

        4. PLOT
            - Altitude of every target vs. UTC time; the transit is marked
              with a dot, and the legend gives the observable dark hours.
              The right axis gives the airmass sec(z).

            - Shaded areas mark the day plus civil twilight, the nautical
              twilight and the astronomical twilight (darker means brighter
              sky); the dashed gray line is the Moon and the dotted red line
              is 'min_alt'. The limits between shaded areas are interpolated
              between grid points, so that consecutive areas touch without
              gaps or overlaps.

            - Sunset, sunrise and the civil, nautical and astronomical
              twilights are marked with dotted vertical lines, labeled on
              top of the plot.

            - Below the legend, a table gives the UTC times of the night
              events (sunset, twilights, midnight, sunrise, moonrise and
              moonset), the length of the night (sunset to sunrise) and the
              dark time (between astronomical twilights), as hh:mm; '--:--'
              means that the event does not happen within the plotted night.

            - Airmass is computed as sec(z), accurate to ~1% for z < 70 deg
              (airmass < 3); closer to the horizon it is overestimated.

    ----------
    References
    ----------
        - Morris, B. M., Tollerud, E., Sipőcz, B. et al. 2018, AJ, 155, 128.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import astropy.units as u
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_bool, _check_min_alt, _check_night, _coord, _event_minute, _format_hhmm, _get_observer, _is_event, _mean_midnight, _observable, _solar_midnight


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coord_keys = {'ra', 'dec', 'right_ascension', 'declination', 'unit', 'frame'}
    if isinstance(targets, dict) and len(targets) > 0 and not {str(k).strip().lower() for k in targets} <= coord_keys:  # {name: coordinates}
        targets_obs = {}  # Target name -> SkyCoord, in the input order
        for name, c in targets.items():
            if not isinstance(name, str):
                raise TypeError(f"Target names (keys of 'targets') must be strings. Got {type(name).__name__}.")
            elif not name.strip():
                raise ValueError(f"Target names (keys of 'targets') must not be empty. Got '{name}'.")
            targets_obs[name] = _coord(c, f"targets['{name}']")  # The target name appears in the error messages
    else:  # Single target
        targets_obs = {'Target' : _coord(targets, "targets")}
    #
    observatory_obs = _get_observer(observatory)
    #
    night_str = _check_night(night, "night", "Observing night")
    #
    min_alt = _check_min_alt(min_alt)  # Scalar Quantity in degrees
    #
    show = _check_bool(show, "show", "Display flag")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'observatory'       : observatory_obs,
                    'night'             : night_str,
                    'times'             : np.nan,
                    'moon_alt'          : np.nan,
                    'moon_illumination' : np.nan,
                    'events'            : {},
                    'night_hours'       : np.nan,
                    'dark_hours'        : np.nan,
                    'targets'           : {},
                    'fig'               : None,
                   }


    # ####################  PARAMETERS  ####################
    # ##########  NIGHT LIMITS  ##########
    midnight = _solar_midnight(observatory_obs, _mean_midnight(observatory_obs, night_str))  # True solar midnight (Sun at lower culmination) of the night starting on 'night_str'

    sunset  = observatory_obs.sun_set_time(midnight, which="previous", horizon=-0.8333 * u.deg)  # Upper limb + refraction
    sunrise = observatory_obs.sun_rise_time(midnight, which="next", horizon=-0.8333 * u.deg)
    if _is_event(sunset) and _is_event(sunrise):
        grid_start, grid_end = sunset, sunrise
    else:  # Polar day or night: fall back to ±12 h around midnight
        grid_start, grid_end = midnight - 12 * u.hour, midnight + 12 * u.hour


    # ##########  TIME GRID  ##########
    step      = 2.  # Minutes
    n_steps   = int(np.ceil((grid_end - grid_start).to(u.min).value / step)) + 1  # Enough points to reach 'grid_end'
    time_grid = grid_start + np.arange(n_steps) * step * u.min  # Array-valued Time (UTC)


    # ##########  SUN & MOON  ##########
    sun_alt       = observatory_obs.sun_altaz(time_grid).alt.deg  # Defines twilights and dark time along the grid
    moon_alt      = observatory_obs.moon_altaz(time_grid).alt.deg  # Plotted as a dashed line
    moon_illum    = float(observatory_obs.moon_illumination(midnight))  # Illuminated fraction (0-1)
    moon_midnight = observatory_obs.moon_altaz(midnight)  # Position of the Moon at midnight (for the Moon separation)


    # ##########  NIGHT EVENTS  ##########
    # Event key -> (astroplan time, label in the event list, short label on top of the plot or None)
    night_events = {
                    'sunset'           : (sunset, "Sunset", "Sunset"),
                    'civil_evening'    : (observatory_obs.twilight_evening_civil(midnight, which="previous"), "Civil twilight (end)", "Civ. twil."),
                    'nautical_evening' : (observatory_obs.twilight_evening_nautical(midnight, which="previous"), "Nautical twilight (end)", "Naut. twil."),
                    'astro_evening'    : (observatory_obs.twilight_evening_astronomical(midnight, which="previous"), "Astro. twilight (end)", "Astro. twil."),
                    'midnight'         : (midnight, "Midnight", None),
                    'astro_morning'    : (observatory_obs.twilight_morning_astronomical(midnight, which="next"), "Astro. twilight (start)", "Astro. twil."),
                    'nautical_morning' : (observatory_obs.twilight_morning_nautical(midnight, which="next"), "Nautical twilight (start)", "Naut. twil."),
                    'civil_morning'    : (observatory_obs.twilight_morning_civil(midnight, which="next"), "Civil twilight (start)", "Civ. twil."),
                    'sunrise'          : (sunrise, "Sunrise", "Sunrise"),
                    'moonrise'         : (observatory_obs.moon_rise_time(grid_start, which="next"), "Moonrise", None),
                    'moonset'          : (observatory_obs.moon_set_time(grid_start, which="next"), "Moonset", None),
                   }
    events = {key: _event_minute(t, grid_start, grid_end) for key, (t, _, _) in night_events.items()}  # None if not within the night


    # ##########  NIGHT LENGTH  ##########
    astro_evening = night_events['astro_evening'][0]
    astro_morning = night_events['astro_morning'][0]
    night_hours   = float((sunrise - sunset).to(u.hour).value) if _is_event(sunset) and _is_event(sunrise) else np.nan  # Sunset to sunrise
    if _is_event(astro_evening) and _is_event(astro_morning):
        dark_hours = float((astro_morning - astro_evening).to(u.hour).value)  # Between astronomical twilights
    else:  # No -18 deg crossing: no dark time, or dark during the whole ±12 h window (polar night)
        dark_hours = 24. if observatory_obs.sun_altaz(midnight).alt.deg < -18. else 0.


    # ####################  COMPUTATION  ####################
    targets_results = {}  # Filled with one entry per target (keys are the target names)
    times_str       = np.array([str(t) for t in time_grid.utc.strftime("%Y-%m-%d %H:%M")])  # Grid times as strings (returned in 'times')
    min_alt_deg     = float(min_alt.to(u.deg).value)  # Plain float, for NumPy comparisons

    for name, c in targets_obs.items():


        # ##########  ALTITUDE & AIRMASS  ##########
        alt     = observatory_obs.altaz(time_grid, c).alt.deg  # Geometric altitude along the grid, in degrees
        airmass = np.where(alt > 0, 1. / np.sin(np.radians(np.clip(alt, 1e-6, 90))), np.nan)  # sec(z), NaN below the horizon


        # ##########  TRANSIT  ##########
        i_max = int(np.argmax(alt))  # Grid point of maximum altitude within the night


        # ##########  OBSERVABLE WINDOW  ##########
        lo, hi, length = _observable(sun_alt, alt, min_alt_deg)  # Interpolated between grid points
        if np.any(length > 0):
            idx        = np.where(length > 0)[0]  # Grid intervals with an observable part
            t_first    = time_grid[idx[0]] + lo[idx[0]] * step * u.min  # Start of the first observable part
            t_last     = time_grid[idx[-1]] + hi[idx[-1]] * step * u.min  # End of the last observable part
            start, end = _event_minute(t_first, grid_start, grid_end), _event_minute(t_last, grid_start, grid_end)
        else:
            start, end = None, None
        hours = float(np.sum(length) * step / 60.)  # Fractions of interval -> hours


        # ##########  MOON SEPARATION  ##########
        moon_sep = float(observatory_obs.altaz(midnight, c).separation(moon_midnight).deg)  # Both positions in the same AltAz frame


        # ##########  RESULTS  ##########
        targets_results[name] = {
                                 'altitude' : alt,
                                 'airmass'  : airmass,
                                 'transit'  : str(times_str[i_max]),
                                 'max_alt'  : float(alt[i_max]),
                                 'start'    : start,
                                 'end'      : end,
                                 'hours'    : hours,
                                 'moon_sep' : moon_sep,
                                }


    # ####################  PLOTTING  ####################
    # ##########  INIT FIGURE  ##########
    fig, ax = plt.subplots(figsize=(12, 6))
    t_plot  = time_grid.to_datetime()  # Python datetimes (UTC) for the Matplotlib date axis


    # ##########  TWILIGHTS  ##########
    # Sky level of each grid point (0: dark, 1: astronomical, 2: nautical, 3: civil twilight or day). The
    # limits between levels are interpolated where the Sun crosses -18, -12 or -6 deg, so that consecutive
    # shaded areas touch exactly (no white gaps between them)
    thresholds = [-18., -12., -6.]  # Sun altitude, in deg
    alphas = {  # Darker shading = brighter sky
              1 : 0.15,
              2 : 0.30,
              3 : 0.45,
             }
    t_num     = mdates.date2num(t_plot)  # Matplotlib date numbers (days), interpolated below
    sky_level = np.digitize(sun_alt, thresholds)  # 0: < -18, 1: [-18, -12), 2: [-12, -6), 3: >= -6 deg
    changes   = np.where(np.diff(sky_level) != 0)[0]  # Grid intervals where the sky level changes

    edges  = [t_num[0]]  # Limits of the constant-level segments
    levels = [sky_level[0]]  # Sky level of each segment
    for i in changes:
        threshold = thresholds[min(sky_level[i], sky_level[i + 1])]  # Threshold crossed between samples i and i + 1
        fraction  = (threshold - sun_alt[i]) / (sun_alt[i + 1] - sun_alt[i])  # Linear interpolation of the crossing
        edges.append(t_num[i] + fraction * (t_num[i + 1] - t_num[i]))
        levels.append(sky_level[i + 1])
    edges.append(t_num[-1])

    for (t_left, t_right), level in zip(zip(edges[:-1], edges[1:]), levels):
        if level > 0:  # Dark time is not shaded
            ax.axvspan(t_left, t_right, color="steelblue", alpha=alphas[level], lw=0, zorder=0)


    # ##########  NIGHT EVENTS (TOP MARKERS)  ##########
    n_top_labels = 0
    for key, (t, _, short_label) in night_events.items():
        if short_label is None or events[key] is None:  # Only Sun events that happen within the night
            continue
        n_top_labels += 1  # Used to leave room above the axes for the labels
        t_event       = t.to_datetime()
        ax.text(
                t_event,
                1.01,
                short_label,
                transform=ax.get_xaxis_transform(),  # x in data coordinates, y in axes coordinates
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=8,
                color="dimgray"
               )


    # ##########  TARGETS  ##########
    for name, res in targets_results.items():
        line, = ax.plot(t_plot, res['altitude'], lw=2, label=f"{name} ({res['hours']:.1f} h)")  # Observable dark hours in the legend
        ax.plot(t_plot[int(np.argmax(res['altitude']))], res['max_alt'], "o", color=line.get_color(), ms=6)  # Transit


    # ##########  MOON & LIMITS  ##########
    ax.plot(t_plot, moon_alt, ls="--", color="gray", lw=1.5, label=f"Moon ({100 * moon_illum:.0f}%)")
    ax.axhline(min_alt_deg, color="red", ls=":", lw=1.2)  # Minimum altitude


    # ##########  AXES  ##########
    # ##### X AXIS #####
    ax.set_xlim(t_plot[0], t_plot[-1])  # Sunset to sunrise (or ±12 h)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))  # UTC hours
    ax.set_xlabel(f"Time [UTC] (night starting {night_str})")


    # ##### Y AXIS #####
    ax.set_ylim(0, 90)
    ax.set_ylabel("Altitude [deg]")


    # ##### AIRMASS AXIS #####
    ax_am         = ax.twinx()  # Right axis sharing the altitude scale
    airmass_ticks = np.array([1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0, 5.0])
    ax_am.set_ylim(0, 90)  # Same limits as the altitude axis, so the ticks line up
    ax_am.set_yticks(90. - np.degrees(np.arccos(1. / airmass_ticks)))  # Altitude of each airmass: alt = 90 - arccos(1 / X)
    ax_am.set_yticklabels([f"{x:.1f}" for x in airmass_ticks])
    ax_am.set_ylabel("Airmass")


    # ##### TITLE & LEGEND #####
    ax.set_title(f"Visibility — {observatory_obs.name} — {night_str}", fontsize=14, pad=70 if n_top_labels else 6)  # Room for the event labels (in points)
    ax.grid(True, ls="--", alpha=0.4)
    legend = ax.legend(loc="upper left", bbox_to_anchor=(1.07, 1.0), fontsize=9, frameon=False)  # Outside the axes, to the right of the airmass axis


    # ##########  NIGHT EVENTS (TABLE BELOW THE LEGEND)  ##########
    event_lines = ["Night events [UTC]"]  # Monospaced table: label padded to 26 characters + time
    for key, (_, label, _) in night_events.items():
        event_lines.append(f"{label:<26}{events[key][-5:] if events[key] is not None else '--:--'}")  # 'HH:MM' or '--:--' if not within the night
    event_lines += [
                    '',
                    f"{'Night length':<26}{_format_hhmm(night_hours)}",
                    f"{'Dark time':<26}{_format_hhmm(dark_hours)}"
                   ]
    events_text = ax.text(
                          1.075,
                          0.95,
                          "\n".join(event_lines),
                          transform=ax.transAxes,
                          ha="left",
                          va="top",
                          fontsize=8,
                          family="monospace"
                         )


    # ##########  SHOW PLOT  ##########
    fig.tight_layout()  # Leaves room for the legend, the event table and the top labels

    fig.canvas.draw()  # The table is placed just below the legend, whose height is only known once drawn
    legend_box = legend.get_window_extent().transformed(ax.transAxes.inverted())  # Legend box in axes coordinates
    events_text.set_position((legend_box.x0 + 0.005, legend_box.y0 - 0.04))  # Just below the legend

    if show:
        plt.show()  # Inline backend: displays the figure; interactive backends: opens a window


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'times'             : times_str,
                         'moon_alt'          : moon_alt,
                         'moon_illumination' : moon_illum * 100,  # Fraction -> percent
                         'events'            : events,
                         'night_hours'       : night_hours,
                         'dark_hours'        : dark_hours,
                         'targets'           : targets_results,
                         'fig'               : fig,
                        })
    return results_dict




# ==============================
def annual_visibility(
                      coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                      observatory : str | dict,
                      start_date  : str | datetime | Time | None = None,
                      n_nights    : int = 365,
                      min_alt     : int | float | u.Quantity = 30,
                      show        : bool = True
    ) -> dict[str, SkyCoord | Observer | list[str] | np.ndarray | str | int | Figure]:
    """
    Compute and plot the altitude and the observable dark hours of a target on each night over a period (e.g. a semester or a year).

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the target. See Notes.

        observatory : str | dict
            Observatory parameters. See Notes.

        start_date : str | datetime | Time | None, optional
            First night ('YYYY-MM-DD'). Default is None (the night starting
            on the current UTC date).

        n_nights : int, optional
            Number of consecutive nights. Must be positive. Default is 365.

        min_alt : int | float | u.Quantity, optional
            Minimum altitude for the target to be considered observable, in
            [0, 90] deg. If int or float, it is read in degrees. Default is
            30.

        show : bool, optional
            If True, display the figure with plt.show(). Default is True.

    -------
    Returns
    -------
        results_dict : dict[str, SkyCoord | Observer | list[str] | np.ndarray | str | int | Figure]
            A dictionary containing the following key-value pairs:
                - 'coords' : SkyCoord
                    The input coordinates converted to a SkyCoord object.

                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'nights' : list[str]
                    Nights ('YYYY-MM-DD').

                - 'hours' : np.ndarray[float]
                    Dark hours (Sun below -18 deg) with the target above
                    'min_alt', per night.

                - 'dark_hours' : np.ndarray[float]
                    Total dark hours per night (Sun below -18 deg; 0.0 if
                    the Sun never reaches -18 deg).

                - 'alt_midnight' : np.ndarray[float]
                    Altitude of the target at the local (true solar)
                    midnight of each night, in degrees.

                - 'max_alt' : np.ndarray[float]
                    Maximum altitude of the target during the dark time of
                    each night (Sun below -18 deg), in degrees (NaN if the
                    Sun never reaches -18 deg).

                - 'moon_illumination' : np.ndarray[float]
                    Illuminated fraction of the Moon at the local midnight
                    of each night, in percent (0-100).

                - 'best_night' : str
                    Night with the largest number of observable hours (the
                    first one in case of ties).

                - 'n_observable' : int
                    Number of nights with at least 1 observable hour.

                - 'fig' : Figure
                    The matplotlib figure, with two panels (altitude and
                    observable hours; see Notes).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the coordinates, the observatory or 'start_date' cannot be
            parsed, 'n_nights' is not positive, or 'min_alt' is not finite or
            lies outside [0, 90] deg.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. COMPUTATION
            - Each night is sampled every 5 min from mean local noon to the
              next mean local noon (from the observatory longitude). The
              twilight and 'min_alt' crossings are interpolated linearly
              between grid points, so that 'hours', 'dark_hours' and
              'max_alt' agree with 'visibility' and 'night_info' to a few
              seconds.

            - Dark time is the time with the Sun below -18 deg, as in
              'night_info', 'visibility', 'periodic_event' and
              'obs_phase_night'.

            - The local midnight is the true solar midnight (Sun at lower
              culmination), as in 'night_info' and 'visibility'. It is
              computed from the apparent hour angle of the Sun at mean local
              midnight; the two differ by the equation of time (< 17 min).

            - Typical use: choosing the semester and dates for a proposal,
              and estimating how many epochs of a spectroscopic binary can
              be secured in a campaign. Use 'visibility' for the detailed
              plan of a given night.

        4. PLOT
            - Upper panel: maximum altitude of the target during the dark
              time of each night (shaded) and altitude at local midnight
              (black line); the dotted red line is 'min_alt'.

            - Lower panel: dark hours with the target above 'min_alt'
              (shaded) and total dark time (black line).

            - Both panels show the Moon illumination at local midnight
              (dotted gray line, right axis) and share the date axis.

            - Each night covers one day on the date axis, from its starting
              date to the next one, so the axis runs from 'start_date' to
              'start_date' + 'n_nights' days (the morning after the last
              night). Both limits are always labeled, e.g. 2026-10-01 and
              2027-10-01 for 365 nights starting on 2026-10-01.

    ----------
    References
    ----------
        - Morris, B. M., Tollerud, E., Sipőcz, B. et al. 2018, AJ, 155, 128.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy.coordinates import AltAz, get_sun
    from astropy.time import Time
    import astropy.units as u
    from datetime import datetime, timedelta
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import numbers
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _above, _check_bool, _check_min_alt, _check_night, _coord, _get_observer, _observable, _solar_midnight


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coords_obs = _coord(coords, "coords")
    #
    observatory_obs = _get_observer(observatory)
    #
    start_str = _check_night(start_date, "start_date", "Start date")
    #
    if not isinstance(n_nights, numbers.Integral) or isinstance(n_nights, bool):
        raise TypeError(f"Number of nights 'n_nights' must be an integer. Got {type(n_nights).__name__}.")
    elif n_nights <= 0:
        raise ValueError(f"Number of nights 'n_nights' must be positive. Got {n_nights}.")
    n_nights = int(n_nights)  # Python int (also for NumPy integers)
    #
    min_alt = _check_min_alt(min_alt)
    #
    show = _check_bool(show, "show", "Display flag")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'coords'            : coords_obs,
                    'observatory'       : observatory_obs,
                    'nights'            : [],
                    'hours'             : np.nan,
                    'dark_hours'        : np.nan,
                    'alt_midnight'      : np.nan,
                    'max_alt'           : np.nan,
                    'moon_illumination' : np.nan,
                    'best_night'        : '',
                    'n_observable'      : np.nan,
                    'fig'               : None,
                   }


    # ####################  PARAMETERS  ####################
    # ##########  TIME GRID  ##########
    step        = 5.  # Minutes
    n_per_night = int(24 * 60 / step)  # Samples from local noon to the next local noon
    lon_hours   = observatory_obs.location.lon.wrap_at(180 * u.deg).deg / 15.  # Longitude in hours (East positive)

    local_noon = Time(start_str, scale="utc") + (12. - lon_hours) * u.hour  # Mean local noon of the first night (UTC)
    offsets    = (np.arange(n_nights)[:, None] * 24. * 60. + np.arange(n_per_night)[None, :] * step) * u.min  # Shape (n_nights, n_per_night)
    time_grid  = local_noon + offsets.ravel()  # Flattened for a single (vectorized) astropy call


    # ####################  COMPUTATION  ####################
    # ##########  ALTITUDES  ##########
    altaz_frame = AltAz(obstime=time_grid, location=observatory_obs.location)  # No pressure: geometric altitudes, as in astroplan with pressure=0
    sun_alt     = get_sun(time_grid).transform_to(altaz_frame).alt.deg.reshape(n_nights, n_per_night)  # One row per night
    target_alt  = coords_obs.transform_to(altaz_frame).alt.deg.reshape(n_nights, n_per_night)


    # ##########  HOURS PER NIGHT  ##########
    # Durations are interpolated linearly between grid points (as in 'visibility'), so that they
    # do not depend on the 5-min grid
    lo_dark, hi_dark = _above(-sun_alt, 18.)  # Sun below -18 deg (astronomical night)
    dark_length      = np.clip(hi_dark - lo_dark, 0., None)  # Dark fraction of each grid interval
    _, _, obs_length = _observable(sun_alt, target_alt, min_alt.to(u.deg).value)  # Dark and above 'min_alt'

    hours      = obs_length.sum(axis=1) * step / 60.  # Observable hours per night
    dark_hours = dark_length.sum(axis=1) * step / 60.  # Dark hours per night


    # ##########  MAXIMUM ALTITUDE IN DARK TIME  ##########
    # Target altitude at both ends of the dark part of every interval
    # (grid points or interpolated twilight crossings)
    alt_0, alt_1                   = target_alt[:, :-1], target_alt[:, 1:]  # Altitude at the start and end of each interval
    alt_edges                      = np.maximum(alt_0 + lo_dark * (alt_1 - alt_0), alt_0 + hi_dark * (alt_1 - alt_0))  # Linear in each interval: the maximum is at one end
    max_alt                        = np.where(dark_length > 0, alt_edges, -np.inf).max(axis=1)  # Only intervals with dark time
    max_alt[~np.isfinite(max_alt)] = np.nan  # Nights without dark time


    # ##########  LOCAL MIDNIGHT  ##########
    mean_midnights = local_noon + (np.arange(n_nights) * 24. + 12.) * u.hour  # Mean local midnights (UTC)
    midnights      = _solar_midnight(observatory_obs, mean_midnights)  # True solar midnights (Sun at lower culmination), as in the other functions

    altaz_midnight = AltAz(obstime=midnights, location=observatory_obs.location)
    alt_midnight   = np.atleast_1d(coords_obs.transform_to(altaz_midnight).alt.deg)  # One value per night


    # ##########  MOON  ##########
    moon_illum = 100. * np.atleast_1d(observatory_obs.moon_illumination(midnights)).astype(float)  # Percent


    # ##########  SUMMARY  ##########
    nights       = [(Time(start_str, scale="utc") + i * u.day).strftime("%Y-%m-%d") for i in range(n_nights)]
    best_night   = nights[int(np.argmax(hours))]  # First night with the maximum (argmax keeps the first on ties)
    n_observable = int(np.sum(hours >= 1.))  # Nights with at least 1 observable hour


    # ####################  PLOTTING  ####################
    # ##########  INIT FIGURE  ##########
    fig, (ax_alt, ax_hours) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    dates      = [datetime.strptime(n, "%Y-%m-%d") for n in nights]  # Start date of each night
    date_end   = dates[-1] + timedelta(days=1)  # Morning after the last night: end of the period
    edges      = dates + [date_end]  # Each night spans [date, date + 1 day)
    centers    = [d + timedelta(hours=12) for d in dates]  # Center of each night interval
    min_alt_lb = f"alt ≥ {min_alt.to(u.deg).value:.0f}°"  # Legend text


    # ##########  PANEL 1: ALTITUDE  ##########
    # ##### PLOT #####
    ax_alt.fill_between(edges, 0, np.append(max_alt, max_alt[-1]), color="steelblue", alpha=0.6, step="post", label="Max. altitude in dark time")  # Step per night; last value repeated to close the last step
    ax_alt.plot(centers, alt_midnight, color="black", lw=1.2, label="Altitude at midnight")
    ax_alt.axhline(min_alt.to(u.deg).value, color="red", ls=":", lw=1.2, label=f"Minimum altitude ({min_alt_lb})")

    moon_alt_ax = ax_alt.twinx()  # Right axis: Moon illumination [%]
    moon_alt_ax.plot(centers, moon_illum, color="gray", lw=0.8, ls=":", label="Moon illumination")
    moon_alt_ax.set_ylim(0, 105)
    moon_alt_ax.set_ylabel("Moon illumination [%]")


    # ##### AXIS #####
    ax_alt.set_ylim(0, 90)
    ax_alt.set_ylabel("Altitude [deg]")
    ax_alt.set_title(f"Annual visibility — {observatory_obs.name}", fontsize=14)
    ax_alt.grid(True, ls="--", alpha=0.4)
    handles_ax, labels_ax     = ax_alt.get_legend_handles_labels()
    handles_moon, labels_moon = moon_alt_ax.get_legend_handles_labels()
    ax_alt.legend(handles_ax + handles_moon, labels_ax + labels_moon, loc="upper right", fontsize=9)  # Single legend for both axes


    # ##########  PANEL 2: OBSERVABLE HOURS  ##########
    # ##### PLOT #####
    ax_hours.fill_between(edges, 0, np.append(hours, hours[-1]), color="steelblue", alpha=0.6, step="post", label=f"Observable ({min_alt_lb})")  # Step per night
    ax_hours.plot(centers, dark_hours, color="black", lw=1.2, label="Dark time")

    moon_hours_ax = ax_hours.twinx()
    moon_hours_ax.plot(centers, moon_illum, color="gray", lw=0.8, ls=":", label="Moon illumination")
    moon_hours_ax.set_ylim(0, 105)
    moon_hours_ax.set_ylabel("Moon illumination [%]")


    # ##### AXIS #####
    ax_hours.set_ylim(0, max(1., float(np.max(dark_hours))) * 1.1)  # 10% margin; at least 1 h if there is no dark time
    ax_hours.set_ylabel("Time per night [h]")
    ax_hours.grid(True, ls="--", alpha=0.4)
    handles_ax, labels_ax     = ax_hours.get_legend_handles_labels()
    handles_moon, labels_moon = moon_hours_ax.get_legend_handles_labels()
    ax_hours.legend(handles_ax + handles_moon, labels_ax + labels_moon, loc="upper right", fontsize=9)


    # ##########  DATE AXIS  ##########
    # Automatic date ticks, plus the first and last dates of the period (always labeled). Automatic
    # ticks closer than 6% of the axis to them are dropped so that labels do not overlap
    x_start, x_end = mdates.date2num(dates[0]), mdates.date2num(date_end)
    auto_ticks     = mdates.AutoDateLocator(minticks=4, maxticks=10).tick_values(dates[0], date_end)
    min_gap        = 0.06 * (x_end - x_start)
    ticks          = [x_start] + [t for t in auto_ticks if x_start + min_gap < t < x_end - min_gap] + [x_end]

    ax_hours.set_xlim(dates[0], date_end)
    ax_hours.set_xticks(ticks)
    ax_hours.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax_hours.set_xlabel("Date")
    fig.autofmt_xdate()  # Rotated date labels


    # ##########  SHOW PLOT  ##########
    fig.tight_layout()
    if show:
        plt.show()  # Inline backend: displays the figure; interactive backends: opens a window


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'nights'            : nights,
                         'hours'             : hours.astype(float),
                         'dark_hours'        : dark_hours.astype(float),
                         'alt_midnight'      : alt_midnight.astype(float),
                         'max_alt'           : max_alt.astype(float),
                         'moon_illumination' : moon_illum,
                         'best_night'        : best_night,
                         'n_observable'      : n_observable,
                         'fig'               : fig,
                        })
    return results_dict




# ==============================================================================
# ============================== PERIODIC OBJECTS ==============================
# ==============================================================================

# ==============================
def periodic_event(
                   coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                   period      : int | float,
                   observatory : str | dict,
                   hjd0        : int | float | None = None,
                   jd0         : int | float | None = None,
                   start_time  : str | datetime | Time | None = None,
                   n_events    : int = 10,
                   min_alt     : int | float | u.Quantity = 30,
                   phase       : int | float = 0.
    ) -> dict[str, SkyCoord | float | Observer | str | list[str | float | bool | int]]:
    """
    Predict the next periodic events of an object (eclipses, maxima, quadratures...) and check whether they are observable from a given observatory.

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the periodic object. See Notes.

        period : int | float
            Period in days. Must be positive.

        observatory : str | dict
            Observatory parameters. See Notes.

        hjd0 : int | float | None, optional
            Reference epoch of the periodic event (phase 0) as a Heliocentric
            Julian Date (UTC). Optional if 'jd0' is provided; if both are
            given, 'hjd0' is used. Default is None.

        jd0 : int | float | None, optional
            Reference epoch of the periodic event (phase 0) as a Julian Date
            (UTC). Optional if 'hjd0' is provided. Default is None.

        start_time : str | datetime | Time | None, optional
            Time from which events are searched (UTC), e.g. 'YYYY-MM-DD' or
            'YYYY-MM-DD HH:MM:SS'. Default is None (current UTC time).

        n_events : int, optional
            Number of upcoming events to compute. Must be positive. Default
            is 10.

        min_alt : int | float | u.Quantity, optional
            Minimum altitude for an event to be considered observable, in
            [0, 90] deg. If int or float, it is read in degrees. Default is
            30.

        phase : int | float, optional
            Phase of the event with respect to the reference epoch, in
            [0, 1). E.g. 0 for primary eclipses, 0.5 for secondary eclipses
            of circular orbits, 0.25 / 0.75 for quadratures. Default is 0.

    -------
    Returns
    -------
        results_dict : dict[str, SkyCoord | float | Observer | str | list[str | float | bool | int]]
            A dictionary containing the following key-value pairs:
                - 'coords' : SkyCoord
                    The input coordinates converted to a SkyCoord object.

                - 'period' : float
                    Period in days (input value).

                - 'hjd0' : float
                    Reference Heliocentric Julian Date (UTC).

                - 'jd0' : float
                    Reference Julian Date (UTC).

                - 'phase' : float
                    Phase of the events (input value).

                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'start_time' : str
                    Start of the search ('YYYY-MM-DD HH:MM:SS' UTC).

                - 'cycles' : list[int]
                    Cycle number E of each event, counted from the reference
                    epoch.

                - 'events' : list[str]
                    UTC times of the events ('YYYY-MM-DD HH:MM:SS').

                - 'events_hjd' : list[float]
                    Heliocentric Julian Dates (UTC) of the events.

                - 'observable' : list[bool]
                    True if the event happens during the astronomical night
                    (Sun below -18 deg) with the target above 'min_alt'.

                - 'altitude' : list[float]
                    Altitude of the target at each event, in degrees.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If any input is not finite or is outside its valid range, no
            reference epoch is given, or the coordinates, the observatory
            or 'start_time' cannot be parsed.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. EPHEMERIS
            - Events follow the linear ephemeris:
                    T(E) = T0 + (phase + E) * P
              in the time scale of the reference epoch.

            - If 'hjd0' is given, the ephemeris is linear in HJD and each
              event is converted to geocentric UTC by removing the
              heliocentric light-travel time at that date (up to ±8.3 min).
              If only 'jd0' is given, the ephemeris is applied directly in JD.

            - HJD values are in UTC (classical convention). The difference
              with BJD_TDB (~1 min) is negligible for planning purposes.

            - The uncertainty of the predicted times grows linearly with the
              number of cycles: σ_T(E) ≈ sqrt(σ_T0^2 + E^2 σ_P^2).

    ----------
    References
    ----------
        - Morris, B. M., Tollerud, E., Sipőcz, B. et al. 2018, AJ, 155, 128.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astroplan import FixedTarget
    from astropy.time import Time
    import astropy.units as u
    import numbers
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_epochs, _check_min_alt, _check_period, _check_time, _coord, _get_observer, _reference_epochs
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coords_obs = _coord(coords, "coords")
    target     = FixedTarget(coord=coords_obs)  # Target object of astroplan, for 'target_is_up'
    #
    period = _check_period(period)
    #
    observatory_obs = _get_observer(observatory)
    #
    hjd0, jd0, hjd_flag = _check_epochs(hjd0, jd0)  # hjd_flag = True: the ephemeris is defined in HJD
    #
    t_start = _check_time(start_time, "start_time", "Start time")  # None -> now
    #
    if not isinstance(n_events, numbers.Integral) or isinstance(n_events, bool):
        raise TypeError(f"Number of periodic events 'n_events' must be an integer. Got {type(n_events).__name__}.")
    elif n_events <= 0:
        raise ValueError(f"Number of periodic events 'n_events' must be positive. Got {n_events}.")
    n_events = int(n_events)
    #
    min_alt = _check_min_alt(min_alt)
    #
    phase = _check_real(phase, "phase", "Event phase", "dimensionless", (0., 1., "[)"))


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'coords'      : coords_obs,
                    'period'      : period,
                    'hjd0'        : np.nan,
                    'jd0'         : np.nan,
                    'phase'       : phase,
                    'observatory' : observatory_obs,
                    'start_time'  : str(t_start.strftime("%Y-%m-%d %H:%M:%S")),
                    'cycles'      : [],
                    'events'      : [],
                    'events_hjd'  : [],
                    'observable'  : [],
                    'altitude'    : [],
                   }


    # ####################  COMPUTATION  ####################
    # ##########  REFERENCE EPOCH  ##########
    location  = observatory_obs.location
    hjd0, jd0 = _reference_epochs(hjd0, jd0, hjd_flag, coords_obs, location)


    # ##########  EVENTS (EPHEMERIS TIME SCALE)  ##########
    if hjd_flag:
        t_ref  = hjd0 + phase * period  # Event epoch in HJD
        t_from = t_start.utc.jd + t_start.light_travel_time(coords_obs, "heliocentric", location=location).to_value(u.day)  # Start time in HJD
    else:
        t_ref  = jd0 + phase * period  # Event epoch in JD
        t_from = t_start.utc.jd

    first_cycle = int(np.ceil((t_from - t_ref) / period))  # First event after the start time
    cycles      = first_cycle + np.arange(n_events)  # Cycle numbers E of the next 'n_events' events
    events_ref  = t_ref + cycles * period  # T(E) = T0 + (phase + E) * P, in the time stamp of the ephemeris


    # ##########  EVENTS (UTC)  ##########
    if hjd_flag:  # Remove the heliocentric light-travel time at each event (two iterations converge to < 1 ms)
        events_jd = events_ref.copy()  # First guess: JD = HJD
        for _ in range(2):
            ltt_helio = Time(events_jd, format="jd", scale="utc").light_travel_time(coords_obs, "heliocentric", location=location)
            events_jd = events_ref - ltt_helio.to_value(u.day)  # JD = HJD - light-travel time (evaluated at the current JD)
        events_hjd = events_ref
    else:
        events_jd  = events_ref  # Ephemeris already in JD
        ltt_helio  = Time(events_jd, format="jd", scale="utc").light_travel_time(coords_obs, "heliocentric", location=location)
        events_hjd = events_jd + ltt_helio.to_value(u.day)  # HJD reported for reference
    event_times = Time(events_jd, format="jd", scale="utc")  # Geocentric UTC times of the events


    # ##########  OBSERVABILITY  ##########
    altitude   = observatory_obs.altaz(event_times, target).alt  # Quantity (deg), compared with 'min_alt' below
    is_up      = np.atleast_1d(observatory_obs.target_is_up(event_times, target))  # Target above the horizon
    is_high    = np.atleast_1d(altitude >= min_alt)  # Target above the minimum altitude
    is_night   = np.atleast_1d(observatory_obs.is_night(event_times, horizon=-18 * u.deg))  # Astronomical night
    observable = is_up & is_high & is_night  # Every condition must hold


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'hjd0'       : hjd0,
                         'jd0'        : jd0,
                         'cycles'     : [int(c) for c in cycles],
                         'events'     : [str(t.strftime("%Y-%m-%d %H:%M:%S")) for t in event_times],
                         'events_hjd' : [float(t) for t in events_hjd],
                         'observable' : [bool(o) for o in observable],
                         'altitude'   : [float(a) for a in np.atleast_1d(altitude.to(u.deg).value)],
                        })
    return results_dict




# ==============================
def obs_phase(
              coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
              period      : int | float,
              observatory : str | dict,
              times       : str | datetime | Time | list[str | datetime | Time],
              hjd0        : int | float | None = None,
              jd0         : int | float | None = None
    ) -> dict[str, SkyCoord | float | Observer | list]:
    """
    Compute the phase of a periodic object at one or several observation times.

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the periodic object. See Notes.

        period : int | float
            Period in days. Must be positive.

        observatory : str | dict
            Observatory parameters. See Notes.

        times : str | datetime | Time | list[str | datetime | Time]
            Observation time(s) in UTC, in any format accepted by
            astropy.time.Time (e.g. 'YYYY-MM-DD HH:MM:SS').

        hjd0 : int | float | None, optional
            Reference epoch (phase 0) as a Heliocentric Julian Date (UTC).
            Optional if 'jd0' is provided; if both are given, 'hjd0' is used.
            Default is None.

        jd0 : int | float | None, optional
            Reference epoch (phase 0) as a Julian Date (UTC). Optional if
            'hjd0' is provided. Default is None.

    -------
    Returns
    -------
        results_dict : dict[str, SkyCoord | float | Observer | list]
            A dictionary containing the following key-value pairs:
                - 'coords' : SkyCoord
                    The input coordinates converted to a SkyCoord object.

                - 'period' : float
                    Period in days (input value).

                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'hjd0' : float
                    Reference Heliocentric Julian Date (UTC).

                - 'jd0' : float
                    Reference Julian Date (UTC).

                - 'times' : list[str]
                    Observation times ('YYYY-MM-DD HH:MM:SS' UTC).

                - 'times_jd' : list[float]
                    Observation times as Julian Dates (UTC).

                - 'times_hjd' : list[float]
                    Observation times as Heliocentric Julian Dates (UTC).

                - 'phases' : list[float]
                    Phases in [0, 1).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If any input is not finite or is outside its valid range, no
            reference epoch is given, or the coordinates, the observatory
            or a time cannot be parsed.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. PHASES
            - Phases are computed as:
                    φ = ((t - T0) mod P) / P
              using HJD for both t and T0 if 'hjd0' is given, and JD
              otherwise.

            - HJD values are in UTC (classical convention). For timing
              analyses use BJD_TDB from 'bary_corr'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import astropy.units as u
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_epochs, _check_period, _check_times, _coord, _get_observer, _phases, _reference_epochs


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coords_obs = _coord(coords, "coords")
    #
    period = _check_period(period)
    #
    observatory_obs = _get_observer(observatory)
    #
    astropy_times = _check_times(times, "times", "Observation time")  # Array-valued Time (UTC)
    #
    hjd0, jd0, hjd_flag = _check_epochs(hjd0, jd0)  # hjd_flag = True: the ephemeris is defined in HJD


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'coords'      : coords_obs,
                    'period'      : period,
                    'observatory' : observatory_obs,
                    'hjd0'        : np.nan,
                    'jd0'         : np.nan,
                    'times'       : [],
                    'times_jd'    : [],
                    'times_hjd'   : [],
                    'phases'      : [],
                   }


    # ####################  COMPUTATION  ####################
    # ##########  REFERENCE EPOCH  ##########
    location  = observatory_obs.location
    hjd0, jd0 = _reference_epochs(hjd0, jd0, hjd_flag, coords_obs, location)


    # ##########  TIME CONVERSIONS  ##########
    ltt_helio = astropy_times.light_travel_time(coords_obs, "heliocentric", location=location)  # Heliocentric light-travel time at each observation
    times_jd  = np.atleast_1d(astropy_times.utc.jd)  # JD (UTC)
    times_hjd = times_jd + np.atleast_1d(ltt_helio.to_value(u.day))  # HJD (UTC)


    # ##########  PHASES  ##########
    phases = _phases(astropy_times, coords_obs, location, period, hjd0, jd0, hjd_flag)  # From HJD if 'hjd0' was given, from JD otherwise


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'hjd0'      : hjd0,
                         'jd0'       : jd0,
                         'times'     : [str(t.strftime("%Y-%m-%d %H:%M:%S")) for t in astropy_times],
                         'times_jd'  : times_jd.tolist(),
                         'times_hjd' : times_hjd.tolist(),
                         'phases'    : phases.tolist(),
                        })
    return results_dict




# ==============================
def obs_phase_night(
                    coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                    period      : int | float,
                    observatory : str | dict,
                    nights      : str | datetime | Time | list[str | datetime | Time],
                    hjd0        : int | float | None = None,
                    jd0         : int | float | None = None,
                    min_alt     : int | float | u.Quantity = 30,
                    show        : bool = True
    ) -> dict[str, SkyCoord | float | Observer | list | Figure]:
    """
    Compute and plot the phases of a periodic object covered on each night of a series of observing nights.

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the periodic object. See Notes.

        period : int | float
            Period in days. Must be positive.

        observatory : str | dict
            Observatory parameters. See Notes.

        nights : str | datetime | Time | list[str | datetime | Time]
            Night(s) of observation, identified by the date on which they
            start ('YYYY-MM-DD'): each date refers to the night starting on
            that date and ending the next day.

        hjd0 : int | float | None, optional
            Reference epoch (phase 0) as a Heliocentric Julian Date (UTC).
            Optional if 'jd0' is provided; if both are given, 'hjd0' is used.
            Default is None.

        jd0 : int | float | None, optional
            Reference epoch (phase 0) as a Julian Date (UTC). Optional if
            'hjd0' is provided. Default is None.

        min_alt : int | float | u.Quantity, optional
            Minimum altitude for the target to be considered observable, in
            [0, 90] deg. If int or float, it is read in degrees. Default is
            30.

        show : bool, optional
            If True, display the figure with plt.show(). Default is True.

    -------
    Returns
    -------
        results_dict : dict[str, SkyCoord | float | Observer | list | Figure]
            A dictionary containing the following key-value pairs:
                - 'coords' : SkyCoord
                    The input coordinates converted to a SkyCoord object.

                - 'period' : float
                    Period in days (input value).

                - 'hjd0' : float
                    Reference Heliocentric Julian Date (UTC).

                - 'jd0' : float
                    Reference Julian Date (UTC).

                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'nights' : list[str]
                    Nights of observation ('YYYY-MM-DD').

                - 'phases' : list[list[float | None]]
                    Phases [start, end] covered on each night while the
                    target is observable ([None, None] if not observable).
                    If start > end, the coverage wraps through phase 0.

                - 'full_coverage' : list[bool]
                    True if the observable window of the night is longer
                    than the period (every phase is covered).

                - 'times' : list[list[str | float]]
                    UTC times [start, end] of the observable window of each
                    night ('YYYY-MM-DD HH:MM:SS'; [NaN, NaN] if not
                    observable).

                - 'fig' : Figure
                    The matplotlib phase-coverage figure.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If any input is not finite or is outside its valid range, no
            reference epoch is given, or the coordinates, the observatory
            or a night cannot be parsed.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. OBSERVABLE WINDOW
            - Each night is sampled every 2 min (as in 'visibility')
              between the evening and morning astronomical twilights
              (Sun at -18 deg); the window is given by the times at
              which the target crosses 'min_alt', interpolated linearly
              between samples. This is the same definition of dark time
              and observable window used by 'night_info', 'visibility',
              'annual_visibility' and 'periodic_event'.

            - Where the Sun never reaches -18 deg (summer at high
              latitudes), a warning is issued and the night is reported as
              not observable. Where it stays below -18 deg during the whole
              day (polar night), the window covers ±12 h around midnight.

            - Phases are computed as in 'obs_phase' (HJD if 'hjd0' is given,
              JD otherwise).

        4. PLOT
            - Upper panel: phase range covered on each night, with the UTC
              start (green) and end (red) times of the window.

            - Lower panel: cumulative phase coverage of all the nights,
              useful to check whether a campaign samples the whole orbit
              (e.g. for SB2 disentangling or light-curve coverage).
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astroplan import FixedTarget
    from astropy.time import Time
    import astropy.units as u
    import matplotlib.pyplot as plt
    import numpy as np
    import warnings


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_bool, _check_epochs, _check_min_alt, _check_period, _check_times, _coord, _draw_phase_range, _get_observer, _is_event, _mean_midnight, _phases, _reference_epochs, _solar_midnight


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coords_obs = _coord(coords, "coords")
    target     = FixedTarget(coord=coords_obs)  # Target object of astroplan, for 'altaz'
    #
    period = _check_period(period)
    #
    observatory_obs = _get_observer(observatory)
    #
    nights_obs = [str(t.strftime("%Y-%m-%d")) for t in _check_times(nights, "nights", "Observing nights")]  # Dates on which the nights start
    #
    hjd0, jd0, hjd_flag = _check_epochs(hjd0, jd0)  # hjd_flag = True: the ephemeris is defined in HJD
    #
    min_alt = _check_min_alt(min_alt)
    #
    show = _check_bool(show, "show", "Display flag")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'coords'        : coords_obs,
                    'period'        : period,
                    'hjd0'          : np.nan,
                    'jd0'           : np.nan,
                    'observatory'   : observatory_obs,
                    'nights'        : nights_obs,
                    'phases'        : [],
                    'full_coverage' : [],
                    'times'         : [],
                    'fig'           : None,
                   }


    # ####################  COMPUTATION  ####################
    # ##########  REFERENCE EPOCH  ##########
    location  = observatory_obs.location
    hjd0, jd0 = _reference_epochs(hjd0, jd0, hjd_flag, coords_obs, location)


    # ##########  NIGHTS  ##########
    phases_per_night     = []  # [min_phase, max_phase] of each night (None if not observable)
    full_per_night       = []  # True if the window is longer than the period
    times_per_night      = []  # UTC start and end of each window, as strings (NaN if not observable)
    times_per_night_plot = []  # The same, as datetimes for the plot labels (None if not observable)
    min_alt_deg          = float(min_alt.to(u.deg).value)
    step                 = 2.  # Minutes (same grid step as 'visibility')

    for night in nights_obs:


        # ##### DARK WINDOW #####
        midnight = _solar_midnight(observatory_obs, _mean_midnight(observatory_obs, night))  # True solar midnight (Sun at lower culmination) of the night starting on 'night'

        night_start = observatory_obs.twilight_evening_astronomical(midnight, which="previous")  # Sun at -18 deg (evening)
        night_end   = observatory_obs.twilight_morning_astronomical(midnight, which="next")  # Sun at -18 deg (morning)
        polar_night = observatory_obs.sun_altaz(midnight).alt.deg < -18.  # Sun below -18 deg at midnight

        if not (_is_event(night_start) and _is_event(night_end)) and polar_night:  # Dark during the whole day: ±12 h around midnight
            night_start, night_end = midnight - 12 * u.hour, midnight + 12 * u.hour
        elif not (_is_event(night_start) and _is_event(night_end)):  # The Sun never reaches -18 deg
            warnings.warn(f"The Sun does not reach -18 deg on the night of {night}; the night is reported as not observable.", stacklevel=2)
            phases_per_night.append([None, None])  # Same placeholders as a night without observable window
            full_per_night.append(False)
            times_per_night_plot.append([None, None])
            times_per_night.append([np.nan, np.nan])
            continue

        n_grid    = max(int(np.ceil((night_end.jd - night_start.jd) * 1440. / step)) + 1, 2)  # Spacing <= 'step', both twilights included
        time_grid = Time(np.linspace(night_start.jd, night_end.jd, n_grid), format="jd", scale="utc")  # Dark window only
        altitudes = observatory_obs.altaz(time_grid, target).alt.deg  # Geometric altitudes, in degrees


        # ##### OBSERVABLE WINDOW #####
        # The 'min_alt' crossings are interpolated linearly between grid points (as in 'visibility')
        grid_jd = time_grid.utc.jd
        idx     = np.where(altitudes >= min_alt_deg)[0]  # Grid points above the minimum altitude

        if len(idx) > 0:
            i_first, i_last   = idx[0], idx[-1]  # First and last points above 'min_alt'
            jd_first, jd_last = grid_jd[i_first], grid_jd[i_last]  # Window limits before interpolation
            if i_first > 0:  # Rises above 'min_alt' within the dark window
                fraction = (min_alt_deg - altitudes[i_first - 1]) / (altitudes[i_first] - altitudes[i_first - 1])
                jd_first = grid_jd[i_first - 1] + fraction * (grid_jd[i_first] - grid_jd[i_first - 1])
            if i_last < len(grid_jd) - 1:  # Sets below 'min_alt' within the dark window
                fraction = (min_alt_deg - altitudes[i_last]) / (altitudes[i_last + 1] - altitudes[i_last])
                jd_last  = grid_jd[i_last] + fraction * (grid_jd[i_last + 1] - grid_jd[i_last])
            window = Time([jd_first, jd_last], format="jd", scale="utc")  # Start and end of the observable window

            min_phase, max_phase = (float(ph) for ph in _phases(window, coords_obs, location, period, hjd0, jd0, hjd_flag))  # Phases at both ends of the window
            times_min, times_max = window[0].to_datetime(), window[1].to_datetime()  # For the time labels of the plot
            full_coverage        = bool((jd_last - jd_first) >= period)  # Window longer than the period
        else:
            min_phase, max_phase = None, None
            times_min, times_max = None, None
            full_coverage        = False

        phases_per_night.append([min_phase, max_phase])
        full_per_night.append(full_coverage)
        times_per_night_plot.append([times_min, times_max])
        times_per_night.append([
                                times_min.strftime("%Y-%m-%d %H:%M:%S") if times_min is not None else np.nan,
                                times_max.strftime("%Y-%m-%d %H:%M:%S") if times_max is not None else np.nan
                               ])


    # ####################  PLOTTING  ####################
    # ##########  INIT FIGURE  ##########
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 8), gridspec_kw={'height_ratios' : [1, 0.05]}, sharex=True)  # Thin lower panel for the cumulative coverage
    ax1.grid(True, linestyle="--", alpha=0.5)

    bbox_props  = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.75)  # Semi-transparent box so that the times are always readable
    x_positions = np.round(np.arange(0, 1.01, 0.1), 1)  # Phase ticks 0.0, 0.1, ..., 1.0


    # ##########  PHASES PER NIGHT  ##########
    # ##### PLOT #####
    for i, (phase_range, time_range, full) in enumerate(zip(phases_per_night, times_per_night_plot, full_per_night)):
        init_phase, final_phase = phase_range  # Phases at the start and end of the window
        init_time, final_time   = time_range  # UTC times at the start and end of the window

        if init_phase is None or final_phase is None:  # Not observable on this night
            continue

        _draw_phase_range(ax1, i, init_phase, final_phase, full, 15)  # One bar per night (y = night index)

        init_ha        = "left" if init_phase < 0.9 else "right"  # Alignment that keeps the labels inside the plot
        final_ha       = "right" if final_phase > 0.1 else "left"
        v_offset_init  = 8  # Points
        v_offset_final = 20 if abs(final_phase - init_phase) < 0.08 else 8  # Stagger the end time if both labels are close

        ax1.annotate(
                     init_time.strftime("%H:%M"),
                     xy=(init_phase, i),
                     xytext=(0, v_offset_init),
                     textcoords="offset points",
                     fontsize=9,
                     color="green",
                     ha=init_ha,
                     va="bottom",
                     bbox=bbox_props
                    )
        ax1.annotate(
                     final_time.strftime("%H:%M"),
                     xy=(final_phase, i),
                     xytext=(0, v_offset_final),
                     textcoords="offset points",
                     fontsize=9,
                     color="red",
                     ha=final_ha,
                     va="bottom",
                     bbox=bbox_props
                    )


    # ##### TITLE #####
    ax1.set_title(f"Observed phases — {observatory_obs.name}", fontsize=16)


    # ##### X AXIS #####
    ax1.set_xlim(0, 1)
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels([f"{x:.1f}" for x in x_positions])


    # ##### Y AXIS #####
    ax1.set_ylabel("Observing night")
    ax1.set_yticks(np.arange(len(nights_obs)))  # One row per night (y = night index)
    ax1.set_yticklabels(nights_obs)
    ax1.set_ylim(-0.5, len(nights_obs) - 0.5)  # Half a row of margin above and below


    # ##########  CUMULATIVE PHASES  ##########
    # ##### PLOT #####
    for phase_range, full in zip(phases_per_night, full_per_night):
        init_phase, final_phase = phase_range
        if init_phase is not None and final_phase is not None:
            _draw_phase_range(ax2, 0., init_phase, final_phase, full, 30)  # Every night on the same bar (y = 0)


    # ##### X AXIS #####
    ax2.set_xlabel("Phase")
    ax2.set_xlim(0, 1)
    ax2.set_xticks(x_positions)
    ax2.set_xticklabels([f"{x:.1f}" for x in x_positions])


    # ##### Y AXIS #####
    ax2.set_ylim(-0.01, 0.01)  # Single bar at y = 0
    ax2.set_yticks([])


    # ##########  SHOW PLOT  ##########
    fig.tight_layout()
    if show:
        plt.show()


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'hjd0'          : hjd0,
                         'jd0'           : jd0,
                         'phases'        : phases_per_night,
                         'full_coverage' : full_per_night,
                         'times'         : times_per_night,
                         'fig'           : fig,
                        })
    return results_dict




# ==========================================================================
# ============================== OBSERVATIONS ==============================
# ==========================================================================

# ==============================
def telescope_pointing(
                       coords      : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
                       observatory : str | dict,
                       obstime     : str | datetime | Time | None = None,
                       show        : bool = True
    ) -> dict[str, SkyCoord | Time | Observer | float | str | Figure]:
    """
    Compute and plot the telescope pointing (azimuth, elevation, airmass, hour and parallactic angles) for a target.

    ----------
    Parameters
    ----------
        coords : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the target. See Notes.

        observatory : str | dict
            Observatory parameters. See Notes.

        obstime : str | datetime | Time | None, optional
            Observation time in UTC (e.g. '2026-06-10T22:30:00'). Default is
            None (current time).

        show : bool, optional
            If True, display the figure with plt.show(). Default is True.

    -------
    Returns
    -------
        results_dict : dict[str, SkyCoord | Time | Observer | float | str | Figure]
            A dictionary containing the following key-value pairs:
                - 'coords' : SkyCoord
                    The input coordinates converted to a SkyCoord object.

                - 'obstime' : Time
                    The observation time (UTC).

                - 'observatory' : Observer
                    The observatory object (location and timezone).

                - 'elevation' : float
                    Elevation (altitude) above the horizon, in degrees.

                - 'azimuth' : float
                    Azimuth in degrees (0 = North, 90 = East).

                - 'direction' : str
                    Cardinal direction of the pointing ('N', 'NE', ...).

                - 'airmass' : float
                    Airmass sec(z) (NaN below the horizon).

                - 'hour_angle' : float
                    Hour angle in hours, in [-12, 12) (negative: East of the
                    meridian, rising).

                - 'lst' : str
                    Local apparent sidereal time ('hh:mm:ss').

                - 'parallactic_angle' : float
                    Parallactic angle in degrees.

                - 'moon_sep' : float
                    Angular separation from the Moon, in degrees.

                - 'moon_illumination' : float
                    Illuminated fraction of the Moon, in percent (0-100).

                - 'fig' : Figure
                    The matplotlib pointing figure (compass and elevation
                    diagrams).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the coordinates, the observatory or 'obstime' cannot be
            parsed, or 'obstime' is not a single time.

    -----
    Notes
    -----
        1. ACCEPTED COORDINATE FORMATS
            - Coordinates are converted into an astropy SkyCoord object (ICRS
              unless another frame is given).

            - Single string with RA and Dec, separated by whitespace:
                >>> "05:35:17.3 -05:23:28"
                >>> "05h35m17.3s -05d23m28s"
                >>> "05 35 17.3 -05 23 28"   (space-separated sexagesimal, as copied from SIMBAD)
                >>> "83.82208 -5.39111"

            - Tuple or list (RA, Dec) of strings, numbers or astropy
              quantities; strings accept the same notations as above:
                >>> ("05:35:17.3", "-05:23:28")
                >>> ("05h35m17.3s", "-05d23m28s")
                >>> ("05 35 17.3", "-05 23 28")
                >>> ("83.82208", "-5.39111")
                >>> [83.82208, -5.39111]
                >>> (5.5881389 * u.hourangle, -5.39111 * u.deg)

            - Dictionary with 'ra' (or 'right_ascension') and 'dec' (or
              'declination') keys, and optional 'unit' and 'frame' keys;
              values accept the same notations as tuples:
                >>> {'ra': "05:35:17.3", 'dec': "-05:23:28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object:
                >>> SkyCoord(83.82208 * u.deg, -5.39111 * u.deg, frame="icrs")

            - RA values containing ':' or 'h' are read in hours; plain
              numbers, decimal strings and 'd'-style strings are read in
              degrees, unless explicit units are given. Dec is always read
              in degrees.

        2. ACCEPTED OBSERVATORY FORMATS
            - Pre-defined codes (case-insensitive):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

            - Custom dictionary ('name' is optional):
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive
                     'latitude': 28.30000000,       # deg, North positive
                     'elevation': 2390.0,           # m above sea level
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

        3. HOUR ANGLE
            - HA = LST - RA, with the local apparent sidereal time and the
              apparent right ascension of date (true equator and equinox,
              TETE frame), so that precession and nutation since J2000 are
              taken into account.

        4. PARALLACTIC ANGLE
            - Atmospheric differential refraction disperses the image along
              the vertical. With a slit spectrograph, align the slit with the
              parallactic angle at large airmass (≳ 1.3) to avoid
              wavelength-dependent slit losses (Filippenko 1982). It is less
              critical for fiber-fed spectrographs, but still affects the
              blue end of the spectrum for small fibers.

        5. AIRMASS
            - Computed as sec(z), accurate to ~1% for z < 70 deg.

    ----------
    References
    ----------
        - Filippenko, A. V. 1982, PASP, 94, 715.
        - Morris, B. M., Tollerud, E., Sipőcz, B. et al. 2018, AJ, 155, 128.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy.coordinates import AltAz, TETE
    import astropy.units as u
    from matplotlib.patches import Arc
    import matplotlib.pyplot as plt
    import numpy as np


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_bool, _check_time, _coord, _get_observer


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    coords_obs = _coord(coords, "coords")
    #
    observatory_obs = _get_observer(observatory)
    #
    obstime = _check_time(obstime, "obstime", "Observation time")  # Scalar Time (UTC); None -> now
    #
    show = _check_bool(show, "show", "Display flag")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'coords'            : coords_obs,
                    'obstime'           : obstime,
                    'observatory'       : observatory_obs,
                    'elevation'         : np.nan,
                    'azimuth'           : np.nan,
                    'direction'         : '',
                    'airmass'           : np.nan,
                    'hour_angle'        : np.nan,
                    'lst'               : '',
                    'parallactic_angle' : np.nan,
                    'moon_sep'          : np.nan,
                    'moon_illumination' : np.nan,
                    'fig'               : None,
                   }


    # ####################  PARAMETERS  ####################
    # ##########  ALT-AZ COORDINATES  ##########
    altaz_frame = AltAz(obstime=obstime, location=observatory_obs.location)  # Geometric (no refraction)
    pointing    = coords_obs.transform_to(altaz_frame)
    elevation   = float(pointing.alt.degree)  # Altitude above the horizon
    azimuth     = float(pointing.az.degree)  # From North through East
    airmass     = float(1. / np.sin(np.radians(elevation))) if elevation > 0 else np.nan  # sec(z)


    # ##########  CARDINAL DIRECTION  ##########
    directions         = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]  # 8-point compass rose, clockwise from North
    cardinal_direction = directions[int((azimuth + 22.5) // 45) % 8]  # 45 deg sectors centered on each direction


    # ##########  HOUR & PARALLACTIC ANGLES  ##########
    lst               = observatory_obs.local_sidereal_time(obstime)  # Apparent sidereal time
    ra_apparent       = coords_obs.transform_to(TETE(obstime=obstime)).ra  # RA of date, consistent with the apparent LST
    hour_angle        = float((lst - ra_apparent).wrap_at(12 * u.hourangle).hour)  # HA = LST - RA, in [-12, 12)
    parallactic_angle = float(observatory_obs.parallactic_angle(obstime, coords_obs).deg)  # Angle, at the target, between the directions to the zenith and to the north celestial pole


    # ##########  MOON  ##########
    moon_altaz = observatory_obs.moon_altaz(obstime)  # Moon in the same AltAz frame as the target
    moon_sep   = float(pointing.separation(moon_altaz).deg)
    moon_illum = float(observatory_obs.moon_illumination(obstime))  # Illuminated fraction (0-1)


    # ####################  PLOTTING  ####################
    # ##########  INIT FIGURE  ##########
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))  # Left: compass (azimuth); right: dome (elevation)


    # ##########  PLOT 1: COMPASS (AZIMUTH)  ##########
    # ##### CONFIG #####
    ax1 = axes[0]
    ax1.set_aspect("equal")  # The compass circle is not distorted
    ax1.axis("off")  # Schematic drawing: no axes


    # ##### COMPASS CIRCLE #####
    ax1.add_artist(plt.Circle((0, 0), 1, color="lightgray", fill=False, lw=2))  # Unit circle
    rad_az = np.deg2rad(azimuth)


    # ##### POINTING ARROW #####
    # 0 deg is North (up): x = sin, y = cos
    ax1.arrow(
              0,
              0,
              np.sin(rad_az) * 0.8,
              np.cos(rad_az) * 0.8,
              head_width=0.08,
              head_length=0.12,
              fc="blue",
              ec="blue",
              linewidth=2.5,
              length_includes_head=True
             )
    ax1.plot(0, 0, "ko", markersize=6)  # Center pivot


    # ##### CARDINAL LABELS #####
    for dir_label, angle in zip(["N 0°", "E 90°", "S 180°", "W 270°"], [0, 90, 180, 270]):
        rad_label = np.deg2rad(angle)
        x, y      = np.sin(rad_label) * 1.15, np.cos(rad_label) * 1.15  # Just outside the circle
        ha        = "center" if angle in [0, 180] else ("left" if angle == 90 else "right")  # Alignment that keeps each label outside the circle
        va        = "bottom" if angle == 0 else ("top" if angle == 180 else "center")
        ax1.text(x, y, dir_label, ha=ha, va=va, fontsize=11, fontweight="bold")


    # ##### TICKS #####
    # Every 10 deg (longer every 30 and 90 deg)
    for angle in np.arange(0, 360, 10):
        rad    = np.deg2rad(angle)
        inner  = 0.94 if angle % 30 == 0 else 0.97  # Inner end of the tick (radius)
        outer  = 1.06 if angle % 90 == 0 else (1.03 if angle % 30 == 0 else 1.01)  # Outer end of the tick (radius)
        lw     = 2.0 if angle % 90 == 0 else (1.2 if angle % 30 == 0 else 0.6)  # Thicker ticks at the cardinal points
        x0, y0 = np.sin(rad) * inner, np.cos(rad) * inner
        x1, y1 = np.sin(rad) * outer, np.cos(rad) * outer
        ax1.plot([x0, x1], [y0, y1], color="black", lw=lw, solid_capstyle="round")


    # ##### TITLE #####
    ax1.set_title(f"Telescope pointing: {cardinal_direction} ({azimuth:.1f}°)", fontsize=14, pad=15)


    # ##### AXIS #####
    ax1.set_xlim(-1.2, 1.2)
    ax1.set_ylim(-1.2, 1.2)


    # ##########  PLOT 2: DOME (ELEVATION)  ##########
    # ##### CONFIG #####
    ax2 = axes[1]
    ax2.set_aspect("equal")  # The dome arc is not distorted
    ax2.axis("off")  # Schematic drawing: no axes


    # ##### ELEVATION COLOR #####
    # Color code: below the horizon / low / good elevation
    is_visible = elevation >= 0
    if elevation < 0:
        el_color = "red"
    elif elevation <= 30:
        el_color = "orange"
    else:
        el_color = "green"
    vis_text = '' if is_visible else " (below horizon)"  # Appended to the title
    y_min    = -0.2 if is_visible else -1.2  # Show the ground below the horizon if the arrow points there


    # ##### HORIZON & GROUND #####
    ax2.hlines(0, -1.2, 1.2, color="saddlebrown", lw=2)  # Horizon
    ax2.fill_between([-1.2, 1.2], y_min, 0, color="sandybrown", alpha=0.3)  # Ground


    # ##### DOME ARC #####
    ax2.add_patch(Arc((0, 0), 2, 2, theta1=0, theta2=180, color="gray", lw=1.5))  # Upper half of the unit circle


    # ##### ELEVATION ARROW #####
    # 0 deg is the horizon: x = cos, y = sin
    rad_el = np.deg2rad(elevation)
    ax2.arrow(
              0,
              0,
              np.cos(rad_el) * 0.9,
              np.sin(rad_el) * 0.9,
              head_width=0.06,
              head_length=0.1,
              fc=el_color,
              ec=el_color,
              linewidth=3,
              length_includes_head=True
             )
    ax2.plot(0, 0, "ko", markersize=6)  # Center pivot


    # ##### ELEVATION REFERENCES #####
    for el_angle in [30, 60, 90]:
        rad_ref      = np.deg2rad(el_angle)
        x_ref, y_ref = np.cos(rad_ref), np.sin(rad_ref)  # Point of the dome at this elevation
        ax2.plot([0, x_ref], [0, y_ref], color="gray", lw=1, linestyle=":")  # Reference line from the center
        ax2.text(x_ref * 1.1, y_ref * 1.1, f"{el_angle}°", ha="center", va="center", fontsize=10, color="dimgray")


    # ##### SUMMARY BOX #####
    summary = (
               f"Airmass: {airmass:.2f}\n"
               f"HA: {hour_angle:+.2f} h\n"
               f"Parallactic angle: {parallactic_angle:.1f}°\n"
               f"Moon: {moon_sep:.0f}° ({100 * moon_illum:.0f}%)"
              )
    ax2.text(-1.15, 1.15, summary, ha="left", va="top", fontsize=10, family="monospace")  # Top-left corner, data coordinates


    # ##### TITLE #####
    ax2.set_title(f"Telescope elevation: {elevation:.1f}°{vis_text}", fontsize=14, color=el_color, pad=15)


    # ##### AXIS #####
    ax2.set_xlim(-1.2, 1.2)
    ax2.set_ylim(y_min, 1.2)


    # ##########  SHOW PLOT  ##########
    fig.tight_layout()
    if show:
        plt.show()  # Inline backend: displays the figure; interactive backends: opens a window


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'elevation'         : elevation,
                         'azimuth'           : azimuth,
                         'direction'         : cardinal_direction,
                         'airmass'           : airmass,
                         'hour_angle'        : hour_angle,
                         'lst'               : str(lst.to_string(unit=u.hourangle, sep=":", precision=0, pad=True)),  # 'hh:mm:ss'
                         'parallactic_angle' : parallactic_angle,
                         'moon_sep'          : moon_sep,
                         'moon_illumination' : moon_illum * 100,
                         'fig'               : fig,
                        })
    return results_dict




# ==============================
def binaralt(
             block : bool = False
    ) -> None:
    """
    Launch the BinarAlt browser interface for visibility analysis, weather monitoring, lunar constraints and observation planning.

    ----------
    Parameters
    ----------
        block : bool, optional
            If True, keep the server in the foreground until interrupted
            (Ctrl+C). If False, run it in a background thread so that the
            Python session can be used meanwhile. Default is False.

    -------
    Returns
    -------
        None
            The interface is opened in the default web browser.

    ------
    Raises
    ------
        TypeError
            If 'block' is not a boolean.

        FileNotFoundError
            If the interface file (binaralt.html) is missing from the
            installation.

    -----
    Notes
    -----
        1. INTERFACE
            - BinarAlt is an interactive, browser-based planning tool for the
              preparation and execution of observing nights: target
              visibility, ephemerides, lunar information, live observatory
              weather, exposure scheduling, an interactive exposure timeline
              and an observation-plan table and logbook.

            - Targets can be imported from and exported to CSV files with
              the columns 'Name', 'RA', 'Dec', 'Period (d)', 'HJD0', 'JD0',
              'Nº exposures', 'Exp. (s)', 'Own overhead (s)',
              'Own min. alt. (°)' and 'Priority (1-5)'.

        2. USER GUIDE
            - The user guide, BinarAlt_UserGuide.pdf, is distributed in the
              'READMEs/' folder of the repository. A copy is installed with
              the package next to the interface ('panorama/support/obs/'),
              from where it is opened with the question-mark icon of the
              interface.

        3. SERVER
            - The interface (binaralt.html) is served on 127.0.0.1, on a
              free port chosen by the operating system, by a lightweight HTTP
              server, and opened in the default web browser.

            - The astronomical calculations run in the browser (Astronomy
              Engine library, bundled in the page). The local server only
              relays the live weather requests to the external observatory
              services, which do not allow direct (cross-origin) requests
              from the page. The weather panel needs internet access.

            - With block=False the server runs in a daemon thread and stops
              when the Python session ends.

        4. REQUIREMENTS
            - The interface requires a modern web browser (e.g. Chrome,
              Edge, Firefox or Safari) with JavaScript enabled. It can be
              launched from a script, the Python console or a local Jupyter
              notebook, but it opens in a separate browser tab, not inline.

            - The server listens on 127.0.0.1, so the interface is only
              reachable from the same machine: it does not work from a
              remote Jupyter server (JupyterHub, Colab, SSH) unless the port
              is forwarded.

            - In a script run with 'python script.py', use block=True: with
              block=False the server stops when the script ends.

            - Without internet access the astronomical calculations still
              work (Astronomy Engine is bundled in the HTML), but the
              weather panel is not available and system fonts are used.

    ----------
    References
    ----------
        - Cross, D., Astronomy Engine (software),
          https://github.com/cosinekitty/astronomy.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import functools
    import http.server
    from pathlib import Path
    import threading
    import webbrowser


    # ##### PANORAMA #####
    from .support.obs.obs_support import _BinarAltRequestHandler, _check_bool, _WeatherService


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    block = _check_bool(block, "block", "Blocking-mode flag")


    # ####################  PARAMETERS  ####################
    # ##########  FILES  ##########
    html_dir  = Path(__file__).parent / "support" / "obs"  # Installed with the package (package-data)
    html_file = html_dir / "binaralt.html"

    if not html_file.is_file():
        raise FileNotFoundError(f"BinarAlt interface file not found: '{html_file}'.")


    # ##########  HTTP SERVER  ##########
    _BinarAltRequestHandler.weather_service = _WeatherService()  # Shared by every request (class attribute)
    handler                                 = functools.partial(_BinarAltRequestHandler, directory=str(html_dir))  # Static files served from 'html_dir'
    server                                  = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)  # Port 0: the OS assigns a free port
    url                                     = f"http://127.0.0.1:{server.server_address[1]}/{html_file.name}"  # Port assigned by the OS


    # ####################  OUTPUT  ####################
    # ##########  OPTION 1: FOREGROUND  ##########
    if block:
        print(f"\tBinarAlt: serving at {url} (Ctrl+C to stop)")
        webbrowser.open(url)  # Default web browser
        try:
            server.serve_forever()  # Runs until interrupted by the user
        except KeyboardInterrupt:
            print("\n\tBinarAlt: shutting down server...")
        finally:
            server.server_close()  # Release the port


    # ##########  OPTION 2: BACKGROUND  ##########
    else:
        print(f"\tBinarAlt: serving at {url} (stops with the Python session)")
        threading.Thread(target=server.serve_forever, daemon=True).start()  # Daemon thread: stops with the Python session
        webbrowser.open(url)
