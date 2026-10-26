"""
Internal helpers of the 'obs' module: coordinate and observatory parsing,
validation of the arguments shared by the public functions, astroplan event
checks and formatting, night, epoch and phase computations, drawing of phase
ranges, FITS header utilities, and the local server of the BinarAlt
interface.

Conventions shared by every function in this module:
    - Functions and classes are private (leading underscore); they are used
      by the public functions of 'panorama.obs' and are not exposed at the
      package level.
    - Coordinates are returned as astropy SkyCoord objects and observatories
      as astroplan Observer objects.
    - Angles are in degrees, elevations in meters and times in UTC.
    - Observing nights are identified by the date on which they start
      ('YYYY-MM-DD').
    - Validation helpers receive the name of the argument of the calling
      function ('name'), so that their error messages refer to that argument
      and the same argument gives the same message in every function
      ('_check_bool' is also reused by 'panorama.prismas').
    - All functions and class constructors validate their inputs before
      computing; nested functions and the HTTP callbacks of the BinarAlt
      server do not raise (they report errors in the JSON response).
"""

from __future__ import annotations

# ==================================================================
# ============================== INIT ==============================
# ==================================================================

# ####################  MODULES & PACKAGES  ####################
"""
- http.server is imported at module level because '_BinarAltRequestHandler'
  inherits from http.server.SimpleHTTPRequestHandler, which must exist when the
  class is defined (it cannot be imported inside a function)
"""
import http.server


# ####################  PRE-DEFINED OBSERVATORIES  ####################
"""
- Geodetic coordinates (longitude East positive and latitude North positive, in
  degrees; elevation in meters above sea level) and IANA timezone of the
  observatories that can be selected by code in every 'obs' function
"""
_OBSERVATORIES = {
                  'OT'       : {
                                'name'      : "Observatorio del Teide (OT)",
                                'longitude' : -16.50972222,
                                'latitude'  : 28.30000000,
                                'elevation' : 2390.0,
                                'timezone'  : "Atlantic/Canary",
                               },
                  'ORM'      : {
                                'name'      : "Observatorio del Roque de los Muchachos (ORM)",
                                'longitude' : -17.89166700,
                                'latitude'  : 28.75611100,
                                'elevation' : 2396.0,
                                'timezone'  : "Atlantic/Canary",
                               },
                  'CAHA'     : {
                                'name'      : "Calar Alto Observatory (CAHA)",
                                'longitude' : -2.54611111,
                                'latitude'  : 37.22361111,
                                'elevation' : 2168.0,
                                'timezone'  : "Europe/Madrid",
                               },
                  'OHP'      : {
                                'name'      : "Observatoire de Haute-Provence (OHP)",
                                'longitude' : 5.71333333,
                                'latitude'  : 43.93083333,
                                'elevation' : 650.0,
                                'timezone'  : "Europe/Paris",
                               },
                  'LASILLA'  : {
                                'name'      : "La Silla Observatory (ESO)",
                                'longitude' : -70.73000000,
                                'latitude'  : -29.25666667,
                                'elevation' : 2400.0,
                                'timezone'  : "America/Santiago",
                               },
                  'PARANAL'  : {
                                'name'      : "Paranal Observatory (ESO)",
                                'longitude' : -70.40416667,
                                'latitude'  : -24.62722222,
                                'elevation' : 2635.0,
                                'timezone'  : "America/Santiago",
                               },
                  'SAAO'     : {
                                'name'      : "South African Astronomical Observatory (SAAO)",
                                'longitude' : 20.81055556,
                                'latitude'  : -32.37944444,
                                'elevation' : 1798.0,
                                'timezone'  : "Africa/Johannesburg",
                               },
                  'MAUNAKEA' : {
                                'name'      : "Maunakea Observatories",
                                'longitude' : -155.47166667,
                                'latitude'  : 19.82666667,
                                'elevation' : 4205.0,
                                'timezone'  : "Pacific/Honolulu",
                               },
                 }


# ####################  INSTRUMENTS  ####################
"""
- Instruments supported by 'panorama.obs.fits_catalog': file pattern of the
  reduced spectra, observatory code, UTC hour at which the observing night
  rolls over to the next calendar day, and FITS keywords of each catalog column
  (aliases tried in order, the first match wins)
"""
_INSTRUMENTS = {
                'HERMES' : {
                            'pattern'          : "*_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fit*",
                            'observatory'      : "ORM",  # Mercator telescope site
                            'night_split_hour' : 9,
                            'header_keywords'  : {
                                                  'telescope'  : ["TELESCOP"],
                                                  'instrument' : ["INSTRUME"],
                                                  'object'     : ["OBJECT"],
                                                  'filter'     : ["FILTER"],
                                                  'exptime_s'  : ["EXPTIME"],
                                                  'date_obs'   : ["DATE-OBS"],
                                                  'ra'         : ["RA"],
                                                  'dec'        : ["DEC"],
                                                  'airmass'    : ["AIRMASS"],
                                                  'altitude'   : ["TELALT"],
                                                 },
                           },
               }




# =========================================================================
# ============================== COORDINATES ==============================
# =========================================================================

# ==============================
def _coord(
           coord_input : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord,
           name        : str = "coords"
    ) -> SkyCoord:
    """
    Convert any of the supported coordinate formats into an astropy SkyCoord object.

    ----------
    Parameters
    ----------
        coord_input : str | tuple[str | float, str | float] | list[str | float] | dict | SkyCoord
            Coordinates of the target. See Notes.

        name : str, optional
            Name of the argument of the calling function, used in the error
            messages. Default is "coords".

    -------
    Returns
    -------
        coords : SkyCoord
            The parsed coordinates (a SkyCoord input is returned unchanged).

    ------
    Raises
    ------
        TypeError
            If 'coord_input' has an unsupported type (including None), or
            'name' is not a string.

        ValueError
            If 'coord_input' is empty or cannot be parsed into valid
            coordinates, or 'name' is empty.

    -----
    Notes
    -----
        1. ACCEPTED FORMATS
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
              'declination') keys (case-insensitive), and optional 'unit'
              (overrides the unit inference) and 'frame' (default 'icrs')
              keys; values accept the same notations as tuples:
                >>> {'ra': "05 35 17.3", 'dec': "-05 23 28"}
                >>> {'ra': 83.82208, 'dec': -5.39111, 'unit': "deg", 'frame': "icrs"}

            - SkyCoord object (returned unchanged).

        2. UNIT INFERENCE
            - The units of RA and Dec are inferred independently: RA is read
              in hours (hourangle) if it contains ':' or 'h', and in degrees
              otherwise (also for 'd'-style strings such as '83d49m19s');
              Dec is always read in degrees.

            - Plain numbers and decimal strings are therefore read as
              degrees. Mixed inputs such as ("83.82208", "-05:23:28") are
              parsed correctly.

            - Space-separated components are joined with ':' only when all
              of them are numbers, so strings with units (e.g. '83.8 deg')
              are passed unchanged to SkyCoord.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(name, str):
        raise TypeError(f"Error-message text 'name' must be a string. Got {type(name).__name__}.")
    elif not name.strip():
        raise ValueError(f"Error-message text 'name' must be a non-empty string. Got '{name}'.")
    #
    if not isinstance(coord_input, (str, tuple, list, dict, SkyCoord)):
        raise TypeError(f"Coordinates '{name}' must be a string, a (RA, Dec) tuple or list, a dictionary or an astropy SkyCoord. Got {type(coord_input).__name__}.")
    elif isinstance(coord_input, str) and not coord_input.strip():
        raise ValueError(f"Coordinates '{name}' cannot be an empty string.")
    elif isinstance(coord_input, (tuple, list, dict)) and len(coord_input) == 0:
        raise ValueError(f"Coordinates '{name}' cannot be an empty {type(coord_input).__name__}.")


    # ##########  FUNCTIONS  ##########
    def _ra_unit(
                 ra_val : str | float | Quantity
        ) -> UnitBase:
        """
        Infer the unit of a right ascension value.

        ----------
        Parameters
        ----------
            ra_val : str | float | Quantity
                Right ascension as given by the user.

        -------
        Returns
        -------
            unit : UnitBase
                u.hourangle for ':' or 'h' strings, u.deg otherwise (also for
                explicit 'd' strings such as '83d49m19s').
        """
        ra_str = str(ra_val).strip().lower()  # Works for strings, numbers and quantities
        if "d" in ra_str:  # Explicit degrees (e.g. '83d49m19s')
            return u.deg
        elif ":" in ra_str or "h" in ra_str:  # Sexagesimal hours
            return u.hourangle
        return u.deg  # Plain numbers and decimal strings
    #
    def _normalize(
                   value : object
        ) -> object:
        """
        Join space-separated sexagesimal components with ':' (e.g. '05 35 17.3' -> '05:35:17.3').

        ----------
        Parameters
        ----------
            value : object
                RA or Dec value as given by the user.

        -------
        Returns
        -------
            value : object
                Normalized string, or the input unchanged if it is not a
                string made of 2 or 3 numeric space-separated components.
        """
        if isinstance(value, str):
            parts = value.strip().split()  # E.g. '05 35 17.3' -> ['05', '35', '17.3']
            if len(parts) in (2, 3) and all(re.fullmatch(r"[+-]?\d+(\.\d*)?", p) for p in parts):  # 'hh mm' or 'hh mm ss.s'
                return ":".join(parts)
            return value.strip()  # Other strings (e.g. '05:35:17.3', '83.8 deg') unchanged
        return value  # Numbers and quantities unchanged
    #
    def _build(
               ra    : object,
               dec   : object,
               frame : str = "icrs",
               unit  : object | None = None
        ) -> SkyCoord:
        """
        Create the SkyCoord object, raising a ValueError with context if it fails.

        ----------
        Parameters
        ----------
            ra : object
                Right ascension (string, number or astropy quantity).

            dec : object
                Declination (string, number or astropy quantity).

            frame : str, optional
                Coordinate frame. Default is "icrs".

            unit : object | None, optional
                Explicit units passed to SkyCoord; None infers them (RA with
                '_ra_unit', Dec in degrees). Default is None.

        -------
        Returns
        -------
            coords : SkyCoord
                The parsed coordinates.

        ------
        Raises
        ------
            ValueError
                If SkyCoord cannot parse the coordinates (the original error
                is chained).
        """
        ra, dec = _normalize(ra), _normalize(dec)  # 'hh mm ss' -> 'hh:mm:ss'
        try:
            if unit is not None:  # Explicit units
                return SkyCoord(ra, dec, unit=unit, frame=frame)
            elif hasattr(ra, "unit") and hasattr(dec, "unit"):  # Astropy quantities
                return SkyCoord(ra, dec, frame=frame)
            return SkyCoord(ra, dec, unit=(_ra_unit(ra), u.deg), frame=frame)  # Inferred units: RA in hours or degrees, Dec in degrees
        except Exception as e:  # Any parsing error of astropy, reported with the argument name
            raise ValueError(f"Coordinates '{name}' could not be parsed (RA = {ra!r}, Dec = {dec!r}): {e}") from e


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: SKYCOORD  ##########
    if isinstance(coord_input, SkyCoord):
        coords = coord_input  # Already parsed (any frame)


    # ##########  OPTION 2: STRING  ##########
    elif isinstance(coord_input, str):
        parts = coord_input.strip().split()  # Whitespace-separated components
        if len(parts) == 2:  # 'RA Dec'
            ra, dec = parts
        elif len(parts) == 6:  # 'hh mm ss.s ±dd mm ss' (space-separated sexagesimal)
            ra, dec = ":".join(parts[:3]), ":".join(parts[3:])
        else:
            raise ValueError(f"Coordinates '{name}' must contain RA and Dec ('RA Dec' or 'hh mm ss ±dd mm ss'). Got '{coord_input}'.")
        coords = _build(ra, dec)


    # ##########  OPTION 3: DICTIONARY  ##########
    elif isinstance(coord_input, dict):
        norm_dict = {str(k).strip().lower(): v for k, v in coord_input.items()}  # Case-insensitive keys
        ra        = norm_dict.get('ra', norm_dict.get('right_ascension'))
        dec       = norm_dict.get('dec', norm_dict.get('declination'))

        if ra is None or dec is None:
            raise ValueError(f"Coordinates '{name}' given as a dictionary must contain 'ra' (or 'right_ascension') and 'dec' (or 'declination') keys. Got {list(coord_input.keys())}.")
        coords = _build(ra, dec, frame=norm_dict.get('frame', "icrs"), unit=norm_dict.get('unit'))


    # ##########  OPTION 4: SEQUENCE (TUPLE / LIST)  ##########
    else:
        if len(coord_input) != 2:
            raise ValueError(f"Coordinates '{name}' given as a tuple or list must have exactly 2 elements (RA, Dec). Got {len(coord_input)}.")
        ra, dec = coord_input  # (RA, Dec) in this order
        coords  = _build(ra, dec)


    # ####################  OUTPUT  ####################
    return coords




# ===========================================================================
# ============================== OBSERVATORIES ==============================
# ===========================================================================

# ==============================
def _get_observer(
                  observatory : str | dict
    ) -> Observer:
    """
    Create an astroplan Observer object from a pre-defined observatory code or a dictionary of parameters.

    ----------
    Parameters
    ----------
        observatory : str | dict
            Observatory code or parameters. See Notes.

    -------
    Returns
    -------
        observer : Observer
            Observer object with the location, timezone and name of the
            observatory.

    ------
    Raises
    ------
        TypeError
            If 'observatory' or any of its values has the wrong type.

        ValueError
            If the code is unknown, a required key is missing, or a value is
            not finite, out of range or not a valid timezone.

    -----
    Notes
    -----
        1. PRE-DEFINED OBSERVATORIES
            - Codes are case-insensitive, and spaces, '-' and '_' are
              ignored (e.g. 'La Silla' -> 'LASILLA'):
                >>> "OT"       : Observatorio del Teide (Tenerife, Spain)
                >>> "ORM"      : Observatorio del Roque de los Muchachos (La Palma, Spain)
                >>> "CAHA"     : Calar Alto Observatory (Almería, Spain)
                >>> "OHP"      : Observatoire de Haute-Provence (France)
                >>> "LASILLA"  : La Silla Observatory (ESO, Chile)
                >>> "PARANAL"  : Paranal Observatory (ESO, Chile)
                >>> "SAAO"     : South African Astronomical Observatory (Sutherland)
                >>> "MAUNAKEA" : Maunakea Observatories (Hawai'i, USA)

        2. CUSTOM DICTIONARY
            - Required keys 'longitude', 'latitude', 'elevation' and
              'timezone'; 'name' is optional:
                >>> {
                     'name': "My Observatory",
                     'longitude': -16.50972222,     # deg, East positive, in [-180, 180]
                     'latitude': 28.30000000,       # deg, North positive, in [-90, 90]
                     'elevation': 2390.0,           # m above sea level, non-negative
                     'timezone': "Atlantic/Canary"  # IANA timezone name
                    }

            - The input dictionary is not modified.

        3. REFRACTION
            - The atmospheric pressure is set to 0, so altitudes are
              geometric (no refraction). Refraction raises objects by
              ~1.7 arcmin at 30 deg and ~35 arcmin at the horizon, which is
              negligible for planning.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astroplan import Observer
    import astropy.units as u
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    # ##### OPTION 1: OBSERVATORY CODE #####
    if isinstance(observatory, str):
        obs_key = observatory.strip().upper().replace(" ", '').replace("-", '').replace("_", '')  # E.g. 'La Silla' -> 'LASILLA'
        if obs_key not in _OBSERVATORIES:
            raise ValueError(f"Observatory code 'observatory' is unknown. Available codes: {list(_OBSERVATORIES.keys())}. Got '{observatory}'.")
        obs_data = dict(_OBSERVATORIES[obs_key])  # Copy: the pre-defined table is never modified


    # ##### OPTION 2: OBSERVATORY DICTIONARY #####
    elif isinstance(observatory, dict):
        obs_data = dict(observatory)  # Copy: never modify the user's dictionary
        #
        missing_keys = [k for k in ('longitude', 'latitude', 'elevation', 'timezone') if obs_data.get(k) is None]
        if missing_keys:
            raise ValueError(f"Observatory dictionary 'observatory' must contain the keys 'longitude', 'latitude', 'elevation' and 'timezone'. Missing: {missing_keys}.")
        #
        for key in ('longitude', 'latitude', 'elevation'):
            value = obs_data[key]
            if not isinstance(value, numbers.Real) or isinstance(value, bool):
                raise TypeError(f"Observatory {key} 'observatory['{key}']' must be a real number. Got {type(value).__name__}.")
            elif not np.isfinite(value):
                raise ValueError(f"Observatory {key} 'observatory['{key}']' must be finite. Got {value}.")
        #
        if not -180. <= obs_data['longitude'] <= 180.:
            raise ValueError(f"Observatory longitude 'observatory['longitude']' must lie in [-180, 180] deg (East positive). Got {obs_data['longitude']}.")
        elif not -90. <= obs_data['latitude'] <= 90.:
            raise ValueError(f"Observatory latitude 'observatory['latitude']' must lie in [-90, 90] deg (North positive). Got {obs_data['latitude']}.")
        elif obs_data['elevation'] < 0:
            raise ValueError(f"Observatory elevation 'observatory['elevation']' must be non-negative (in m). Got {obs_data['elevation']}.")
        #
        if not isinstance(obs_data['timezone'], str):
            raise TypeError(f"Observatory timezone 'observatory['timezone']' must be an IANA timezone name (e.g. 'Atlantic/Canary'). Got {type(obs_data['timezone']).__name__}.")
        #
        if obs_data.get('name') is None:  # 'name' is optional
            obs_data['name'] = "Unknown Observatory"
        elif not isinstance(obs_data['name'], str):
            raise TypeError(f"Observatory name 'observatory['name']' must be a string. Got {type(obs_data['name']).__name__}.")


    # ##### OPTION 3: UNSUPPORTED TYPE #####
    else:
        raise TypeError(f"Observatory 'observatory' must be a pre-defined code (str) or a dictionary of parameters. Got {type(observatory).__name__}.")


    # ####################  COMPUTATION  ####################
    try:
        observer = Observer(
                            longitude=float(obs_data['longitude']) * u.deg,  # East positive
                            latitude=float(obs_data['latitude']) * u.deg,  # North positive
                            elevation=float(obs_data['elevation']) * u.m,  # Above sea level
                            timezone=obs_data['timezone'],  # IANA name, used for local times
                            name=obs_data['name'],
                            pressure=0 * u.mbar  # No refraction (geometric altitudes)
                           )
    except Exception as e:  # Typically an unknown timezone name
        raise ValueError(f"Observatory 'observatory' could not be created (check the IANA 'timezone' name, e.g. 'Atlantic/Canary'). Got '{obs_data['timezone']}'. Error: {e}") from e


    # ####################  OUTPUT  ####################
    return observer




# ========================================================================
# ============================== VALIDATION ==============================
# ========================================================================

# ==============================
def _check_bool(
                value : object,
                name  : str,
                label : str
    ) -> bool:
    """
    Validate a boolean flag argument.

    ----------
    Parameters
    ----------
        value : object
            Argument to check.

        name : str
            Name of the argument in the error messages (e.g. 'show').

        label : str
            Role of the argument in the error messages (e.g. 'Display flag').

    -------
    Returns
    -------
        value : bool
            The validated flag.

    ------
    Raises
    ------
        TypeError
            If 'value' is not a boolean, or 'name' or 'label' are not
            strings.

        ValueError
            If 'name' or 'label' are empty.

    -----
    Notes
    -----
        - Shared by the boolean flags of 'panorama.obs' and
          'panorama.prismas' (and of the 'prismas' helpers), so that the same
          argument gives the same error messages.
    """

    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")


    # ####################  COMPUTATION  ####################
    if not isinstance(value, bool):
        raise TypeError(f"{label} '{name}' must be a boolean. Got {type(value).__name__}.")


    # ####################  OUTPUT  ####################
    return value




# ==============================
def _check_period(
                  period : object
    ) -> float:
    """
    Validate the period of a periodic object (in days).

    ----------
    Parameters
    ----------
        period : object
            Period in days.

    -------
    Returns
    -------
        period : float
            The validated period, as a Python float.

    ------
    Raises
    ------
        TypeError
            If 'period' is not a real number.

        ValueError
            If 'period' is not finite or not positive.

    -----
    Notes
    -----
        - Shared by 'periodic_event', 'obs_phase' and 'obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  OUTPUT  ####################
    return _check_real(period, "period", "Period", "in days", "positive")




# ==============================
def _check_epochs(
                  hjd0 : object,
                  jd0  : object
    ) -> tuple[float | None, float | None, bool]:
    """
    Validate the reference epochs (phase 0) of a periodic object, given as HJD and / or JD.

    ----------
    Parameters
    ----------
        hjd0 : object
            Reference epoch as a Heliocentric Julian Date (UTC), or None.

        jd0 : object
            Reference epoch as a Julian Date (UTC), or None.

    -------
    Returns
    -------
        epochs : tuple[float | None, float | None, bool]
            Tuple (hjd0, jd0, hjd_flag): the validated epochs (Python floats
            or None) and True if the ephemeris is defined in HJD ('hjd0'
            given), False if it is defined in JD.

    ------
    Raises
    ------
        TypeError
            If an epoch is neither a real number nor None.

        ValueError
            If an epoch is not finite or not positive, or both are None.

    -----
    Notes
    -----
        - If both epochs are given, 'hjd0' is used and a warning is issued.

        - Shared by 'periodic_event', 'obs_phase' and 'obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import warnings


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    hjd0 = _check_real(hjd0, "hjd0", "Reference Heliocentric Julian Date", None, "positive", allow_none=True)
    #
    jd0 = _check_real(jd0, "jd0", "Reference Julian Date", None, "positive", allow_none=True)


    # ####################  COMPUTATION  ####################
    if hjd0 is None and jd0 is None:
        raise ValueError("Either the reference Heliocentric Julian Date 'hjd0' or the reference Julian Date 'jd0' must be provided. Got None for both.")
    elif hjd0 is not None and jd0 is not None:
        warnings.warn("Both the reference Heliocentric Julian Date 'hjd0' and the reference Julian Date 'jd0' were given; 'hjd0' is used and 'jd0' is ignored.", stacklevel=3)


    # ####################  OUTPUT  ####################
    return hjd0, jd0, hjd0 is not None




# ==============================
def _check_min_alt(
                   min_alt : object
    ) -> Quantity:
    """
    Validate a minimum altitude given in degrees or as an angular astropy Quantity.

    ----------
    Parameters
    ----------
        min_alt : object
            Minimum altitude: a real number (read in degrees) or a scalar
            astropy Quantity with angular units.

    -------
    Returns
    -------
        min_alt : Quantity
            The minimum altitude as a scalar Quantity in degrees.

    ------
    Raises
    ------
        TypeError
            If 'min_alt' is neither a real number nor an astropy Quantity.

        ValueError
            If 'min_alt' is not a single value, has no angular units, is not
            finite or lies outside [0, 90] deg.

    -----
    Notes
    -----
        - Shared by 'visibility', 'annual_visibility', 'periodic_event' and
          'obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    import astropy.units as u
    import numbers
    import numpy as np


    # ####################  COMPUTATION  ####################
    # ##########  TYPE AND UNITS  ##########
    if isinstance(min_alt, u.Quantity):
        if not min_alt.isscalar:
            raise ValueError(f"Minimum altitude 'min_alt' must be a single value. Got an array of shape {min_alt.shape}.")
        elif not min_alt.unit.is_equivalent(u.deg):
            raise ValueError(f"Minimum altitude 'min_alt' must have angular units (e.g. u.deg). Got {min_alt.unit}.")
        min_alt = min_alt.to(u.deg)  # E.g. arcmin or rad -> deg
    elif isinstance(min_alt, numbers.Real) and not isinstance(min_alt, bool):
        min_alt = float(min_alt) * u.deg  # Plain numbers are read in degrees
    else:
        raise TypeError(f"Minimum altitude 'min_alt' must be a real number (in deg) or an astropy Quantity with angular units. Got {type(min_alt).__name__}.")


    # ##########  VALUE  ##########
    if not np.isfinite(min_alt.value):
        raise ValueError(f"Minimum altitude 'min_alt' must be finite. Got {min_alt}.")
    elif not 0. <= min_alt.value <= 90.:
        raise ValueError(f"Minimum altitude 'min_alt' must lie in [0, 90] deg. Got {min_alt}.")


    # ####################  OUTPUT  ####################
    return min_alt




# ==============================
def _check_night(
                 night : object,
                 name  : str,
                 label : str
    ) -> str:
    """
    Validate a single night (or date) and return it as a 'YYYY-MM-DD' string; None means the current UTC date.

    ----------
    Parameters
    ----------
        night : object
            Night identified by the date on which it starts: a string
            ('YYYY-MM-DD'), a datetime, an astropy Time or None.

        name : str
            Name of the argument in the error messages (e.g. 'night').

        label : str
            Role of the argument in the error messages (e.g. 'Observing
            night').

    -------
    Returns
    -------
        night_str : str
            The night as 'YYYY-MM-DD' (UTC date).

    ------
    Raises
    ------
        TypeError
            If 'night' has an unsupported type, or 'name' or 'label' are not
            strings.

        ValueError
            If 'night' cannot be parsed or is not a single date, or 'name' or
            'label' are empty.

    -----
    Notes
    -----
        - Shared by 'night_info', 'visibility' and 'annual_visibility'.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.time import Time
    from datetime import datetime


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")


    # ####################  COMPUTATION  ####################
    if night is None:  # Night starting on the current UTC date
        return str(Time.now().utc.strftime("%Y-%m-%d"))
    elif not isinstance(night, (str, datetime, Time)):
        raise TypeError(f"{label} '{name}' must be a string ('YYYY-MM-DD'), a datetime, an astropy Time or None. Got {type(night).__name__}.")
    try:
        night_time = Time(night)  # Any format accepted by astropy (ISO string, datetime, Time)
    except Exception as e:
        raise ValueError(f"{label} '{name}' could not be parsed with astropy.time.Time (expected 'YYYY-MM-DD'). Got '{night}'. Error: {e}") from e
    if not night_time.isscalar:
        raise ValueError(f"{label} '{name}' must be a single date. Got {len(night_time)} dates.")


    # ####################  OUTPUT  ####################
    return str(night_time.utc.strftime("%Y-%m-%d"))




# ==============================
def _check_time(
                time  : object,
                name  : str,
                label : str
    ) -> Time:
    """
    Validate a single time (UTC) and return it as an astropy Time; None means now.

    ----------
    Parameters
    ----------
        time : object
            Time in UTC: a string (e.g. 'YYYY-MM-DD HH:MM:SS'), a datetime,
            an astropy Time or None.

        name : str
            Name of the argument in the error messages (e.g. 'obstime').

        label : str
            Role of the argument in the error messages (e.g. 'Observation
            time').

    -------
    Returns
    -------
        time : Time
            The time as a scalar astropy Time in the UTC scale.

    ------
    Raises
    ------
        TypeError
            If 'time' has an unsupported type, or 'name' or 'label' are not
            strings.

        ValueError
            If 'time' cannot be parsed or is not a single time, or 'name' or
            'label' are empty.

    -----
    Notes
    -----
        - Shared by 'periodic_event' and 'telescope_pointing'.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.time import Time
    from datetime import datetime


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")


    # ####################  COMPUTATION  ####################
    if time is None:  # Current time
        return Time.now().utc
    elif not isinstance(time, (str, datetime, Time)):
        raise TypeError(f"{label} '{name}' must be a string, a datetime, an astropy Time or None (in UTC). Got {type(time).__name__}.")
    try:
        time_obs = Time(time).utc  # Converted to the UTC scale
    except Exception as e:
        raise ValueError(f"{label} '{name}' could not be parsed with astropy.time.Time (e.g. 'YYYY-MM-DD HH:MM:SS', in UTC). Got '{time}'. Error: {e}") from e
    if not time_obs.isscalar:
        raise ValueError(f"{label} '{name}' must be a single time. Got {len(time_obs)} times.")


    # ####################  OUTPUT  ####################
    return time_obs




# ==============================
def _check_times(
                 times : object,
                 name  : str,
                 label : str
    ) -> Time:
    """
    Validate one or several times (UTC) and return them as an array-valued astropy Time.

    ----------
    Parameters
    ----------
        times : object
            Time or times in UTC: a string, a datetime, an astropy Time
            (scalar or array) or a list / tuple of them.

        name : str
            Name of the argument in the error messages (e.g. 'obstime').

        label : str
            Role of the argument in the error messages (e.g. 'Observation
            time').

    -------
    Returns
    -------
        times_obs : Time
            The times as a 1D astropy Time in the UTC scale.

    ------
    Raises
    ------
        TypeError
            If 'times' or any of its elements has an unsupported type, or
            'name' or 'label' are not strings.

        ValueError
            If 'times' is empty or any time cannot be parsed, or 'name' or
            'label' are empty.

    -----
    Notes
    -----
        - A single time is accepted and converted to a one-element list.

        - Shared by 'bary_corr', 'obs_phase' and 'obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.time import Time
    from datetime import datetime


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")


    # ####################  COMPUTATION  ####################
    # ##########  CONTAINER  ##########
    if isinstance(times, Time) and not times.isscalar:  # Array-valued Time object
        times = list(times)
    elif isinstance(times, (str, datetime, Time)):  # Single time
        times = [times]
    if not isinstance(times, (list, tuple)):
        raise TypeError(f"{label} '{name}' must be a string, a datetime, an astropy Time or a list of these (in UTC). Got {type(times).__name__}.")
    elif len(times) == 0:
        raise ValueError(f"{label} '{name}' must contain at least one time. Got an empty list.")


    # ##########  ELEMENTS  ##########
    times_list = []  # Scalar Time objects in UTC, in the input order
    for t in times:
        if not isinstance(t, (str, datetime, Time)):
            raise TypeError(f"Every element of {label.lower()} '{name}' must be a string, a datetime or an astropy Time (in UTC). Got {type(t).__name__}.")
        try:
            times_list.append(Time(t).utc)
        except Exception as e:
            raise ValueError(f"{label} '{name}' contains a time that could not be parsed with astropy.time.Time (in UTC). Got '{t}'. Error: {e}") from e


    # ####################  OUTPUT  ####################
    return Time(times_list)  # Array-valued Time




# ==============================
def _check_per_observation(
                           value      : object,
                           n_values   : int,
                           name       : str,
                           label      : str,
                           unit       : str,
                           constraint : str | None = None
    ) -> np.ndarray:
    """
    Validate a real number or a list with one real number per observation, and return it as an array.

    ----------
    Parameters
    ----------
        value : object
            A single real number (applied to every observation) or a list /
            tuple with one real number per observation.

        n_values : int
            Number of observations. Must be positive.

        name : str
            Name of the argument in the error messages (e.g. 'exptime').

        label : str
            Physical quantity of the argument in the error messages (e.g.
            'Exposure time').

        unit : str
            Unit shown in the type error (e.g. 'in s').

        constraint : str | None, optional
            Constraint of every value: 'positive', 'non-negative' or None
            (any finite value). Default is None.

    -------
    Returns
    -------
        values : np.ndarray
            Float array with 'n_values' elements.

    ------
    Raises
    ------
        TypeError
            If 'value' is neither a real number nor a list / tuple of real
            numbers, or 'n_values' is not an integer.

        ValueError
            If the list length differs from 'n_values', a value is not
            finite or does not satisfy 'constraint', or 'n_values' is not
            positive.

    -----
    Notes
    -----
        - Shared by the 'exptime' and 'rv' arguments of 'bary_corr'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numbers
    import numpy as np


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(n_values, numbers.Integral) or isinstance(n_values, bool):
        raise TypeError(f"Number of observations 'n_values' must be an integer. Got {type(n_values).__name__}.")
    elif n_values <= 0:
        raise ValueError(f"Number of observations 'n_values' must be positive. Got {n_values}.")


    # ####################  COMPUTATION  ####################
    if isinstance(value, numbers.Real) and not isinstance(value, bool):  # Single value for every observation
        values = [value] * n_values
    elif isinstance(value, (list, tuple)):
        if len(value) != n_values:
            raise ValueError(f"{label} '{name}' must have one value per observation time ({n_values}). Got {len(value)}.")
        values = list(value)  # Copy: the input list is never modified
    else:
        raise TypeError(f"{label} '{name}' must be a real number ({unit}), a list of real numbers or None. Got {type(value).__name__}.")
    values = [_check_real(v, name, label, unit, constraint) for v in values]  # Type, finiteness and constraint of every value


    # ####################  OUTPUT  ####################
    return np.array(values, dtype=float)




# ==================================================================
# ============================== TIME ==============================
# ==================================================================

# ==============================
def _is_event(
              t : Time | None
    ) -> bool:
    """
    Check whether an astroplan event time (rise, set, twilight...) actually occurs.

    ----------
    Parameters
    ----------
        t : Time | None
            Event time returned by astroplan, possibly masked.

    -------
    Returns
    -------
        is_event : bool
            False if 't' is None, masked or not finite (the event does not
            happen, e.g. no astronomical twilight in summer at high
            latitudes), True otherwise.

    ------
    Raises
    ------
        TypeError
            If 't' is neither an astropy Time nor None.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.time import Time
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if t is not None and not isinstance(t, Time):
        raise TypeError(f"Event time 't' must be an astropy Time or None. Got {type(t).__name__}.")


    # ####################  COMPUTATION  ####################
    if t is None:
        return False
    elif getattr(t, "masked", False) and np.any(t.mask):  # Masked times: the event does not occur
        return False
    jd = np.asarray(getattr(t.jd, "unmasked", t.jd), dtype=float)  # Underlying values if 't.jd' is a masked array


    # ####################  OUTPUT  ####################
    return bool(np.all(np.isfinite(jd)))




# ==============================
def _event_iso(
               t : Time | None
    ) -> str | None:
    """
    Format an astroplan event time as a UTC string with seconds, or return None if the event does not occur.

    ----------
    Parameters
    ----------
        t : Time | None
            Event time returned by astroplan (possibly masked).

    -------
    Returns
    -------
        iso : str | None
            Time as 'YYYY-MM-DD HH:MM:SS' UTC, or None if the event does not
            occur (see '_is_event').

    ------
    Raises
    ------
        TypeError
            If 't' is neither an astropy Time nor None.

    -----
    Notes
    -----
        - Used by 'panorama.obs.night_info'.
    """

    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    occurs = _is_event(t)  # Also checks the type of 't'


    # ####################  OUTPUT  ####################
    return str(t.utc.strftime("%Y-%m-%d %H:%M:%S")) if occurs else None




# ==============================
def _event_minute(
                  t       : Time | None,
                  t_start : Time,
                  t_end   : Time
    ) -> str | None:
    """
    Format an astroplan event time as a UTC string rounded to the minute, or return None if it does not occur within a time interval.

    ----------
    Parameters
    ----------
        t : Time | None
            Event time returned by astroplan (possibly masked or None).

        t_start : Time
            Start of the interval (e.g. start of the plotted night).

        t_end : Time
            End of the interval (e.g. end of the plotted night).

    -------
    Returns
    -------
        time_str : str | None
            Time as 'YYYY-MM-DD HH:MM' UTC (rounded to the nearest minute),
            or None if the event does not occur or falls outside
            [t_start, t_end].

    ------
    Raises
    ------
        TypeError
            If 't' is neither an astropy Time nor None, or 't_start' or
            't_end' is not an astropy Time.

    -----
    Notes
    -----
        - Used by 'panorama.obs.visibility'.
    """

    # ####################  MODULES & PACKAGES  ####################
    import astropy.units as u


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    t_night = _event_within(t, t_start, t_end)  # Also checks the types of the arguments


    # ####################  OUTPUT  ####################
    if t_night is None:
        return None
    return str((t_night + 30 * u.s).utc.strftime("%Y-%m-%d %H:%M"))  # Rounded to the nearest minute (strftime truncates the seconds)




# ==============================
def _event_within(
                  t       : Time | None,
                  t_start : Time,
                  t_end   : Time
    ) -> Time | None:
    """
    Keep an astroplan event time only if it occurs within a time interval (e.g. the night).

    ----------
    Parameters
    ----------
        t : Time | None
            Event time returned by astroplan (possibly masked or None).

        t_start : Time
            Start of the interval (e.g. sunset, or midnight - 12 h).

        t_end : Time
            End of the interval (e.g. sunrise, or midnight + 12 h).

    -------
    Returns
    -------
        t_night : Time | None
            The event time, or None if it does not occur or falls outside
            [t_start, t_end].

    ------
    Raises
    ------
        TypeError
            If 't' is neither an astropy Time nor None, or 't_start' or
            't_end' is not an astropy Time.

    -----
    Notes
    -----
        - Used by 'panorama.obs.night_info' (moonrise and moonset) and by
          '_event_minute'.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.time import Time


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    occurs = _is_event(t)  # Also checks the type of 't'
    #
    if not isinstance(t_start, Time):
        raise TypeError(f"Start of the interval 't_start' must be an astropy Time. Got {type(t_start).__name__}.")
    #
    if not isinstance(t_end, Time):
        raise TypeError(f"End of the interval 't_end' must be an astropy Time. Got {type(t_end).__name__}.")


    # ####################  OUTPUT  ####################
    if not occurs or not (t_start.jd <= t.jd <= t_end.jd):  # Compared as Julian Dates (any time scale offset is negligible here)
        return None
    return t




# ==============================
def _event_hours(
                 t1 : Time | None,
                 t2 : Time | None
    ) -> float:
    """
    Compute the time between two astroplan events, in hours.

    ----------
    Parameters
    ----------
        t1 : Time | None
            First event (possibly masked or None).

        t2 : Time | None
            Second event (possibly masked or None).

    -------
    Returns
    -------
        hours : float
            t2 - t1 in hours, or NaN if any of the events does not occur.

    ------
    Raises
    ------
        TypeError
            If 't1' or 't2' is neither an astropy Time nor None.

    -----
    Notes
    -----
        - Used by 'panorama.obs.night_info'.
    """

    # ####################  MODULES & PACKAGES  ####################
    import astropy.units as u
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    occur = _is_event(t1) and _is_event(t2)  # Also checks the types of 't1' and 't2'


    # ####################  OUTPUT  ####################
    if not occur:
        return np.nan
    return float((t2 - t1).to(u.hour).value)




# ==============================
def _format_hhmm(
                 hours : int | float
    ) -> str:
    """
    Format a duration in hours as 'hh:mm'.

    ----------
    Parameters
    ----------
        hours : int | float
            Duration in hours (NaN if not defined).

    -------
    Returns
    -------
        duration : str
            Duration as 'hh:mm', or '--:--' if 'hours' is not finite.

    ------
    Raises
    ------
        TypeError
            If 'hours' is not a real number.

    -----
    Notes
    -----
        - Used by 'panorama.obs.visibility' for the table of night events.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(hours, numbers.Real) or isinstance(hours, bool):
        raise TypeError(f"Duration 'hours' (in h) must be a real number. Got {type(hours).__name__}.")


    # ####################  COMPUTATION  ####################
    if not np.isfinite(hours):  # Undefined duration (e.g. the Sun does not set)
        return "--:--"
    minutes = int(round(hours * 60.))  # Rounded to the minute before splitting, to avoid '60' minutes


    # ####################  OUTPUT  ####################
    return f"{minutes // 60:02d}:{minutes % 60:02d}"




# ==============================
def _above(
           y         : np.ndarray[float],
           threshold : float
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Find the part of each grid interval where a linearly interpolated quantity is above a threshold.

    ----------
    Parameters
    ----------
        y : np.ndarray[float]
            Quantity sampled on a regular grid along the last axis (e.g.
            an altitude, in degrees). At least two grid points.

        threshold : float
            Threshold, in the same units as 'y'.

    -------
    Returns
    -------
        limits : tuple[np.ndarray[float], np.ndarray[float]]
            Tuple (lo, hi) with the start and end of the part of each
            interval [i, i + 1] with y >= threshold, as fractions of the
            interval (0-1); lo >= hi if no part is above the threshold.

    ------
    Raises
    ------
        TypeError
            If 'y' is not a real NumPy array or 'threshold' is not a real
            number.

        ValueError
            If 'y' has fewer than two grid points along its last axis, or
            'threshold' is not finite.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(y, np.ndarray) or not np.issubdtype(y.dtype, np.number) or y.dtype == bool:
        raise TypeError(f"Sampled quantity 'y' must be a real NumPy array. Got {type(y).__name__}.")
    elif y.ndim == 0 or y.shape[-1] < 2:
        raise ValueError(f"Sampled quantity 'y' must have at least 2 grid points along its last axis. Got shape {y.shape}.")
    #
    if not isinstance(threshold, numbers.Real) or isinstance(threshold, bool):
        raise TypeError(f"Threshold 'threshold' must be a real number. Got {type(threshold).__name__}.")
    elif not np.isfinite(threshold):
        raise ValueError(f"Threshold 'threshold' must be finite. Got {threshold}.")


    # ####################  COMPUTATION  ####################
    y0, y1 = y[..., :-1], y[..., 1:]  # Values at the start and end of each interval
    with np.errstate(divide="ignore", invalid="ignore"):  # y0 == y1 only when both are on the same side (fraction unused)
        fraction = np.clip((threshold - y0) / (y1 - y0), 0., 1.)  # Linear interpolation of the crossing
    lo = np.where(y0 >= threshold, 0., np.where(y1 >= threshold, fraction, 1.))  # Starts above (0), crosses upward (fraction) or never above (1)
    hi = np.where(y1 >= threshold, 1., np.where(y0 >= threshold, fraction, 0.))  # Ends above (1), crosses downward (fraction) or never above (0)


    # ####################  OUTPUT  ####################
    return lo, hi




# ==============================
def _observable(
                sun_alt     : np.ndarray[float],
                target_alt  : np.ndarray[float],
                min_alt_deg : float
    ) -> tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]]:
    """
    Find the part of each grid interval with the Sun below -18 deg and the target above the minimum altitude.

    ----------
    Parameters
    ----------
        sun_alt : np.ndarray[float]
            Altitude of the Sun along the grid, in degrees.

        target_alt : np.ndarray[float]
            Altitude of the target along the grid, in degrees (same last
            axis as 'sun_alt'; leading axes are broadcast, e.g. one row per
            target).

        min_alt_deg : float
            Minimum altitude of the target, in degrees, in [-90, 90].

    -------
    Returns
    -------
        intervals : tuple[np.ndarray[float], np.ndarray[float], np.ndarray[float]]
            Tuple (lo, hi, length) with the start, end and length of the
            observable part of each interval, as fractions of the
            interval (length is 0 where the target is not observable).

    ------
    Raises
    ------
        TypeError
            If the altitudes are not real NumPy arrays or 'min_alt_deg' is
            not a real number.

        ValueError
            If the arrays have fewer than two grid points or different grid
            lengths, or 'min_alt_deg' is not finite or outside [-90, 90].
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(sun_alt, np.ndarray) or not np.issubdtype(sun_alt.dtype, np.number) or sun_alt.dtype == bool:
        raise TypeError(f"Altitude of the Sun 'sun_alt' (in deg) must be a real NumPy array. Got {type(sun_alt).__name__}.")
    elif sun_alt.ndim == 0 or sun_alt.shape[-1] < 2:
        raise ValueError(f"Altitude of the Sun 'sun_alt' must have at least 2 grid points. Got shape {sun_alt.shape}.")
    #
    if not isinstance(target_alt, np.ndarray) or not np.issubdtype(target_alt.dtype, np.number) or target_alt.dtype == bool:
        raise TypeError(f"Altitude of the target 'target_alt' (in deg) must be a real NumPy array. Got {type(target_alt).__name__}.")
    elif target_alt.ndim == 0 or target_alt.shape[-1] != sun_alt.shape[-1]:
        raise ValueError(f"Altitude of the target 'target_alt' must have the same number of grid points as 'sun_alt' ({sun_alt.shape[-1]}). Got shape {target_alt.shape}.")
    #
    if not isinstance(min_alt_deg, numbers.Real) or isinstance(min_alt_deg, bool):
        raise TypeError(f"Minimum altitude 'min_alt_deg' (in deg) must be a real number. Got {type(min_alt_deg).__name__}.")
    elif not np.isfinite(min_alt_deg):
        raise ValueError(f"Minimum altitude 'min_alt_deg' (in deg) must be finite. Got {min_alt_deg}.")
    elif not -90. <= min_alt_deg <= 90.:
        raise ValueError(f"Minimum altitude 'min_alt_deg' must lie in [-90, 90] deg. Got {min_alt_deg}.")


    # ####################  COMPUTATION  ####################
    lo_sun, hi_sun       = _above(-sun_alt, 18.)  # Sun below -18 deg (astronomical night)
    lo_target, hi_target = _above(target_alt, min_alt_deg)  # Target above the minimum altitude
    lo, hi               = np.maximum(lo_sun, lo_target), np.minimum(hi_sun, hi_target)  # Intersection of both conditions


    # ####################  OUTPUT  ####################
    return lo, hi, np.clip(hi - lo, 0., None)  # Negative lengths (empty intersections) set to 0




# ==============================
def _mean_midnight(
                   observer : Observer,
                   night    : str
    ) -> Time:
    """
    Compute the mean local midnight (UTC) of the observing night that starts on a given date.

    ----------
    Parameters
    ----------
        observer : Observer
            Observatory (see '_get_observer').

        night : str
            Date on which the night starts ('YYYY-MM-DD').

    -------
    Returns
    -------
        mean_midnight : Time
            Mean local midnight: 00:00 UTC of the next day minus longitude /
            15 h (East positive).

    ------
    Raises
    ------
        TypeError
            If 'observer' is not an astroplan Observer or 'night' is not a
            string.

        ValueError
            If 'night' is not a valid 'YYYY-MM-DD' date.

    -----
    Notes
    -----
        - See '_solar_midnight' for the true (apparent) solar midnight.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astroplan import Observer
    from astropy.time import Time
    import astropy.units as u
    from datetime import datetime


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(observer, Observer):
        raise TypeError(f"Observatory 'observer' must be an astroplan Observer. Got {type(observer).__name__}.")
    #
    if not isinstance(night, str):
        raise TypeError(f"Night date 'night' must be a 'YYYY-MM-DD' string. Got {type(night).__name__}.")
    try:
        datetime.strptime(night.strip(), "%Y-%m-%d")
    except ValueError as e:
        raise ValueError(f"Night date 'night' must be a valid 'YYYY-MM-DD' date. Got '{night}'.") from e


    # ####################  COMPUTATION  ####################
    lon_deg = observer.location.lon.wrap_at(180 * u.deg).deg  # East positive, in [-180, 180)


    # ####################  OUTPUT  ####################
    return Time(night.strip(), scale="utc") + 1 * u.day - (lon_deg / 360.) * u.day  # 00:00 UTC of the next day minus longitude / 15 h




# ==============================
def _solar_midnight(
                    observer      : Observer,
                    mean_midnight : Time
    ) -> Time:
    """
    Compute the true (apparent) solar midnight, with the Sun at lower culmination, from the mean local midnight.

    ----------
    Parameters
    ----------
        observer : Observer
            Observatory (see '_get_observer').

        mean_midnight : Time
            Mean local midnight (UTC), scalar or array (see '_mean_midnight').

    -------
    Returns
    -------
        midnight : Time
            True solar midnight (UTC), with the same shape as
            'mean_midnight'.

    ------
    Raises
    ------
        TypeError
            If 'observer' is not an astroplan Observer or 'mean_midnight' is
            not an astropy Time.

    -----
    Notes
    -----
        - The true solar midnight (hour angle of the Sun of 12 h) differs
          from the mean local midnight by the equation of time (< 17 min).
          It is computed from the apparent hour angle of the Sun at mean
          midnight, while 'Observer.midnight' interpolates on a coarse grid
          and can be off by ~2-3 min.

        - The hour angle of the Sun grows by 1 h per hour of time: the
          sidereal rate (1.00274) minus the apparent motion of the Sun in
          right ascension (~0.00274, i.e. ~1 deg/day). The hour-angle
          offset is therefore used directly as a time offset; the residual
          error, due to the variation of the solar motion in right
          ascension along the year, is below 0.01 s.

        - Shared by 'night_info', 'visibility', 'annual_visibility' and
          'obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astroplan import Observer
    from astropy.coordinates import get_sun, TETE
    from astropy.time import Time
    import astropy.units as u


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(observer, Observer):
        raise TypeError(f"Observatory 'observer' must be an astroplan Observer. Got {type(observer).__name__}.")
    #
    if not isinstance(mean_midnight, Time):
        raise TypeError(f"Mean local midnight 'mean_midnight' must be an astropy Time. Got {type(mean_midnight).__name__}.")


    # ####################  COMPUTATION  ####################
    sun_ra = get_sun(mean_midnight).transform_to(TETE(obstime=mean_midnight)).ra  # Apparent RA of the Sun
    sun_ha = (observer.local_sidereal_time(mean_midnight) - sun_ra).wrap_at(24 * u.hourangle).hour  # Close to 12 h


    # ####################  OUTPUT  ####################
    return mean_midnight + (12. - sun_ha) * u.hour  # The hour angle of the Sun grows by 1 h per hour (see Notes)




# ==============================
def _reference_epochs(
                      hjd0     : float | None,
                      jd0      : float | None,
                      hjd_flag : bool,
                      coords   : SkyCoord,
                      location : EarthLocation
    ) -> tuple[float, float]:
    """
    Convert the reference epoch of a periodic object between HJD and JD (both in UTC).

    ----------
    Parameters
    ----------
        hjd0 : float | None
            Reference epoch as a Heliocentric Julian Date (UTC), or None.

        jd0 : float | None
            Reference epoch as a Julian Date (UTC), or None.

        hjd_flag : bool
            True if the ephemeris is defined in HJD ('hjd0' is used), False
            if it is defined in JD ('jd0' is used).

        coords : SkyCoord
            Coordinates of the object.

        location : EarthLocation
            Location of the observatory.

    -------
    Returns
    -------
        epochs : tuple[float, float]
            Tuple (hjd0, jd0) with the reference epoch in both time stamps.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the epoch selected by 'hjd_flag' is None.

    -----
    Notes
    -----
        - The heliocentric light-travel time is evaluated at the given epoch
          (for HJD -> JD the error is < 0.01 s).

        - Shared by 'periodic_event', 'obs_phase' and 'obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    import astropy.units as u
    import numbers


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for epoch_name, epoch in (("hjd0", hjd0), ("jd0", jd0)):
        if epoch is not None and (not isinstance(epoch, numbers.Real) or isinstance(epoch, bool)):
            raise TypeError(f"Reference epoch '{epoch_name}' must be a real number or None. Got {type(epoch).__name__}.")
    #
    if not isinstance(hjd_flag, bool):
        raise TypeError(f"HJD flag 'hjd_flag' must be a boolean. Got {type(hjd_flag).__name__}.")
    elif (hjd0 if hjd_flag else jd0) is None:
        raise ValueError(f"Reference epoch '{'hjd0' if hjd_flag else 'jd0'}' selected by 'hjd_flag' must not be None.")
    #
    if not isinstance(coords, SkyCoord):
        raise TypeError(f"Coordinates 'coords' must be an astropy SkyCoord. Got {type(coords).__name__}.")
    #
    if not isinstance(location, EarthLocation):
        raise TypeError(f"Observatory location 'location' must be an astropy EarthLocation. Got {type(location).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: HJD -> JD  ##########
    if hjd_flag:
        ltt_helio = Time(hjd0, format="jd", scale="utc").light_travel_time(coords, "heliocentric", location=location)  # Evaluated at the HJD (error < 0.01 s)
        jd0       = hjd0 - ltt_helio.to_value(u.day)  # JD = HJD - light-travel time


    # ##########  OPTION 2: JD -> HJD  ##########
    else:
        ltt_helio = Time(jd0, format="jd", scale="utc").light_travel_time(coords, "heliocentric", location=location)
        hjd0      = jd0 + ltt_helio.to_value(u.day)  # HJD = JD + light-travel time


    # ####################  OUTPUT  ####################
    return float(hjd0), float(jd0)




# ==============================
def _phases(
            times    : Time,
            coords   : SkyCoord,
            location : EarthLocation,
            period   : float,
            hjd0     : float,
            jd0      : float,
            hjd_flag : bool
    ) -> np.ndarray[float]:
    """
    Compute the phases of a periodic object at the given UTC times, from its reference epoch and period.

    ----------
    Parameters
    ----------
        times : Time
            Observation times (UTC), scalar or array.

        coords : SkyCoord
            Coordinates of the object.

        location : EarthLocation
            Location of the observatory.

        period : float
            Period in days. Must be positive.

        hjd0 : float
            Reference epoch (phase 0) as a Heliocentric Julian Date (UTC).

        jd0 : float
            Reference epoch (phase 0) as a Julian Date (UTC).

        hjd_flag : bool
            True if the ephemeris is defined in HJD (phases from HJD and
            'hjd0'), False if it is defined in JD (phases from JD and 'jd0').

    -------
    Returns
    -------
        phases : np.ndarray[float]
            Phases in [0, 1), as a 1D array.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If 'period' or the reference epoch selected by 'hjd_flag' is not
            finite or not positive.

    -----
    Notes
    -----
        - Phases are computed as:
                φ = ((t - T0) mod P) / P
          using HJD for both t and T0 if 'hjd_flag' is True, and JD
          otherwise. HJD values are in UTC (classical convention).

        - Shared by 'panorama.obs.obs_phase' and 'panorama.obs.obs_phase_night'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    import astropy.units as u
    import numpy as np


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(times, Time):
        raise TypeError(f"Observation times 'times' must be an astropy Time. Got {type(times).__name__}.")
    #
    if not isinstance(coords, SkyCoord):
        raise TypeError(f"Coordinates 'coords' must be an astropy SkyCoord. Got {type(coords).__name__}.")
    #
    if not isinstance(location, EarthLocation):
        raise TypeError(f"Observatory location 'location' must be an astropy EarthLocation. Got {type(location).__name__}.")
    #
    period = _check_real(period, "period", "Period", "in days", "positive")
    #
    if not isinstance(hjd_flag, bool):
        raise TypeError(f"HJD flag 'hjd_flag' must be a boolean. Got {type(hjd_flag).__name__}.")
    #
    if hjd_flag:  # Only the epoch used by the ephemeris is needed
        hjd0 = _check_real(hjd0, "hjd0", "Reference Heliocentric Julian Date", None, "positive")
    else:
        jd0 = _check_real(jd0, "jd0", "Reference Julian Date", None, "positive")


    # ####################  COMPUTATION  ####################
    times_jd = np.atleast_1d(times.utc.jd)  # JD (UTC)
    if hjd_flag:
        ltt_helio = times.light_travel_time(coords, "heliocentric", location=location)  # Heliocentric light-travel time at each time
        phases    = ((times_jd + np.atleast_1d(ltt_helio.to_value(u.day)) - hjd0) % period) / period  # From HJD (UTC)
    else:
        phases = ((times_jd - jd0) % period) / period  # From JD (UTC)


    # ####################  OUTPUT  ####################
    return phases




# ======================================================================
# ============================== PLOTTING ==============================
# ======================================================================

# ==============================
def _draw_phase_range(
                      ax          : Axes,
                      y           : int | float,
                      init_phase  : float,
                      final_phase : float,
                      full        : bool,
                      linewidth   : int | float
    ) -> None:
    """
    Draw a phase range as horizontal bar(s), splitting it when it wraps through phase 0.

    ----------
    Parameters
    ----------
        ax : Axes
            matplotlib.axes.Axes where the range is drawn.

        y : int | float
            Vertical position of the bar(s), in data coordinates.

        init_phase : float
            Phase at the start of the observable window, in [0, 1].

        final_phase : float
            Phase at the end of the observable window, in [0, 1]. If smaller
            than 'init_phase', the range wraps through phase 0.

        full : bool
            If True, the whole phase interval [0, 1] is drawn.

        linewidth : int | float
            Line width of the bar(s), in points. Must be positive.

    -------
    Returns
    -------
        None
            The bar(s) are drawn on 'ax'.

    ------
    Raises
    ------
        TypeError
            If 'ax' is not a matplotlib Axes, 'full' is not a boolean, or
            any other input is not a real number.

        ValueError
            If any value is not finite, a phase lies outside [0, 1] or
            'linewidth' is not positive.

    -----
    Notes
    -----
        - Used by 'panorama.obs.obs_phase_night' for the phases covered on
          each night and for the cumulative coverage.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from matplotlib.axes import Axes


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(ax, Axes):
        raise TypeError(f"Axes 'ax' must be a matplotlib Axes. Got {type(ax).__name__}.")
    #
    y = _check_real(y, "y", "Vertical position of the bar")
    #
    init_phase = _check_real(init_phase, "init_phase", "Initial phase", "dimensionless", (0., 1., "[]"))
    #
    final_phase = _check_real(final_phase, "final_phase", "Final phase", "dimensionless", (0., 1., "[]"))
    #
    if not isinstance(full, bool):
        raise TypeError(f"Full coverage flag 'full' must be a boolean. Got {type(full).__name__}.")
    #
    linewidth = _check_real(linewidth, "linewidth", "Line width", "in points", "positive")


    # ####################  PLOTTING  ####################
    # ##########  OPTION 1: FULL PERIOD  ##########
    if full:  # Window longer than the period: every phase is covered
        ax.hlines(y=y, xmin=0., xmax=1., color="black", linewidth=linewidth)


    # ##########  OPTION 2: SINGLE RANGE  ##########
    elif init_phase <= final_phase:
        ax.hlines(y=y, xmin=init_phase, xmax=final_phase, color="black", linewidth=linewidth)


    # ##########  OPTION 3: RANGE WRAPPING THROUGH PHASE 0  ##########
    else:
        ax.hlines(y=y, xmin=init_phase, xmax=1., color="black", linewidth=linewidth)  # From the start to phase 1
        ax.hlines(y=y, xmin=0., xmax=final_phase, color="black", linewidth=linewidth)  # From phase 0 to the end




# ==================================================================
# ============================== FITS ==============================
# ==================================================================

# ==============================
def _compute_obs_night(
                       date_obs         : str | None,
                       night_split_hour : int
    ) -> list[str]:
    """
    Compute the observing night and the observation time from a FITS DATE-OBS value.

    ----------
    Parameters
    ----------
        date_obs : str | None
            Observation date and time in ISO format ('YYYY-MM-DDTHH:MM:SS',
            with or without fractional seconds or a trailing 'Z', or
            'YYYY-MM-DD').

        night_split_hour : int
            UTC hour (0-23) before which an observation is assigned to the
            night that started on the previous calendar date.

    -------
    Returns
    -------
        night_time : list[str]
            Two-element list with the observing night ('YYYY-MM-DD', date on
            which the night started) and the UTC observation time
            ('HH:MM:SS'). A date without time is taken as the night itself
            (time 'N/A'), and ['N/A', 'N/A'] is returned if 'date_obs' is
            missing or cannot be parsed.

    ------
    Raises
    ------
        TypeError
            If 'night_split_hour' is not an integer.

        ValueError
            If 'night_split_hour' is outside [0, 23].

    -----
    Notes
    -----
        - Example with night_split_hour = 9 (HERMES): an exposure at
          '2026-11-02T03:10:00' belongs to the night of '2026-11-01', while
          one at '2026-11-01T22:00:00' belongs to the same date.

        - The observation time is the UTC time of DATE-OBS, not shifted.

        - A missing or unreadable DATE-OBS does not raise an error, so that
          one bad header never stops the cataloging of a whole folder.
    """

    # ####################  MODULES & PACKAGES  ####################
    from datetime import datetime, timedelta
    import numbers
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(night_split_hour, numbers.Integral) or isinstance(night_split_hour, bool):
        raise TypeError(f"Night split hour 'night_split_hour' must be an integer (UTC hour). Got {type(night_split_hour).__name__}.")
    elif not 0 <= night_split_hour <= 23:
        raise ValueError(f"Night split hour 'night_split_hour' must lie in [0, 23] (UTC hour). Got {night_split_hour}.")
    #
    if not isinstance(date_obs, str) or not date_obs.strip():  # Missing DATE-OBS: no error (see Notes)
        return ["N/A", "N/A"]


    # ####################  COMPUTATION  ####################
    # ##########  PARSING  ##########
    date_str = date_obs.strip()
    if date_str.endswith("Z"):  # UTC designator, not understood by fromisoformat in Python < 3.11
        date_str = date_str[:-1]

    try:
        dt = datetime.fromisoformat(date_str)  # 'YYYY-MM-DD', 'YYYY-MM-DDTHH:MM:SS[.ffffff]'
    except ValueError:
        dt = None
        for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):  # Fallback for non-standard fractional seconds
            try:
                dt = datetime.strptime(date_str, fmt)
                break
            except ValueError:
                pass

    if dt is None:  # Unreadable DATE-OBS: no error (see Notes)
        return ["N/A", "N/A"]
    elif "T" not in date_str and " " not in date_str:  # Date without time: it already identifies the night
        warnings.warn(f"DATE-OBS '{date_obs}' has no time component. The night is taken as '{dt.strftime('%Y-%m-%d')}' and the time as 'N/A'.", UserWarning, stacklevel=2)
        return [dt.strftime("%Y-%m-%d"), "N/A"]


    # ##########  OBSERVING NIGHT  ##########
    time_str = dt.strftime("%H:%M:%S")  # UTC time of the exposure, not shifted
    if dt.hour < night_split_hour:  # Early-morning exposure: the night started on the previous date
        dt = dt - timedelta(days=1)


    # ####################  OUTPUT  ####################
    return [dt.strftime("%Y-%m-%d"), time_str]




# ==============================
def _get_header_value(
                      header  : Header,
                      aliases : list[str] | tuple[str, ...],
                      default : object
    ) -> object:
    """
    Return the value of the first keyword of a list that is present in a FITS header.

    ----------
    Parameters
    ----------
        header : Header
            astropy.io.fits.Header to search.

        aliases : list[str] | tuple[str, ...]
            Possible keywords, tried in order; the first one present wins.

        default : object
            Value returned if none of the keywords is present.

    -------
    Returns
    -------
        value : object
            Value of the first matching keyword, or 'default'.

    ------
    Raises
    ------
        TypeError
            If 'header' is not a FITS header, or 'aliases' is not a list or
            tuple of strings.

    -----
    Notes
    -----
        - Keywords present but empty (e.g. OBJECT = '') are skipped, so that
          the next alias or the default is used.

        - String values are returned without surrounding whitespace.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.io.fits import Header


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(header, Header):
        raise TypeError(f"FITS header 'header' must be an astropy.io.fits.Header. Got {type(header).__name__}.")
    #
    if not isinstance(aliases, (list, tuple)) or not all(isinstance(key, str) for key in aliases):
        raise TypeError(f"Header keywords 'aliases' must be a list or tuple of strings. Got {aliases!r}.")


    # ####################  COMPUTATION  ####################
    for key in aliases:
        if key in header:
            value = header[key]
            if isinstance(value, str):
                value = value.strip()  # FITS strings are padded with spaces
                if not value:  # Empty keyword: try the next alias
                    continue
            return value  # First present, non-empty keyword


    # ####################  OUTPUT  ####################
    return default  # None of the keywords is present




# ======================================================================
# ============================== BINARALT ==============================
# ======================================================================

# ==============================
class _WeatherService:
    """
    Retrieve and normalize the live weather data of the observatories shown in BinarAlt.

    BinarAlt's weather panel needs live data from each observatory's own
    weather page, but none of those sites send CORS headers, so a fetch()
    from the page itself (served from this local server's own origin) would
    always be blocked by the browser. Fetching server-side has no such
    restriction, so the page calls this server's own '/weather' endpoint and
    this proxy does the upstream request.

    None of these sites expose a documented, stable API: each fetcher
    targets the data source used by that site's own front-end (a JSON API
    for the OT stations and the TNG, or a small regular expression over the
    plain HTML for the others), found by inspecting each site's network
    requests. They can break if a site changes its markup. Each fetcher
    fails independently, falling back to '-' per field, and never raises
    into the HTTP response.

    ----------
    Parameters
    ----------
        cache_ttl : int | float, optional
            Time, in seconds, during which a reading is reused before
            fetching it again. Must be finite and non-negative. Default is
            20.

    -----
    Notes
    -----
        1. STATIONS
            - Observatorio del Teide: 'song', 'stella', 'speculoos'.

            - Observatorio del Roque de los Muchachos: 'not', 'tng', 'gtc',
              'int', 'magic', 'liverpool'.

        2. RESULTS
            - Every result has the same keys: 'ok', 'source_url',
              'temperature', 'humidity', 'wind_speed', 'wind_direction',
              'precipitation', 'dust' and 'seeing' (formatted strings, '-'
              if not available; see '_blank').

            - Wind directions follow the meteorological convention: the
              direction the wind blows FROM (0 deg / N = from the north).

        3. THREADS
            - The cache is protected by a lock, since the server handles
              requests in several threads.

        4. METHODS LAYOUT
            - Methods are grouped by level-1 headers inside the class, so
              their own body sections use level-2 and level-3 headers.
    """

    # ####################  CONFIGURATION  ####################
    _OCAN_URL = "https://ocan.iac.es/meteo/api/v1/station/{station}/last"  # OT stations JSON API
    _TNG_URL  = "https://tngweb.tng.iac.es/api/meteo/weather/{field}"  # TNG time-series JSON API
    _OCAN_FIELDS = {  # Fields shown for each OT station
                    'song'      : ('temperature', 'humidity', 'wind_speed', 'wind_direction', 'rain', 'seeing'),
                    'stella'    : ('temperature', 'humidity', 'wind_speed', 'wind_direction', 'dust'),
                    'speculoos' : ('temperature', 'humidity', 'wind_speed'),
                   }


    # ####################  INIT  ####################
    def __init__(
                 self,
                 cache_ttl : int | float = 20
        ) -> None:
        """
        Initialize the cache and the table of weather fetchers.

        ----------
        Parameters
        ----------
            cache_ttl : int | float, optional
                Time, in seconds, during which a reading is reused. Must be
                finite and non-negative. Default is 20.

        -------
        Returns
        -------
            None
                The cache, its lock and the table of fetchers are stored as
                attributes.

        ------
        Raises
        ------
            TypeError
                If 'cache_ttl' is not a real number.

            ValueError
                If 'cache_ttl' is not finite or is negative.
        """

        # ##########  MODULES & PACKAGES  ##########
        # ##### PYTHON #####
        import threading


        # ##### PANORAMA #####
        from ..units.units_support import _check_real


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        cache_ttl = _check_real(cache_ttl, "cache_ttl", "Cache lifetime", "in s", "non-negative")


        # ##########  PARAMETERS  ##########
        # ##### CACHE #####
        self.cache     = {}  # Station -> (monotonic time of the reading, payload)
        self.cache_ttl = cache_ttl
        self._lock     = threading.Lock()


        # ##### FETCHERS #####
        self.fetchers = {
                         'song'      : lambda: self._fetch_ocan('song'),
                         'stella'    : lambda: self._fetch_ocan('stella'),
                         'speculoos' : lambda: self._fetch_ocan('speculoos'),
                         'not'       : self._fetch_not,
                         'tng'       : self._fetch_tng,
                         'gtc'       : self._fetch_gtc,
                         'int'       : self._fetch_int,
                         'magic'     : self._fetch_magic,
                         'liverpool' : self._fetch_liverpool,
                        }


    # ####################  HELPERS  ####################
    def _fetch(
               self,
               url     : str,
               timeout : int | float = 8
        ) -> str:
        """
        Fetch a URL and return its response as UTF-8 text.

        ----------
        Parameters
        ----------
            url : str
                URL to request.

            timeout : int | float, optional
                Timeout of the request in seconds. Default is 8.

        -------
        Returns
        -------
            text : str
                Response body decoded as UTF-8 (invalid bytes replaced).

        ------
        Raises
        ------
            urllib.error.URLError
                If the request fails or times out (caught by the fetchers).
        """

        # ##########  MODULES & PACKAGES  ##########
        import urllib.request


        # ##########  PARAMETERS  ##########
        request = urllib.request.Request(
                                         url,
                                         headers={'User-Agent' : "Mozilla/5.0 (BinarAlt weather proxy)"}  # Some sites reject requests without a browser-like agent
                                        )


        # ##########  OUTPUT  ##########
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")  # Invalid bytes replaced instead of raising
    #
    def _blank(
               self,
               source_url : str
        ) -> dict[str, bool | str]:
        """
        Return an empty, normalized weather result.

        ----------
        Parameters
        ----------
            source_url : str
                Public weather page of the station.

        -------
        Returns
        -------
            result : dict[str, bool | str]
                A dictionary containing the following key-value pairs:
                    - 'ok' : bool
                        True once data have been read (False here).

                    - 'source_url' : str
                        Public weather page of the station (input value).

                    - 'temperature' : str
                        Air temperature with units, e.g. '8.5 °C' ('-' if
                        not available).

                    - 'humidity' : str
                        Relative humidity with units, e.g. '35.0 %' ('-' if
                        not available).

                    - 'wind_speed' : str
                        Wind speed in m/s, e.g. '6.2 m/s' ('-' if not
                        available).

                    - 'wind_direction' : str
                        Direction the wind blows FROM, e.g. '225° (SW)' ('-'
                        if not available).

                    - 'precipitation' : str
                        Rain / precipitation state or amount ('-' if not
                        available).

                    - 'dust' : str
                        Dust concentration in µg/m³ ('-' if not available).

                    - 'seeing' : str
                        Seeing in arcsec, e.g. '0.85"' ('-' if not
                        available).
        """

        # ##########  OUTPUT  ##########
        return {
                'ok'             : False,
                'source_url'     : source_url,
                'temperature'    : "-",
                'humidity'       : "-",
                'wind_speed'     : "-",
                'wind_direction' : "-",
                'precipitation'  : "-",
                'dust'           : "-",
                'seeing'         : "-",
               }
    #
    def _compass(
                 self,
                 deg : object
        ) -> str:
        """
        Convert a numerical wind direction into meteorological notation, e.g. '225° (SW)'.

        ----------
        Parameters
        ----------
            deg : object
                Wind direction in degrees (number or numeric string); 0 deg / N
                means wind blowing FROM the north.

        -------
        Returns
        -------
            direction : str
                '{deg}° ({16-point compass label})', or '-' if 'deg' is not a
                finite number.
        """

        # ##########  MODULES & PACKAGES  ##########
        import math


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        try:
            deg = float(deg)
        except (TypeError, ValueError):
            return "-"
        if not math.isfinite(deg):
            return "-"


        # ##########  COMPUTATION  ##########
        directions = [
                      "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                      "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"
                     ]
        index      = int((deg % 360) / 22.5 + 0.5) % 16  # 16-point compass rose, 22.5 deg per sector


        # ##########  OUTPUT  ##########
        return f"{deg % 360:.0f}° ({directions[index]})"
    #
    @staticmethod
    def _is_number(
                   value : object
        ) -> bool:
        """
        Check whether a value is a finite int or float (bool excluded).

        ----------
        Parameters
        ----------
            value : object
                Value to check.

        -------
        Returns
        -------
            is_number : bool
                True for finite int or float values.
        """

        # ##########  MODULES & PACKAGES  ##########
        import math


        # ##########  OUTPUT  ##########
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
    #
    @staticmethod
    def _search(
                page    : str,
                pattern : str
        ) -> str | None:
        """
        Extract the first group of a regular expression from a web page.

        ----------
        Parameters
        ----------
            page : str
                Text of the web page.

            pattern : str
                Regular expression with one capturing group.

        -------
        Returns
        -------
            value : str | None
                Stripped content of the group, or None if there is no match.
        """

        # ##########  MODULES & PACKAGES  ##########
        import re


        # ##########  OUTPUT  ##########
        match = re.search(pattern, page)
        return match.group(1).strip() if match else None
    #
    @staticmethod
    def _feed_number(
                     states : list[str],
                     index  : int
        ) -> float | None:
        """
        Read a number of a whitespace-separated status feed.

        ----------
        Parameters
        ----------
            states : list[str]
                Values of the feed.

            index : int
                Position of the value in the feed.

        -------
        Returns
        -------
            value : float | None
                Finite value, or None if missing or not numeric.
        """

        # ##########  MODULES & PACKAGES  ##########
        import math


        # ##########  OUTPUT  ##########
        try:
            value = float(states[index])
        except (IndexError, ValueError):
            return None
        return value if math.isfinite(value) else None
    #
    def _tng_last(
                  self,
                  field : str
        ) -> float | None:
        """
        Return the last numeric value of a TNG weather time series.

        ----------
        Parameters
        ----------
            field : str
                API field name (e.g. 'temperature', 'windspeed').

        -------
        Returns
        -------
            value : float | None
                Last finite value of the series, or None if unavailable.
        """

        # ##########  MODULES & PACKAGES  ##########
        import json


        # ##########  OUTPUT  ##########
        try:
            data = json.loads(self._fetch(self._TNG_URL.format(field=field)))
            for _, value in reversed(data['series'][0]['data']):  # (timestamp, value) pairs, latest last
                if self._is_number(value):
                    return value
        except Exception:  # Offline station or changed API: the field is left empty
            pass
        return None


    # ####################  QUERY  ####################
    def _get(
             self,
             station : str,
             force   : bool = False
        ) -> dict[str, bool | str]:
        """
        Get the weather information of a station.

        ----------
        Parameters
        ----------
            station : str
                Station code (case-insensitive): 'song', 'stella', 'speculoos',
                'not', 'tng', 'gtc', 'int', 'magic' or 'liverpool'.

            force : bool, optional
                If True (manual refresh button of BinarAlt's WEATHER panel),
                the cache is bypassed and a fresh reading is fetched.
                Default is False.

        -------
        Returns
        -------
            payload : dict[str, bool | str]
                Normalized result (see '_blank'), plus an 'error' key (str)
                when the station is unknown or the fetcher fails.

        -----
        Notes
        -----
            - Results are cached for 'cache_ttl' seconds.

            - Invalid inputs never raise: the error is reported in the
              payload, so that the browser request is always answered.
        """

        # ##########  MODULES & PACKAGES  ##########
        import time


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        station = str(station or '').strip().lower()  # None or empty -> '' (unknown station)
        fetcher = self.fetchers.get(station)  # None for unknown stations
        force   = bool(force)

        if fetcher is None:
            return {
                    'ok'    : False,
                    'error' : f"Unknown weather station '{station}'. Available: {list(self.fetchers.keys())}.",
                   }


        # ##########  COMPUTATION  ##########
        # ##### CACHE CHECK #####
        now = time.monotonic()  # Monotonic clock: unaffected by changes of the system time

        if not force:
            with self._lock:  # The cache is shared by the server threads
                cached = self.cache.get(station)
            if cached is not None:
                cached_time, cached_payload = cached
                if now - cached_time < self.cache_ttl:  # Recent reading: reused
                    return cached_payload


        # ##### FETCH #####
        try:
            payload = fetcher()
        except Exception as e:  # A fetcher must never break the HTTP response
            payload          = self._blank('')  # Every field '-'
            payload['error'] = str(e)


        # ##### UPDATE CACHE #####
        with self._lock:
            self.cache[station] = (now, payload)  # Failed readings are also cached, so a broken site is not requested continuously


        # ##########  OUTPUT  ##########
        return payload


    # ####################  OBSERVATORY-SPECIFIC FETCHERS  ####################
    def _fetch_ocan(
                    self,
                    station : str
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the SONG, STELLA or SPECULOOS stations (Observatorio del Teide).

        ----------
        Parameters
        ----------
            station : str
                'song', 'stella' or 'speculoos'.

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - The three stations share the same JSON API (ocan.iac.es); only
              the fields shown for each of them differ (see '_OCAN_FIELDS').

            - The API gives dust in mg/m³, converted here to µg/m³, and a
              seeing of 0 when there is no measurement.
        """

        # ##########  MODULES & PACKAGES  ##########
        import json


        # ##########  COMPUTATION  ##########
        # ##### READ #####
        source_url = "https://ocan.iac.es/meteo/OT/"
        result     = self._blank(source_url)

        try:
            data = json.loads(self._fetch(self._OCAN_URL.format(station=station)))
        except Exception:  # Unreachable station: empty result
            return result

        result['ok'] = True
        fields       = self._OCAN_FIELDS[station]  # Only the fields shown on the station page


        # ##### NORMALIZATION #####
        if 'temperature' in fields and self._is_number(data.get('temperature')):
            result['temperature'] = f"{data['temperature']:.1f} °C"

        if 'humidity' in fields and self._is_number(data.get('humidity')):
            result['humidity'] = f"{data['humidity']:.1f} %"

        if 'wind_speed' in fields and self._is_number(data.get('wind_speed')):
            result['wind_speed'] = f"{data['wind_speed']:.1f} m/s"

        if 'wind_direction' in fields and self._is_number(data.get('wind_direction')):
            result['wind_direction'] = self._compass(data['wind_direction'])

        if 'rain' in fields and self._is_number(data.get('rain')):
            result['precipitation'] = f"{data['rain']:.1f} mm"

        if 'seeing' in fields and self._is_number(data.get('seeing')) and data['seeing'] != 0:  # 0 = no measurement
            result['seeing'] = f"{data['seeing']:.2f}\""

        if 'dust' in fields and self._is_number(data.get('dust')):
            result['dust'] = f"{data['dust'] * 1000:.1f} µg/m³"  # mg/m³ -> µg/m³


        # ##########  OUTPUT  ##########
        return result
    #
    def _fetch_not(
                   self
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the Nordic Optical Telescope (ORM).

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - The public page is an iframe shell; the actual readings are in
              latest_not.php. The seeing is the ING DIMM and the dust the TNG
              aerosol monitor, as shown on the NOT page.
        """

        # ##########  COMPUTATION  ##########
        # ##### READ #####
        source_url = "https://www.not.iac.es/weather/"
        result     = self._blank(source_url)

        try:
            page = self._fetch("https://www.not.iac.es/weather/latest_not.php?c=")
        except Exception:  # Unreachable station: empty result
            return result

        result['ok'] = True


        # ##### NORMALIZATION #####
        value = self._search(page, r"Temperature:</td>\s*<td[^>]*><b>([\d.\-]+)</b>")  # Number in bold after the label cell
        if value:
            result['temperature'] = f"{value} °C"

        value = self._search(page, r"Humidity:</td>\s*<td[^>]*><b>([\d.\-]+)</b>")
        if value:
            result['humidity'] = f"{value} %"

        value = self._search(page, r"Wind Speed:</td>\s*<td[^>]*><b>([\d.\-]+)</b>")
        if value:
            result['wind_speed'] = f"{value} m/s"

        value = self._search(page, r"Wind Dir:</td>\s*<td[^>]*><b>([\d.\-]+)</b>")
        if value:
            result['wind_direction'] = self._compass(value)

        value = self._search(page, r"Rain:</td>\s*<td[^>]*><b>([^<]+)</b>")
        if value:
            result['precipitation'] = "Wet" if value == "True" else "Dry"  # The page gives the rain sensor state as True / False

        value = self._search(page, r"DIMM \(<a[^>]*>ING</a>\)</td>\s*<td[^>]*><b>([\d.\-]+)</b>")
        if value:
            result['seeing'] = f"{value}\""

        value = self._search(page, r"Aerosols \(<a[^>]*>TNG</a>\)</td>\s*<td[^>]*><b>([\d.\-]+)</b>")
        if value:
            result['dust'] = f"{value} µg/m³"


        # ##########  OUTPUT  ##########
        return result
    #
    def _fetch_tng(
                   self
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the Telescopio Nazionale Galileo (ORM).

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - Each field is obtained independently from its time-series
              endpoint (see '_tng_last'), because the summary endpoint can
              fail when the station is offline.
        """

        # ##########  PARAMETERS  ##########
        source_url = "https://tngweb.tng.iac.es/weather/"
        result     = self._blank(source_url)
        formats = {  # Result key -> (API field, formatter)
                   'temperature'    : ("temperature", lambda v: f"{v:.1f} °C"),
                   'humidity'       : ("humidity", lambda v: f"{v:.1f} %"),
                   'wind_speed'     : ("windspeed", lambda v: f"{v:.1f} m/s"),
                   'wind_direction' : ("winddir", self._compass),
                   'dust'           : ("dust", lambda v: f"{v:.1f} µg/m³"),
                   'seeing'         : ("seeing", lambda v: f"{v:.2f}\""),
                  }


        # ##########  COMPUTATION  ##########
        got_any = False  # True if at least one field could be read
        for key, (field, formatter) in formats.items():
            value = self._tng_last(field)  # One request per field
            if value is not None:
                result[key] = formatter(value)
                got_any     = True

        result['ok'] = got_any  # The station counts as reachable if any field was read


        # ##########  OUTPUT  ##########
        return result
    #
    def _fetch_gtc(
                   self
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the Gran Telescopio Canarias (ORM).

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - Weather values are contained directly in the server-rendered
              HTML, already formatted with units.

            - The wind direction is read from the CSS class of the wind
              indicator (e.g. 'NE'), and the seeing from the ING DIMM row.
        """

        # ##########  MODULES & PACKAGES  ##########
        import re


        # ##########  PARAMETERS  ##########
        source_url = "https://atmosportal.gtc.iac.es/"
        result     = self._blank(source_url)
        row        = r"<td>{label}</td><td></td><td[^>]*>(?:<a[^>]*>)?([^<]+)</td>"  # Value cell of a labeled table row


        # ##########  COMPUTATION  ##########
        # ##### READ #####
        try:
            page = self._fetch(source_url)
        except Exception:  # Unreachable station: empty result
            return result

        result['ok'] = True


        # ##### NORMALIZATION #####
        value = self._search(page, row.format(label=re.escape("Outdoor Temp")))
        if value:
            result['temperature'] = value

        value = self._search(page, row.format(label=re.escape("Outdoor Humidity")))
        if value:
            result['humidity'] = value

        value = self._search(page, row.format(label=re.escape("Wind Mean")))
        if value:
            result['wind_speed'] = value

        value = self._search(page, row.format(label=re.escape("Rain")))
        if value:
            result['precipitation'] = value

        directions = re.findall(r"id='windDir' class='([A-Za-z]+)'", page)  # Compass label in the CSS class of the wind indicator
        if directions:
            result['wind_direction'] = directions[-1]  # Last indicator: current reading

        value = self._search(page, row.format(label=re.escape("Internal PM10 ")))  # The label has a trailing space on the page
        if value:
            result['dust'] = value.replace("ug/m3", "µg/m³")

        value = self._search(page, r"ING DIMM</a></td>\s*<td[^>]*>[^<]*</td>\s*<td[^>]*>\s*([\d.\-]+)")
        if value:
            result['seeing'] = f"{value}\""


        # ##########  OUTPUT  ##########
        return result
    #
    def _fetch_int(
                   self
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the Isaac Newton Telescope / ING (ORM).

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - The weather page is server-rendered and contains the current
              measurements, already formatted with units, directly in the
              HTML. The wind speed is converted from km/h to m/s.
        """

        # ##########  MODULES & PACKAGES  ##########
        import html
        import re


        # ##########  PARAMETERS  ##########
        source_url = "https://astro.ing.iac.es/weather/"
        result     = self._blank(source_url)
        row        = r"<B>{label}</B>.{{0,20}}?<div class=\"data\">([^<]+)"  # Value that follows a bold label


        # ##########  COMPUTATION  ##########
        # ##### READ #####
        try:
            page = self._fetch(source_url)
        except Exception:  # Unreachable station: empty result
            return result

        result['ok'] = True
        values       = {}
        for label in ('Temperature', 'Humidity', 'Wind Speed', 'Wind Direction', 'Rain', 'DIMM Seeing'):  # Keys of 'values'
            value         = self._search(page, row.format(label=re.escape(label)))
            values[label] = html.unescape(value).strip() if value else None  # E.g. '&deg;' -> '°'


        # ##### NORMALIZATION #####
        if values['Temperature']:
            result['temperature'] = values['Temperature']

        if values['Humidity']:
            result['humidity'] = values['Humidity']

        if values['Wind Speed']:
            match                = re.match(r"([\d.\-]+)\s*km/h", values['Wind Speed'])
            result['wind_speed'] = f"{float(match.group(1)) / 3.6:.1f} m/s" if match else values['Wind Speed']  # km/h -> m/s

        if values['Wind Direction']:
            result['wind_direction'] = values['Wind Direction']

        if values['Rain']:
            result['precipitation'] = values['Rain']

        if values['DIMM Seeing']:
            result['seeing'] = values['DIMM Seeing']


        # ##########  OUTPUT  ##########
        return result
    #
    def _fetch_magic(
                     self
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the MAGIC telescopes (ORM).

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - The weather page is server-rendered; fields shown as 'n/a' are
              treated as missing. The wind speed is converted from km/h to
              m/s.

            - The wind direction is shown as an image / CSS element rather
              than a number, so it is not reported.
        """

        # ##########  MODULES & PACKAGES  ##########
        import re


        # ##########  PARAMETERS  ##########
        source_url = "http://www.magic.iac.es/site/weather/index.html"
        result     = self._blank(source_url)
        row        = r">{label}:</a></td>\s*<td[^>]*><b>([^<]+)</b>"  # Value that follows a linked label


        # ##########  COMPUTATION  ##########
        # ##### READ #####
        try:
            page = self._fetch(source_url)
        except Exception:  # Unreachable station: empty result
            return result

        result['ok'] = True
        values       = {}
        for label in ('Temperature', 'Humidity', 'Wind Speed', 'Rain', 'TNG Dust'):  # Keys of 'values'
            value         = self._search(page, row.format(label=re.escape(label)))
            values[label] = None if not value or value == "n/a" else value  # 'n/a' on the page: missing reading


        # ##### NORMALIZATION #####
        if values['Temperature']:
            result['temperature'] = f"{values['Temperature']} °C"

        if values['Humidity']:
            result['humidity'] = f"{values['Humidity']} %"

        if values['Wind Speed']:
            try:
                result['wind_speed'] = f"{float(values['Wind Speed']) / 3.6:.1f} m/s"  # km/h -> m/s
            except ValueError:  # Non-numeric reading: the field is left empty
                pass

        if values['Rain']:
            result['precipitation'] = values['Rain']

        if values['TNG Dust']:
            result['dust'] = f"{values['TNG Dust']} µg/m³"


        # ##########  OUTPUT  ##########
        return result
    #
    def _fetch_liverpool(
                         self
        ) -> dict[str, bool | str]:
        """
        Fetch the weather of the Liverpool Telescope (ORM).

        -------
        Returns
        -------
            result : dict[str, bool | str]
                Normalized result (see '_blank').

        -----
        Notes
        -----
            - The dashboard obtains its data from the V2/met_now status
              feed, a whitespace-separated list of values: rain state (2),
              wind speed in m/s (6), wind bearing (7), temperature in °C (8)
              and relative humidity as a fraction (10).

            - The feed gives the direction the wind blows TO, so 180 deg are
              added to follow the meteorological convention (FROM).
        """

        # ##########  COMPUTATION  ##########
        # ##### READ #####
        source_url = "https://telescope.ljmu.ac.uk/Live/v2_index.php"
        result     = self._blank(source_url)

        try:
            states = self._fetch("https://telescope.ljmu.ac.uk/Live/V2/met_now").split()
        except Exception:  # Unreachable station: empty result
            return result

        result['ok'] = True


        # ##### NORMALIZATION #####
        if len(states) > 2 and states[2]:  # Rain state (e.g. 'dry', 'wet')
            result['precipitation'] = states[2].capitalize()

        wind_speed = self._feed_number(states, 6)
        if wind_speed is not None:
            result['wind_speed'] = f"{wind_speed:.1f} m/s"

        bearing = self._feed_number(states, 7)
        if bearing is not None:
            result['wind_direction'] = self._compass((bearing + 180) % 360)  # TO -> FROM

        temperature = self._feed_number(states, 8)
        if temperature is not None:
            result['temperature'] = f"{temperature:.1f} °C"

        humidity = self._feed_number(states, 10)
        if humidity is not None:
            result['humidity'] = f"{humidity * 100:.0f} %"  # Fraction (0-1) -> percent


        # ##########  OUTPUT  ##########
        return result




# ==============================
class _BinarAltRequestHandler(http.server.SimpleHTTPRequestHandler):
    """
    HTTP request handler of the BinarAlt local server.

    Serves the static files of the BinarAlt interface (binaralt.html and its
    user guide, BinarAlt_UserGuide.pdf) and a '/weather' JSON endpoint backed
    by the _WeatherService instance stored in 'weather_service' (set as a
    class attribute before the server starts).

    ----------
    Parameters
    ----------
        *args, **kwargs
            Arguments of http.server.SimpleHTTPRequestHandler (request,
            client address, server and the 'directory' served), passed by
            the HTTP server for every request.

    -----
    Notes
    -----
        - Endpoint: GET /weather?station=<station>[&force=1]; 'force=1'
          bypasses the cache (manual refresh button).

        - The server only listens on 127.0.0.1, so it is not reachable from
          other machines.

        - The default per-request access log is silenced.
    """

    # ####################  CONFIGURATION  ####################
    weather_service : _WeatherService | None = None


    # ####################  HTTP HANDLING  ####################
    # Overrides BaseHTTPRequestHandler.log_message; name fixed by the library
    def log_message(
                    self,
                    format : str,
                    *args  : object
        ) -> None:
        """
        Silence the default per-request access log to stderr.

        ----------
        Parameters
        ----------
            format : str
                Format string of the log message (name imposed by the base
                class).

            *args : object
                Values of the log message.

        -------
        Returns
        -------
            None
                Nothing is logged.
        """
        pass
    #
    # Overrides BaseHTTPRequestHandler.do_GET; name fixed by the library
    def do_GET(
               self
        ) -> None:
        """
        Route '/weather' requests to the weather service and everything else to the static file server.

        -------
        Returns
        -------
            None
                The response is sent to the browser.
        """

        # ##########  MODULES & PACKAGES  ##########
        import urllib.parse


        # ##########  COMPUTATION  ##########
        parsed = urllib.parse.urlsplit(self.path)  # Path and query string of the request

        if parsed.path == "/weather":
            self._serve_weather(urllib.parse.parse_qs(parsed.query))  # E.g. {'station': ['tng'], 'force': ['1']}
        else:
            super().do_GET()  # Static files (binaralt.html, user guide...)
    #
    def _send_json(
                   self,
                   payload : dict,
                   status  : int = 200
        ) -> None:
        """
        Send a JSON response with no-cache headers.

        ----------
        Parameters
        ----------
            payload : dict
                JSON-serializable content of the response.

            status : int, optional
                HTTP status code. Default is 200.

        -------
        Returns
        -------
            None
                The JSON response is written to the client.
        """

        # ##########  MODULES & PACKAGES  ##########
        import json


        # ##########  OUTPUT  ##########
        body = json.dumps(payload).encode("utf-8")  # Bytes of the response body

        self.send_response(status)  # Status line (e.g. 200 OK)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))  # Length in bytes, not characters
        self.send_header("Cache-Control", "no-store")  # Every refresh must reach the weather service
        self.end_headers()  # Blank line between headers and body
        self.wfile.write(body)
    #
    def _serve_weather(
                       self,
                       query : dict[str, list[str]]
        ) -> None:
        """
        Serve /weather?station=<station>[&force=1] as JSON.

        ----------
        Parameters
        ----------
            query : dict[str, list[str]]
                Parsed query string (urllib.parse.parse_qs).

        -------
        Returns
        -------
            None
                The weather payload (or an error) is sent as JSON.

        -----
        Notes
        -----
            - 'force=1' (or 'force=true') bypasses the cache (see
              '_WeatherService._get'), as used by BinarAlt's manual refresh
              button.
        """

        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        if self.weather_service is None:  # 'binaralt' did not set the service before starting the server
            self._send_json(
                            {
                             'ok'    : False,
                             'error' : "Weather service not configured.",
                            },
                            status=503
                           )
            return


        # ##########  PARAMETERS  ##########
        station = (query.get('station') or [''])[0].strip().lower()  # First value of the 'station' parameter ('' if absent)
        force   = (query.get('force') or [''])[0].strip().lower() in ("1", "true")  # Manual refresh: bypass the cache


        # ##########  OUTPUT  ##########
        try:
            payload = self.weather_service._get(station, force=force)
        except Exception as e:  # Never leave the browser request hanging
            self._send_json(
                            {
                             'ok'    : False,
                             'error' : str(e),
                            },
                            status=500
                           )
            return

        self._send_json(payload)  # 200 OK, also when a station could not be read (see 'ok' in the payload)
