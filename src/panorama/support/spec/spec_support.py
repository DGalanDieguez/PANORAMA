"""
Internal helpers of the 'spec' module: input validation shared by the public
functions, spectrum file formats and readers ('readspec', 'writespec'), the
check buttons of 'plotspec', the interactive spectrum cleaner ('rmvspec') and
the continuum fitting machinery ('normspec').

Conventions shared by every function in this module:
    - Functions, classes and constants are private (leading underscore);
      they are used by the public functions of 'panorama.spec' (and
      '_check_array' by the 'prismas' helpers) and are not exposed at the
      package level.
    - Spectrum files are read and written in the formats given by their
      extension ('_ASCII_EXTENSIONS', '_FITS_EXTENSIONS'), optionally
      followed by '.gz'.
    - Spectra are pairs (wave, flux) of 1D float NumPy arrays with the same
      length. Wavelengths are in the units of the input spectrum (usually Å)
      and are never converted.
    - File arguments are paths (path + filename), given either as strings
      or as pathlib.Path objects. Output directories are never created.
    - Output spectra are written with 'panorama.spec.writespec'.
    - Scalar real arguments are validated with '_check_real'
      ('support/units/units_support.py'), and the file arguments of
      'panorama.spec' with '_check_input_file' / '_check_output_file'
      ('support/tex/tex_support.py'): these helpers are shared with the
      other modules instead of being duplicated here.
    - Imports are done inside each function, so the type hints are not
      evaluated when the functions are defined.
    - Methods of the classes are grouped by level-1 headers, so their own
      body sections use level-2 and level-3 headers.
    - Interactive tools print their messages to the terminal and never raise
      exceptions inside Matplotlib callbacks.
    - All functions, class constructors and the methods called from outside
      their class validate their inputs before computing; nested functions
      and Matplotlib callbacks do not (see the PANORAMA conventions).
"""

from __future__ import annotations

# ==================================================================
# ============================== INIT ==============================
# ==================================================================

# ####################  SPECTRUM FILES  ####################
"""
- Extensions (lowercase, without '.gz') of the ASCII and FITS spectrum files
  read by 'readspec' and written by 'writespec'
"""
_ASCII_EXTENSIONS = (".txt", ".dat", ".asc", ".ascii")
_FITS_EXTENSIONS  = (".fits", ".fit", ".fts")


# ####################  FITS SPECTRA  ####################
"""
- Names of the wavelength and flux columns of FITS tables (case-insensitive),
  and upper limit of ln(lambda): any ln(lambda) is below it, and any lambda in
  Å or nm is above it
"""
_FITS_WAVE_COLS = ("wave", "wavelength", "lambda", "lam")
_FITS_FLUX_COLS = ("flux", "intensity", "flx", "spec", "flambda", "counts", "cts")
_LN_WAVE_MAX    = 20.


# ####################  MATPLOTLIB  ####################
"""
- Backends without a window: they never receive mouse / keyboard events
"""
_NON_INTERACTIVE_BACKENDS = ("agg", "cairo", "pdf", "pgf", "ps", "svg", "template")


# ####################  RMVSPEC  ####################
"""
- Keyboard help of 'rmvspec', printed when the window opens and with 'h'
"""
_RMVSPEC_HELP = """\
      \nrmvspec -- keyboard shortcuts (they take priority over Matplotlib's)
            mouse drag : select a span (red)
            r : remove selected spans          x : discard last selection
            c : clear all selections           a : restore all removed points
            u : undo                           U : redo
            g : toggle grid                    l : show / hide OB-type lines
            h : print this help
            s : save cleaned spectrum          S : save spectrum + PNG figure
            d : export removed spans           D : export pending selected spans
            Other Matplotlib keys still work: o (zoom), p (pan), q (close), arrows...
        """


# ####################  NORMSPEC  ####################
"""
- Continuum functions and output spectra accepted by 'normspec', and maximum
  number of sample points drawn (keeps the plot responsive)
"""
_FUNCTIONS    = ("spline3", "spline1", "legendre", "chebyshev")
_OUTPUT_TYPES = ("ratio", "difference", "fit", "data")
_MAX_MARKERS  = 20000


# ####################  NORMSPEC KEYS  ####################
"""
- Keys of the interactive window of 'normspec' that change the continuum
  function
"""
_FUNCTION_KEYS = {
                  '1' : "spline1",
                  '3' : "spline3",
                  'l' : "legendre",
                  'c' : "chebyshev",
                 }


# ####################  NORMSPEC HELP  ####################
"""
- Controls of 'normspec', printed when the window opens and with '?'
"""
_NORMSPEC_HELP = """\
        \nnormspec -- controls (Matplotlib's keyboard shortcuts are disabled)
            Automatic continuum:
                f       : fit the continuum with the current parameters
                + / -   : increase / decrease the order
                1 / 3   : linear / cubic spline
                l / c   : Legendre / Chebyshev polynomial
                s s     : mark a sample region (press 's' at both ends)
                t       : reset the sample to the whole spectrum
                d / u   : delete / undelete the nearest sample point
            Manual continuum:
                left-click           : add a continuum point
                right-click on point : remove it
                enter                : cubic spline through the manual points
            General:
                r : reset    w : save    q : quit    ? : this help
        """




# ====================================================================
# ============================== COMMON ==============================
# ====================================================================

# ==============================
def _check_array(
                 value : object,
                 name  : str,
                 label : str,
                 kind  : str = "float",
                 size  : int | None = None
    ) -> None:
    """
    Check that an argument is a 1D NumPy array of the expected kind (and, optionally, size).

    ----------
    Parameters
    ----------
        value : object
            Argument to check.

        name : str
            Name of the argument in the error messages (e.g. 'wave').

        label : str
            Physical meaning of the argument in the error messages (e.g.
            'Wavelength array').

        kind : str, optional
            Expected kind of the elements: 'float' (any real numeric type,
            booleans excluded), 'int' (integers) or 'bool'. Default is
            "float".

        size : int | None, optional
            Expected number of elements, or None (any size). Default is None.

    -------
    Returns
    -------
        None
            Nothing is returned; an exception is raised if the check fails.

    ------
    Raises
    ------
        TypeError
            If 'value' is not a NumPy array, its elements are not of the
            expected kind, or 'name', 'label', 'kind' or 'size' have the
            wrong type.

        ValueError
            If 'value' is not 1D or does not have 'size' elements, 'name' or
            'label' are empty, 'kind' is unknown or 'size' is negative.

    -----
    Notes
    -----
        - Shared by the 'spec' helpers that receive spectra and by the
          'prismas' helpers that receive whole catalogs ('_passes',
          '_plot_needed_models'), so that the same argument gives the same
          error messages.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(name, str):
        raise TypeError(f"Error-message text 'name' must be a string. Got {type(name).__name__}.")
    elif not name.strip():
        raise ValueError(f"Error-message text 'name' must be a non-empty string. Got '{name}'.")
    #
    if not isinstance(label, str):
        raise TypeError(f"Error-message text 'label' must be a string. Got {type(label).__name__}.")
    elif not label.strip():
        raise ValueError(f"Error-message text 'label' must be a non-empty string. Got '{label}'.")
    #
    if not isinstance(kind, str):
        raise TypeError(f"Array kind 'kind' must be a string. Got {type(kind).__name__}.")
    elif kind not in ("float", "int", "bool"):
        raise ValueError(f"Array kind 'kind' must be 'float', 'int' or 'bool'. Got '{kind}'.")
    #
    if size is not None:
        if not isinstance(size, numbers.Integral) or isinstance(size, bool):
            raise TypeError(f"Array size 'size' must be an integer or None. Got {type(size).__name__}.")
        elif size < 0:
            raise ValueError(f"Array size 'size' must be >= 0. Got {size}.")


    # ####################  COMPUTATION  ####################
    if not isinstance(value, np.ndarray):
        raise TypeError(f"{label} '{name}' must be a 1D NumPy array. Got {type(value).__name__}.")
    elif value.ndim != 1:
        raise ValueError(f"{label} '{name}' must be 1D. Got shape {value.shape}.")
    #
    if kind == "bool":
        valid_kind = value.dtype == bool
    elif kind == "int":
        valid_kind = np.issubdtype(value.dtype, np.integer)  # Any signed or unsigned NumPy integer
    else:
        valid_kind = np.issubdtype(value.dtype, np.number) and value.dtype != bool and not np.issubdtype(value.dtype, np.complexfloating)  # Real numbers only
    if not valid_kind:
        raise TypeError(f"{label} '{name}' must contain {kind} values. Got dtype {value.dtype}.")
    #
    if size is not None and value.size != size:
        raise ValueError(f"{label} '{name}' must contain {size} elements. Got {value.size}.")




# ==============================
def _check_wave_lim(
                    wave_lim : tuple[float, float] | None
    ) -> None:
    """
    Validate the optional wavelength range (min, max) of the public 'spec' functions.

    ----------
    Parameters
    ----------
        wave_lim : tuple[float, float] | None
            Wavelength range (min, max), in the units of the spectrum, or
            None (no range).

    -------
    Returns
    -------
        None
            Nothing is returned; an exception is raised if the check fails.

    ------
    Raises
    ------
        TypeError
            If 'wave_lim' is not a tuple / list or None, or its elements are
            not real numbers.

        ValueError
            If 'wave_lim' does not have two finite elements with min < max.

    -----
    Notes
    -----
        - Shared by 'plotspec', 'rmvspec' and 'normspec', so that the same
          argument gives the same error messages.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if wave_lim is None:
        return
    elif not isinstance(wave_lim, (tuple, list)):
        raise TypeError(f"Wavelength range 'wave_lim' must be a (min, max) tuple or list, or None. Got {type(wave_lim).__name__}.")
    elif len(wave_lim) != 2:
        raise ValueError(f"Wavelength range 'wave_lim' must contain exactly two elements (min, max). Got {len(wave_lim)}.")
    elif not all(isinstance(w, numbers.Real) and not isinstance(w, bool) for w in wave_lim):
        raise TypeError(f"Wavelength range 'wave_lim' must contain real numbers. Got {tuple(type(w).__name__ for w in wave_lim)}.")
    elif not all(np.isfinite(w) for w in wave_lim):
        raise ValueError(f"Wavelength range 'wave_lim' must contain finite values. Got {tuple(wave_lim)}.")
    elif wave_lim[0] >= wave_lim[1]:
        raise ValueError(f"Wavelength range 'wave_lim' must satisfy min < max. Got {tuple(wave_lim)}.")




# ==============================
def _prepare_spectrum(
                      wave             : np.ndarray[float],
                      flux             : np.ndarray[float],
                      wave_lim         : tuple[float, float] | None,
                      keep_flux_nonfin : bool
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Prepare a spectrum for an interactive tool: discard non-finite points, sort it by wavelength and cut it to the wavelength range.

    ----------
    Parameters
    ----------
        wave : np.ndarray[float]
            Wavelength array (1D).

        flux : np.ndarray[float]
            Flux array (1D, same length as 'wave').

        wave_lim : tuple[float, float] | None
            Wavelength range (min, max) kept, already validated with
            '_check_wave_lim', or None (whole spectrum).

        keep_flux_nonfin : bool
            If True, points with non-finite fluxes are kept (only non-finite
            wavelengths are discarded, as in 'normspec'); if False, they are
            discarded too (as in 'rmvspec').

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelengths of the kept points, sorted.

        flux : np.ndarray[float]
            Fluxes of the kept points.

    ------
    Raises
    ------
        TypeError
            If 'wave' or 'flux' are not 1D float arrays, 'wave_lim' is not a
            valid range type or 'keep_flux_nonfin' is not a boolean.

        ValueError
            If the arrays have different lengths, 'wave_lim' is not a valid
            range, or fewer than two points remain.

    -----
    Notes
    -----
        - Discarded non-finite points issue a RuntimeWarning that points to
          the caller of the public function.

        - The sort is stable: points with equal wavelengths keep their order.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(wave, "wave", "Wavelength array")
    #
    _check_array(flux, "flux", "Flux array", size=wave.size)
    #
    _check_wave_lim(wave_lim)
    #
    if not isinstance(keep_flux_nonfin, bool):
        raise TypeError(f"Non-finite flux flag 'keep_flux_nonfin' must be a boolean. Got {type(keep_flux_nonfin).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  NON-FINITE POINTS  ##########
    keep = np.isfinite(wave) if keep_flux_nonfin else np.isfinite(wave) & np.isfinite(flux)  # Non-finite fluxes can be neither drawn nor selected in 'rmvspec'
    if not keep.all():
        warnings.warn(f"{int((~keep).sum())} non-finite point(s) discarded from the spectrum 'spec'.", RuntimeWarning, stacklevel=3)  # stacklevel=3: points to the caller of the public function
        wave, flux = wave[keep], flux[keep]


    # ##########  WAVELENGTH ORDER  ##########
    if np.any(np.diff(wave) < 0):  # Sorted only if needed
        sort_idx   = np.argsort(wave, kind="stable")  # Stable: points with equal wavelengths keep their order
        wave, flux = wave[sort_idx], flux[sort_idx]


    # ##########  WAVELENGTH RANGE  ##########
    if wave_lim is not None:
        sel        = (wave >= wave_lim[0]) & (wave <= wave_lim[1])  # Inclusive limits
        wave, flux = wave[sel], flux[sel]
    if wave.size < 2:  # Too few points to draw, clean or fit
        points = "points" if keep_flux_nonfin else "finite points"
        raise ValueError(f"Spectrum 'spec' must contain at least 2 {points} inside the wavelength range 'wave_lim'. Got {wave.size}.")


    # ####################  OUTPUT  ####################
    return wave, flux




# ==============================
def _load_spectrum(
                   spec : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
                   name : str = "spec"
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Return the (wave, flux) arrays of a spectrum given as a pair of arrays or as a spectrum file.

    ----------
    Parameters
    ----------
        spec : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path
            Spectrum: a (wave, flux) pair of 1D array-like objects with the
            same length, or the path to a spectrum file readable by
            'readspec'.

        name : str, optional
            Name of the argument in the error messages (e.g. 'spec[2]').
            Default is "spec".

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelength array (1D, float). Strictly increasing if read from a
            file; as given otherwise.

        flux : np.ndarray[float]
            Flux array (1D, float, same length as 'wave').

    ------
    Raises
    ------
        TypeError
            If 'spec' is not a pair or a path, its arrays are not numeric, or
            'name' is not a string.

        ValueError
            If the pair does not have two elements, its arrays are not 1D,
            have different lengths or are empty, or 'name' is empty.

    -----
    Notes
    -----
        - Spectrum files are read with 'readspec', which raises its own
          exceptions (see 'readspec').

        - Non-finite values are kept: each calling function decides how to
          handle them (discarding them or not).
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numpy as np
    import os


    # ##### PANORAMA #####
    from ...spec import readspec


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(name, str):
        raise TypeError(f"Error-message text 'name' must be a string. Got {type(name).__name__}.")
    elif not name.strip():
        raise ValueError(f"Error-message text 'name' must be a non-empty string. Got '{name}'.")
    #
    if not isinstance(spec, (tuple, list, str, os.PathLike)):
        raise TypeError(f"Spectrum '{name}' must be a (wave, flux) pair of 1D arrays or a path (str or pathlib.Path) to a spectrum file. Got {type(spec).__name__}.")
    elif isinstance(spec, (tuple, list)) and len(spec) != 2:
        raise ValueError(f"Spectrum '{name}' must be a (wave, flux) pair. Got {len(spec)} elements.")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: SPECTRUM FILE  ##########
    if isinstance(spec, (str, os.PathLike)):
        wave, flux = readspec(spec)  # Validated, finite and sorted wavelengths


    # ##########  OPTION 2: (WAVE, FLUX) PAIR  ##########
    else:
        try:
            wave = np.asarray(spec[0], dtype=float)  # Lists, tuples or arrays of any numeric type
            flux = np.asarray(spec[1], dtype=float)
        except (TypeError, ValueError) as e:  # Non-numeric or ragged sequences
            raise TypeError(f"Wavelength and flux arrays of the spectrum '{name}' must be numeric array-like objects. Got {type(spec[0]).__name__} and {type(spec[1]).__name__}.") from e
        #
        if wave.ndim != 1 or flux.ndim != 1:
            raise ValueError(f"Wavelength and flux arrays of the spectrum '{name}' must be 1D. Got shapes {wave.shape} and {flux.shape}.")
        elif wave.shape != flux.shape:
            raise ValueError(f"Wavelength and flux arrays of the spectrum '{name}' must have the same length. Got {wave.size} and {flux.size}.")
        elif wave.size == 0:
            raise ValueError(f"Spectrum '{name}' contains no data points.")


    # ####################  OUTPUT  ####################
    return wave, flux




# ==============================
def _merge_ranges(
                  pairs : list[tuple[float, float]]
    ) -> list[tuple[float, float]]:
    """
    Sort the (min, max) wavelength ranges and merge the overlapping ones.

    ----------
    Parameters
    ----------
        pairs : list[tuple[float, float]]
            Wavelength ranges (min, max). Each pair may be given in any order
            (it is sorted first).

    -------
    Returns
    -------
        merged : list[tuple[float, float]]
            Sorted, non-overlapping wavelength ranges (min, max).

    ------
    Raises
    ------
        TypeError
            If 'pairs' is not a list / tuple of (min, max) pairs of real
            numbers.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(pairs, (list, tuple)):
        raise TypeError(f"Wavelength ranges 'pairs' must be a list or tuple of (min, max) pairs. Got {type(pairs).__name__}.")
    for pair in pairs:
        if not (isinstance(pair, (tuple, list)) and len(pair) == 2 and all(isinstance(b, numbers.Real) and not isinstance(b, bool) for b in pair)):
            raise TypeError(f"Each wavelength range of 'pairs' must be a (min, max) pair of real numbers. Got {pair}.")


    # ####################  COMPUTATION  ####################
    merged = []
    for w1, w2 in sorted((min(a, b), max(a, b)) for a, b in pairs):  # Each pair as (min, max), sorted by its start
        if merged and w1 <= merged[-1][1]:  # Overlaps (or touches) the last merged range: extend it
            merged[-1] = (merged[-1][0], max(w2, merged[-1][1]))
        else:  # Disjoint: new range
            merged.append((w1, w2))


    # ####################  OUTPUT  ####################
    return merged




# ==============================
def _take_keyboard(
                   fig  : Figure,
                   keys : tuple[str, ...]
    ) -> tuple[bool, dict[str, list[str]]]:
    """
    Give the keys of an interactive tool priority over Matplotlib's shortcuts, for one figure only.

    ----------
    Parameters
    ----------
        fig : Figure
            Figure of the interactive tool.

        keys : tuple[str, ...]
            Keys used by the tool.

    -------
    Returns
    -------
        disconnected : bool
            True if Matplotlib's default key handler of the figure was
            disconnected (classic 'toolbar2' backends). The tool may then
            forward the keys it does not use with
            'matplotlib.backend_bases.key_press_handler'.

        keymaps : dict[str, list[str]]
            A dictionary keyed by toolmanager tool name; each value is its
            original keymap, to be restored with '_release_keyboard'. Empty if
            the backend does not use 'toolmanager'.

    ------
    Raises
    ------
        TypeError
            If 'fig' is not a matplotlib Figure or 'keys' is not a tuple of
            strings.

    -----
    Notes
    -----
        - No global Matplotlib setting (rcParams) is modified.

        - Backends with 'toolmanager' register their shortcuts per tool: the
          tool keys are removed from their keymaps until '_release_keyboard'
          is called.
    """

    # ####################  MODULES & PACKAGES  ####################
    from matplotlib.figure import Figure


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(fig, Figure):
        raise TypeError(f"Figure 'fig' must be a matplotlib Figure. Got {type(fig).__name__}.")
    #
    if not isinstance(keys, tuple) or not all(isinstance(k, str) for k in keys):
        raise TypeError(f"Tool keys 'keys' must be a tuple of strings. Got {keys!r}.")


    # ####################  COMPUTATION  ####################
    # ##########  CLASSIC TOOLBAR  ##########
    canvas       = fig.canvas
    manager      = getattr(canvas, "manager", None)  # None for figures without a window (e.g. Agg)
    disconnected = False
    handler_id   = getattr(manager, "key_press_handler_id", None)  # Connection id of Matplotlib's default key handler
    if handler_id is not None:
        canvas.mpl_disconnect(handler_id)  # Matplotlib no longer receives the keys directly
        disconnected = True


    # ##########  TOOLMANAGER  ##########
    keymaps     = {}
    toolmanager = getattr(manager, "toolmanager", None)  # Only set when rcParams['toolbar'] == 'toolmanager'
    if toolmanager is not None:
        for tool_name in list(toolmanager.tools):
            tool_keys = list(toolmanager.get_tool_keymap(tool_name) or [])  # Shortcuts of the tool
            kept      = [k for k in tool_keys if k not in keys]  # Shortcuts that do not clash with the tool keys
            if kept != tool_keys:
                keymaps[tool_name] = tool_keys  # Original keymap, restored by '_release_keyboard'
                toolmanager.update_keymap(tool_name, kept)


    # ####################  OUTPUT  ####################
    return disconnected, keymaps




# ==============================
def _release_keyboard(
                      fig     : Figure,
                      keymaps : dict[str, list[str]]
    ) -> None:
    """
    Restore the toolmanager keymaps modified by '_take_keyboard'.

    ----------
    Parameters
    ----------
        fig : Figure
            Figure of the interactive tool.

        keymaps : dict[str, list[str]]
            A dictionary keyed by toolmanager tool name; each value is the
            original keymap to restore (as returned by '_take_keyboard').

    -------
    Returns
    -------
        None
            The keymaps of the figure's toolmanager are restored in place.

    ------
    Raises
    ------
        TypeError
            If 'fig' is not a matplotlib Figure or 'keymaps' is not a
            dictionary.
    """

    # ####################  MODULES & PACKAGES  ####################
    from matplotlib.figure import Figure


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(fig, Figure):
        raise TypeError(f"Figure 'fig' must be a matplotlib Figure. Got {type(fig).__name__}.")
    #
    if not isinstance(keymaps, dict):
        raise TypeError(f"Original keymaps 'keymaps' must be a dictionary. Got {type(keymaps).__name__}.")


    # ####################  COMPUTATION  ####################
    manager     = getattr(fig.canvas, "manager", None)
    toolmanager = getattr(manager, "toolmanager", None)
    if toolmanager is not None:
        for tool_name, tool_keys in keymaps.items():
            try:
                toolmanager.update_keymap(tool_name, tool_keys)  # Original shortcuts back
            except Exception:
                pass  # The tool may no longer exist while the window closes




# ==============================
def _check_interactive_backend(
                               tool : str
    ) -> bool:
    """
    Warn if the current Matplotlib backend cannot receive mouse / keyboard events.

    ----------
    Parameters
    ----------
        tool : str
            Name of the interactive tool, used in the warning (e.g. 'rmvspec').

    -------
    Returns
    -------
        interactive : bool
            True if the backend is interactive; False (and a RuntimeWarning)
            otherwise.

    ------
    Raises
    ------
        TypeError
            If 'tool' is not a string.

        ValueError
            If 'tool' is empty.
    """

    # ####################  MODULES & PACKAGES  ####################
    import matplotlib
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(tool, str):
        raise TypeError(f"Tool name 'tool' must be a string. Got {type(tool).__name__}.")
    elif not tool.strip():
        raise ValueError(f"Tool name 'tool' must be a non-empty string. Got '{tool}'.")


    # ####################  COMPUTATION  ####################
    backend     = matplotlib.get_backend().lower()  # E.g. 'qtagg', 'tkagg', 'agg', 'module://matplotlib_inline.backend_inline'
    interactive = "inline" not in backend and backend not in _NON_INTERACTIVE_BACKENDS  # Inline and file backends never receive events
    if not interactive:
        warnings.warn(
                      f"Matplotlib backend '{backend}' is not interactive: {tool} cannot receive mouse / keyboard events. "
                      "Use e.g. '%matplotlib qt' in Jupyter or the QtAgg / TkAgg backend.",
                      RuntimeWarning,
                      stacklevel=4
                     )


    # ####################  OUTPUT  ####################
    return interactive




# ==============================
def _draw_spectral_lines(
                         ax       : Axes,
                         wave_min : float,
                         wave_max : float
    ) -> list[Artist]:
    """
    Draw the OB-type star diagnostic lines inside a wavelength range, labeled at the bottom of the axes.

    ----------
    Parameters
    ----------
        ax : Axes
            Axes where the lines are drawn.

        wave_min : float
            Minimum wavelength of the data, in Å. Lines below it are not drawn.

        wave_max : float
            Maximum wavelength of the data, in Å. Lines above it are not drawn.

    -------
    Returns
    -------
        line_artists : list[Artist]
            Vertical lines and labels drawn (two artists per line), e.g. to
            toggle their visibility. Empty if no line is inside the range.

    ------
    Raises
    ------
        TypeError
            If 'ax' is not a matplotlib Axes or the wavelengths are not real
            numbers.

        ValueError
            If the wavelengths are not finite or wave_min > wave_max.

    -----
    Notes
    -----
        - The lines and their colors come from '_OB_LINES' and
          '_LINE_COLORS' ('spec_OBlines'): rest-frame air wavelengths in Å.

        - Lines outside [wave_min, wave_max] are skipped, since they would
          widen the wavelength axis.

        - Each label uses a blended transform: its x coordinate is in
          wavelength units and its y coordinate in axes units, so it follows
          its line when zooming / panning in wavelength and always stays at
          the bottom of the axes when zooming in flux. Labels of lines outside
          the visible range are clipped.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from matplotlib.axes import Axes


    # ##### PANORAMA #####
    from .spec_OBlines import _LINE_COLORS, _OB_LINES
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(ax, Axes):
        raise TypeError(f"Axes 'ax' must be a matplotlib Axes. Got {type(ax).__name__}.")
    #
    _check_real(wave_min, "wave_min", "Minimum wavelength")
    #
    _check_real(wave_max, "wave_max", "Maximum wavelength")
    if wave_min > wave_max:
        raise ValueError(f"Wavelength limits 'wave_min' and 'wave_max' must satisfy wave_min <= wave_max. Got {wave_min} and {wave_max}.")


    # ####################  PLOTTING  ####################
    line_artists     = []  # Two artists per line: the vertical line and its label
    bottom_transform = ax.get_xaxis_transform()  # Wavelength units in x, axes units in y (0 = bottom, 1 = top)
    for name, (line_wave, species) in _OB_LINES.items():
        if not wave_min <= line_wave <= wave_max:  # Outside the data: it would widen the wavelength axis
            continue
        color = _LINE_COLORS.get(species, _LINE_COLORS.get(species.split()[0], "0.4"))  # Species first (H, He), then element (metals)
        line_artists.append(ax.axvline(line_wave, color=color, linestyle=":", linewidth=0.8, alpha=0.7, zorder=0))  # Behind the spectra (zorder=0)
        label_artist = ax.text(
                               line_wave,
                               0.02,  # 2% above the bottom of the axes (axes units)
                               name,
                               transform=bottom_transform,
                               rotation=90,
                               ha="center",
                               va="bottom",
                               fontsize=7,
                               color=color,
                               clip_on=True,  # Hidden when the line is outside the visible range
                               bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1.)
                              )
        line_artists.append(label_artist)


    # ####################  OUTPUT  ####################
    return line_artists




# ============================================================================
# ============================== SPECTRUM FILES ==============================
# ============================================================================

# ==============================
def _split_spectrum_name(
                         file : str | Path
    ) -> tuple[str, bool, str]:
    """
    Split a spectrum file name into its basename, gzip flag and extension.

    ----------
    Parameters
    ----------
        file : str | Path
            Path and filename of the spectrum file.

    -------
    Returns
    -------
        basename : str
            File name without its directory (e.g. 'HD886.fits.gz').

        compressed : bool
            True if the file name ends in '.gz' (any case).

        ext : str
            Lowercase extension, without the '.gz' suffix (e.g. '.fits' for
            'HD886.fits.gz'); '' if the name has no extension.

    ------
    Raises
    ------
        TypeError
            If 'file' is not a path.
    """

    # ####################  MODULES & PACKAGES  ####################
    import os
    from pathlib import Path


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(file, (str, os.PathLike)):
        raise TypeError(f"Spectrum file 'file' must be a path (str or pathlib.Path) to the spectrum file. Got {type(file).__name__}.")


    # ####################  COMPUTATION  ####################
    basename   = Path(file).name  # Without the directory
    compressed = basename.lower().endswith(".gz")
    ext        = Path(basename[:-3] if compressed else basename).suffix.lower()  # Extension without '.gz' (e.g. 'HD886.fits.gz' -> '.fits')


    # ####################  OUTPUT  ####################
    return basename, compressed, ext




# ==============================
def _spectrum_file_format(
                          file  : str | Path,
                          name  : str = "output_file",
                          label : str = "Output file"
    ) -> tuple[str, Path]:
    """
    Return the format of a spectrum file ('ascii' or 'fits'), chosen from its extension, and its path without the extension.

    ----------
    Parameters
    ----------
        file : str | Path
            Path and filename of the spectrum file (e.g.
            'results/HD886.fits.gz').

        name : str, optional
            Name of the argument in the error messages (e.g. 'file').
            Default is "output_file".

        label : str, optional
            Role of the argument in the error messages (e.g. 'Spectrum
            file'). Default is "Output file".

    -------
    Returns
    -------
        file_format : str
            'ascii' for the extensions in '_ASCII_EXTENSIONS' (.txt, .dat,
            .asc, .ascii) and 'fits' for those in '_FITS_EXTENSIONS' (.fits,
            .fit, .fts).

        base_path : Path
            'file' without its extension and the '.gz' suffix (e.g.
            'results/HD886'), used to name the files written next to it.

    ------
    Raises
    ------
        TypeError
            If 'file' is not a path, or 'name' or 'label' are not strings.

        ValueError
            If the extension of 'file' (without '.gz') is not supported, or
            'name' or 'label' are empty.

    -----
    Notes
    -----
        - Shared by 'readspec', 'writespec', 'rmvspec' and 'normspec', so that
          the files read and written have the same formats and unsupported
          extensions give the same error message (they are rejected before
          any interactive window opens).
    """

    # ####################  MODULES & PACKAGES  ####################
    from pathlib import Path


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    basename, compressed, ext = _split_spectrum_name(file)  # Also checks that 'file' is a path
    #
    if not isinstance(name, str):
        raise TypeError(f"Error-message text 'name' must be a string. Got {type(name).__name__}.")
    elif not name.strip():
        raise ValueError(f"Error-message text 'name' must be a non-empty string. Got '{name}'.")
    #
    if not isinstance(label, str):
        raise TypeError(f"Error-message text 'label' must be a string. Got {type(label).__name__}.")
    elif not label.strip():
        raise ValueError(f"Error-message text 'label' must be a non-empty string. Got '{label}'.")


    # ####################  COMPUTATION  ####################
    # ##########  FORMAT  ##########
    if ext in _ASCII_EXTENSIONS:
        file_format = "ascii"
    elif ext in _FITS_EXTENSIONS:
        file_format = "fits"
    else:
        raise ValueError(f"{label} '{name}' must have an ASCII {_ASCII_EXTENSIONS} or FITS {_FITS_EXTENSIONS} extension (optionally followed by .gz). Got '{basename}'.")


    # ##########  BASE PATH  ##########
    stem      = basename[:-3] if compressed else basename  # Without '.gz'
    base_path = Path(file).with_name(stem[:-len(ext)])  # Same directory, without the extension


    # ####################  OUTPUT  ####################
    return file_format, base_path




# ==============================
def _read_ascii(
                file : str | Path
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Read the wavelength and flux columns of an ASCII spectrum file (optionally gzip-compressed).

    ----------
    Parameters
    ----------
        file : str | Path
            Path and filename of the ASCII spectrum file (already validated
            by 'readspec').

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelength column (first column), as stored in the file.

        flux : np.ndarray[float]
            Flux column (second column), as stored in the file.

    ------
    Raises
    ------
        TypeError
            If 'file' is not a path.

        ValueError
            If a data line has fewer than two columns, the number of columns
            changes, the wavelength / flux of a data line are not numeric or
            there are no data.

    -----
    Notes
    -----
        - Lines starting with '#' and empty lines are ignored. Commas and
          semicolons are accepted as separators, besides whitespace.

        - Non-numeric lines before the first data line (column names
          without '#', e.g. 'wavelength flux') are skipped with a warning.
          After the first data line, a non-numeric wavelength or flux raises
          an error.

        - Only the first two columns must be numeric: extra columns (e.g.
          quality flags) may contain text and trigger a warning.

        - Every data line must have the same number of columns.
    """

    # ####################  MODULES & PACKAGES  ####################
    import gzip
    import numpy as np
    from pathlib import Path
    import re
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    compressed = _split_spectrum_name(file)[1]  # Also checks that 'file' is a path
    file       = Path(file)


    # ####################  PARAMETERS  ####################
    rows     = []  # (wavelength, flux) of every data line
    ncols    = None  # Set by the first data line
    n_header = 0  # Non-numeric lines skipped before the first data line
    opener   = gzip.open if compressed else open  # Same text interface for plain and gzip-compressed files


    # ####################  COMPUTATION  ####################
    # ##########  READ  ##########
    with opener(file, "rt", encoding="utf-8", errors="replace") as f:  # Undecodable bytes in comments do not prevent reading
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):  # Empty lines and comments
                continue
            parts = re.sub(r"[,;]", " ", line).split()  # Commas and semicolons treated as whitespace
            #
            if ncols is None:
                try:
                    [float(p) for p in parts[:2]]  # Only the wavelength and flux columns must be numeric
                except ValueError:
                    n_header += 1  # Column names before the first data line (e.g. 'wavelength flux')
                    continue
            #
            if len(parts) < 2:
                raise ValueError(f"Spectrum file '{file}' must contain at least 2 columns (wavelength, flux). Problematic line:\n{line}")
            elif ncols is None:
                ncols = len(parts)  # Reference number of columns, fixed by the first data line
            elif len(parts) != ncols:
                raise ValueError(f"Inconsistent number of columns in spectrum file '{file}': expected {ncols}, found {len(parts)}. Problematic line:\n{line}")
            #
            try:
                rows.append((float(parts[0]), float(parts[1])))  # Only columns 1 and 2 are used
            except ValueError as e:
                raise ValueError(f"Non-numeric wavelength or flux in spectrum file '{file}'. Problematic line:\n{line}") from e


    # ##########  CHECKS  ##########
    if not rows:
        raise ValueError(f"Spectrum file '{file}' contains no usable numeric data.")
    if n_header:
        warnings.warn(f"Spectrum file '{file}': {n_header} non-numeric line(s) before the first data line skipped as column names.", stacklevel=3)
    if ncols > 2:
        warnings.warn(f"Spectrum file '{file}': contains {ncols} columns; columns 1 and 2 are used as wavelength and flux.", stacklevel=3)


    # ####################  OUTPUT  ####################
    data = np.array(rows, dtype=float)  # Shape (n_rows, 2)
    return data[:, 0], data[:, 1]




# ==============================
def _read_fits(
               file : str | Path
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Read a FITS spectrum file stored as a binary table or as an image with WCS keywords (optionally gzip-compressed).

    ----------
    Parameters
    ----------
        file : str | Path
            Path and filename of the FITS spectrum file (already validated by
            'readspec').

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelength array, as stored (table) or built from the WCS keywords
            (image).

        flux : np.ndarray[float]
            Flux array.

    ------
    Raises
    ------
        TypeError
            If 'file' is not a path.

        ValueError
            If the file contains neither a table with wavelength / flux
            columns nor an image with data, or the WCS keywords of the image
            are missing.

    -----
    Notes
    -----
        - The binary table is tried first ('_read_fits_table'), then the
          image ('_read_fits_image').

        - Instrument-specific layouts (e.g. a pipeline that stores the
          spectrum in another way) are added here, as another reading
          function tried before the generic ones.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.io import fits
    from pathlib import Path


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _split_spectrum_name(file)  # Checks that 'file' is a path
    file = Path(file)


    # ####################  COMPUTATION  ####################
    with fits.open(file) as hdul:  # Gzip-compressed files are opened transparently by astropy
        spectrum = _read_fits_table(hdul, file)  # Binary table first
        if spectrum is None:
            spectrum = _read_fits_image(hdul, file)  # Then image with WCS keywords
    if spectrum is None:
        raise ValueError(f"FITS file '{file}' contains no valid spectrum (table with wavelength / flux columns or image with WCS keywords).")


    # ####################  OUTPUT  ####################
    return spectrum




# ==============================
def _read_fits_table(
                     hdul : HDUList,
                     file : Path
    ) -> tuple[np.ndarray[float], np.ndarray[float]] | None:
    """
    Read the spectrum from the first binary table with wavelength and flux columns.

    ----------
    Parameters
    ----------
        hdul : HDUList
            Open FITS file.

        file : Path
            Path and filename of the FITS file (used in the warnings).

    -------
    Returns
    -------
        spectrum : tuple[np.ndarray[float], np.ndarray[float]] | None
            Wavelength and flux arrays (wave, flux), or None if no table has
            recognizable columns.

    ------
    Raises
    ------
        TypeError
            If 'hdul' is not an astropy HDUList or 'file' is not a path.

    -----
    Notes
    -----
        - Columns are searched by name ('_FITS_WAVE_COLS' / '_FITS_FLUX_COLS',
          case-insensitive) across every extension (e.g. HERMES / IACOB
          spectra, files written by 'writespec').

        - Tables without recognizable columns (observing logs, order
          lists...) are skipped: they may be auxiliary tables next to an
          image.

        - If the columns store several rows, only the first one is used
          (with a warning).
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.io import fits
    import numpy as np
    import os
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(hdul, fits.HDUList):
        raise TypeError(f"Open FITS file 'hdul' must be an astropy HDUList. Got {type(hdul).__name__}.")
    #
    if not isinstance(file, (str, os.PathLike)):
        raise TypeError(f"Spectrum file 'file' must be a path (str or pathlib.Path) to the spectrum file. Got {type(file).__name__}.")


    # ####################  COMPUTATION  ####################
    for hdu in hdul[1:]:  # Tables are always extensions (never the primary HDU)
        # The FITS_rec check also excludes compressed images, stored internally as tables
        if not isinstance(hdu, fits.BinTableHDU) or not isinstance(hdu.data, fits.FITS_rec):
            continue
        #
        colmap   = {col.lower(): col for col in hdu.data.columns.names}  # Lowercase name -> actual column name (case-insensitive search)
        wave_col = next((colmap[c] for c in _FITS_WAVE_COLS if c in colmap), None)  # First recognized wavelength column
        flux_col = next((colmap[c] for c in _FITS_FLUX_COLS if c in colmap), None)  # First recognized flux column
        if wave_col is None or flux_col is None:  # Auxiliary table (observing log, order list...): skipped
            continue
        #
        wave_data = np.asarray(hdu.data[wave_col])  # 1D (one value per row) or 2D (one array per row)
        flux_data = np.asarray(hdu.data[flux_col])
        if wave_data.ndim > 1 and wave_data.shape[0] > 1:  # Several spectra stored in the table
            warnings.warn(f"Spectrum file '{file}': wavelength / flux columns store {wave_data.shape[0]} rows; only the first row is used.", stacklevel=4)
        #
        wave = np.asarray(wave_data[0] if wave_data.ndim > 1 else wave_data, dtype=float)  # First row of array columns
        flux = np.asarray(flux_data[0] if flux_data.ndim > 1 else flux_data, dtype=float)
        return wave, flux  # First table with both columns


    # ####################  OUTPUT  ####################
    return None  # No table with wavelength and flux columns




# ==============================
def _read_fits_image(
                     hdul : HDUList,
                     file : Path
    ) -> tuple[np.ndarray[float], np.ndarray[float]] | None:
    """
    Read the spectrum from the first image with data, building the wavelengths from its WCS keywords.

    ----------
    Parameters
    ----------
        hdul : HDUList
            Open FITS file.

        file : Path
            Path and filename of the FITS file (used in the messages and to
            detect ln(lambda) scales from 'log' in its name).

    -------
    Returns
    -------
        spectrum : tuple[np.ndarray[float], np.ndarray[float]] | None
            Wavelength and flux arrays (wave, flux), or None if the file has
            no image with data.

    ------
    Raises
    ------
        TypeError
            If 'hdul' is not an astropy HDUList or 'file' is not a path.

        ValueError
            If CRVAL1 or CDELT1 / CD1_1 are missing.

    -----
    Notes
    -----
        1. IMAGE
            - Normally the PrimaryHDU, but some pipelines leave it empty and
              put the spectrum (with its own WCS header) in an extension.
              Keywords are searched in the image header first, then in the
              primary header.

            - Multi-dimensional images (e.g. flux + error, IRAF multispec):
              the first 1D spectrum (first row) is used, with a warning.

            - Missing CRPIX1: pixel 1 is assumed (FEROS / IRAF convention).

        2. WAVELENGTH SCALE
            - DC-FLAG = 1: log10(lambda) (IRAF).

            - CTYPE1 ending in '-LOG': logarithmic scale of the FITS standard
              (Greisen et al. 2006).

            - Values below '_LN_WAVE_MAX' with '-LOG' in CTYPE1 or 'log' in
              the file name: ln(lambda) (e.g. HERMES '*_log_merged_cf.fits'
              spectra, whose headers do not flag the scale).

    ----------
    References
    ----------
        - Greisen, E. W., Calabretta, M. R., Valdes, F. G. et al. 2006,
          A&A, 446, 747.
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.io import fits
    import numpy as np
    import os
    from pathlib import Path
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(hdul, fits.HDUList):
        raise TypeError(f"Open FITS file 'hdul' must be an astropy HDUList. Got {type(hdul).__name__}.")
    #
    if not isinstance(file, (str, os.PathLike)):
        raise TypeError(f"Spectrum file 'file' must be a path (str or pathlib.Path) to the spectrum file. Got {type(file).__name__}.")


    # ##########  FUNCTIONS  ##########
    def _header_value(
                      key     : str,
                      default : object = None
        ) -> object:
        """
        Return the value of a keyword from the image header or, if missing, from the primary header.

        ----------
        Parameters
        ----------
            key : str
                Keyword name (e.g. 'CRVAL1').

            default : object, optional
                Value returned if no header contains the keyword. Default is
                None.

        -------
        Returns
        -------
            value : object
                Keyword value, or 'default'.
        """
        return next((h[key] for h in headers if key in h), default)


    # ####################  COMPUTATION  ####################
    # ##########  IMAGE  ##########
    images = [h for h in hdul if isinstance(h, (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)) and h.data is not None and h.data.size > 0]  # Image HDUs with data
    if not images:
        return None  # No image: not a valid spectrum (handled by '_read_fits')
    image_hdu = images[0]  # Usually the primary HDU; an extension if the primary one is empty
    headers   = (image_hdu.header, hdul[0].header)  # Searched in this order by '_header_value'


    # ##########  FLUX  ##########
    data = np.squeeze(np.asarray(image_hdu.data, dtype=float))  # Removes axes of length 1 (e.g. shape (1, 1, N) -> (N,))
    if data.ndim > 1:
        warnings.warn(f"Spectrum file '{file}': image data has shape {data.shape}; only the first 1D spectrum (first row) is used.", stacklevel=4)
        data = data.reshape(-1, data.shape[-1])[0]  # NAXIS1 is the last NumPy axis
    flux = np.atleast_1d(data)  # A single-pixel image is still a 1D array


    # ##########  WCS KEYWORDS  ##########
    crval = _header_value("CRVAL1")
    cdelt = _header_value("CDELT1", _header_value("CD1_1"))
    crpix = _header_value("CRPIX1", 1.)  # FEROS / IRAF convention when missing
    if crval is None:
        raise ValueError(f"Missing WCS keyword 'CRVAL1' in FITS file '{file}'.")
    elif cdelt is None:
        raise ValueError(f"Missing WCS keyword 'CDELT1' (or 'CD1_1') in FITS file '{file}'.")
    crval, cdelt, crpix = float(crval), float(cdelt), float(crpix)  # Header values may be stored as strings or integers


    # ##########  WAVELENGTH SCALE  ##########
    x             = crval + (np.arange(flux.size) - (crpix - 1.)) * cdelt  # FITS pixels are 1-based
    ctype         = str(_header_value("CTYPE1", '')).strip().upper()
    is_log_values = np.nanmax(x) < _LN_WAVE_MAX  # 'x' holds log(lambda), not lambda
    #
    if _header_value("DC-FLAG", 0) == 1:
        wave = 10. ** x  # IRAF: log10(lambda)
    elif ctype.endswith("-LOG") and not is_log_values:
        wave = crval * np.exp((x - crval) / crval)  # FITS standard (Greisen et al. 2006)
    elif is_log_values and (ctype.endswith("-LOG") or "log" in Path(file).name.lower()):
        wave = np.exp(x)  # Natural log ln(lambda): e.g. HERMES '*_log_merged_cf' files, non-standard -LOG headers
    else:
        wave = x  # Linear scale


    # ####################  OUTPUT  ####################
    return wave, flux




# ======================================================================
# ============================== PLOTSPEC ==============================
# ======================================================================

# ==============================
class _SpectrumToggle:
    """
    Callback of the check buttons of 'plotspec' that shows / hides each spectrum and the OB-type star lines.

    ----------
    Parameters
    ----------
        fig : Figure
            Figure of 'plotspec'.

        spectra_plots : list[Line2D]
            Curves of the spectra, labeled with the text of their check
            button.

        line_artists : list[Artist]
            Vertical lines and labels of the OB-type star diagnostic lines
            (see '_draw_spectral_lines'); empty if they are not drawn.

        lines_label : str
            Text of the check button of the spectral lines.

    -----
    Notes
    -----
        - Matplotlib may keep only weak references to bound-method callbacks,
          so 'plotspec' stores this object in the figure to keep the check
          buttons working while the figure is open.
    """

    # ####################  INIT  ####################
    def __init__(
                 self,
                 fig           : Figure,
                 spectra_plots : list[Line2D],
                 line_artists  : list[Artist],
                 lines_label   : str
        ) -> None:
        """
        Store the figure and the artists whose visibility is toggled.

        ----------
        Parameters
        ----------
            fig : Figure
                Figure of 'plotspec'.

            spectra_plots : list[Line2D]
                Curves of the spectra, labeled with the text of their check
                button.

            line_artists : list[Artist]
                Vertical lines and labels of the OB-type star diagnostic
                lines; empty if they are not drawn.

            lines_label : str
                Text of the check button of the spectral lines.

        -------
        Returns
        -------
            None
                The figure and the artists are stored as attributes.

        ------
        Raises
        ------
            TypeError
                If any input has the wrong type.

            ValueError
                If 'lines_label' is empty.
        """

        # ##########  MODULES & PACKAGES  ##########
        from matplotlib.artist import Artist
        from matplotlib.figure import Figure
        from matplotlib.lines import Line2D


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        if not isinstance(fig, Figure):
            raise TypeError(f"Figure 'fig' must be a matplotlib Figure. Got {type(fig).__name__}.")
        #
        if not isinstance(spectra_plots, list) or not all(isinstance(p, Line2D) for p in spectra_plots):
            raise TypeError(f"Spectrum curves 'spectra_plots' must be a list of matplotlib Line2D. Got {type(spectra_plots).__name__}.")
        #
        if not isinstance(line_artists, list) or not all(isinstance(a, Artist) for a in line_artists):
            raise TypeError(f"Spectral line artists 'line_artists' must be a list of matplotlib artists. Got {type(line_artists).__name__}.")
        #
        if not isinstance(lines_label, str):
            raise TypeError(f"Spectral lines button label 'lines_label' must be a string. Got {type(lines_label).__name__}.")
        elif not lines_label.strip():
            raise ValueError(f"Spectral lines button label 'lines_label' must be a non-empty string. Got '{lines_label}'.")


        # ##########  PARAMETERS  ##########
        self.fig           = fig
        self.spectra_plots = spectra_plots
        self.line_artists  = line_artists
        self.lines_label   = lines_label


    # ####################  EVENTS  ####################
    def _on_clicked(
                    self,
                    label : str
        ) -> None:
        """
        Toggle the visibility of the spectrum with the given label, or of all the spectral lines; callback of the check buttons.

        ----------
        Parameters
        ----------
            label : str
                Text of the check button clicked.

        -------
        Returns
        -------
            None
                The visibility of the spectrum (or of the spectral lines and
                their labels) is toggled and the figure redrawn.
        """

        # ##########  COMPUTATION  ##########
        # ##### OPTION 1: SPECTRAL LINES #####
        if label == self.lines_label:
            for artist in self.line_artists:  # Vertical lines and their labels
                artist.set_visible(not artist.get_visible())


        # ##### OPTION 2: SPECTRUM #####
        else:
            for spectrum_plot in self.spectra_plots:
                if spectrum_plot.get_label() == label:  # Labels are unique (see 'plotspec')
                    spectrum_plot.set_visible(not spectrum_plot.get_visible())
                    break


        # ##########  PLOTTING  ##########
        self.fig.canvas.draw_idle()  # Redraw when the event loop is idle




# =====================================================================
# ============================== RMVSPEC ==============================
# =====================================================================

# ==============================
class _InteractiveSpectrumCleaner:
    """
    Interactive window of 'rmvspec' to clean a spectrum by removing selected wavelength spans.

    ----------
    Parameters
    ----------
        wave : np.ndarray[float]
            Wavelength array (finite, sorted, at least two points).

        flux : np.ndarray[float]
            Flux array (finite, same length as 'wave').

        output_file : Path | None
            File where the cleaned spectrum is saved with 'writespec', or None
            (nothing is saved).

        show_lines : bool, optional
            Draw the OB-type star diagnostic lines. Default is True.

        autosave : bool, optional
            Save the cleaned spectrum when the window is closed. Default is
            True.

    -----
    Notes
    -----
        1. USE
            - Spans are selected by dragging the left mouse button (red) and
              removed with 'r'. Removed points stay visible in light gray;
              every removal can be undone / redone.

            - The cleaned spectrum is returned by '_run' when the window is
              closed.

        2. KEYBOARD
            - Matplotlib's default key handler is disconnected for this
              figure only (no global rcParams are modified). Every key goes to
              '_on_key': the keys of the tool (see 'KEYS') run the tool
              action, and any other key is forwarded to Matplotlib's default
              handler, so zoom ('o'), pan ('p'), close ('q'), arrows... keep
              working.

            - Backends with 'toolmanager' register their shortcuts per tool:
              the keys of the tool are removed from their keymaps while the
              window is open and restored when it is closed (see
              '_take_keyboard').

            - 'shift+<letter>' (sent by some backends) is normalized to the
              uppercase letter, so 'U', 'S' and 'D' work everywhere.

            - 'l' shows / hides the OB-type star lines, replacing Matplotlib's
              log-scale toggle of the flux axis in this window.

        3. OUTPUT FILES
            - The cleaned spectrum is written with 'writespec' (format chosen
              from the extension of 'output_file'). Figures and exported spans
              are named after 'base_path': 'output_file' without its
              extension (and '.gz').
    """

    # ####################  CONFIGURATION  ####################
    KEYS = ("r", "x", "c", "a", "u", "U", "g", "l", "h", "s", "S", "d", "D")  # They always take priority over Matplotlib's shortcuts


    # ####################  INIT  ####################
    def __init__(
                 self,
                 wave        : np.ndarray[float],
                 flux        : np.ndarray[float],
                 output_file : Path | None,
                 show_lines  : bool = True,
                 autosave    : bool = True
        ) -> None:
        """
        Initialize the state of the tool and build the figure.

        ----------
        Parameters
        ----------
            wave : np.ndarray[float]
                Wavelength array (finite, sorted, at least two points).

            flux : np.ndarray[float]
                Flux array (finite, same length as 'wave').

            output_file : Path | None
                File where the cleaned spectrum is saved with 'writespec' ('s',
                'S' and, if 'autosave', on closing), in an existing directory
                and with an extension supported by 'writespec'. None disables
                saving and exporting.

            show_lines : bool, optional
                Draw the OB-type star diagnostic lines inside the spectrum
                ('_OB_LINES'), labeled at the bottom of the axes. They can be
                hidden / shown with 'l'. Default is True.

            autosave : bool, optional
                Save the cleaned spectrum to 'output_file' when the window is
                closed, unless it has not changed since the last save. Default
                is True.

        -------
        Returns
        -------
            None
                The state of the tool is stored and its figure is created.

        ------
        Raises
        ------
            TypeError
                If any input has the wrong type.

            ValueError
                If the arrays have different lengths, fewer than two points or
                non-finite / unsorted values, or the extension of
                'output_file' is not supported by 'writespec'.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np
        from pathlib import Path


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        _check_array(wave, "wave", "Wavelength array")
        if wave.size < 2:
            raise ValueError(f"Wavelength array 'wave' must contain at least 2 points. Got {wave.size}.")
        elif not np.all(np.isfinite(wave)):
            raise ValueError("Wavelength array 'wave' must contain only finite values.")
        elif np.any(np.diff(wave) < 0):
            raise ValueError("Wavelength array 'wave' must be sorted in increasing order.")
        #
        _check_array(flux, "flux", "Flux array", size=wave.size)
        if not np.all(np.isfinite(flux)):
            raise ValueError("Flux array 'flux' must contain only finite values.")
        #
        if output_file is not None and not isinstance(output_file, Path):
            raise TypeError(f"Output file 'output_file' must be a pathlib.Path (already validated by 'rmvspec') or None. Got {type(output_file).__name__}.")
        base_path = _spectrum_file_format(output_file)[1] if output_file is not None else None  # Also rejects unsupported extensions
        #
        if not isinstance(show_lines, bool):
            raise TypeError(f"Spectral lines flag 'show_lines' must be a boolean. Got {type(show_lines).__name__}.")
        #
        if not isinstance(autosave, bool):
            raise TypeError(f"Autosave flag 'autosave' must be a boolean. Got {type(autosave).__name__}.")


        # ##########  PARAMETERS  ##########
        # ##### DATA #####
        self.x           = wave
        self.y           = flux
        self.output_file = output_file
        self.base_path   = base_path  # 'output_file' without extension: base name of the figures and exported spans
        self.show_lines  = show_lines
        self.autosave    = autosave


        # ##### SELECTIONS #####
        self.mask            = np.ones(self.x.size, dtype=bool)  # True: point kept
        self.pending         = []  # Selected spans (xmin, xmax), not removed yet
        self.pending_patches = []  # Red patches of the pending spans
        self.removed         = []  # Spans effectively removed


        # ##### UNDO / REDO #####
        self.undo_stack = []  # States (mask, removed)
        self.redo_stack = []


        # ##### MISCELLANEOUS #####
        self.show_grid     = False
        self._saved_mask   = None  # Mask at the last save (autosave only if it changed)
        self._forward_keys = False  # Forward non-tool keys to Matplotlib's default handler
        self._keymaps      = {}  # Original keymaps of the toolmanager tools (restored on close)


        # ##########  COMPUTATION  ##########
        self._build_figure()
        self._forward_keys, self._keymaps = _take_keyboard(self.fig, self.KEYS)  # The tool keys take priority over Matplotlib's shortcuts


    # ####################  FIGURE  ####################
    def _build_figure(
                      self
        ) -> None:
        """
        Create the figure, the curves, the reference lines and the interactions.

        -------
        Returns
        -------
            None
                The figure and its artists are stored as attributes.
        """

        # ##########  MODULES & PACKAGES  ##########
        import matplotlib.pyplot as plt
        from matplotlib.widgets import SpanSelector
        import numpy as np


        # ##########  PLOTTING  ##########
        # ##### FIGURE #####
        self.fig, self.ax = plt.subplots(figsize=(9, 5))
        self.fig.subplots_adjust(bottom=0.18)  # Room for the status line


        # ##### CURVES #####
        # Removed points in light gray (visual reference); gaps are drawn as gaps
        self.ghost_line, = self.ax.plot(self.x, np.full_like(self.y, np.nan), "-", lw=0.8, color="0.75", zorder=1)
        self.line,       = self.ax.plot(self.x, self.y, "-", lw=0.8, color="black", zorder=2)


        # ##### SPECTRAL LINES #####
        self.xlim_orig    = (self.x[0], self.x[-1])  # Whole spectrum (wavelengths already sorted)
        self.line_artists = _draw_spectral_lines(self.ax, self.x[0], self.x[-1]) if self.show_lines else []


        # ##### AXES #####
        self.ax.set_xlabel("Wavelength [Å]")
        self.ax.set_ylabel("Flux")
        self.ax.set_xlim(self.xlim_orig)
        self._autoscale_y()


        # ##### STATUS #####
        self.status_text = self.fig.text(0.01, 0.015, '', fontsize=7.5, family="monospace", color="0.35")  # Bottom-left corner, figure coordinates
        self._update_status()


        # ##### SPAN SELECTOR #####
        step = np.median(np.diff(self.x))  # Typical pixel size, in wavelength units
        self.span_selector = SpanSelector(
                                          self.ax,
                                          self._on_select,
                                          "horizontal",  # Spans in wavelength
                                          useblit=True,  # Fast redraw while dragging
                                          button=1,  # Left mouse button
                                          minspan=0.5 * step,  # A single click does not create an empty span
                                          props=dict(facecolor="red", alpha=0.25)
                                         )


        # ##### EVENTS #####
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)  # Tool keys (and forwarding of the others)
        self.fig.canvas.mpl_connect("close_event", self._on_close)  # Autosave and keymap restoration
    #
    def _autoscale_y(
                     self
        ) -> None:
        """
        Fit the y axis to the kept points inside the visible wavelength range.

        -------
        Returns
        -------
            None
                The y limits of the axes are updated (unchanged if no kept
                point is visible).
        """

        # ##########  PLOTTING  ##########
        x0, x1 = sorted(self.ax.get_xlim())  # Visible wavelength range (sorted in case the axis is inverted)
        sel    = self.mask & (self.x >= x0) & (self.x <= x1)  # Kept points inside it
        if not sel.any():
            return
        ymin, ymax = self.y[sel].min(), self.y[sel].max()
        pad        = 0.05 * (ymax - ymin) if ymax > ymin else (0.05 * abs(ymax) or 1.)  # 5% margin; constant or zero flux also handled
        self.ax.set_ylim(ymin - pad, ymax + pad)
    #
    def _update_status(
                       self
        ) -> None:
        """
        Update the status line (pending spans, removed points, undo / redo).

        -------
        Returns
        -------
            None
                The status text of the figure is updated.
        """

        # ##########  PLOTTING  ##########
        n_removed = int((~self.mask).sum())  # Points removed so far
        percent   = 100. * n_removed / self.mask.size
        self.status_text.set_text(
                                  f"pending: {len(self.pending)} | removed: {n_removed} pts ({percent:.1f}%) | "
                                  f"undo: {len(self.undo_stack)}  redo: {len(self.redo_stack)}"
                                 )
    #
    def _redraw(
                self
        ) -> None:
        """
        Update the curves keeping the current zoom of the user.

        -------
        Returns
        -------
            None
                The figure is redrawn.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  PLOTTING  ##########
        self.line.set_ydata(np.where(self.mask, self.y, np.nan))  # Kept points (NaN breaks the line at removed points)
        self.ghost_line.set_ydata(np.where(self.mask, np.nan, self.y))  # Removed points, in light gray
        self._autoscale_y()
        self._update_status()
        self.fig.canvas.draw_idle()  # Redraw when the event loop is idle
    #
    @staticmethod
    def _msg(
             text : str
        ) -> None:
        """
        Print a message of the tool.

        ----------
        Parameters
        ----------
            text : str
                Message to print, after the '[rmvspec]' prefix.

        -------
        Returns
        -------
            None
                The message is printed to the terminal.
        """

        # ##########  OUTPUT  ##########
        print(f"[rmvspec] {text}")


    # ####################  EDITING  ####################
    def _on_select(
                   self,
                   xmin : float,
                   xmax : float
        ) -> None:
        """
        Store a span selected with the mouse as pending (red); callback of the SpanSelector.

        ----------
        Parameters
        ----------
            xmin : float
                One end of the span, in wavelength units.

            xmax : float
                Other end of the span, in wavelength units.

        -------
        Returns
        -------
            None
                The span is added to 'pending' and drawn.
        """

        # ##########  COMPUTATION  ##########
        xmin, xmax = sorted((float(xmin), float(xmax)))  # Dragging to the left gives xmin > xmax
        self.pending.append((xmin, xmax))
        self.pending_patches.append(self.ax.axvspan(xmin, xmax, color="red", alpha=0.25))  # Red patch until removed or discarded
        self._update_status()
        self.fig.canvas.draw_idle()
    #
    def _discard_last(
                      self
        ) -> None:
        """
        Discard the last pending span ('x').

        -------
        Returns
        -------
            None
                The last pending span and its patch are removed.
        """

        # ##########  COMPUTATION  ##########
        if not self.pending:
            self._msg("No pending selections.")
            return
        self.pending.pop()
        self.pending_patches.pop().remove()  # Also removes the patch from the axes
        self._update_status()
        self.fig.canvas.draw_idle()
    #
    def _clear_pending(
                       self
        ) -> None:
        """
        Discard every pending span ('c').

        -------
        Returns
        -------
            None
                All pending spans and their patches are removed.
        """

        # ##########  COMPUTATION  ##########
        for patch in self.pending_patches:
            patch.remove()  # Removed from the axes
        self.pending_patches.clear()
        self.pending.clear()
        self._update_status()
        self.fig.canvas.draw_idle()
    #
    def _snapshot(
                  self
        ) -> tuple[np.ndarray[bool], list[tuple[float, float]]]:
        """
        Return the current editing state.

        -------
        Returns
        -------
            state : tuple[np.ndarray[bool], list[tuple[float, float]]]
                Copies of the mask of kept points and of the removed spans,
                in this order.
        """

        # ##########  OUTPUT  ##########
        return self.mask.copy(), list(self.removed)  # Copies, so that later edits do not modify the stored state
    #
    def _push_undo(
                   self
        ) -> None:
        """
        Store the current state before an edit; a new edit clears the redo stack.

        -------
        Returns
        -------
            None
                The undo and redo stacks are updated.
        """

        # ##########  COMPUTATION  ##########
        self.undo_stack.append(self._snapshot())
        self.redo_stack.clear()  # A new edit invalidates the undone states
    #
    def _remove_pending(
                        self
        ) -> None:
        """
        Remove the points inside every pending span ('r').

        -------
        Returns
        -------
            None
                The mask, the removed spans and the figure are updated.
        """

        # ##########  COMPUTATION  ##########
        # ##### PENDING SPANS #####
        if not self.pending:
            self._msg("No spans selected.")
            return


        # ##### NEW MASK #####
        new_mask  = self.mask.copy()
        effective = []  # Spans that still contain kept points
        for xmin, xmax in self.pending:
            span = (self.x >= xmin) & (self.x <= xmax)  # Points inside the span
            if (span & new_mask).any():  # Ignore spans where every point was already removed
                effective.append((xmin, xmax))
                new_mask &= ~span
        #
        if not effective:
            self._msg("Selected spans contain no remaining data points.")
            self._clear_pending()
            return


        # ##### UPDATE #####
        self._push_undo()  # State before the removal
        n_before  = int(self.mask.sum())
        self.mask = new_mask
        self.removed.extend(effective)
        self._clear_pending()  # Pending spans become removed spans
        self._redraw()
        self._msg(f"Removed {n_before - int(self.mask.sum())} points in {len(effective)} span(s).")
    #
    def _restore_all(
                     self
        ) -> None:
        """
        Restore every removed point ('a'); the action can be undone.

        -------
        Returns
        -------
            None
                The mask and the removed spans are reset.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  COMPUTATION  ##########
        if self.mask.all():
            self._msg("No removed points to restore.")
            return
        self._push_undo()  # Restoring can be undone
        self.mask    = np.ones_like(self.mask)  # Every point kept again
        self.removed = []
        self._redraw()
        self._msg("All points restored.")
    #
    def _undo(
              self
        ) -> None:
        """
        Undo the last edit ('u').

        -------
        Returns
        -------
            None
                The previous state is restored.
        """

        # ##########  COMPUTATION  ##########
        if not self.undo_stack:
            self._msg("Nothing to undo.")
            return
        self.redo_stack.append(self._snapshot())  # Current state, recoverable with 'U'
        self.mask, self.removed = self.undo_stack.pop()  # Previous state
        self._redraw()
    #
    def _redo(
              self
        ) -> None:
        """
        Redo the last undone edit ('U').

        -------
        Returns
        -------
            None
                The undone state is restored.
        """

        # ##########  COMPUTATION  ##########
        if not self.redo_stack:
            self._msg("Nothing to redo.")
            return
        self.undo_stack.append(self._snapshot())  # Current state, recoverable with 'u'
        self.mask, self.removed = self.redo_stack.pop()  # Last undone state
        self._redraw()
    #
    def _toggle_grid(
                     self
        ) -> None:
        """
        Show / hide the grid ('g').

        -------
        Returns
        -------
            None
                The grid of the axes is toggled.
        """

        # ##########  PLOTTING  ##########
        self.show_grid = not self.show_grid
        self.ax.grid(self.show_grid)
        self.fig.canvas.draw_idle()
    #
    def _toggle_lines(
                      self
        ) -> None:
        """
        Show / hide the OB-type star diagnostic lines and their labels ('l').

        -------
        Returns
        -------
            None
                The visibility of the lines is toggled (a message is printed
                if no line was drawn).
        """

        # ##########  PLOTTING  ##########
        if not self.line_artists:
            self._msg("No OB-type lines drawn (outside the spectrum, or show_lines=False).")
            return
        visible = not self.line_artists[0].get_visible()  # Every line shares the same visibility
        for artist in self.line_artists:
            artist.set_visible(visible)
        self.fig.canvas.draw_idle()


    # ####################  SAVING  ####################
    def _has_output(
                    self
        ) -> bool:
        """
        Check that an output file was given, printing a message otherwise.

        -------
        Returns
        -------
            has_output : bool
                True if 'output_file' is not None.
        """

        # ##########  OUTPUT  ##########
        if self.output_file is None:
            self._msg("No output file given. Call rmvspec(..., output_file='file.txt') to save.")
            return False
        return True
    #
    def _save_spectrum(
                       self,
                       png : bool = False
        ) -> None:
        """
        Save the cleaned spectrum ('s') and, optionally, a PNG of the figure and its pickle copy ('S').

        ----------
        Parameters
        ----------
            png : bool, optional
                If True, also save the figure as '<stem>.png' and '<stem>.pkl'
                next to 'output_file' ('<stem>' = 'base_path'), with
                'savefig_pkl'. Default is False.

        -------
        Returns
        -------
            None
                The cleaned spectrum is written to 'output_file' with
                'writespec' (ASCII or FITS, from its extension), replacing the
                previous file.

        -----
        Notes
        -----
            - Wavelengths and fluxes are written with the precision of
              'writespec' ('%.12g' and '%.10g' in ASCII files, float64 in FITS
              files), so absolute fluxes (~1e-15) keep their precision.

            - Pending (red) selections and the status line are hidden in the
              saved figure, so the tight crop of 'savefig_pkl' fits the axes
              without blank margins. They are shown again afterwards and the
              window stays open.

            - Errors while writing (e.g. permissions, or repeated wavelengths
              rejected by 'writespec') are printed, not raised, so the window
              stays usable.
        """

        # ##########  MODULES & PACKAGES  ##########
        # ##### PANORAMA #####
        from ...spec import writespec
        from ...tex import savefig_pkl


        # ##########  OUTPUT  ##########
        # ##### SPECTRUM #####
        if not self._has_output():
            return
        try:
            writespec((self.x[self.mask], self.y[self.mask]), self.output_file, overwrite=True)  # Replaced at every save
        except Exception as e:  # Never raise inside a Matplotlib callback
            self._msg(f"Could not save the cleaned spectrum: {e}")
            return
        self._saved_mask = self.mask.copy()  # Autosave on closing only if the mask changes after this save
        self._msg(f"Cleaned spectrum saved as '{self.output_file}'.")


        # ##### FIGURE #####
        if png:
            # The status line sits at the bottom-left corner of the figure: with the tight crop of 'savefig_pkl' it would add blank margins on the left and at the bottom
            hidden = [artist for artist in self.pending_patches + [self.status_text] if artist.get_visible()]
            for artist in hidden:
                artist.set_visible(False)  # Not shown in the saved figure
            try:
                written_files = savefig_pkl(self.fig, self.base_path, formats=("png",), close=False)  # The window stays open
                self._msg(f"Figure saved as {', '.join(repr(str(f)) for f in written_files)}.")
            except Exception as e:  # Never raise inside a Matplotlib callback
                self._msg(f"Could not save the figure: {e}")
            finally:
                for artist in hidden:
                    artist.set_visible(True)  # Shown again in the window
                self.fig.canvas.draw_idle()
    #
    def _export_spans(
                      self,
                      spans  : list[tuple[float, float]],
                      suffix : str
        ) -> None:
        """
        Export spans, merged, to '<stem>_<suffix>.txt' ('d': removed, 'D': pending).

        ----------
        Parameters
        ----------
            spans : list[tuple[float, float]]
                Wavelength spans (min, max) to export.

            suffix : str
                Suffix of the exported file name (e.g. 'rmvspans').

        -------
        Returns
        -------
            None
                The spans are written as two tab-separated columns (wave_min,
                wave_max) with a '#' header, next to 'output_file'
                ('<stem>' = 'base_path').
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  OUTPUT  ##########
        if not spans:
            self._msg(f"No {suffix} spans available to export.")
            return
        if not self._has_output():
            return
        export_path = self.base_path.with_name(f"{self.base_path.name}_{suffix}.txt")  # E.g. 'HD886_rmvspans.txt'
        try:
            np.savetxt(export_path, np.asarray(_merge_ranges(spans)), fmt="%.6f", delimiter="\t", header="wave_min\twave_max", encoding="utf-8")  # Overlapping spans merged
        except OSError as e:
            self._msg(f"Could not export the spans: {e}")
            return
        self._msg(f"Spans exported to '{export_path}'.")


    # ####################  EVENTS  ####################
    @staticmethod
    def _normalize_key(
                       key : str | None
        ) -> str | None:
        """
        Normalize 'shift+<letter>' (sent by some backends) to the uppercase letter.

        ----------
        Parameters
        ----------
            key : str | None
                Key of the event (e.g. 'shift+u').

        -------
        Returns
        -------
            key : str | None
                Normalized key (e.g. 'U'); other keys are returned unchanged.
        """

        # ##########  OUTPUT  ##########
        if key and key.lower().startswith("shift+") and len(key) == 7:  # 'shift+' followed by a single character
            return key[-1].upper()
        return key
    #
    def _on_key(
                self,
                event : KeyEvent
        ) -> None:
        """
        Run the action of a tool key, or forward any other key to Matplotlib's default handler.

        ----------
        Parameters
        ----------
            event : KeyEvent
                Key press event of the figure.

        -------
        Returns
        -------
            None
                The action of the key is executed.
        """

        # ##########  MODULES & PACKAGES  ##########
        from matplotlib.backend_bases import key_press_handler


        # ##########  PARAMETERS  ##########
        actions = {
                   'r' : self._remove_pending,
                   'x' : self._discard_last,
                   'c' : self._clear_pending,
                   'a' : self._restore_all,
                   'u' : self._undo,
                   'U' : self._redo,
                   'g' : self._toggle_grid,
                   'l' : self._toggle_lines,
                   'h' : lambda: print(_RMVSPEC_HELP),
                   's' : lambda: self._save_spectrum(png=False),
                   'S' : lambda: self._save_spectrum(png=True),
                   'd' : lambda: self._export_spans(self.removed, "rmvspans"),
                   'D' : lambda: self._export_spans(self.pending, "selectedspans"),
                  }


        # ##########  COMPUTATION  ##########
        action = actions.get(self._normalize_key(event.key))  # None for keys not used by the tool
        if action is not None:
            action()  # Tool key: always the tool action
        elif self._forward_keys:
            key_press_handler(event, self.fig.canvas, self.fig.canvas.toolbar)  # Other keys: Matplotlib's behavior
    #
    def _on_close(
                  self,
                  event : CloseEvent
        ) -> None:
        """
        Autosave the cleaned spectrum (if enabled and changed) and restore the keymaps when the window closes.

        ----------
        Parameters
        ----------
            event : CloseEvent
                Close event of the figure (unused).

        -------
        Returns
        -------
            None
                The spectrum may be written to 'output_file'.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  OUTPUT  ##########
        changed = self._saved_mask is None or not np.array_equal(self._saved_mask, self.mask)  # Never saved, or edited since the last save
        if self.autosave and self.output_file is not None and changed:
            self._save_spectrum(png=False)
        _release_keyboard(self.fig, self._keymaps)
        self._keymaps = {}  # Restored only once


    # ####################  RUN  ####################
    def _run(
             self
        ) -> tuple[np.ndarray[float], np.ndarray[float]]:
        """
        Show the window and return the cleaned spectrum when it is closed.

        -------
        Returns
        -------
            wave_clean : np.ndarray[float]
                Wavelengths of the kept points.

            flux_clean : np.ndarray[float]
                Fluxes of the kept points.

        -----
        Notes
        -----
            - With a non-interactive backend (e.g. Agg, inline) a
              RuntimeWarning is issued, the figure is closed without
              autosaving and the spectrum is returned unmodified.
        """

        # ##########  MODULES & PACKAGES  ##########
        import matplotlib.pyplot as plt


        # ##########  COMPUTATION  ##########
        # ##### OPTION 1: NON-INTERACTIVE BACKEND #####
        if not _check_interactive_backend("rmvspec"):
            _release_keyboard(self.fig, self._keymaps)
            self._keymaps = {}
            plt.close(self.fig)  # 'plt.show' returns at once and the figure would stay open in pyplot
            return self.x.copy(), self.y.copy()


        # ##### OPTION 2: INTERACTIVE BACKEND #####
        print(_RMVSPEC_HELP)  # Keyboard shortcuts printed when the window opens
        try:
            plt.show(block=True)  # Returns when the window is closed
        finally:
            _release_keyboard(self.fig, self._keymaps)
            self._keymaps = {}


        # ##########  OUTPUT  ##########
        return self.x[self.mask].copy(), self.y[self.mask].copy()  # Kept points only




# ======================================================================
# ============================== NORMSPEC ==============================
# ======================================================================

# ==============================
def _parse_sample(
                  sample : str | list[tuple[float, float]] | tuple[float, float] | None,
                  wmin   : float,
                  wmax   : float
    ) -> tuple[list[tuple[float, float]], bool]:
    """
    Parse the sample regions of the continuum fit: '*' / None, an IRAF range string or (min, max) pairs.

    ----------
    Parameters
    ----------
        sample : str | list[tuple[float, float]] | tuple[float, float] | None
            Wavelength ranges used in the fit:
                - '*' or None: the whole spectrum.
                - An IRAF range string 'w1:w2,w3:w4,...' (e.g.
                  '4000:4090,4150:4300').
                - A list of (w1, w2) pairs, or a single (w1, w2) pair.

        wmin : float
            Minimum wavelength of the spectrum.

        wmax : float
            Maximum wavelength of the spectrum.

    -------
    Returns
    -------
        regions : list[tuple[float, float]]
            Sorted, merged wavelength ranges (min, max).

        sample_all : bool
            True if the whole spectrum is sampled ('*' or None).

    ------
    Raises
    ------
        TypeError
            If 'sample' is not a string, a list / tuple of pairs or None, or a
            pair does not contain real numbers.

        ValueError
            If a range of the string is malformed, a pair has non-finite
            values, 'sample' is empty, no range overlaps [wmin, wmax], or
            'wmin' / 'wmax' are NaN or wmin > wmax.

    -----
    Notes
    -----
        - 'wmin' and 'wmax' may be -inf / +inf, to check only the format of
          'sample' (see '_check_fit_parameters').
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numbers
    import numpy as np


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if sample is not None and not isinstance(sample, (str, tuple, list)):
        raise TypeError(f"Sample regions 'sample' must be '*', an IRAF range string ('w1:w2,w3:w4') or a list of (min, max) pairs. Got {type(sample).__name__}.")
    #
    _check_real(wmin, "wmin", "Minimum wavelength of the spectrum", allow_inf=True)
    #
    _check_real(wmax, "wmax", "Maximum wavelength of the spectrum", allow_inf=True)
    if wmin > wmax:
        raise ValueError(f"Wavelength limits 'wmin' and 'wmax' must satisfy wmin <= wmax. Got {wmin} and {wmax}.")


    # ##########  FUNCTIONS  ##########
    def _is_real(
                 value : object
        ) -> bool:
        """
        Check whether a value is a real number (booleans excluded).

        ----------
        Parameters
        ----------
            value : object
                Value to check.

        -------
        Returns
        -------
            is_real : bool
                True if 'value' is a real number and not a boolean.
        """
        return isinstance(value, numbers.Real) and not isinstance(value, bool)


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: WHOLE SPECTRUM  ##########
    if sample is None or (isinstance(sample, str) and sample.strip() == "*"):
        return [(wmin, wmax)], True


    # ##########  OPTION 2: IRAF STRING  ##########
    pairs = []
    if isinstance(sample, str):
        for piece in sample.split(","):  # 'w1:w2,w3:w4' -> ['w1:w2', 'w3:w4']
            bounds = piece.strip().split(":")
            try:
                if len(bounds) != 2:
                    raise ValueError  # Caught below to give a single error message
                pairs.append((float(bounds[0]), float(bounds[1])))
            except ValueError as e:
                raise ValueError(f"Sample regions 'sample' contain an invalid range '{piece}'. Expected 'w1:w2' (e.g. '4000:4090,4150:4300').") from e


    # ##########  OPTION 3: (MIN, MAX) PAIRS  ##########
    else:
        if len(sample) == 2 and all(_is_real(b) for b in sample):
            sample = [sample]  # Single (w1, w2) pair
        for pair in sample:
            if not (isinstance(pair, (tuple, list)) and len(pair) == 2 and all(_is_real(b) for b in pair)):
                raise TypeError(f"Each range of the sample regions 'sample' must be a (min, max) pair of real numbers. Got {pair}.")
            pairs.append((float(pair[0]), float(pair[1])))  # Python floats


    # ##########  CHECKS  ##########
    if not pairs:
        raise ValueError("Sample regions 'sample' must contain at least one (min, max) range. Got an empty list.")
    elif not all(np.isfinite(b) for pair in pairs for b in pair):
        raise ValueError(f"Sample regions 'sample' must contain finite wavelengths. Got {pairs}.")
    regions = _merge_ranges(pairs)  # Sorted and without overlaps
    if not any(w2 >= wmin and w1 <= wmax for w1, w2 in regions):  # At least one range must overlap [wmin, wmax]
        raise ValueError(f"No range of the sample regions 'sample' overlaps the spectrum ({wmin:.3f}-{wmax:.3f}).")


    # ####################  OUTPUT  ####################
    return regions, False




# ==============================
def _check_fit_parameters(
                          function    : str,
                          order       : int,
                          sample      : str | list[tuple[float, float]] | tuple[float, float] | None,
                          naverage    : int,
                          low_reject  : float,
                          high_reject : float,
                          niterate    : int,
                          grow        : float,
                          output_type : str
    ) -> dict[str, str | int | float]:
    """
    Validate the continuum fitting parameters of 'normspec' and return them normalized.

    ----------
    Parameters
    ----------
        function : str
            Continuum function: 'spline3', 'spline1', 'legendre' or
            'chebyshev'.

        order : int
            Number of spline pieces, or number of polynomial terms (>= 1).

        sample : str | list[tuple[float, float]] | tuple[float, float] | None
            Sample regions (see '_parse_sample'). Only their format is
            checked here.

        naverage : int
            Pixels combined into each sample point (non-zero).

        low_reject : float
            Rejection limit below the fit, in rms units (>= 0).

        high_reject : float
            Rejection limit above the fit, in rms units (>= 0).

        niterate : int
            Maximum number of rejection iterations (>= 1).

        grow : float
            Rejection growing radius, in sample points (>= 0).

        output_type : str
            Output spectrum ('ratio', 'difference', 'fit' or 'data') or a
            unique abbreviation.

    -------
    Returns
    -------
        params : dict[str, str | int | float]
            A dictionary containing the following key-value pairs:
                - 'function' : str
                    Continuum function (input value).

                - 'order' : int
                    Order of the function, as a Python int.

                - 'naverage' : int
                    Pixels per sample point, as a Python int.

                - 'low_reject' : float
                    Low rejection limit (rms units), as a Python float.

                - 'high_reject' : float
                    High rejection limit (rms units), as a Python float.

                - 'niterate' : int
                    Maximum number of rejection iterations, as a Python int.

                - 'grow' : float
                    Rejection growing radius (sample points), as a Python
                    float.

                - 'output_type' : str
                    Full name of the output spectrum (abbreviations
                    resolved).

    ------
    Raises
    ------
        TypeError
            If any parameter has the wrong type.

        ValueError
            If any parameter has an invalid value.

    -----
    Notes
    -----
        - Shared by 'normspec' and '_ContinuumFit', so that both give the same
          error messages. Resolving 'output_type' twice is harmless.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(function, str):
        raise TypeError(f"Continuum function 'function' must be a string. Got {type(function).__name__}.")
    elif function not in _FUNCTIONS:
        raise ValueError(f"Continuum function 'function' must be one of {_FUNCTIONS}. Got '{function}'.")
    #
    if not isinstance(order, numbers.Integral) or isinstance(order, bool):
        raise TypeError(f"Continuum order 'order' must be an integer. Got {type(order).__name__}.")
    elif order < 1:
        raise ValueError(f"Continuum order 'order' must be >= 1. Got {order}.")
    #
    _parse_sample(sample, -np.inf, np.inf)  # Format only; the overlap with the spectrum is checked by '_ContinuumFit'
    #
    if not isinstance(naverage, numbers.Integral) or isinstance(naverage, bool):
        raise TypeError(f"Number of pixels per sample point 'naverage' must be an integer. Got {type(naverage).__name__}.")
    elif naverage == 0:
        raise ValueError(f"Number of pixels per sample point 'naverage' must be non-zero. Got {naverage}.")
    #
    if not isinstance(low_reject, numbers.Real) or isinstance(low_reject, bool):
        raise TypeError(f"Low rejection limit 'low_reject' (in rms units) must be a real number. Got {type(low_reject).__name__}.")
    elif not np.isfinite(low_reject):
        raise ValueError(f"Low rejection limit 'low_reject' (in rms units) must be finite. Got {low_reject}.")
    elif low_reject < 0:
        raise ValueError(f"Low rejection limit 'low_reject' (in rms units) must be >= 0. Got {low_reject}.")
    #
    if not isinstance(high_reject, numbers.Real) or isinstance(high_reject, bool):
        raise TypeError(f"High rejection limit 'high_reject' (in rms units) must be a real number. Got {type(high_reject).__name__}.")
    elif not np.isfinite(high_reject):
        raise ValueError(f"High rejection limit 'high_reject' (in rms units) must be finite. Got {high_reject}.")
    elif high_reject < 0:
        raise ValueError(f"High rejection limit 'high_reject' (in rms units) must be >= 0. Got {high_reject}.")
    #
    if not isinstance(niterate, numbers.Integral) or isinstance(niterate, bool):
        raise TypeError(f"Maximum number of rejection iterations 'niterate' must be an integer. Got {type(niterate).__name__}.")
    elif niterate < 1:
        raise ValueError(f"Maximum number of rejection iterations 'niterate' must be >= 1. Got {niterate}.")
    #
    if not isinstance(grow, numbers.Real) or isinstance(grow, bool):
        raise TypeError(f"Rejection growing radius 'grow' (in sample points) must be a real number. Got {type(grow).__name__}.")
    elif not np.isfinite(grow):
        raise ValueError(f"Rejection growing radius 'grow' (in sample points) must be finite. Got {grow}.")
    elif grow < 0:
        raise ValueError(f"Rejection growing radius 'grow' (in sample points) must be >= 0. Got {grow}.")
    #
    if not isinstance(output_type, str):
        raise TypeError(f"Output spectrum type 'output_type' must be a string. Got {type(output_type).__name__}.")
    matches = [t for t in _OUTPUT_TYPES if output_type and t.startswith(output_type)]  # Full names that start with the given abbreviation
    if len(matches) != 1:  # Empty, unknown or ambiguous abbreviation
        raise ValueError(f"Output spectrum type 'output_type' must be one of {_OUTPUT_TYPES} or a unique abbreviation. Got '{output_type}'.")


    # ####################  OUTPUT  ####################
    return {
            'function'    : function,
            'order'       : int(order),
            'naverage'    : int(naverage),
            'low_reject'  : float(low_reject),
            'high_reject' : float(high_reject),
            'niterate'    : int(niterate),
            'grow'        : float(grow),
            'output_type' : matches[0],
           }




# ==============================
def _sample_points(
                   wave     : np.ndarray[float],
                   flux     : np.ndarray[float],
                   regions  : list[tuple[float, float]],
                   naverage : int
    ) -> tuple[np.ndarray[float], np.ndarray[float], np.ndarray[int]]:
    """
    Build the sample points of the continuum fit by combining groups of |naverage| consecutive pixels in each region.

    ----------
    Parameters
    ----------
        wave : np.ndarray[float]
            Wavelength array (finite, sorted).

        flux : np.ndarray[float]
            Flux array (same length as 'wave'). Non-finite fluxes are skipped.

        regions : list[tuple[float, float]]
            Sorted, merged wavelength ranges (min, max) to sample.

        naverage : int
            Number of consecutive pixels combined into each sample point:
            their average if positive, their median if negative. Non-zero.

    -------
    Returns
    -------
        x_s : np.ndarray[float]
            Wavelengths of the sample points (mean wavelength of each group),
            sorted.

        y_s : np.ndarray[float]
            Fluxes of the sample points.

        pix_group : np.ndarray[int]
            Index of the sample point of each pixel (-1 = pixel not sampled).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the arrays have different lengths or 'naverage' is zero.

    -----
    Notes
    -----
        - Groups never cross region borders: the last group of a region may
          have fewer than |naverage| pixels (IRAF 'sample' + 'naverage').
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(wave, "wave", "Wavelength array")
    #
    _check_array(flux, "flux", "Flux array", size=wave.size)
    #
    _merge_ranges(regions)  # Type check of the (min, max) ranges only
    #
    if not isinstance(naverage, numbers.Integral) or isinstance(naverage, bool):
        raise TypeError(f"Number of pixels per sample point 'naverage' must be an integer. Got {type(naverage).__name__}.")
    elif naverage == 0:
        raise ValueError(f"Number of pixels per sample point 'naverage' must be non-zero. Got {naverage}.")


    # ####################  PARAMETERS  ####################
    m              = abs(naverage)  # Pixels per sample point
    finite         = np.isfinite(flux)  # Non-finite fluxes are never sampled
    pix_group      = np.full(wave.size, -1)  # -1: pixel not sampled
    x_list, y_list = [], []  # Sample points of each region
    n_groups       = 0  # Sample points created so far (offset of the group indices of the next region)


    # ####################  COMPUTATION  ####################
    for w1, w2 in regions:
        idx = np.flatnonzero((wave >= w1) & (wave <= w2) & finite)  # Usable pixels of the region
        if idx.size == 0:
            continue
        gid            = np.arange(idx.size) // m  # Group of each pixel within the region
        pix_group[idx] = n_groups + gid  # Global index of the sample point of each pixel
        n_groups      += gid[-1] + 1
        #
        counts = np.bincount(gid)  # Pixels in each group (the last one may have fewer than m)
        x_list.append(np.bincount(gid, weights=wave[idx]) / counts)  # Mean wavelength of each group
        if m == 1:  # One pixel per sample point
            y_list.append(flux[idx])
        elif naverage > 0:  # Mean flux of each group
            y_list.append(np.bincount(gid, weights=flux[idx]) / counts)
        else:  # Median flux of each group
            y_list.append(np.array([np.median(c) for c in np.split(flux[idx], np.arange(m, idx.size, m))]))


    # ####################  OUTPUT  ####################
    if not x_list:  # No usable pixel in any region
        return np.empty(0), np.empty(0), pix_group
    return np.concatenate(x_list), np.concatenate(y_list), pix_group




# ==============================
def _fit_function(
                  x        : np.ndarray[float],
                  y        : np.ndarray[float],
                  function : str,
                  order    : int,
                  wmin     : float,
                  wmax     : float,
                  quiet    : bool = False
    ) -> Callable[[np.ndarray[float]], np.ndarray[float]]:
    """
    Fit the continuum function to the sample points by least squares (IRAF 'function' and 'order').

    ----------
    Parameters
    ----------
        x : np.ndarray[float]
            Sample wavelengths (finite, sorted).

        y : np.ndarray[float]
            Sample fluxes (finite, same length as 'x').

        function : str
            Continuum function: 'spline3', 'spline1', 'legendre' or
            'chebyshev'.

        order : int
            Number of spline pieces of equal length over [wmin, wmax], or
            number of polynomial terms (degree + 1).

        wmin : float
            Minimum wavelength of the spectrum (polynomial domain and spline
            interval).

        wmax : float
            Maximum wavelength of the spectrum (polynomial domain and spline
            interval).

        quiet : bool, optional
            If True, do not warn when the spline knots are moved to sample
            quantiles. Default is False.

    -------
    Returns
    -------
        fit : Callable[[np.ndarray[float]], np.ndarray[float]]
            Fitted function: fit(wave) returns the continuum.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the arrays have different lengths, 'function' is unknown,
            'order' < 1, or the wavelength limits are not finite with
            wmin < wmax.

    -----
    Notes
    -----
        - If some spline pieces have no sample points (e.g. gaps between
          sample regions), the knots are placed at quantiles of the points
          instead, so the Schoenberg-Whitney conditions hold (Schoenberg &
          Whitney 1953).

    ----------
    References
    ----------
        - Schoenberg, I. J. & Whitney, A. 1953,
          Trans. Amer. Math. Soc., 74, 246.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numbers
    import numpy as np
    from scipy.interpolate import LSQUnivariateSpline
    import warnings


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(x, "x", "Sample wavelengths")
    #
    _check_array(y, "y", "Sample fluxes", size=x.size)
    #
    if not isinstance(function, str):
        raise TypeError(f"Continuum function 'function' must be a string. Got {type(function).__name__}.")
    elif function not in _FUNCTIONS:
        raise ValueError(f"Continuum function 'function' must be one of {_FUNCTIONS}. Got '{function}'.")
    #
    if not isinstance(order, numbers.Integral) or isinstance(order, bool):
        raise TypeError(f"Continuum order 'order' must be an integer. Got {type(order).__name__}.")
    elif order < 1:
        raise ValueError(f"Continuum order 'order' must be >= 1. Got {order}.")
    #
    _check_real(wmin, "wmin", "Minimum wavelength of the spectrum")
    #
    _check_real(wmax, "wmax", "Maximum wavelength of the spectrum")
    if wmin >= wmax:
        raise ValueError(f"Wavelength limits 'wmin' and 'wmax' must satisfy wmin < wmax. Got {wmin} and {wmax}.")
    #
    if not isinstance(quiet, bool):
        raise TypeError(f"Quiet flag 'quiet' must be a boolean. Got {type(quiet).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: POLYNOMIALS  ##########
    if function in ("legendre", "chebyshev"):
        poly = np.polynomial.Legendre if function == "legendre" else np.polynomial.Chebyshev
        return poly.fit(x, y, deg=order - 1, domain=[wmin, wmax])  # 'order' terms = degree + 1; domain mapped to [-1, 1]


    # ##########  OPTION 2: SPLINES  ##########
    k     = 3 if function == "spline3" else 1  # Spline degree
    knots = np.linspace(wmin, wmax, order + 1)[1:-1]  # Interior knots of 'order' equal pieces
    try:
        return LSQUnivariateSpline(x, y, knots, bbox=[wmin, wmax], k=k)  # Least-squares spline with fixed knots
    except ValueError:  # Schoenberg-Whitney conditions not satisfied (pieces without sample points)
        if not quiet:
            warnings.warn(f"Some of the {order} equal spline pieces have no sample points; knots placed at sample quantiles instead.")
        knots = np.unique(np.quantile(x, np.linspace(0., 1., order + 1)[1:-1]))  # Same number of points in every piece
        knots = knots[(knots > x[0]) & (knots < x[-1])]  # Interior knots only
        return LSQUnivariateSpline(x, y, knots, bbox=[wmin, wmax], k=k)




# ==============================
def _fit_continuum(
                   wave        : np.ndarray[float],
                   x_s         : np.ndarray[float],
                   y_s         : np.ndarray[float],
                   deleted     : np.ndarray[bool],
                   function    : str,
                   order       : int,
                   low_reject  : float,
                   high_reject : float,
                   niterate    : int,
                   grow        : float
    ) -> tuple[np.ndarray[float], np.ndarray[bool], float]:
    """
    Fit the continuum to the sample points with iterative sigma rejection.

    ----------
    Parameters
    ----------
        wave : np.ndarray[float]
            Wavelength array (sorted) where the continuum is evaluated.

        x_s : np.ndarray[float]
            Wavelengths of the sample points (sorted).

        y_s : np.ndarray[float]
            Fluxes of the sample points.

        deleted : np.ndarray[bool]
            Sample points deleted by the user (never used in the fit).

        function : str
            Continuum function (see '_fit_function').

        order : int
            Order of the function (see '_fit_function').

        low_reject : float
            Rejection limit below the fit, in units of the rms of the
            residuals (0 disables it).

        high_reject : float
            Rejection limit above the fit, in units of the rms of the
            residuals (0 disables it).

        niterate : int
            Maximum number of rejection iterations.

        grow : float
            Rejection growing radius, in sample points (only its integer part
            is used).

    -------
    Returns
    -------
        continuum : np.ndarray[float]
            Continuum evaluated over 'wave'.

        rejected : np.ndarray[bool]
            Sample points rejected by the iterations.

        rms : float
            Rms of the residuals of the final fit (flux units).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If any input has an invalid value (arrays of different lengths,
            unknown 'function', negative rejection limits or 'grow',
            'order' or 'niterate' < 1), or there are fewer usable sample
            points than free coefficients.

    -----
    Notes
    -----
        - Each iteration rejects the points with residuals below
          -low_reject * rms or above high_reject * rms, plus the points within
          'grow' sample points of them, and refits. It stops after 'niterate'
          iterations, when no new points are rejected, or when too few points
          would remain (with a warning). Rejected points stay rejected.

        - The rms uses N - 1 degrees of freedom.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(wave, "wave", "Wavelength array")
    if wave.size < 2:
        raise ValueError(f"Wavelength array 'wave' must contain at least 2 points. Got {wave.size}.")
    #
    _check_array(x_s, "x_s", "Sample wavelengths")
    #
    _check_array(y_s, "y_s", "Sample fluxes", size=x_s.size)
    #
    _check_array(deleted, "deleted", "Mask of deleted sample points", kind="bool", size=x_s.size)
    #
    _check_fit_parameters(function, order, None, 1, low_reject, high_reject, niterate, grow, "ratio")  # 'sample', 'naverage' and 'output_type' are not used here


    # ##########  SAMPLE POINTS  ##########
    spline_degrees = {
                      'spline3' : 3,
                      'spline1' : 1,
                     }
    k     = spline_degrees.get(function, 0)  # 0 for polynomials
    ncoef = order + k if k else order  # Number of free coefficients
    use   = ~deleted  # Sample points available to the fit
    if use.sum() < ncoef:
        raise ValueError(f"Too few sample points ({int(use.sum())}) for a '{function}' continuum of order {order} ({ncoef} needed).")


    # ####################  PARAMETERS  ####################
    wmin, wmax = float(wave[0]), float(wave[-1])  # Domain of the function (whole spectrum)
    rejected   = np.zeros(x_s.size, dtype=bool)  # Accumulated over the iterations


    # ####################  COMPUTATION  ####################
    # ##########  FIRST FIT  ##########
    fit = _fit_function(x_s[use], y_s[use], function, order, wmin, wmax)


    # ##########  REJECTION ITERATIONS  ##########
    for _ in range(niterate):
        if low_reject <= 0 and high_reject <= 0:  # Rejection disabled
            break
        good  = use & ~rejected  # Points used in the current fit
        resid = y_s - fit(x_s)
        sigma = np.sqrt(np.sum(resid[good]**2) / max(good.sum() - 1, 1))  # Rms of the residuals (N - 1 degrees of freedom)
        new   = good & (((low_reject > 0) & (resid < -low_reject * sigma)) | ((high_reject > 0) & (resid > high_reject * sigma)))  # Newly rejected points
        #
        if grow >= 1 and new.any():  # Also reject the sample points within 'grow' points of each rejected one
            idx   = np.flatnonzero(new)
            cover = np.zeros(x_s.size + 1, dtype=int)  # Difference array: +1 at the start and -1 after the end of each grown interval
            np.add.at(cover, np.clip(idx - int(grow), 0, x_s.size), 1)
            np.add.at(cover, np.clip(idx + int(grow) + 1, 0, x_s.size), -1)
            new |= good & (np.cumsum(cover)[:-1] > 0)  # Points covered by at least one grown interval
        #
        if not new.any():
            break  # Converged
        elif (good & ~new).sum() < ncoef:
            warnings.warn("Rejection stopped: too few sample points would remain for the fit.")
            break
        #
        rejected |= new  # Rejected points stay rejected
        keep      = use & ~rejected
        fit       = _fit_function(x_s[keep], y_s[keep], function, order, wmin, wmax, quiet=True)  # Refit without the rejected points (knot warning already issued)


    # ##########  FINAL RESIDUALS  ##########
    good  = use & ~rejected
    resid = y_s[good] - fit(x_s[good])
    rms   = float(np.sqrt(np.sum(resid**2) / max(good.sum() - 1, 1)))


    # ####################  OUTPUT  ####################
    return np.asarray(fit(wave), dtype=float), rejected, rms  # Continuum evaluated at every pixel




# ==============================
def _pixel_mask(
                mask_s    : np.ndarray[bool],
                pix_group : np.ndarray[int]
    ) -> np.ndarray[bool]:
    """
    Map a mask of sample points onto the spectrum pixels.

    ----------
    Parameters
    ----------
        mask_s : np.ndarray[bool]
            Mask of the sample points.

        pix_group : np.ndarray[int]
            Index of the sample point of each pixel (-1 = pixel not sampled).

    -------
    Returns
    -------
        mask : np.ndarray[bool]
            Mask of the pixels: the value of their sample point, False for
            pixels not sampled.

    ------
    Raises
    ------
        TypeError
            If 'mask_s' is not a boolean array or 'pix_group' is not an
            integer array.

        ValueError
            If 'pix_group' points to sample points outside 'mask_s'.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(mask_s, "mask_s", "Mask of the sample points", kind="bool")
    #
    _check_array(pix_group, "pix_group", "Sample point of each pixel", kind="int")
    if pix_group.size and (pix_group.min() < -1 or pix_group.max() >= mask_s.size):
        raise ValueError(f"Sample point of each pixel 'pix_group' must lie in [-1, {mask_s.size - 1}]. Got [{pix_group.min()}, {pix_group.max()}].")


    # ####################  COMPUTATION  ####################
    mask          = np.zeros(pix_group.size, dtype=bool)  # Pixels not sampled stay False
    sampled       = pix_group >= 0
    mask[sampled] = mask_s[pix_group[sampled]]  # Each sampled pixel takes the value of its sample point


    # ####################  OUTPUT  ####################
    return mask




# ==============================
def _output_spectrum(
                     flux        : np.ndarray[float],
                     continuum   : np.ndarray[float],
                     pix_rej     : np.ndarray[bool] | None,
                     output_type : str,
                     replace     : bool
    ) -> np.ndarray[float]:
    """
    Compute the output spectrum ('ratio', 'difference', 'fit' or 'data') from the continuum.

    ----------
    Parameters
    ----------
        flux : np.ndarray[float]
            Flux array.

        continuum : np.ndarray[float]
            Continuum array (same length as 'flux').

        pix_rej : np.ndarray[bool] | None
            Mask of the rejected / deleted pixels (True = rejected), or None.

        output_type : str
            Output spectrum: 'ratio' (flux / continuum), 'difference'
            (flux - continuum), 'fit' (continuum) or 'data' (flux).

        replace : bool
            If True, the rejected pixels are replaced by the continuum before
            computing the 'ratio', 'difference' and 'data' outputs.

    -------
    Returns
    -------
        output : np.ndarray[float]
            Output spectrum (same length as 'flux').

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the arrays have different lengths or 'output_type' is unknown.

    -----
    Notes
    -----
        - Ratios where the continuum is <= 0 are set to NaN (with a warning).
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(flux, "flux", "Flux array")
    #
    _check_array(continuum, "continuum", "Continuum array", size=flux.size)
    #
    if pix_rej is not None:
        _check_array(pix_rej, "pix_rej", "Mask of the rejected pixels", kind="bool", size=flux.size)
    #
    if not isinstance(output_type, str):
        raise TypeError(f"Output spectrum type 'output_type' must be a string. Got {type(output_type).__name__}.")
    elif output_type not in _OUTPUT_TYPES:
        raise ValueError(f"Output spectrum type 'output_type' must be one of {_OUTPUT_TYPES}. Got '{output_type}'.")
    #
    if not isinstance(replace, bool):
        raise TypeError(f"Replace flag 'replace' must be a boolean. Got {type(replace).__name__}.")


    # ####################  COMPUTATION  ####################
    data = flux.copy()  # The input flux is never modified
    if replace and pix_rej is not None:
        data[pix_rej] = continuum[pix_rej]  # Rejected / deleted pixels replaced by the fit (IRAF 'replace')
    #
    if output_type == "fit":
        return continuum.copy()
    elif output_type == "data":
        return data
    elif output_type == "difference":
        return data - continuum
    #
    bad = continuum <= 0  # The ratio is undefined or meaningless there
    if bad.any():
        warnings.warn(f"Continuum <= 0 at {int(bad.sum())} pixel(s); ratio set to NaN there.")
    with np.errstate(divide="ignore", invalid="ignore"):  # Divisions by zero are set to NaN just below
        output = data / continuum  # Normalized flux
    output[bad] = np.nan


    # ####################  OUTPUT  ####################
    return output




# ==============================
class _ContinuumFit:
    """
    Spectrum, fitting parameters and current continuum of 'normspec', shared by its non-interactive and interactive modes.

    ----------
    Parameters
    ----------
        wave : np.ndarray[float]
            Wavelength array (finite, sorted).

        flux : np.ndarray[float]
            Flux array (same length as 'wave').

        function : str
            Continuum function (keyword-only, as the following arguments).

        order : int
            Number of spline pieces, or number of polynomial terms.

        sample : str | list[tuple[float, float]] | tuple[float, float] | None
            Sample regions (see '_parse_sample').

        naverage : int
            Pixels combined into each sample point.

        low_reject : float
            Rejection limit below the fit, in rms units.

        high_reject : float
            Rejection limit above the fit, in rms units.

        niterate : int
            Maximum number of rejection iterations.

        grow : float
            Rejection growing radius, in sample points.

        output_type : str
            Output spectrum.

        replace : bool
            Replace the rejected pixels by the fit before computing the
            output.

    -----
    Notes
    -----
        - 'mode' is 'auto' (fit to the sample points) or 'manual' (cubic
          spline through clicked points); 'continuum' and 'result' always
          come from the last successful fit of either mode.

        - The parameters of the last successful automatic fit are stored in
          'fit_info', so the status line and the header of the output file
          describe the continuum actually shown / saved even if a later fit
          fails.

        - The output spectrum is written with 'writespec' (see '_save').
    """

    # ####################  INIT  ####################
    def __init__(
                 self,
                 wave        : np.ndarray[float],
                 flux        : np.ndarray[float],
                 *,
                 function    : str,
                 order       : int,
                 sample      : str | list[tuple[float, float]] | tuple[float, float] | None,
                 naverage    : int,
                 low_reject  : float,
                 high_reject : float,
                 niterate    : int,
                 grow        : float,
                 output_type : str,
                 replace     : bool
        ) -> None:
        """
        Store the spectrum and the fixed parameters, parse the sample and set the resettable parameters.

        ----------
        Parameters
        ----------
            wave : np.ndarray[float]
                Wavelength array (finite, sorted).

            flux : np.ndarray[float]
                Flux array (same length as 'wave').

            function : str
                Continuum function: 'spline3', 'spline1', 'legendre' or
                'chebyshev'.

            order : int
                Number of spline pieces, or number of polynomial terms.

            sample : str | list[tuple[float, float]] | tuple[float, float] | None
                Sample regions (see '_parse_sample').

            naverage : int
                Pixels combined into each sample point (see '_sample_points').

            low_reject : float
                Rejection limit below the fit, in rms units.

            high_reject : float
                Rejection limit above the fit, in rms units.

            niterate : int
                Maximum number of rejection iterations.

            grow : float
                Rejection growing radius, in sample points.

            output_type : str
                Output spectrum (see '_output_spectrum').

            replace : bool
                Replace the rejected pixels by the fit before computing the
                output.

        -------
        Returns
        -------
            None
                The state is initialized with no continuum.

        ------
        Raises
        ------
            TypeError
                If any input has the wrong type.

            ValueError
                If any input has an invalid value (see
                '_check_fit_parameters' and '_parse_sample'), or 'wave' has
                fewer than two finite, sorted points.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        _check_array(wave, "wave", "Wavelength array")
        if wave.size < 2:
            raise ValueError(f"Wavelength array 'wave' must contain at least 2 points. Got {wave.size}.")
        elif not np.all(np.isfinite(wave)):
            raise ValueError("Wavelength array 'wave' must contain only finite values.")
        elif np.any(np.diff(wave) < 0):
            raise ValueError("Wavelength array 'wave' must be sorted in increasing order.")
        #
        _check_array(flux, "flux", "Flux array", size=wave.size)
        #
        params = _check_fit_parameters(function, order, sample, naverage, low_reject, high_reject, niterate, grow, output_type)
        #
        if not isinstance(replace, bool):
            raise TypeError(f"Replace flag 'replace' must be a boolean. Got {type(replace).__name__}.")


        # ##########  PARAMETERS  ##########
        # ##### DATA & FIXED PARAMETERS #####
        self.wave, self.flux              = wave, flux
        self.wmin, self.wmax              = float(wave[0]), float(wave[-1])
        self.naverage                     = params['naverage']
        self.niterate                     = params['niterate']
        self.grow                         = params['grow']
        self.low_reject, self.high_reject = params['low_reject'], params['high_reject']
        self.output_type, self.replace    = params['output_type'], replace


        # ##### INITIAL (RESETTABLE) PARAMETERS #####
        regions, sample_all = _parse_sample(sample, self.wmin, self.wmax)
        self._initial       = {  # Parameters restored by '_reset' ('r' in the interactive window)
                               'function'   : params['function'],
                               'order'      : params['order'],
                               'regions'    : regions,
                               'sample_all' : sample_all,
                              }
        self._reset()


    # ####################  SAMPLE  ####################
    def _reset(
               self
        ) -> None:
        """
        Restore the initial parameters and remove the continuum.

        -------
        Returns
        -------
            None
                The parameters, sample points and fit state are reset.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  COMPUTATION  ##########
        self.function   = self._initial['function']
        self.order      = self._initial['order']
        self.regions    = list(self._initial['regions'])
        self.sample_all = self._initial['sample_all']
        self.mode       = "auto"  # 'auto' (fit to the sample points) or 'manual' (spline through clicked points)
        self.continuum  = None  # No continuum until the first fit
        self.result     = None  # Output spectrum of the last successful fit
        self.rms        = np.nan
        self.n_manual   = 0  # Number of manual continuum points
        self.fit_info   = None  # Parameters of the last successful automatic fit
        self._resample()
    #
    def _resample(
                  self
        ) -> None:
        """
        Recompute the sample points from the current regions; deleted points are forgotten.

        -------
        Returns
        -------
            None
                The sample points and their deleted / rejected masks are reset.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  COMPUTATION  ##########
        self.x_s, self.y_s, self.pix_group = _sample_points(self.wave, self.flux, self.regions, self.naverage)
        self.deleted                       = np.zeros(self.x_s.size, dtype=bool)  # Sample points deleted by the user ('d')
        self.rejected                      = np.zeros(self.x_s.size, dtype=bool)  # Sample points rejected by the sigma clipping
    #
    def _sample_string(
                       self
        ) -> str:
        """
        Return the current sample regions as an IRAF range string.

        -------
        Returns
        -------
            sample : str
                '*' for the whole spectrum, or 'w1:w2,w3:w4,...'.
        """

        # ##########  OUTPUT  ##########
        return "*" if self.sample_all else ",".join(f"{a:.3f}:{b:.3f}" for a, b in self.regions)
    #
    def _add_sample_region(
                           self,
                           w1 : float,
                           w2 : float
        ) -> None:
        """
        Add a sample region; the first one replaces the whole-spectrum sample ('*').

        ----------
        Parameters
        ----------
            w1 : float
                One end of the region, in wavelength units.

            w2 : float
                Other end of the region, in wavelength units.

        -------
        Returns
        -------
            None
                The regions and the sample points are updated.

        ------
        Raises
        ------
            TypeError
                If 'w1' or 'w2' is not a real number.

            ValueError
                If 'w1' or 'w2' is not finite.
        """

        # ##########  MODULES & PACKAGES  ##########
        # ##### PANORAMA #####
        from ..units.units_support import _check_real


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        w1 = _check_real(w1, "w1", "Sample region end")
        #
        w2 = _check_real(w2, "w2", "Sample region end")


        # ##########  COMPUTATION  ##########
        regions                       = [(w1, w2)] if self.sample_all else self.regions + [(w1, w2)]  # The first region replaces the whole-spectrum sample
        self.regions, self.sample_all = _merge_ranges(regions), False
        self._resample()
    #
    def _sample_whole_spectrum(
                               self
        ) -> None:
        """
        Reset the sample to the whole spectrum ('*').

        -------
        Returns
        -------
            None
                The regions and the sample points are updated.
        """

        # ##########  COMPUTATION  ##########
        self.regions, self.sample_all = [(self.wmin, self.wmax)], True
        self._resample()


    # ####################  FIT  ####################
    def _fit(
             self
        ) -> None:
        """
        Fit the continuum automatically with the current parameters.

        -------
        Returns
        -------
            None
                'continuum', 'result', 'rejected', 'rms' and 'fit_info' are
                updated.

        ------
        Raises
        ------
            ValueError
                If the fit fails (e.g. too few sample points); the previous
                continuum is kept.
        """

        # ##########  COMPUTATION  ##########
        continuum, rejected, rms = _fit_continuum(self.wave, self.x_s, self.y_s, self.deleted, self.function, self.order, self.low_reject, self.high_reject, self.niterate, self.grow)  # Raises ValueError on failure: the previous state is kept
        pix_rej                  = _pixel_mask(rejected | self.deleted, self.pix_group)  # Rejected / deleted sample points mapped onto pixels


        # ##########  OUTPUT  ##########
        self.mode      = "auto"
        self.continuum = continuum
        self.rejected  = rejected
        self.rms       = rms
        self.result    = _output_spectrum(self.flux, continuum, pix_rej, self.output_type, self.replace)
        self.fit_info  = {  # Describes the continuum shown, even if a later fit fails
                          'function'   : self.function,
                          'order'      : self.order,
                          'sample'     : self._sample_string(),
                          'rms'        : rms,
                          'n_points'   : int(self.x_s.size),
                          'n_rejected' : int(rejected.sum()),
                          'n_deleted'  : int(self.deleted.sum()),
                         }
    #
    def _fit_manual(
                    self,
                    x : np.ndarray[float],
                    y : np.ndarray[float]
        ) -> None:
        """
        Set a manual continuum: cubic spline through the points (x, y).

        ----------
        Parameters
        ----------
            x : np.ndarray[float]
                Wavelengths of the manual continuum points (any order).

            y : np.ndarray[float]
                Fluxes of the manual continuum points (same length as 'x').

        -------
        Returns
        -------
            None
                'mode', 'continuum', 'result' and 'n_manual' are updated.

        ------
        Raises
        ------
            TypeError
                If 'x' or 'y' is not a 1D numeric NumPy array.

            ValueError
                If the arrays have different lengths, contain non-finite
                values, or there are fewer than two points with different
                wavelengths.

        -----
        Notes
        -----
            - Points with duplicated wavelengths: only the first one is used.

            - Outside the manual points the spline is extrapolated.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np
        from scipy.interpolate import CubicSpline


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        _check_array(x, "x", "Manual continuum wavelengths")
        if not np.all(np.isfinite(x)):
            raise ValueError("Manual continuum wavelengths 'x' must contain only finite values.")
        #
        _check_array(y, "y", "Manual continuum fluxes", size=x.size)
        if not np.all(np.isfinite(y)):
            raise ValueError("Manual continuum fluxes 'y' must contain only finite values.")


        # ##########  COMPUTATION  ##########
        x_unique, idx = np.unique(x.astype(float), return_index=True)  # Sorted, without duplicated wavelengths
        if x_unique.size < 2:
            raise ValueError(f"Manual continuum points must include at least two different wavelengths. Got {x_unique.size}.")


        # ##########  OUTPUT  ##########
        self.mode, self.n_manual = "manual", int(x_unique.size)
        self.continuum           = CubicSpline(x_unique, y.astype(float)[idx])(self.wave)  # Interpolating cubic spline, extrapolated outside the points
        self.result              = _output_spectrum(self.flux, self.continuum, None, self.output_type, self.replace)


    # ####################  REPORT & SAVE  ####################
    def _status(
                self
        ) -> str:
        """
        Describe the current continuum in one line.

        -------
        Returns
        -------
            status : str
                Description of the continuum shown (function, parameters,
                rms, rejected / deleted points).
        """

        # ##########  OUTPUT  ##########
        if self.continuum is None:
            return "No continuum fitted"
        elif self.mode == "manual":
            return f"Manual continuum: cubic spline through {self.n_manual} points"
        info = self.fit_info
        return (f"{info['function']}, order={info['order']} | naverage={self.naverage} | "
                f"reject -{self.low_reject}/+{self.high_reject} sigma, niterate={self.niterate}, grow={self.grow} | "
                f"rms={info['rms']:.3g} | rejected {info['n_rejected']}, deleted {info['n_deleted']} of {info['n_points']} points"
               )
    #
    def _header(
                self
        ) -> str:
        """
        Build the header of the output file, describing the continuum.

        -------
        Returns
        -------
            header : str
                Header lines (without the '#' prefix): output type and
                continuum fit. The column names are added by 'writespec'.
        """

        # ##########  COMPUTATION  ##########
        if self.mode == "manual":
            fit_text = f"manual cubic spline through {self.n_manual} points"
        else:
            info     = self.fit_info
            fit_text = (f"function={info['function']}, order={info['order']}, sample={info['sample']}, "
                        f"naverage={self.naverage}, low_reject={self.low_reject}, high_reject={self.high_reject}, "
                        f"niterate={self.niterate}, grow={self.grow}"
                       )


        # ##########  OUTPUT  ##########
        return f"normspec: output_type={self.output_type}, replace={self.replace}\ncontinuum: {fit_text}"
    #
    def _save(
              self,
              output_file : Path | None
        ) -> None:
        """
        Write the output spectrum (wavelength, output flux) with 'writespec', with a header describing the fit.

        ----------
        Parameters
        ----------
            output_file : Path | None
                Output file (in an existing directory, with an extension
                supported by 'writespec'). If None, nothing is written (a
                message is printed).

        -------
        Returns
        -------
            None
                The spectrum is written to 'output_file' (ASCII or FITS, from
                its extension), replacing the previous file.

        ------
        Raises
        ------
            TypeError
                If 'output_file' is neither a pathlib.Path nor None.

            ValueError
                If the extension of 'output_file' is not supported, or the
                spectrum cannot be written (e.g. repeated wavelengths).

            RuntimeError
                If the file cannot be written for any other reason (e.g. no
                write permission).
        """

        # ##########  MODULES & PACKAGES  ##########
        # ##### PYTHON #####
        from pathlib import Path


        # ##### PANORAMA #####
        from ...spec import writespec


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        if output_file is not None:
            if not isinstance(output_file, Path):
                raise TypeError(f"Output file 'output_file' must be a pathlib.Path (already validated by 'normspec') or None. Got {type(output_file).__name__}.")
            _spectrum_file_format(output_file)


        # ##########  OUTPUT  ##########
        if output_file is None:
            print("No 'output_file' given: nothing saved.")
        elif self.result is None:
            print("No continuum fitted yet: nothing saved.")
        else:
            writespec((self.wave, self.result), output_file, overwrite=True, header=self._header())  # Replaced at every save; the header describes the fit
            print(f"Saved {self.output_type} spectrum to '{output_file}'.")




# ==============================
class _InteractiveNormalization:
    """
    Interactive window of 'normspec': spectrum, sample regions and points, rejected / deleted points and continuum (upper panel), and output spectrum (lower panel).

    ----------
    Parameters
    ----------
        fit : _ContinuumFit
            Continuum fit to edit (modified in place).

        output_file : Path | None
            File written with 'w', or None (saving disabled).

        markrej : bool, optional
            Mark the rejected sample points. Default is True.

    -----
    Notes
    -----
        - The controls are listed in '_NORMSPEC_HELP' (printed with '?').

        - Matplotlib's default keyboard shortcuts clash with the fitting keys
          ('s' save, 'f' fullscreen, 'l' log scale, 'r' home...), so they
          are disabled for this figure only; the toolbar buttons still work.
    """

    # ####################  CONFIGURATION  ####################
    KEYS = ("enter", "f", "+", "=", "-", "1", "3", "l", "c", "s", "t", "d", "u", "r", "w", "q", "?")


    # ####################  INIT  ####################
    def __init__(
                 self,
                 fit         : _ContinuumFit,
                 output_file : Path | None,
                 markrej     : bool = True
        ) -> None:
        """
        Store the continuum fit and the plotting options.

        ----------
        Parameters
        ----------
            fit : _ContinuumFit
                Continuum fit to edit (modified in place).

            output_file : Path | None
                File written with 'w'. None disables saving.

            markrej : bool, optional
                Mark the rejected sample points. Default is True.

        -------
        Returns
        -------
            None
                The state of the window is initialized.

        ------
        Raises
        ------
            TypeError
                If any input has the wrong type.
        """

        # ##########  MODULES & PACKAGES  ##########
        from pathlib import Path


        # ##########  INPUTS VALIDATION  ##########
        # ##### VALUES #####
        if not isinstance(fit, _ContinuumFit):
            raise TypeError(f"Continuum fit 'fit' must be a _ContinuumFit. Got {type(fit).__name__}.")
        #
        if output_file is not None and not isinstance(output_file, Path):
            raise TypeError(f"Output file 'output_file' must be a pathlib.Path (already validated by 'normspec') or None. Got {type(output_file).__name__}.")
        #
        if not isinstance(markrej, bool):
            raise TypeError(f"Mark rejected points flag 'markrej' must be a boolean. Got {type(markrej).__name__}.")


        # ##########  PARAMETERS  ##########
        self.fit            = fit
        self.output_file    = output_file
        self.markrej        = markrej
        self.manual_artists = []  # Manual continuum points (Line2D)
        self.span_artists   = []  # Shaded sample regions
        self.sample_start   = None  # Wavelength of the first 's' press of a sample region
        self._keymaps       = {}  # Original keymaps of the toolmanager tools (restored on close)


    # ####################  RUN  ####################
    def _run(
             self
        ) -> None:
        """
        Build the figure, do the initial fit and show the window.

        -------
        Returns
        -------
            None
                The continuum fit is edited in place.

        -----
        Notes
        -----
            - With blocking backends, it returns when the window is closed.
              With non-blocking ones (e.g. Jupyter widgets) it returns at once
              and the controls keep working: the figure holds a reference to
              this object.

            - With a non-interactive backend (e.g. Agg, inline) a
              RuntimeWarning is issued, the initial automatic fit is shown
              (if the backend displays figures) and the figure is closed.
        """

        # ##########  MODULES & PACKAGES  ##########
        import matplotlib.pyplot as plt


        # ##########  COMPUTATION  ##########
        interactive = _check_interactive_backend("normspec")  # Warns if the backend cannot receive events
        self._build_figure()
        if interactive:
            print(_NORMSPEC_HELP)  # Controls printed when the window opens
        self._refit()  # Initial automatic fit
        if interactive:
            plt.show(block=True)  # Returns when the window is closed, also in IPython / Jupyter, where a plain 'plt.show' does not block
        else:
            plt.show()  # Draws the figure (inline in Jupyter) and returns at once
        if not interactive:  # No close event will restore the keymaps
            _release_keyboard(self.fig, self._keymaps)
            self._keymaps = {}
            plt.close(self.fig)  # 'plt.show' returns at once and the figure would stay open in pyplot


    # ####################  FIGURE  ####################
    def _build_figure(
                      self
        ) -> None:
        """
        Create the two panels, their artists and the event connections.

        -------
        Returns
        -------
            None
                The figure and its artists are stored as attributes.
        """

        # ##########  MODULES & PACKAGES  ##########
        import matplotlib.pyplot as plt
        import numpy as np


        # ##########  PARAMETERS  ##########
        fit = self.fit
        reference_levels = {
                            'ratio'      : 1.,
                            'difference' : 0.,
                           }
        output_labels    = {
                            'ratio'      : "Normalized flux",
                            'difference' : "Flux - continuum",
                            'fit'        : "Continuum",
                            'data'       : "Flux",
                           }


        # ##########  PLOTTING  ##########
        # ##### PANELS #####
        self.fig, (self.ax, self.axn) = plt.subplots(2, 1, sharex=True, figsize=(11, 7), gridspec_kw={'height_ratios' : [3, 1]})
        self.fig.subplots_adjust(hspace=0.05)  # Panels almost touching (shared wavelength axis)


        # ##### UPPER PANEL: SPECTRUM & CONTINUUM #####
        self.ax.plot(fit.wave, fit.flux, "k-", lw=0.8, label="spectrum")  # Static curve; the other artists are updated by '_update_plot'
        self.good_pts,  = self.ax.plot([], [], "o", color="tab:blue", ms=3, ls="none", label="sample points")
        self.rej_pts,   = self.ax.plot([], [], "x", color="tab:red", ms=6, ls="none", label="rejected")
        self.del_pts,   = self.ax.plot([], [], "+", color="tab:orange", ms=8, ls="none", label="deleted")
        self.cont_line, = self.ax.plot([], [], "r-", lw=1.5, label="continuum")
        self.ax.set_ylabel("Flux")
        self.ax.legend(loc="best", fontsize=8)


        # ##### LOWER PANEL: OUTPUT SPECTRUM #####
        self.out_line, = self.axn.plot(fit.wave, np.full(fit.wave.size, np.nan), "k-", lw=0.8)  # Empty (NaN) until the first fit
        if fit.output_type in reference_levels:
            self.axn.axhline(reference_levels[fit.output_type], color="tab:red", ls=":", lw=1)  # 1 for 'ratio', 0 for 'difference'
        self.axn.set_xlabel("Wavelength [Å]")
        self.axn.set_ylabel(output_labels[fit.output_type])


        # ##### EVENTS #####
        # Matplotlib keeps only weak references to bound-method callbacks: the figure holds this object
        self.fig._normspec_session = self
        self._keymaps              = _take_keyboard(self.fig, self.KEYS)[1]  # Default shortcuts disabled, not forwarded
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)  # Left-click: add a manual point
        self.fig.canvas.mpl_connect("pick_event", self._on_pick)  # Right-click on a manual point: remove it
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)  # Fitting keys
        self.fig.canvas.mpl_connect("close_event", self._on_close)  # Keymap restoration
    #
    def _refit(
               self
        ) -> None:
        """
        Fit the continuum automatically (errors are printed, not raised) and redraw.

        -------
        Returns
        -------
            None
                The fit and the figure are updated.
        """

        # ##########  COMPUTATION  ##########
        try:
            self.fit._fit()
        except ValueError as e:  # Never raise inside a Matplotlib callback: the previous continuum is kept
            print(f"Continuum fit failed: {e}")
        self._update_plot()
    #
    def _update_plot(
                     self
        ) -> None:
        """
        Redraw the sample regions and points, the continuum and the output spectrum, keeping the zoom.

        -------
        Returns
        -------
            None
                The artists of the figure are updated.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  PLOTTING  ##########
        fit = self.fit


        # ##### SAMPLE REGIONS #####
        for span in self.span_artists:
            span.remove()  # Shaded regions of the previous sample
        self.span_artists.clear()
        if not fit.sample_all:  # The whole-spectrum sample is not shaded
            for w1, w2 in fit.regions:
                self.span_artists.append(self.ax.axvspan(w1, w2, color="tab:blue", alpha=0.08, lw=0))


        # ##### SAMPLE POINTS #####
        auto = fit.mode == "auto"  # Sample points are only meaningful for the automatic fit
        none = np.zeros(fit.x_s.size, dtype=bool)  # Mask that hides every point
        good = (~fit.rejected & ~fit.deleted) if auto else none
        if good.sum() > _MAX_MARKERS:  # Too many markers would make the window unresponsive
            good = none
        rejected = fit.rejected if (auto and self.markrej) else none
        deleted  = fit.deleted if auto else none
        self.good_pts.set_data(fit.x_s[good], fit.y_s[good])
        self.rej_pts.set_data(fit.x_s[rejected], fit.y_s[rejected])
        self.del_pts.set_data(fit.x_s[deleted], fit.y_s[deleted])


        # ##### CONTINUUM & OUTPUT #####
        if fit.continuum is not None:
            self.cont_line.set_data(fit.wave, fit.continuum)
            self.out_line.set_ydata(fit.result)
            finite = fit.result[np.isfinite(fit.result)]
            if finite.size:  # Robust y limits: outliers do not squash the plot
                lo, hi = np.percentile(finite, [0.5, 99.5])
                pad    = 0.1 * (hi - lo) if hi > lo else 0.1 * max(abs(hi), 1.)
                self.axn.set_ylim(lo - pad, hi + pad)
        else:
            self.cont_line.set_data([], [])
            self.out_line.set_ydata(np.full(fit.wave.size, np.nan))
        #
        self.ax.set_title(fit._status(), fontsize=9)  # Status line: function, parameters, rms, rejected / deleted points
        self.fig.canvas.draw_idle()


    # ####################  EVENTS  ####################
    def _on_click(
                  self,
                  event : MouseEvent
        ) -> None:
        """
        Add a manual continuum point with a left-click in the upper panel (not while zooming / panning).

        ----------
        Parameters
        ----------
            event : MouseEvent
                Button press event of the figure.

        -------
        Returns
        -------
            None
                A manual point is drawn and stored.
        """

        # ##########  COMPUTATION  ##########
        toolbar = getattr(self.fig.canvas.manager, "toolbar", None)
        if event.button == 1 and event.inaxes is self.ax and getattr(toolbar, "mode", '') == '':  # Left-click in the upper panel, not zooming / panning
            artist, = self.ax.plot(event.xdata, event.ydata, "rs", ms=5, picker=5)  # Red square, pickable within 5 points
            self.manual_artists.append(artist)
            self.fig.canvas.draw_idle()
    #
    def _on_pick(
                 self,
                 event : PickEvent
        ) -> None:
        """
        Remove a manual continuum point with a right-click on it.

        ----------
        Parameters
        ----------
            event : PickEvent
                Pick event of the figure.

        -------
        Returns
        -------
            None
                The manual point is removed.
        """

        # ##########  COMPUTATION  ##########
        if event.mouseevent.button == 3 and event.artist in self.manual_artists:  # Right-click on a manual point
            self.manual_artists.remove(event.artist)
            event.artist.remove()
            self.fig.canvas.draw_idle()
    #
    def _on_key(
                self,
                event : KeyEvent
        ) -> None:
        """
        Run the action of a key (see '_NORMSPEC_HELP').

        ----------
        Parameters
        ----------
            event : KeyEvent
                Key press event of the figure.

        -------
        Returns
        -------
            None
                The action of the key is executed.
        """

        # ##########  MODULES & PACKAGES  ##########
        import matplotlib.pyplot as plt
        import numpy as np


        # ##########  COMPUTATION  ##########
        fit, key = self.fit, event.key


        # ##### MANUAL CONTINUUM #####
        if key == "enter":
            if not self.manual_artists:
                print("No manual continuum points selected.")
                return
            points = np.array([a.get_data() for a in self.manual_artists], dtype=float).reshape(-1, 2)  # One (x, y) row per manual point
            try:
                fit._fit_manual(points[:, 0], points[:, 1])
            except ValueError as e:  # E.g. fewer than two different wavelengths
                print(e)
                return
            self._update_plot()


        # ##### AUTOMATIC FIT #####
        elif key == "f":
            self._refit()
        #
        elif key in ("+", "="):  # '=' is '+' without shift on many keyboards
            fit.order += 1
            self._refit()
        #
        elif key == "-":
            if fit.order > 1:
                fit.order -= 1
                self._refit()
            else:
                print("Order is already 1.")
        #
        elif key in _FUNCTION_KEYS:
            fit.function = _FUNCTION_KEYS[key]
            self._refit()


        # ##### SAMPLE REGIONS #####
        elif key == "s":
            if event.inaxes not in (self.ax, self.axn):
                print("Place the cursor inside the plot to mark a sample region.")
            elif self.sample_start is None:  # First press: start of the region
                self.sample_start = event.xdata
                print(f"Sample region starts at {event.xdata:.3f}: press 's' again at the other end.")
            else:  # Second press: end of the region
                fit._add_sample_region(float(self.sample_start), float(event.xdata))
                self.sample_start = None
                self._refit()
        #
        elif key == "t":
            self.sample_start = None
            fit._sample_whole_spectrum()
            self._refit()


        # ##### DELETE / UNDELETE SAMPLE POINTS #####
        elif key in ("d", "u"):
            self._toggle_nearest(event, delete=(key == "d"))


        # ##### RESET #####
        elif key == "r":
            for artist in self.manual_artists:
                artist.remove()
            self.manual_artists.clear()
            self.sample_start = None
            fit._reset()
            self._refit()


        # ##### SAVE / QUIT / HELP #####
        elif key == "w":
            try:
                fit._save(self.output_file)
            except Exception as e:  # Never raise inside a Matplotlib callback
                print(f"Could not save the output spectrum: {e}")
        #
        elif key == "q":
            plt.close(self.fig)  # Matplotlib's 'q' shortcut is disabled in this window
        #
        elif key == "?":
            print(_NORMSPEC_HELP)
    #
    def _on_close(
                  self,
                  event : CloseEvent
        ) -> None:
        """
        Restore the toolmanager keymaps when the window closes.

        ----------
        Parameters
        ----------
            event : CloseEvent
                Close event of the figure (unused).

        -------
        Returns
        -------
            None
                The keymaps are restored.
        """

        # ##########  COMPUTATION  ##########
        _release_keyboard(self.fig, self._keymaps)
        self._keymaps = {}
    #
    def _toggle_nearest(
                        self,
                        event  : KeyEvent,
                        delete : bool
        ) -> None:
        """
        Delete (or undelete) the sample point nearest to the cursor, in display coordinates.

        ----------
        Parameters
        ----------
            event : KeyEvent
                Key press event ('d' or 'u') with the cursor position.

            delete : bool
                True to delete the nearest kept point, False to undelete the
                nearest deleted point.

        -------
        Returns
        -------
            None
                The deleted mask is updated and the continuum refitted.
        """

        # ##########  MODULES & PACKAGES  ##########
        import numpy as np


        # ##########  COMPUTATION  ##########
        fit = self.fit
        if fit.mode != "auto" or event.inaxes is not self.ax:
            print("Place the cursor in the upper panel (automatic fit) to delete / undelete points.")
            return
        candidates = np.flatnonzero(~fit.deleted if delete else fit.deleted)  # Kept points to delete, or deleted points to restore
        if candidates.size == 0:
            print("No points to delete." if delete else "No deleted points.")
            return
        disp                 = self.ax.transData.transform(np.column_stack([fit.x_s[candidates], fit.y_s[candidates]]))  # Display (pixel) coordinates: wavelength and flux scales are comparable
        nearest              = candidates[np.argmin(np.hypot(disp[:, 0] - event.x, disp[:, 1] - event.y))]  # Closest to the cursor
        fit.deleted[nearest] = delete
        self._refit()
