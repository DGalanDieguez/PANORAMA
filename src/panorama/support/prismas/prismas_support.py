"""
Internal helpers of the 'prismas' module: registry of the supported model
atmosphere grids, validation of the arguments shared by the public functions,
reading of the model catalogs distributed with PANORAMA, decoding of the
stellar parameters from the model, SPAMMS and archive names, selection of the
models ('find_models'), SPAMMS grid check and figure ('needed_models'), model
file search, and parallel renaming of the PRISMAS archives ('rename_tars').

Conventions shared by every function in this module:
    - Functions and constants are private (leading underscore); they are used
      by the public functions of 'panorama.prismas' and are not exposed at
      the package level.
    - Every grid is described by one entry of '_grid_registry' (catalogs,
      model and SPAMMS naming patterns, decoders, criteria and file
      extensions); adding a grid only requires its decoders, a registry entry
      and its catalog files in 'support/prismas/models_grids/'.
    - Model names are given without file extension; the files of a model are
      the model name followed by each extension of its grid.
    - Stellar parameters are Python floats (NumPy arrays for whole catalogs):
      'teff' in K, 'logg' as log10(g / [cm s^-2]), 'met' as [M/H] in dex
      (Kurucz) or Z / Zsun (TLUSTY), 'alpha' and 'carbon' in dex; parameters
      not defined for a grid are NaN. Microturbulences are in km/s.
    - Scalar real arguments are validated with '_check_real'
      ('support/units/units_support.py'), boolean flags with '_check_bool'
      ('support/obs/obs_support.py') and files with '_check_input_file'
      ('support/tex/tex_support.py'): these helpers are shared with the other
      modules instead of being duplicated here.
    - All functions validate their inputs before computing.
"""

from __future__ import annotations

# ==================================================================
# ============================== INIT ==============================
# ==================================================================

# ####################  MODULES & PACKAGES  ####################
"""
- functools is imported at module level because '_km_value' is decorated with
  functools.lru_cache, which must exist when the function is defined (it
  cannot be imported inside a function)
"""
import functools


# ####################  CATALOG CACHE  ####################
"""
- Decoded catalogs, keyed by (grid code, catalog file names) and filled by
  '_read_catalog' on the first request, so that a catalog file (about 850000
  names for 'KM') is read and decoded only once per Python session
"""
_CATALOG_CACHE = {}


# ####################  TOLERANCE  ####################
"""
- Absolute tolerance used to compare stellar parameters (in the units of each
  parameter), so that the fractional TLUSTY metallicities (0.032, 0.020, ...)
  and the values decoded from the names are matched despite rounding
"""
_PARAM_TOL = 1e-6


# ####################  SPAMMS  ####################
"""
- SPAMMS script of each input file, keyed by the ending of the input file name
  without extension (case-insensitive; tried in this order, '' matches any
  other file), and sentence that precedes the list of missing grid points
  in the output of 'spamms.py -c'
"""
_SPAMMS_SCRIPTS = {
                   '_disk'    : "spamms_disk.py",
                   'noreflec' : "spamms_noreflec.py",
                   ''         : "spamms.py",
                  }
_SPAMMS_PHRASE  = "The missing grid points are:"


# ####################  DUMMY GRID ENTRIES  ####################
"""
- Effective temperatures in K of the dummy grid entries of SPAMMS (disk
  patches), which are neither plotted nor matched with a model atmosphere
"""
_DUMMY_TEFF = (0.0, 1000.0)


# ####################  ARCHIVES  ####################
"""
- Size in bytes of the I/O blocks used to stream the PRISMAS archives
  ('_rewrite_tar')
"""
_TAR_BUFFER = 4 * 1024 * 1024  # 4 MB




# ===================================================================
# ============================== GRIDS ==============================
# ===================================================================

# ==============================
def _grid_registry(
    ) -> dict[str, dict]:
    """
    Return the registry of the model atmosphere grids supported by 'find_models', 'needed_models' and 'rename_tars'.

    -------
    Returns
    -------
        registry : dict[str, dict]
            A dictionary keyed by grid code (same code as '--mode' in
            'prismas.py', e.g. 'KC'); each value is a dictionary containing
            the following key-value pairs:
                - 'name' : str
                    Human-readable name of the grid, used in the messages.

                - 'catalogs' : tuple[str, ...]
                    File names of the model catalogs of the grid, stored in
                    'support/prismas/models_grids/'; they are read in this
                    order and merged.

                - 'pattern' : str
                    Regular expression that a model name (without extension)
                    must match (used with re.search); its named groups are
                    read by 'decoder'.

                - 'decoder' : Callable[[re.Match], dict | None]
                    Function that converts the match into the stellar
                    parameters (see '_decode_kc'), or returns None if the
                    name is not a model of the grid.

                - 'params' : tuple[str, ...]
                    Selection criteria defined for the grid (subset of
                    'teff', 'logg', 'met', 'alpha', 'carbon').

                - 'met_label' : str
                    Short label of the metallicity column in the summary
                    tables (e.g. '[M/H]', 'Z/Zsun').

                - 'met_unit' : str
                    Meaning and unit of 'met', used in the error messages.

                - 'met_constraint' : str | None
                    Constraint on 'met' passed to '_check_real'
                    ('non-negative' for Z / Zsun; None for [M/H]).

                - 'extensions' : tuple[str, ...]
                    File extensions of a model, in the order in which its
                    files are listed: the files of a model are its name
                    followed by each extension. The first one is the file
                    passed to '--model' in 'prismas.py'. They are also
                    removed from the catalog names that include them.

                - 'spamms_pattern' : str
                    Regular expression of the SPAMMS grid-entry names of the
                    grid (used with re.search); its named groups are read by
                    'spamms_decoder'.

                - 'spamms_decoder' : Callable[[re.Match], dict]
                    Function that converts the match into the stellar
                    parameters and microturbulence of the entry (see
                    '_spamms_kc').

                - 'subgrid_order' : Callable[[float], tuple[str, ...]] | None
                    Function that gives the subgrids to try for an effective
                    temperature, in order of preference (see
                    '_tlusty_subgrid_order'), or None for grids without
                    subgrids.

                - 'archive_name' : str
                    Template of the name of the PRISMAS archives of modes
                    0 / 1 (without '.tar.gz'), filled by '_archive_name' with
                    the fields 'teff', 'logg', 'met', 'alpha', 'carbon'
                    (signed, 'm' / 'p'), 'z' (unsigned metallicity) and
                    'vmic'.

                - 'subgrid_tags' : dict[str, str] | None
                    A dictionary keyed by subgrid name (e.g. 'OSTAR2002');
                    each value is the tag appended by 'rename_tars' to the
                    name of an archive (e.g. '_OSTAR02') when archives of
                    different subgrids would get the same name. None for
                    grids without subgrids.

                - 'plot_catalogs' : dict[str, str] | None
                    A dictionary keyed by legend label; each value is the
                    file name of a catalog (in 'support/prismas/models_grids/')
                    whose nodes are drawn as one series by
                    '_plot_needed_models'. Only used for the figure. None to
                    draw the nodes of 'catalogs', one series per subgrid.

    -----
    Notes
    -----
        1. ADDING A GRID
            - Write the decoders of the model and SPAMMS names (same
              signatures and keys as '_decode_kc' and '_spamms_kc'), add an
              entry to the dictionary below and store its catalog files (one
              model name per line, without extension) in
              'support/prismas/models_grids/'. 'find_models', 'needed_models'
              and 'rename_tars' need no other change.

        2. SUPPORTED GRIDS
            - 'KC': Kurucz-Castelli ATLAS9 ODFNEW models
              ('m0.50t03500g0.0am05ak2odfnew', file '.dat'; SPAMMS entries
              'T15000_G3.00_Mm0.500_V2', 'T15000_G3.00_Mam0.500_V2' or
              'T15000_G3.00_Mm0.500_ap0.400_V2').

            - 'KM': Kurucz-Mészáros ATLAS9 APOGEE models
              ('amm03cm03om03t10000g20v20', file '.mod'; SPAMMS entries
              'T5000_G2.50_Mm0.500_am0.250_cp0.000_V2').

            - 'T': TLUSTY OSTAR2002 and BSTAR2006 models ('G30000g400v10',
              'BG20000g300v2', files '.7', '.5' and '.nst'; SPAMMS entries
              'T30000_G4.00_Z1.000_V10').

        3. REGISTRY
            - The registry is built at every call (it refers to the decoders,
              defined below in this module); callers that need it in a loop
              read it once before the loop.
    """

    # ####################  OUTPUT  ####################
    registry = {
                'KC' : {
                        'name'           : "Kurucz-Castelli ATLAS9 ODFNEW",
                        'catalogs'       : ("KC_models.dat",),
                        'pattern'        : r"^(?P<sign>[mp])(?P<met>\d+\.\d+)t(?P<teff>\d+)g(?P<logg>\d+\.\d+)a[mp]\d{2}(?P<alpha>a?)k(?P<vturb>\d+)odfnew$",
                        'decoder'        : _decode_kc,
                        'params'         : ('teff', 'logg', 'met', 'alpha'),
                        'met_label'      : "[M/H]",
                        'met_unit'       : "[M/H] in dex",
                        'met_constraint' : None,
                        'extensions'     : (".dat",),
                        'spamms_pattern' : r"^T(?P<teff>\d+)_G(?P<logg>\d+(?:\.\d+)?)_M(?P<flag>a?)(?P<met>[mp]\d+(?:\.\d+)?)(?:_a(?P<alpha>[mp]\d+(?:\.\d+)?))?_V(?P<vmic>\d+(?:\.\d+)?)$",
                        'spamms_decoder' : _spamms_kc,
                        'subgrid_order'  : None,
                        'archive_name'   : "T{teff}_G{logg}_M{met}_a{alpha}_V{vmic}",
                        'subgrid_tags'   : None,
                        'plot_catalogs'  : None,
                       },
                'KM' : {
                        'name'           : "Kurucz-Meszaros ATLAS9 APOGEE",
                        'catalogs'       : ("KM_models.dat",),
                        'pattern'        : r"^am(?P<mh>[mp]\d{2})c(?P<cm>[mp]\d{2})o(?P<am>[mp]\d{2})t(?P<teff>\d+)g(?P<logg>\d+)v(?P<vturb>\d+)$",
                        'decoder'        : _decode_km,
                        'params'         : ('teff', 'logg', 'met', 'alpha', 'carbon'),
                        'met_label'      : "[M/H]",
                        'met_unit'       : "[M/H] in dex",
                        'met_constraint' : None,
                        'extensions'     : (".mod",),
                        'spamms_pattern' : r"^T(?P<teff>\d+)_G(?P<logg>\d+(?:\.\d+)?)_M(?P<met>[mp]\d+(?:\.\d+)?)_a(?P<alpha>[mp]\d+(?:\.\d+)?)_c(?P<carbon>[mp]\d+(?:\.\d+)?)_V(?P<vmic>\d+(?:\.\d+)?)$",
                        'spamms_decoder' : _spamms_km,
                        'subgrid_order'  : None,
                        'archive_name'   : "T{teff}_G{logg}_M{met}_a{alpha}_c{carbon}_V{vmic}",
                        'subgrid_tags'   : None,
                        'plot_catalogs'  : None,
                       },
                'T'  : {
                        'name'           : "TLUSTY OSTAR2002/BSTAR2006",
                        'catalogs'       : ("T_OSTAR02+BSTAR06models.dat",),
                        'pattern'        : r"^(?P<sub>B?)(?P<code>[A-Z])(?P<teff>\d+)g(?P<logg>\d+)v(?P<vturb>\d+)$",
                        'decoder'        : _decode_tlusty,
                        'params'         : ('teff', 'logg', 'met'),
                        'met_label'      : "Z/Zsun",
                        'met_unit'       : "Z / Zsun",
                        'met_constraint' : "non-negative",
                        'extensions'     : (".7", ".5", ".nst"),
                        'spamms_pattern' : r"^T(?P<teff>\d+)_G(?P<logg>\d+(?:\.\d+)?)_Z(?P<met>\d+(?:\.\d+)?)_V(?P<vmic>\d+(?:\.\d+)?)$",
                        'spamms_decoder' : _spamms_tlusty,
                        'subgrid_order'  : _tlusty_subgrid_order,
                        'archive_name'   : "T{teff}_G{logg}_Z{z}_V{vmic}",
                        'subgrid_tags'   : {
                                            'OSTAR2002' : "OSTAR02",
                                            'BSTAR2006' : "BSTAR06",
                                           },
                        'plot_catalogs'  : {
                                            'OSTAR2002' : "T_OSTAR02_models.dat",
                                            'BSTAR2006' : "T_BSTAR06_models.dat",
                                           },
                       },
               }
    return registry




# ========================================================================
# ============================== VALIDATION ==============================
# ========================================================================

# ==============================
def _check_grid(
                grid : object
    ) -> str:
    """
    Validate a model grid code (case-insensitive) against '_grid_registry' and return it in upper case.

    ----------
    Parameters
    ----------
        grid : object
            Model grid code given by the user (same code as '--mode' in
            'prismas.py', e.g. 'KC', 't').

    -------
    Returns
    -------
        grid : str
            The grid code without surrounding spaces and in upper case, a key
            of '_grid_registry'.

    ------
    Raises
    ------
        TypeError
            If 'grid' is not a string.

        ValueError
            If 'grid' is not a supported grid.

    -----
    Notes
    -----
        - Shared by every function of 'prismas' that receives a grid, so that
          the same argument gives the same error messages.
    """

    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(grid, str):
        raise TypeError(f"Model grid 'grid' must be a string (e.g. 'KC', 'KM' or 'T'). Got {type(grid).__name__}.")


    # ####################  COMPUTATION  ####################
    registry  = _grid_registry()
    grid_code = grid.strip().upper()  # Case-insensitive: 't' -> 'T'
    if grid_code not in registry:
        raise ValueError(f"Model grid 'grid' must be one of {', '.join(repr(k) for k in registry)} (case-insensitive). Got '{grid}'.")


    # ####################  OUTPUT  ####################
    return grid_code




# ==============================
def _check_directory(
                     directory  : object,
                     name       : str,
                     label      : str,
                     kind       : str,
                     readable   : bool = True,
                     writable   : bool = False,
                     allow_none : bool = False
    ) -> Path | None:
    """
    Validate a directory argument (path, existence, type and permissions) and convert it to a Path.

    ----------
    Parameters
    ----------
        directory : object
            Directory given by the user (str or pathlib.Path), or None if
            'allow_none'.

        name : str
            Name of the argument in the error messages (e.g. 'modeldir').

        label : str
            Role of the argument in the error messages (e.g. 'Model
            directory').

        kind : str
            Kind of directory in the error messages (e.g. 'model directory').

        readable : bool, optional
            If True, the directory must be readable. Default is True.

        writable : bool, optional
            If True, the directory must be writable. Default is False.

        allow_none : bool, optional
            If True, None is accepted and returned unchanged (optional
            arguments). Default is False.

    -------
    Returns
    -------
        directory : Path | None
            The directory as a Path, or None if it was None and 'allow_none'
            is True.

    ------
    Raises
    ------
        TypeError
            If 'directory' is not a path (nor None when 'allow_none'), or
            'name', 'label', 'kind', 'readable', 'writable' or 'allow_none'
            have the wrong type.

        ValueError
            If 'name', 'label' or 'kind' are empty.

        FileNotFoundError
            If the directory does not exist.

        NotADirectoryError
            If the path exists but is not a directory.

        PermissionError
            If the directory cannot be read ('readable') or written
            ('writable').

    -----
    Notes
    -----
        - Shared by 'find_models' ('modeldir', 'output_dir'), 'needed_models'
          ('spamms_dir'), 'rename_tars' ('tardir', 'output_dir'),
          '_locate_files' and '_run_spamms', so that the same argument gives
          the same error messages. Directories are never created.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import os
    from pathlib import Path


    # ##### PANORAMA #####
    from ..obs.obs_support import _check_bool


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label), ("kind", kind)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")
    #
    readable = _check_bool(readable, "readable", "Read permission flag")
    #
    writable = _check_bool(writable, "writable", "Write permission flag")
    #
    allow_none = _check_bool(allow_none, "allow_none", "None flag")


    # ####################  COMPUTATION  ####################
    if directory is None and allow_none:  # Optional argument not given
        return None
    elif not isinstance(directory, (str, os.PathLike)):
        none_text = " or None" if allow_none else ''
        raise TypeError(f"{label} '{name}' must be a path (str or pathlib.Path) to a directory{none_text}. Got {type(directory).__name__}.")
    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError(f"{label} not found: '{directory}'.")
    elif not directory.is_dir():
        article = "an" if kind.strip()[0].lower() in "aeiou" else "a"  # 'an archive directory', 'a model directory'
        raise NotADirectoryError(f"Expected {article} {kind} but found a file: '{directory}'.")
    elif readable and not os.access(directory, os.R_OK):
        raise PermissionError(f"Permission denied to read {kind}: '{directory}'.")
    elif writable and not os.access(directory, os.W_OK):
        raise PermissionError(f"Permission denied to write to {kind}: '{directory}'.")


    # ####################  OUTPUT  ####################
    return directory




# ==============================
def _check_ncores(
                  ncores : object
    ) -> int | None:
    """
    Validate the number of parallel workers (threads) of 'rename_tars'.

    ----------
    Parameters
    ----------
        ncores : object
            Number of parallel workers (>= 1), or None (self-tuned).

    -------
    Returns
    -------
        ncores : int | None
            The number of workers as a Python int, or None.

    ------
    Raises
    ------
        TypeError
            If 'ncores' is neither an integer nor None.

        ValueError
            If 'ncores' is smaller than 1.

    -----
    Notes
    -----
        - Shared by 'rename_tars' and '_run_parallel', so that the same
          argument gives the same error messages.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if ncores is None:  # Self-tuned
        return None
    elif not isinstance(ncores, numbers.Integral) or isinstance(ncores, bool):
        raise TypeError(f"Number of parallel workers 'ncores' must be an integer or None (self-tuned). Got {type(ncores).__name__}.")
    elif ncores < 1:
        raise ValueError(f"Number of parallel workers 'ncores' must be >= 1. Got {ncores}.")


    # ####################  OUTPUT  ####################
    return int(ncores)  # Python int (also for NumPy integers)




# ======================================================================
# ============================== CATALOGS ==============================
# ======================================================================

# ==============================
def _read_catalog(
                  grid     : str,
                  catalogs : tuple[str, ...] | None = None
    ) -> dict[str, tuple | np.ndarray]:
    """
    Read the model catalogs of a grid and decode the stellar parameters of every model, as one column per parameter.

    ----------
    Parameters
    ----------
        grid : str
            Model grid code, a key of '_grid_registry' (e.g. 'KC', 'KM',
            'T'; case-insensitive).

        catalogs : tuple[str, ...] | None, optional
            File names of the catalogs to read, stored in
            'support/prismas/models_grids/' (e.g. the 'plot_catalogs' of the
            grid). Default is None (the 'catalogs' of the grid in
            '_grid_registry').

    -------
    Returns
    -------
        catalog : dict[str, tuple | np.ndarray]
            Decoded models, in the order of the catalog files; every column
            has one element per model.
            A dictionary containing the following key-value pairs:
                - 'model' : tuple[str, ...]
                    Model names, without file extension.

                - 'subgrid' : tuple[str, ...]
                    'OSTAR2002' or 'BSTAR2006' for 'T'; '' for Kurucz.

                - 'teff' : np.ndarray
                    Effective temperatures in K (read-only float array).

                - 'logg' : np.ndarray
                    Surface gravities, log10(g / [cm s^-2]) (read-only float
                    array).

                - 'met' : np.ndarray
                    [M/H] in dex for 'KC' / 'KM'; Z / Zsun for 'T' (read-only
                    float array).

                - 'alpha' : np.ndarray
                    Alpha enhancements [alpha/M] in dex (read-only float
                    array; NaN for 'T').

                - 'carbon' : np.ndarray
                    Carbon abundances [C/M] in dex (read-only float array;
                    NaN except for 'KM').

    ------
    Raises
    ------
        TypeError
            If 'grid' is not a string, or 'catalogs' is not a tuple of
            strings or None.

        ValueError
            If 'grid' is not a supported grid, 'catalogs' is empty or
            contains an empty name, or a catalog file contains no model
            name.

        FileNotFoundError
            If a catalog file of the grid is not installed.

        IsADirectoryError
            If a catalog path is a directory.

        PermissionError
            If a catalog file cannot be read.

        RuntimeError
            If the catalog files contain no valid model name of the grid.

    -----
    Notes
    -----
        1. READING
            - The names are read with '_read_object_names' (support of
              'tex'): lines starting with '#' and empty lines are ignored,
              and only the first whitespace-separated column is read. A file
              extension of the grid at the end of a name (e.g.
              'm0.50t03500g0.0am05k2odfnew.dat') is removed, so catalogs with
              or without extensions are accepted.

            - Names that do not follow the naming of the grid, and repeated
              names (also across the catalog files of a grid), are skipped
              with a single warning each.

            - The names are decoded with the pattern and decoder of the grid
              in '_grid_registry', as in '_decode_model', but with the
              pattern compiled and the registry read only once (the 'KM'
              catalog has about 850000 names).

        2. CACHE
            - The decoded catalog is kept in '_CATALOG_CACHE' (keyed by grid
              and catalog files), so each catalog file is read once per
              Python session; the arrays are read-only so that the cache
              cannot be modified by mistake. Restart the session after
              changing a catalog file.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numpy as np
    from pathlib import Path
    import re
    import warnings


    # ##### PANORAMA #####
    from ..tex.tex_support import _check_input_file, _read_object_names


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]
    #
    if catalogs is None:
        catalogs = grid_info['catalogs']
    elif not isinstance(catalogs, tuple):
        raise TypeError(f"Catalog file names 'catalogs' must be a tuple of strings or None. Got {type(catalogs).__name__}.")
    elif len(catalogs) == 0:
        raise ValueError("Catalog file names 'catalogs' must contain at least one file name. Got an empty tuple.")
    for file_name in catalogs:
        if not isinstance(file_name, str):
            raise TypeError(f"Catalog file names 'catalogs' must contain only strings. Got {type(file_name).__name__}.")
        elif not file_name.strip():
            raise ValueError(f"Catalog file names 'catalogs' must not contain empty names. Got '{file_name}'.")


    # ##########  CACHE  ##########
    cache_key = (grid, catalogs)
    if cache_key in _CATALOG_CACHE:
        return _CATALOG_CACHE[cache_key]


    # ##########  CATALOG FILES  ##########
    catalog_dir   = Path(__file__).resolve().parent / "models_grids"  # support/prismas/models_grids/
    catalog_files = [_check_input_file(catalog_dir / file_name, "catalogs", f"{grid} model catalog", "catalog file") for file_name in catalogs]


    # ####################  PARAMETERS  ####################
    extensions = grid_info['extensions']
    pattern    = re.compile(grid_info['pattern'])  # Compiled once for the whole catalog
    decoder    = grid_info['decoder']
    columns    = ('subgrid', 'teff', 'logg', 'met', 'alpha', 'carbon')


    # ####################  COMPUTATION  ####################
    # ##########  READING  ##########
    names     = []
    values    = {key: [] for key in columns}
    seen      = set()
    undecoded = []
    repeated  = []
    raw_names = [name for catalog_file in catalog_files for name in _read_object_names(catalog_file)]  # First column; empty lines and comments skipped
    for model_name in raw_names:
        for ext in extensions:
            if model_name.endswith(ext):  # Names are kept without file extension
                model_name = model_name[:-len(ext)]
                break
        if model_name in seen:
            repeated.append(model_name)
            continue
        seen.add(model_name)


        # ##### DECODING #####
        match  = pattern.search(model_name)
        params = None if match is None else decoder(match)  # None: not a model of the grid
        if params is None:
            undecoded.append(model_name)
            continue
        names.append(model_name)
        for key in columns:
            values[key].append(params[key])


    # ##########  CHECKS  ##########
    if undecoded:
        warnings.warn(f"{len(undecoded)} name(s) of the {grid} catalog do not follow the naming of the grid and were skipped (e.g. '{undecoded[0]}').", stacklevel=3)
    if repeated:
        warnings.warn(f"{len(repeated)} repeated name(s) of the {grid} catalog were skipped (e.g. '{repeated[0]}').", stacklevel=3)
    if not names:
        raise RuntimeError(f"{grid} model catalogs contain no valid model name: {', '.join(repr(str(f)) for f in catalog_files)}.")


    # ##########  COLUMNS  ##########
    catalog = {
               'model'   : tuple(names),
               'subgrid' : tuple(values['subgrid']),
              }
    for key in columns[1:]:
        column                 = np.array(values[key], dtype=float)
        column.flags.writeable = False  # Protect the cached catalog
        catalog[key]           = column
    _CATALOG_CACHE[cache_key] = catalog


    # ####################  OUTPUT  ####################
    return catalog




# ======================================================================
# ============================== DECODING ==============================
# ======================================================================

# ==============================
def _decode_model(
                  model_name : str,
                  grid       : str
    ) -> dict | None:
    """
    Decode the stellar parameters of a model atmosphere from its name, using the naming pattern and decoder of its grid.

    ----------
    Parameters
    ----------
        model_name : str
            Model name without file extension (e.g.
            'm0.50t03500g0.0am05ak2odfnew', 'G30000g400v10'); also the name
            of the PRISMAS archives of HPC mode without '.tar.gz'.

        grid : str
            Model grid code, a key of '_grid_registry' (case-insensitive).

    -------
    Returns
    -------
        params : dict | None
            Stellar parameters with the keys of '_decode_kc', or None if the
            name does not follow the naming of the grid.

    ------
    Raises
    ------
        TypeError
            If 'model_name' or 'grid' is not a string.

        ValueError
            If 'grid' is not a supported grid.

    -----
    Notes
    -----
        - Used by 'rename_tars' (archive names). '_read_catalog' applies the
          same pattern and decoder to every name of a catalog, so that models
          and archives are decoded in the same way.
    """

    # ####################  MODULES & PACKAGES  ####################
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(model_name, str):
        raise TypeError(f"Model name 'model_name' must be a string. Got {type(model_name).__name__}.")
    #
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]


    # ####################  COMPUTATION  ####################
    match = re.search(grid_info['pattern'], model_name)  # re caches the compiled patterns


    # ####################  OUTPUT  ####################
    return None if match is None else grid_info['decoder'](match)




# ==============================
def _decode_kc(
               match : re.Match
    ) -> dict:
    """
    Convert the match of a Kurucz-Castelli (ODFNEW) model name into its stellar parameters.

    ----------
    Parameters
    ----------
        match : re.Match
            Match of the 'KC' pattern of '_grid_registry', with the groups
            'sign', 'met', 'teff', 'logg' and 'alpha'.

    -------
    Returns
    -------
        params : dict
            A dictionary containing the following key-value pairs:
                - 'subgrid' : str
                    Subgrid of the model ('OSTAR2002' or 'BSTAR2006' for 'T';
                    '' for Kurucz).

                - 'teff' : float
                    Effective temperature in K.

                - 'logg' : float
                    Surface gravity, log10(g / [cm s^-2]).

                - 'met' : float
                    [M/H] in dex ('KC', 'KM') or Z / Zsun ('T').

                - 'alpha' : float
                    Alpha enhancement [alpha/M] in dex (NaN for 'T').

                - 'carbon' : float
                    Carbon abundance [C/M] in dex (NaN except for 'KM').

    ------
    Raises
    ------
        TypeError
            If 'match' is not a re.Match object.

    -----
    Notes
    -----
        - Names have the form 'm0.50t03500g0.0am05ak2odfnew': '[mp]' and the
          number after it give [M/H] ('m' negative, 'p' positive), 't' Teff in
          K and 'g' log g. The token 'a[mp]dd' repeats [M/H]; an 'a' right
          after it (before 'k2', microturbulence of 2 km/s) marks the
          alpha-enhanced models (AODFNEW, [alpha/M] = +0.4); otherwise
          [alpha/M] = 0.

        - The keys of 'params' are shared by every model-name decoder
          ('_decode_km', '_decode_tlusty').

    ----------
    References
    ----------
        - Castelli, F. & Kurucz, R. L. 2003, IAUS, 210, A20.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(match, re.Match):
        raise TypeError(f"Model-name match 'match' must be a re.Match object. Got {type(match).__name__}.")


    # ####################  OUTPUT  ####################
    params = {
              'subgrid' : '',
              'teff'    : float(match['teff']),
              'logg'    : float(match['logg']),
              'met'     : _signed(match['sign'] + match['met']),
              'alpha'   : 0.4 if match['alpha'] else 0.0,  # AODFNEW: [alpha/M] = +0.4
              'carbon'  : np.nan,
             }
    return params




# ==============================
def _decode_km(
               match : re.Match
    ) -> dict:
    """
    Convert the match of a Kurucz-Mészáros (APOGEE) model name into its stellar parameters.

    ----------
    Parameters
    ----------
        match : re.Match
            Match of the 'KM' pattern of '_grid_registry', with the groups
            'mh', 'cm', 'am', 'teff' and 'logg'.

    -------
    Returns
    -------
        params : dict
            Stellar parameters with the keys of '_decode_kc'.

    ------
    Raises
    ------
        TypeError
            If 'match' is not a re.Match object.

    -----
    Notes
    -----
        - Names have the form 'am<mh>c<cm>o<am>t<teff>g<logg>v<vturb>' (e.g.
          'amm03cm03om03t10000g20v20'): [M/H], [C/M] and [alpha/M] are
          3-character codes decoded by '_km_value', Teff is in K, 'g' is
          10 log g and 'v' is 10 times the microturbulence in km/s.

    ----------
    References
    ----------
        - Mészáros, Sz., Allende Prieto, C., Edvardsson, B. et al. 2012,
          AJ, 144, 120.
    """

    # ####################  MODULES & PACKAGES  ####################
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(match, re.Match):
        raise TypeError(f"Model-name match 'match' must be a re.Match object. Got {type(match).__name__}.")


    # ####################  OUTPUT  ####################
    params = {
              'subgrid' : '',
              'teff'    : float(match['teff']),
              'logg'    : float(match['logg']) / 10.0,  # 'g25' -> log g = 2.5
              'met'     : _km_value(match['mh']),
              'alpha'   : _km_value(match['am']),
              'carbon'  : _km_value(match['cm']),
             }
    return params




# ==============================
def _decode_tlusty(
                   match : re.Match
    ) -> dict | None:
    """
    Convert the match of a TLUSTY (OSTAR2002 / BSTAR2006) model name into its stellar parameters.

    ----------
    Parameters
    ----------
        match : re.Match
            Match of the 'T' pattern of '_grid_registry', with the groups
            'sub', 'code', 'teff' and 'logg'.

    -------
    Returns
    -------
        params : dict | None
            Stellar parameters with the keys of '_decode_kc' ('alpha' and
            'carbon' are NaN), or None if the composition code is unknown.

    ------
    Raises
    ------
        TypeError
            If 'match' is not a re.Match object.

    -----
    Notes
    -----
        - Names have the form '[B]<code><teff>g<100 log g>v<vturb>' (e.g.
          'G30000g400v10', 'BG20000g300v2'): a leading 'B' marks BSTAR2006
          (otherwise OSTAR2002), and the composition code gives Z / Zsun
          (see '_tlusty_metallicity').

    ----------
    References
    ----------
        - Lanz, T. & Hubeny, I. 2003, ApJS, 146, 417.
        - Lanz, T. & Hubeny, I. 2007, ApJS, 169, 83.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(match, re.Match):
        raise TypeError(f"Model-name match 'match' must be a re.Match object. Got {type(match).__name__}.")


    # ####################  COMPUTATION  ####################
    met = _tlusty_metallicity(match['code'])
    if met is None:  # Unknown composition code: not a model of the grid
        return None


    # ####################  OUTPUT  ####################
    params = {
              'subgrid' : "BSTAR2006" if match['sub'] else "OSTAR2002",
              'teff'    : float(match['teff']),
              'logg'    : float(match['logg']) / 100.0,  # 'g400' -> log g = 4.00
              'met'     : met,
              'alpha'   : np.nan,
              'carbon'  : np.nan,
             }
    return params




# ==============================
def _signed(
            text : str
    ) -> float:
    """
    Convert a signed abundance code with an 'm' (minus) or 'p' (plus) prefix into a float.

    ----------
    Parameters
    ----------
        text : str
            Code such as 'm0.50' or 'p0.20' (Kurucz-Castelli [M/H], SPAMMS
            and PRISMAS abundances).

    -------
    Returns
    -------
        value : float
            The signed value in dex (e.g. -0.5, +0.2).

    ------
    Raises
    ------
        TypeError
            If 'text' is not a string.

        ValueError
            If 'text' is not an 'm' / 'p' prefix followed by a number.
    """

    # ####################  MODULES & PACKAGES  ####################
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(text, str):
        raise TypeError(f"Signed abundance code 'text' must be a string. Got {type(text).__name__}.")
    elif re.fullmatch(r"[mp]\d+(\.\d+)?", text) is None:
        raise ValueError(f"Signed abundance code 'text' must be 'm' or 'p' followed by a number (e.g. 'm0.50'). Got '{text}'.")


    # ####################  OUTPUT  ####################
    return -float(text[1:]) if text[0] == "m" else float(text[1:])




# ==============================
@functools.lru_cache(maxsize=None)
def _km_value(
              code : str
    ) -> float:
    """
    Convert a 3-character Kurucz-Mészáros abundance code (e.g. 'm05', 'p03') into a value in dex.

    ----------
    Parameters
    ----------
        code : str
            Sign ('m' minus, 'p' plus) followed by two digits giving 10 times
            the absolute value, rounded up for the quarter steps (e.g. 'm05'
            -> -0.5, 'p03' -> +0.25, 'm08' -> -0.75).

    -------
    Returns
    -------
        value : float
            The abundance in dex.

    ------
    Raises
    ------
        TypeError
            If 'code' is not a string.

        ValueError
            If 'code' is not 'm' or 'p' followed by two digits.

    -----
    Notes
    -----
        - Codes ending in 3 or 8 are quarter steps (x.25, x.75), since the
          grid steps are multiples of 0.25 dex. With this convention the 'KM'
          catalog gives [M/H] from -5.0 to +1.5 dex and [C/M] and [alpha/M]
          from -1.5 to +1.0 dex in steps of 0.25 dex (0.5 dex below
          [M/H] = -3.0), as expected for the grid.

        - The function is called three times per 'KM' model (about 850000),
          but there are only a few tens of different codes: the results are
          cached ('functools.lru_cache', defined in the 'INIT' section) and
          the format is checked without regular expressions. Invalid codes
          raise and are never cached.
    """

    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(code, str):
        raise TypeError(f"Kurucz-Meszaros abundance code 'code' must be a string. Got {type(code).__name__}.")
    elif len(code) != 3 or code[0] not in ("m", "p") or not code[1:].isdigit():
        raise ValueError(f"Kurucz-Meszaros abundance code 'code' must be 'm' or 'p' followed by two digits (e.g. 'm05'). Got '{code}'.")


    # ####################  COMPUTATION  ####################
    tenths = int(code[1:])
    value  = (tenths - 0.5) / 10.0 if tenths % 5 == 3 else tenths / 10.0  # 03 -> 0.25, 08 -> 0.75, 13 -> 1.25


    # ####################  OUTPUT  ####################
    return -value if code[0] == "m" else value




# ==============================
def _tlusty_metallicity(
                        code : str
    ) -> float | None:
    """
    Return the metallicity Z / Zsun of a TLUSTY composition code.

    ----------
    Parameters
    ----------
        code : str
            One-letter composition code of the model name (e.g. 'G').

    -------
    Returns
    -------
        met : float | None
            Z / Zsun of the code, or None if the code is unknown.

    ------
    Raises
    ------
        TypeError
            If 'code' is not a string.

    -----
    Notes
    -----
        - Codes of OSTAR2002 (Lanz & Hubeny 2003): C = 2, G = 1, L = 1/2,
          S = 1/5, T = 1/10, V = 1/30, W = 1/50, X = 1/100, Y = 1/1000 and
          Z = 0 (metal-free). The values are those of 'prismas.py', which
          writes V as 0.032 (not 0.033) in the output names; they are used so
          that the SPAMMS grid entries ('..._Z0.032_...') match.

        - BSTAR2006 (Lanz & Hubeny 2007) uses the same letters after a 'B'
          (BC, BG, BL, BS, BT, BV, BW, BX and BZ in the catalogs distributed
          with PANORAMA), with the same values.

    ----------
    References
    ----------
        - Lanz, T. & Hubeny, I. 2003, ApJS, 146, 417.
        - Lanz, T. & Hubeny, I. 2007, ApJS, 169, 83.
    """

    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(code, str):
        raise TypeError(f"TLUSTY composition code 'code' must be a string. Got {type(code).__name__}.")


    # ####################  PARAMETERS  ####################
    met_codes = {
                 'C' : 2.0,
                 'G' : 1.0,
                 'L' : 1 / 2,
                 'S' : 1 / 5,
                 'T' : 1 / 10,
                 'V' : 0.032,  # Nominal 1/30, written 0.032 by prismas.py
                 'W' : 1 / 50,
                 'X' : 1 / 100,
                 'Y' : 1 / 1000,
                 'Z' : 0.0,
                }


    # ####################  OUTPUT  ####################
    return met_codes.get(code)




# ==============================
def _spamms_kc(
               match : re.Match
    ) -> dict:
    """
    Convert the match of a Kurucz-Castelli SPAMMS grid-entry name into its stellar parameters and microturbulence.

    ----------
    Parameters
    ----------
        match : re.Match
            Match of the 'KC' SPAMMS pattern of '_grid_registry', with the
            groups 'teff', 'logg', 'flag', 'met', 'alpha' (optional) and
            'vmic'.

    -------
    Returns
    -------
        params : dict
            A dictionary containing the following key-value pairs:
                - 'teff' : float
                    Effective temperature in K.

                - 'logg' : float
                    Surface gravity, log10(g / [cm s^-2]).

                - 'met' : float
                    [M/H] in dex ('KC', 'KM') or Z / Zsun ('T').

                - 'alpha' : float
                    Alpha enhancement [alpha/M] in dex (NaN for 'T').

                - 'carbon' : float
                    Carbon abundance [C/M] in dex (NaN except for 'KM').

                - 'vmic' : float
                    Microturbulence of the synthesis in km/s.

    ------
    Raises
    ------
        TypeError
            If 'match' is not a re.Match object.

    -----
    Notes
    -----
        - Accepts the SPAMMS naming 'T15000_G3.00_M[a]m0.500_V2', where 'a'
          marks the alpha-enhanced models ([alpha/M] = +0.4; 0 otherwise),
          and the PRISMAS naming 'T15000_G3.00_Mm0.500_ap0.400_V2'.

        - The keys of 'params' are shared by every SPAMMS-name decoder
          ('_spamms_km', '_spamms_tlusty').
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(match, re.Match):
        raise TypeError(f"SPAMMS-name match 'match' must be a re.Match object. Got {type(match).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: PRISMAS NAMING  ##########
    if match['alpha']:  # E.g. '_ap0.400'
        alpha = _signed(match['alpha'])


    # ##########  OPTION 2: SPAMMS NAMING  ##########
    else:  # 'Ma...' = alpha-enhanced
        alpha = 0.4 if match['flag'] else 0.0


    # ####################  OUTPUT  ####################
    params = {
              'teff'   : float(match['teff']),
              'logg'   : float(match['logg']),
              'met'    : _signed(match['met']),
              'alpha'  : alpha,
              'carbon' : np.nan,
              'vmic'   : float(match['vmic']),
             }
    return params




# ==============================
def _spamms_km(
               match : re.Match
    ) -> dict:
    """
    Convert the match of a Kurucz-Mészáros SPAMMS grid-entry name into its stellar parameters and microturbulence.

    ----------
    Parameters
    ----------
        match : re.Match
            Match of the 'KM' SPAMMS pattern of '_grid_registry', with the
            groups 'teff', 'logg', 'met', 'alpha', 'carbon' and 'vmic' (e.g.
            'T5000_G2.50_Mm0.500_am0.250_cp0.000_V2').

    -------
    Returns
    -------
        params : dict
            Stellar parameters and microturbulence with the keys of
            '_spamms_kc'.

    ------
    Raises
    ------
        TypeError
            If 'match' is not a re.Match object.
    """

    # ####################  MODULES & PACKAGES  ####################
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(match, re.Match):
        raise TypeError(f"SPAMMS-name match 'match' must be a re.Match object. Got {type(match).__name__}.")


    # ####################  OUTPUT  ####################
    params = {
              'teff'   : float(match['teff']),
              'logg'   : float(match['logg']),
              'met'    : _signed(match['met']),
              'alpha'  : _signed(match['alpha']),
              'carbon' : _signed(match['carbon']),
              'vmic'   : float(match['vmic']),
             }
    return params




# ==============================
def _spamms_tlusty(
                   match : re.Match
    ) -> dict:
    """
    Convert the match of a TLUSTY SPAMMS grid-entry name into its stellar parameters and microturbulence.

    ----------
    Parameters
    ----------
        match : re.Match
            Match of the 'T' SPAMMS pattern of '_grid_registry', with the
            groups 'teff', 'logg', 'met' (Z / Zsun) and 'vmic' (e.g.
            'T30000_G4.00_Z1.000_V10').

    -------
    Returns
    -------
        params : dict
            Stellar parameters and microturbulence with the keys of
            '_spamms_kc' ('alpha' and 'carbon' are NaN).

    ------
    Raises
    ------
        TypeError
            If 'match' is not a re.Match object.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(match, re.Match):
        raise TypeError(f"SPAMMS-name match 'match' must be a re.Match object. Got {type(match).__name__}.")


    # ####################  OUTPUT  ####################
    params = {
              'teff'   : float(match['teff']),
              'logg'   : float(match['logg']),
              'met'    : float(match['met']),
              'alpha'  : np.nan,
              'carbon' : np.nan,
              'vmic'   : float(match['vmic']),
             }
    return params




# ==============================
def _decode_spamms(
                   name : str,
                   grid : str
    ) -> dict | None:
    """
    Decode a SPAMMS grid-entry name into its stellar parameters and microturbulence, using the SPAMMS pattern and decoder of its grid.

    ----------
    Parameters
    ----------
        name : str
            SPAMMS grid entry (e.g. 'T30000_G4.00_Z1.000_V10'); also the name
            of the PRISMAS archives of modes 0 / 1 without '.tar.gz'.

        grid : str
            Model grid code, a key of '_grid_registry' (case-insensitive).

    -------
    Returns
    -------
        params : dict | None
            Stellar parameters and microturbulence with the keys of
            '_spamms_kc', or None if the name does not follow the SPAMMS
            naming of the grid.

    ------
    Raises
    ------
        TypeError
            If 'name' or 'grid' is not a string.

        ValueError
            If 'grid' is not a supported grid.

    -----
    Notes
    -----
        - Used by 'needed_models' (SPAMMS grid entries) and 'rename_tars'
          (archives already named as in modes 0 / 1, which follow the same
          format).
    """

    # ####################  MODULES & PACKAGES  ####################
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(name, str):
        raise TypeError(f"SPAMMS grid entry 'name' must be a string. Got {type(name).__name__}.")
    #
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]


    # ####################  COMPUTATION  ####################
    match = re.search(grid_info['spamms_pattern'], name.strip())  # re caches the compiled patterns


    # ####################  OUTPUT  ####################
    return None if match is None else grid_info['spamms_decoder'](match)




# ==============================
def _tlusty_subgrid_order(
                          teff : float
    ) -> tuple[str, str]:
    """
    Return the TLUSTY subgrids to try for an effective temperature, in order of preference.

    ----------
    Parameters
    ----------
        teff : float
            Effective temperature in K (non-negative).

    -------
    Returns
    -------
        subgrids : tuple[str, str]
            ('OSTAR2002', 'BSTAR2006') for Teff >= 27500 K, except 28000 K,
            and ('BSTAR2006', 'OSTAR2002') otherwise (also at 28000 K).

    ------
    Raises
    ------
        TypeError
            If 'teff' is not a real number.

        ValueError
            If 'teff' is not finite or is negative.

    -----
    Notes
    -----
        - OSTAR2002 is sampled every 2500 K from 27500 to 55000 K and
          BSTAR2006 every 1000 K from 15000 to 30000 K. In their overlap
          (27500-30000 K) the OSTAR2002 models are recommended, so OSTAR2002
          is preferred from 27500 K.

        - Exception at 28000 K: a Teff between 27000 and 27500 K (e.g.
          27400 K) is interpolated between the BSTAR2006 models at
          Ti = 27000 K and Ts = 28000 K, so that models of different subgrids
          are never mixed in the interpolation. The BSTAR2006 model is
          therefore used at 28000 K (OSTAR2002 has no model there).

        - The selection catalog of 'T' ('T_OSTAR02+BSTAR06models.dat') lists
          BSTAR2006 from 15000 to 28000 K and OSTAR2002 from 27500 to
          55000 K. The other BSTAR2006 models of the overlap (29000 and
          30000 K) are never needed: from 27500 K, SPAMMS interpolates
          between the OSTAR2002 nodes (e.g. 28500 K between 27500 and
          30000 K), so it never requests 29000 K. The second subgrid is a
          fallback when the model only exists there.

        - The temperatures are compared exactly, since Teff is decoded from
          integer names.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    teff = _check_real(teff, "teff", "Effective temperature", unit="in K", constraint="non-negative")


    # ####################  OUTPUT  ####################
    if teff >= 27500 and teff != 28000:  # 28000 K: upper BSTAR2006 node for 27000 K < Teff < 27500 K
        return ("OSTAR2002", "BSTAR2006")
    return ("BSTAR2006", "OSTAR2002")




# =======================================================================
# ============================== SELECTION ==============================
# =======================================================================

# ==============================
def _check_criterion(
                     value      : object,
                     name       : str,
                     label      : str,
                     unit       : str | None = None,
                     constraint : str | None = None
    ) -> float | tuple[float | None, float | None] | list[float] | None:
    """
    Validate a selection criterion of 'find_models' and convert its values to floats.

    ----------
    Parameters
    ----------
        value : object
            Criterion given by the user: None, a number, a (min, max) tuple
            (either limit can be None) or a non-empty list of numbers.

        name : str
            Name of the argument in the error messages (e.g. 'teff').

        label : str
            Physical quantity of the argument in the error messages (e.g.
            'Effective temperature').

        unit : str | None, optional
            Unit or meaning shown in the type errors (e.g. 'in K'). Default
            is None (nothing shown).

        constraint : str | None, optional
            Constraint on every value, as in '_check_real' ('positive' or
            'non-negative'). Default is None (any finite value).

    -------
    Returns
    -------
        criterion : float | tuple[float | None, float | None] | list[float] | None
            The validated criterion, with the same structure as 'value'.

    ------
    Raises
    ------
        TypeError
            If the criterion or any of its values has an unsupported type,
            or 'name', 'label', 'unit' or 'constraint' have the wrong type.

        ValueError
            If any value is not finite or violates 'constraint', a range
            does not have two elements, has no limit or has min > max, a list
            is empty, or 'name', 'label' or 'constraint' are not valid.

    -----
    Notes
    -----
        - Every value is checked with '_check_real', so the error messages
          are the same as in the rest of the package; the position of the
          value is appended to the argument name (e.g. 'teff[0]' for the
          lower limit of a range, 'teff[2]' for the third element of a list).
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numbers


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")
    #
    if unit is not None and not isinstance(unit, str):
        raise TypeError(f"Unit text 'unit' must be a string or None. Got {type(unit).__name__}.")
    #
    if constraint is not None:
        if not isinstance(constraint, str):
            raise TypeError(f"Constraint 'constraint' must be a string or None. Got {type(constraint).__name__}.")
        elif constraint not in ("positive", "non-negative"):
            raise ValueError(f"Constraint 'constraint' must be 'positive', 'non-negative' or None. Got '{constraint}'.")
    #
    is_number = isinstance(value, numbers.Real) and not isinstance(value, bool)
    if value is not None and not is_number and not isinstance(value, (tuple, list)):
        raise TypeError(f"{label} '{name}' must be None, a number, a (min, max) tuple or a list of numbers. Got {type(value).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: NO FILTER  ##########
    if value is None:
        criterion = None


    # ##########  OPTION 2: SINGLE VALUE  ##########
    elif is_number:
        criterion = _check_real(value, name, label, unit=unit, constraint=constraint)


    # ##########  OPTION 3: RANGE  ##########
    elif isinstance(value, tuple):
        if len(value) != 2:
            raise ValueError(f"{label} range '{name}' must be a (min, max) tuple. Got {len(value)} element(s).")
        elif value[0] is None and value[1] is None:
            raise ValueError(f"{label} range '{name}' must have at least one limit (use None for no filter). Got (None, None).")
        criterion = tuple(_check_real(lim, f"{name}[{i}]", label, unit=unit, constraint=constraint, allow_none=True) for i, lim in enumerate(value))
        if None not in criterion and criterion[0] > criterion[1]:
            raise ValueError(f"{label} range '{name}' must satisfy min <= max. Got ({criterion[0]}, {criterion[1]}).")


    # ##########  OPTION 4: DISCRETE VALUES  ##########
    else:
        if len(value) == 0:
            raise ValueError(f"{label} list '{name}' must contain at least one value. Got an empty list.")
        criterion = [_check_real(v, f"{name}[{i}]", label, unit=unit, constraint=constraint) for i, v in enumerate(value)]


    # ####################  OUTPUT  ####################
    return criterion




# ==============================
def _passes(
            values : np.ndarray,
            crit   : float | tuple[float | None, float | None] | list[float] | None,
            tol    : float = _PARAM_TOL
    ) -> np.ndarray:
    """
    Check which values of a stellar parameter satisfy a selection criterion.

    ----------
    Parameters
    ----------
        values : np.ndarray
            1-D float array with the parameter of every model (NaN if not
            defined for the grid).

        crit : float | tuple[float | None, float | None] | list[float] | None
            Criterion validated by '_check_criterion': a single value, a
            (min, max) range (None = open limit), a list of values, or None
            (no filter).

        tol : float, optional
            Absolute tolerance of the comparisons, in the units of 'values'
            (non-negative). Default is _PARAM_TOL (1e-6).

    -------
    Returns
    -------
        mask : np.ndarray
            1-D boolean array, True where 'values' satisfies 'crit' (all
            True if 'crit' is None; False where 'values' is NaN otherwise).

    ------
    Raises
    ------
        TypeError
            If 'values' is not a 1-D NumPy array of real numbers, or 'crit'
            or 'tol' has the wrong type.

        ValueError
            If 'values' is not 1-D, or 'tol' is not finite or is negative.

    -----
    Notes
    -----
        - Comparisons with NaN are False, so a parameter not defined for the
          grid never satisfies a criterion other than None.

        - 'values' is checked with '_check_array' (support of 'spec'), and
          'tol' with '_check_real' (support of 'units').
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numbers
    import numpy as np


    # ##### PANORAMA #####
    from ..spec.spec_support import _check_array
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(values, "values", "Model parameter values")
    #
    if crit is not None and (not isinstance(crit, (numbers.Real, tuple, list)) or isinstance(crit, bool)):
        raise TypeError(f"Selection criterion 'crit' must be None, a number, a (min, max) tuple or a list of numbers. Got {type(crit).__name__}.")
    #
    tol = _check_real(tol, "tol", "Comparison tolerance", constraint="non-negative")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: NO FILTER  ##########
    if crit is None:
        mask = np.ones(values.shape, dtype=bool)


    # ##########  OPTION 2: RANGE  ##########
    elif isinstance(crit, tuple):
        low, high = crit
        mask      = ~np.isnan(values)
        if low is not None:
            mask &= values >= low - tol
        if high is not None:
            mask &= values <= high + tol


    # ##########  OPTION 3: DISCRETE VALUES  ##########
    elif isinstance(crit, list):
        mask = np.zeros(values.shape, dtype=bool)
        for c in crit:
            mask |= np.abs(values - c) <= tol


    # ##########  OPTION 4: SINGLE VALUE  ##########
    else:
        mask = np.abs(values - crit) <= tol


    # ####################  OUTPUT  ####################
    return mask




# ==============================
def _match_key(
               params : dict,
               keys   : tuple[str, ...]
    ) -> tuple[float, ...]:
    """
    Return the matching key of a grid entry or model: its stellar parameters defined for the grid, rounded to 3 decimals.

    ----------
    Parameters
    ----------
        params : dict
            Grid entry or model, with a real value for each name of 'keys'.

        keys : tuple[str, ...]
            Stellar parameters defined for the grid, in order (the 'params'
            of the grid in '_grid_registry', e.g. ('teff', 'logg', 'met')).

    -------
    Returns
    -------
        key : tuple[float, ...]
            The values of 'keys' in 'params', rounded to 3 decimals (the
            precision of the SPAMMS names).

    ------
    Raises
    ------
        TypeError
            If 'params' is not a dictionary, 'keys' is not a tuple of
            strings, or a value is not a real number.

        ValueError
            If a name of 'keys' is not in 'params'.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(params, dict):
        raise TypeError(f"Stellar parameters 'params' must be a dictionary. Got {type(params).__name__}.")
    #
    if not isinstance(keys, tuple) or not all(isinstance(k, str) for k in keys):
        raise TypeError(f"Parameter names 'keys' must be a tuple of strings. Got {type(keys).__name__}.")
    for k in keys:
        if k not in params:
            raise ValueError(f"Stellar parameters 'params' must contain the key '{k}'. Got keys {list(params)}.")
        elif not isinstance(params[k], numbers.Real) or isinstance(params[k], bool):
            raise TypeError(f"Stellar parameter '{k}' of 'params' must be a real number. Got {type(params[k]).__name__}.")


    # ####################  OUTPUT  ####################
    return tuple(round(float(params[k]), 3) for k in keys)




# ==============================
def _sort_key(
              entry : dict,
              keys  : tuple[str, ...] = ('teff', 'logg', 'met', 'alpha', 'carbon')
    ) -> tuple[float, ...]:
    """
    Return the sorting key of a model or archive entry, placing missing or NaN values last.

    ----------
    Parameters
    ----------
        entry : dict
            Model or archive entry (missing, None or NaN values allowed).

        keys : tuple[str, ...], optional
            Keys of the sorting key, in order. Default is ('teff', 'logg',
            'met', 'alpha', 'carbon'), the order of 'find_models'.

    -------
    Returns
    -------
        key : tuple[float, ...]
            The values of 'keys', with missing, None or NaN values replaced
            by +inf.

    ------
    Raises
    ------
        TypeError
            If 'entry' is not a dictionary, 'keys' is not a tuple of strings,
            or a value is neither a real number nor None.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(entry, dict):
        raise TypeError(f"Model entry 'entry' must be a dictionary. Got {type(entry).__name__}.")
    #
    if not isinstance(keys, tuple) or not all(isinstance(k, str) for k in keys):
        raise TypeError(f"Sorting keys 'keys' must be a tuple of strings. Got {type(keys).__name__}.")
    for k in keys:
        x = entry.get(k)
        if x is not None and (not isinstance(x, numbers.Real) or isinstance(x, bool)):
            raise TypeError(f"Value '{k}' of the model entry 'entry' must be a real number or None. Got {type(x).__name__}.")


    # ####################  COMPUTATION  ####################
    key = []
    for k in keys:
        x = entry.get(k)
        key.append(float("inf") if x is None or x != x else float(x))  # x != x only for NaN


    # ####################  OUTPUT  ####################
    return tuple(key)




# ====================================================================
# ============================== SPAMMS ==============================
# ====================================================================

# ==============================
def _run_spamms(
                input_file  : str | Path,
                script_file : str | Path,
                spamms_dir  : str | Path
    ) -> dict:
    """
    Run the grid check of SPAMMS ('-c') for one input file and read the missing grid points, the warnings and the errors from its output.

    ----------
    Parameters
    ----------
        input_file : str | Path
            Path and filename of the SPAMMS input file.

        script_file : str | Path
            Path and filename of the SPAMMS script ('spamms.py' or one of the
            modified versions of '_SPAMMS_SCRIPTS').

        spamms_dir : str | Path
            Directory where SPAMMS is run (its own directory).

    -------
    Returns
    -------
        result : dict
            A dictionary containing the following key-value pairs:
                - 'entries' : list[str] | None
                    Missing grid points printed by SPAMMS ([] if none is
                    missing; None if the run failed or the list could not be
                    read).

                - 'error' : str
                    Reason of the failure ('' if the run succeeded).

                - 'warnings' : list[str]
                    Warnings printed by SPAMMS, as 'Category: message',
                    without duplicates and without the paths of the Python
                    files that issued them ([] if none).

    ------
    Raises
    ------
        TypeError
            If any argument is not a path.

        FileNotFoundError
            If 'input_file', 'script_file' or 'spamms_dir' does not exist.

        IsADirectoryError
            If 'input_file' or 'script_file' is a directory.

        NotADirectoryError
            If 'spamms_dir' is not a directory.

        PermissionError
            If 'input_file', 'script_file' or 'spamms_dir' cannot be read.

    -----
    Notes
    -----
        1. RUN
            - SPAMMS is run with the Python interpreter of the current
              session as 'python <script> -i <input> -c -n <n_cores>' from
              'spamms_dir', with every path absolute so that relative paths
              given by the user are not resolved twice. 'n_cores' is the
              number of cores allocated to the process ('_allocated_cores'):
              every core of the computer, or only those assigned to the job
              by HTCondor / Slurm.

            - If the input file is not in 'spamms_dir', it is copied there
              for the run and always removed afterwards: also if the copy is
              partial or fails, SPAMMS fails, or the run raises or is
              interrupted (e.g. Ctrl+C). Only a killed Python process (e.g.
              SIGKILL) can leave it. A copy that cannot be removed issues a
              warning, and a copy that cannot be made is reported in 'error'.

            - An input file already in 'spamms_dir' (the same file, or one
              with identical content) is used and never removed. A different
              file (or a directory) with the same name already there is never
              overwritten: the run is skipped and reported in 'error'.

        2. OUTPUT
            - The missing grid points are the Python list printed after
              '_SPAMMS_PHRASE', read with ast.literal_eval (no code is
              executed).

            - A run with a non-zero exit status and no list is a failure, not
              "nothing missing"; 'error' gives the exit status and the last
              line of the standard error (usually the exception), or of the
              standard output if the standard error is empty.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import ast
    import filecmp
    import re
    import shutil
    import subprocess
    import sys
    import warnings


    # ##### PANORAMA #####
    from ..tex.tex_support import _check_input_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    input_file = _check_input_file(input_file, "input_file", "SPAMMS input file", "SPAMMS input file").resolve()
    #
    script_file = _check_input_file(script_file, "script_file", "SPAMMS script", "Python script").resolve()
    #
    spamms_dir = _check_directory(spamms_dir, "spamms_dir", "SPAMMS directory", "SPAMMS directory").resolve()


    # ##########  DICTIONARY  ##########
    result = {
              'entries'  : None,
              'error'    : '',
              'warnings' : [],
             }


    # ####################  PARAMETERS  ####################
    warning_pattern = re.compile(r"^(?:.*?:\d+:\s*)?(?P<category>[A-Za-z]*Warning)\s*:\s*(?P<message>.+)$")  # '/path/file.py:12: UserWarning: text' or 'Warning: text'


    # ####################  COMPUTATION  ####################
    # ##########  INPUT FILE IN THE SPAMMS DIRECTORY  ##########
    run_input = spamms_dir / input_file.name
    copied    = False  # True once a temporary copy may exist in 'spamms_dir' (also a partial one)
    if run_input.exists():
        if not run_input.is_file():
            result['error'] = f"a directory with the same name as the input file exists in the SPAMMS directory ('{input_file.name}')"
            return result
        elif not (run_input.samefile(input_file) or filecmp.cmp(run_input, input_file, shallow=False)):
            result['error'] = f"a different file with the same name already exists in the SPAMMS directory ('{input_file.name}')"
            return result


    # ##########  RUN  ##########
    # The copy is made inside the 'try' block, so that it is removed whatever happens after it starts: failed or partial copy,
    # SPAMMS error, exception or interruption (e.g. Ctrl+C) during the run. An input file already in 'spamms_dir' is never removed
    try:
        if not run_input.exists():
            copied = True  # Set before copying: a partial copy is also removed
            try:
                shutil.copy2(input_file, run_input)
            except OSError as e:  # E.g. no write permission in 'spamms_dir' or disk full: this input file fails, the others go on
                result['error'] = f"the input file could not be copied to the SPAMMS directory ({e})"
                return result
        process = subprocess.run(
                                 [sys.executable, str(script_file), "-i", str(run_input), "-c", "-n", str(_allocated_cores())],  # Every allocated core
                                 capture_output=True,
                                 text=True,
                                 cwd=spamms_dir  # SPAMMS must run from its own directory
                                )
    finally:
        if copied:
            try:
                run_input.unlink(missing_ok=True)  # Temporary copy removed in every case
            except OSError as e:  # Never hide the original error: report the leftover copy instead
                warnings.warn(f"Temporary copy of the SPAMMS input file could not be removed: '{run_input}' ({e}).", stacklevel=3)
    output = process.stdout + "\n" + process.stderr


    # ##########  WARNINGS  ##########
    for line in output.splitlines():
        match = warning_pattern.match(line.strip())
        if match:
            text = f"{match['category']}: {match['message'].strip()}"
            if text not in result['warnings']:
                result['warnings'].append(text)


    # ##########  MISSING GRID POINTS  ##########
    # ##### OPTION 1: LIST PRINTED #####
    if _SPAMMS_PHRASE in output:
        raw   = output.split(_SPAMMS_PHRASE, 1)[1]
        start = raw.find("[")
        end   = raw.find("]", start + 1) if start >= 0 else -1
        if start < 0 or end < 0:
            result['error'] = f"no Python list found after '{_SPAMMS_PHRASE}' in the SPAMMS output"
        else:
            try:
                entries = ast.literal_eval(raw[start:end + 1])  # Safe parsing of the printed Python list
            except (ValueError, SyntaxError) as e:
                result['error'] = f"the list of missing grid points could not be read ({e})"
            else:
                if isinstance(entries, list) and all(isinstance(m, str) for m in entries):
                    result['entries'] = entries
                else:
                    result['error'] = "the list of missing grid points is not a list of names"


    # ##### OPTION 2: NOTHING MISSING #####
    elif process.returncode == 0:
        result['entries'] = []


    # ##### OPTION 3: FAILED RUN #####
    else:
        stderr_lines    = [line.strip() for line in process.stderr.splitlines() if line.strip()]
        stdout_lines    = [line.strip() for line in process.stdout.splitlines() if line.strip()]
        last_line       = (stderr_lines or stdout_lines or ["no error message"])[-1]  # Usually the exception; SPAMMS may also report errors on the standard output
        result['error'] = f"SPAMMS exited with status {process.returncode} ({last_line})"


    # ####################  OUTPUT  ####################
    return result




# =====================================================================
# ============================== FIGURES ==============================
# =====================================================================

# ==============================
def _plot_needed_models(
                        teff    : np.ndarray,
                        logg    : np.ndarray,
                        in_grid : np.ndarray,
                        grid    : str
    ) -> Figure:
    """
    Plot the grid entries needed by SPAMMS on top of the nodes of the model atmosphere grid (Kiel diagram).

    ----------
    Parameters
    ----------
        teff : np.ndarray
            1-D float array with the effective temperatures of the needed
            entries, in K.

        logg : np.ndarray
            1-D float array with the surface gravities of the needed entries,
            log10(g / [cm s^-2]); same size as 'teff'.

        in_grid : np.ndarray
            1-D boolean array, True for the needed entries whose model is in
            the catalogs of the grid; same size as 'teff'.

        grid : str
            Model grid code, a key of '_grid_registry' (case-insensitive).

    -------
    Returns
    -------
        fig : Figure
            The figure (not shown, saved nor closed).

    ------
    Raises
    ------
        TypeError
            If 'teff' or 'logg' is not a 1-D float NumPy array, 'in_grid' is
            not a 1-D boolean NumPy array, or 'grid' is not a string.

        ValueError
            If 'teff', 'logg' and 'in_grid' have different sizes, or 'grid'
            is not a supported grid.

    -----
    Notes
    -----
        1. GRID NODES
            - The grid nodes are the distinct (Teff, log g) pairs of the
              'plot_catalogs' of the grid in '_grid_registry', one series per
              catalog ('T_OSTAR02_models.dat' and 'T_BSTAR06_models.dat' for
              'T', so that the full OSTAR2002 and BSTAR2006 grids are drawn,
              also in their overlap), or of its 'catalogs', one series per
              subgrid, if 'plot_catalogs' is None.

            - The series are drawn with decreasing marker sizes, so that
              overlapping subgrids remain visible.

        2. NEEDED ENTRIES
            - Black crosses if their model is in the grid, and red crosses
              surrounded by a red circle otherwise (outside the grid), drawn
              larger than the grid nodes and on top of them.

            - The legend is placed above the axes, so that it never hides
              grid nodes or needed entries.

        3. AXES
            - The axes are reversed (hot / low gravity to the left / top)
              and span the grid nodes and the needed entries.

            - The figure is two A&A columns wide (17 cm = 6.69 in) and
              4.8 in high.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from matplotlib.legend_handler import HandlerTuple
    import matplotlib.pyplot as plt
    import numpy as np


    # ##### PANORAMA #####
    from ..spec.spec_support import _check_array


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    _check_array(teff, "teff", "Effective temperatures of the needed entries")
    #
    _check_array(logg, "logg", "Surface gravities of the needed entries", size=teff.size)
    #
    _check_array(in_grid, "in_grid", "In-grid flags of the needed entries", kind="bool", size=teff.size)
    #
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]


    # ####################  PARAMETERS  ####################
    markers = ["D", "o", "s", "^", "v", "P"]  # One per series of grid nodes (cycled)


    # ####################  COMPUTATION  ####################
    # ##########  GRID NODES  ##########
    nodes = {}  # Legend label -> distinct (Teff, log g) pairs


    # ##### OPTION 1: CATALOGS OF THE FIGURE #####
    if grid_info['plot_catalogs'] is not None:
        for label, file_name in grid_info['plot_catalogs'].items():
            catalog      = _read_catalog(grid, (file_name,))
            nodes[label] = np.unique(np.column_stack((catalog['teff'], catalog['logg'])), axis=0)


    # ##### OPTION 2: CATALOGS OF THE GRID, ONE SERIES PER SUBGRID #####
    else:
        catalog  = _read_catalog(grid)
        subgrids = np.array(catalog['subgrid'])
        for subgrid in dict.fromkeys(catalog['subgrid']):  # Subgrids in catalog order
            mask                                = subgrids == subgrid
            nodes[subgrid or grid_info['name']] = np.unique(np.column_stack((catalog['teff'][mask], catalog['logg'][mask])), axis=0)


    # ####################  PLOTTING  ####################
    # ##########  FIGURE  ##########
    fig, ax = plt.subplots(figsize=(6.69, 4.8))  # Two A&A columns (17 cm) wide


    # ##########  GRID NODES  ##########
    handles = []  # Legend symbols, filled together with 'labels' so that both lists always stay paired
    labels  = []
    for i, (label, xy) in enumerate(nodes.items()):
        handles.append(ax.scatter(
                                  xy[:, 0] / 1e3,  # K -> kK
                                  xy[:, 1],
                                  s=max(30 - 12 * i, 10),  # Decreasing sizes (30, 18, ... points^2): overlapping subgrids remain visible
                                  marker=markers[i % len(markers)],
                                  color=f"C{i % 10}",
                                  linewidths=0,
                                  zorder=1
                                 ))
        labels.append(f"{label} ({len(xy)})")


    # ##########  NEEDED ENTRIES  ##########
    # ##### INSIDE THE GRID #####
    if in_grid.any():
        handles.append(ax.scatter(
                                  teff[in_grid] / 1e3,
                                  logg[in_grid],
                                  s=80,  # Larger than the grid nodes
                                  marker="x",
                                  color="black",
                                  linewidths=1.5,
                                  zorder=3
                                 ))
        labels.append(f"Needed, in the grid ({int(in_grid.sum())})")


    # ##### OUTSIDE THE GRID #####
    if (~in_grid).any():
        cross  = ax.scatter(
                            teff[~in_grid] / 1e3,
                            logg[~in_grid],
                            s=80,
                            marker="x",
                            color="red",
                            linewidths=1.5,
                            zorder=3
                           )
        circle = ax.scatter(
                            teff[~in_grid] / 1e3,
                            logg[~in_grid],
                            s=260,  # Circle around the cross
                            marker="o",
                            facecolors="none",
                            edgecolors="red",
                            linewidths=1.2,
                            zorder=3
                           )
        handles.append((circle, cross))  # Drawn as one symbol in the legend (HandlerTuple)
        labels.append(f"Needed, outside the grid ({int((~in_grid).sum())})")


    # ##########  AXES  ##########
    all_teff = np.concatenate([xy[:, 0] for xy in nodes.values()] + [teff]) / 1e3  # Grid nodes and needed entries, in kK
    all_logg = np.concatenate([xy[:, 1] for xy in nodes.values()] + [logg])
    ax.set_xlim(all_teff.max() + 1.5, all_teff.min() - 1.5)  # Reversed: hot to the left; margin for the circled crosses
    ax.set_ylim(all_logg.max() + 0.25, all_logg.min() - 0.25)  # Reversed: low gravity at the top
    ax.set_xlabel(r"$T_{\rm eff}$ [kK]")
    ax.set_ylabel(r"$\log g$ [dex]")
    ax.tick_params(which="both", direction="out", top=True, right=True)


    # ##########  LEGEND  ##########
    # Above the axes, so that it never hides grid nodes or needed entries
    ax.legend(
              handles,
              labels,
              handler_map={tuple: HandlerTuple(ndivide=1)},  # Circle and cross overlaid in a single symbol
              loc="lower center",
              bbox_to_anchor=(0.5, 1.01),
              ncol=len(handles),
              fontsize=8,
              frameon=False,
              handletextpad=0.3,
              columnspacing=1.2
             )
    fig.tight_layout()


    # ####################  OUTPUT  ####################
    return fig




# ===================================================================
# ============================== FILES ==============================
# ===================================================================

# ==============================
def _locate_files(
                  model_names : list[str],
                  modeldir    : str | Path,
                  recursive   : bool,
                  extensions  : tuple[str, ...]
    ) -> dict[str, list[Path]]:
    """
    Locate in a directory the files of a list of models (model name followed by each extension of the grid).

    ----------
    Parameters
    ----------
        model_names : list[str]
            Model names, without file extension.

        modeldir : str | Path
            Directory where the model files are searched.

        recursive : bool
            If True, also search the subdirectories of 'modeldir'.

        extensions : tuple[str, ...]
            File extensions of a model, including the dot (e.g. ('.dat',),
            ('.7', '.5', '.nst')).

    -------
    Returns
    -------
        model_files : dict[str, list[Path]]
            A dictionary keyed by model name (every name of 'model_names');
            each value is the list of the files of the model found in
            'modeldir', in the order of 'extensions' ([] if none is found).

    ------
    Raises
    ------
        TypeError
            If 'model_names' is not a list of strings, 'modeldir' is not a
            path, 'recursive' is not a boolean, or 'extensions' is not a
            tuple of strings.

        ValueError
            If 'extensions' is empty or contains an extension that does not
            start with a dot.

        FileNotFoundError
            If 'modeldir' does not exist.

        NotADirectoryError
            If 'modeldir' is not a directory.

        PermissionError
            If 'modeldir' cannot be read.

    -----
    Notes
    -----
        - The directory is listed once (os.scandir, or os.walk with
          'recursive'), so the cost does not depend on the number of models.

        - With 'recursive', subdirectories are visited in alphabetical order;
          a file name found in more than one directory is taken from the
          first one, and a single warning is issued, since all the files are
          copied to the same output directory by 'find_models'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import os
    from pathlib import Path
    import warnings


    # ##### PANORAMA #####
    from ..obs.obs_support import _check_bool


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(model_names, list):
        raise TypeError(f"Model names 'model_names' must be a list of strings. Got {type(model_names).__name__}.")
    for model_name in model_names:
        if not isinstance(model_name, str):
            raise TypeError(f"Model names 'model_names' must contain only strings. Got {type(model_name).__name__}.")
    #
    modeldir = _check_directory(modeldir, "modeldir", "Model directory", "model directory")
    #
    recursive = _check_bool(recursive, "recursive", "Recursive search flag")
    #
    if not isinstance(extensions, tuple):
        raise TypeError(f"Model file extensions 'extensions' must be a tuple of strings. Got {type(extensions).__name__}.")
    elif len(extensions) == 0:
        raise ValueError("Model file extensions 'extensions' must contain at least one extension. Got an empty tuple.")
    for ext in extensions:
        if not isinstance(ext, str):
            raise TypeError(f"Model file extensions 'extensions' must contain only strings. Got {type(ext).__name__}.")
        elif not ext.startswith(".") or len(ext) < 2:
            raise ValueError(f"Model file extensions 'extensions' must start with a dot (e.g. '.dat'). Got '{ext}'.")


    # ####################  COMPUTATION  ####################
    # ##########  WANTED FILES  ##########
    wanted   = {model_name + ext for model_name in model_names for ext in extensions}
    found    = {}  # File name -> path of the first file with that name
    repeated = []


    # ##########  OPTION 1: RECURSIVE SEARCH  ##########
    if recursive:
        for root, dirs, files in os.walk(modeldir):
            dirs.sort()  # Visit the subdirectories in alphabetical order
            for file_name in sorted(files):
                if file_name in wanted:
                    if file_name in found:
                        repeated.append(Path(root) / file_name)
                    else:
                        found[file_name] = Path(root) / file_name


    # ##########  OPTION 2: TOP-LEVEL SEARCH  ##########
    else:
        with os.scandir(modeldir) as entries:
            for entry in entries:
                if entry.name in wanted and entry.is_file():
                    found[entry.name] = Path(entry.path)


    # ##########  CHECKS  ##########
    if repeated:
        warnings.warn(f"{len(repeated)} model file(s) of 'modeldir' share their name with a file of another subdirectory and were ignored (e.g. '{repeated[0]}').", stacklevel=3)


    # ####################  OUTPUT  ####################
    return {model_name: [found[model_name + ext] for ext in extensions if model_name + ext in found] for model_name in model_names}




# ======================================================================
# ============================== ARCHIVES ==============================
# ======================================================================

# ==============================
def _archive_name(
                  params : dict,
                  grid   : str,
                  vmic   : int
    ) -> str:
    """
    Build the name of a PRISMAS archive of modes 0 / 1 (without '.tar.gz') from the stellar parameters of its model, as 'prismas.py' does.

    ----------
    Parameters
    ----------
        params : dict
            Stellar parameters of the model, with the keys of '_decode_kc'
            ('teff', 'logg', 'met', 'alpha', 'carbon'; NaN if not defined
            for the grid).

        grid : str
            Model grid code, a key of '_grid_registry' (case-insensitive).

        vmic : int
            Microturbulence of the synthesis in km/s (integer, as in the
            names of 'prismas.py'; non-negative).

    -------
    Returns
    -------
        name : str
            The archive name, from the 'archive_name' template of the grid
            (e.g. 'T49000_G4.50_Mm0.500_ap0.400_V2', 'T30000_G4.00_Z1.000_V10').

    ------
    Raises
    ------
        TypeError
            If 'params' is not a dictionary, 'grid' is not a string, a
            parameter of the grid is not a real number, or 'vmic' is not an
            integer.

        ValueError
            If 'grid' is not a supported grid, 'vmic' is negative, or a
            parameter of the grid is missing or not finite.

    -----
    Notes
    -----
        - Teff is written without decimals, log g with 2 decimals, and the
          abundances with 3 decimals, signed with 'm' (-) / 'p' (+) for
          Kurucz ('Mm0.500', 'ap0.400'; 'p0.000' for 0) and unsigned for the
          TLUSTY Z / Zsun ('Z0.032'), as in the outputs of 'prismas.py'.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numbers
    import numpy as np


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(params, dict):
        raise TypeError(f"Stellar parameters 'params' must be a dictionary. Got {type(params).__name__}.")
    #
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]
    for key in grid_info['params']:
        if key not in params:
            raise ValueError(f"Stellar parameters 'params' must contain the key '{key}' for the '{grid}' grid. Got keys {list(params)}.")
        _check_real(params[key], key, "Stellar parameter")
    #
    if not isinstance(vmic, numbers.Integral) or isinstance(vmic, bool):
        raise TypeError(f"Microturbulence 'vmic' must be an integer (in km/s). Got {type(vmic).__name__}.")
    elif vmic < 0:
        raise ValueError(f"Microturbulence 'vmic' must be non-negative. Got {vmic}.")


    # ####################  COMPUTATION  ####################
    fields = {
              'teff' : f"{params['teff']:.0f}",
              'logg' : f"{params['logg']:.2f}",
              'z'    : f"{abs(params['met']):.3f}",  # TLUSTY Z / Zsun (unsigned)
              'vmic' : f"{int(vmic)}",
             }
    for key in ('met', 'alpha', 'carbon'):
        value       = params.get(key, np.nan)
        fields[key] = '' if np.isnan(value) else f"{'m' if value < 0 else 'p'}{abs(value):.3f}"  # 'p0.000' for 0


    # ####################  OUTPUT  ####################
    return grid_info['archive_name'].format(**fields)




# ==============================
def _allocated_cores(
    ) -> int:
    """
    Return the number of cores actually allocated to the process, not those of the whole node.

    -------
    Returns
    -------
        n_cores : int
            Number of allocated cores (>= 1).

    -----
    Notes
    -----
        - CPU affinity (set by HTCondor / Slurm), further limited by
          'SLURM_CPUS_PER_TASK' or 'OMP_NUM_THREADS' (set by HTCondor to
          'request_cpus') when they are defined.
    """

    # ####################  MODULES & PACKAGES  ####################
    import os


    # ####################  COMPUTATION  ####################
    # ##########  AFFINITY  ##########
    if hasattr(os, "process_cpu_count"):  # Python >= 3.13
        n_cores = os.process_cpu_count()
    elif hasattr(os, "sched_getaffinity"):  # Linux
        n_cores = len(os.sched_getaffinity(0))
    else:  # macOS, Windows
        n_cores = os.cpu_count()
    n_cores = n_cores or 1  # None if the count cannot be determined


    # ##########  SCHEDULER LIMITS  ##########
    for var in ("SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS"):
        try:
            value = int(os.environ.get(var, ''))
        except ValueError:  # Not defined or not an integer
            continue
        if value > 0:
            n_cores = min(n_cores, value)


    # ####################  OUTPUT  ####################
    return max(1, n_cores)




# ==============================
def _run_parallel(
                  jobs     : list,
                  func     : Callable[[object], object],
                  ncores   : int | None,
                  callback : Callable[[], object] | None = None,
                  gain     : float = 1.10
    ) -> tuple[int, bool]:
    """
    Run 'func(job)' for every job in threads, with a fixed or self-tuned number of workers.

    ----------
    Parameters
    ----------
        jobs : list
            Items passed one by one to 'func'.

        func : Callable[[object], object]
            Function run for each job (in a worker thread).

        ncores : int | None
            Fixed number of workers (>= 1), or None for self-tuning.

        callback : Callable[[], object] | None, optional
            Function without arguments called in the main thread after each
            completed job (e.g. the 'update' method of a progress bar).
            Default is None.

        gain : float, optional
            Minimum throughput ratio (> 1) to keep doubling the number of
            workers. Default is 1.10.

    -------
    Returns
    -------
        workers : int
            Number of workers in use at the end (0 if 'jobs' is empty).

        tuned : bool
            True if the number of workers was self-tuned.

    ------
    Raises
    ------
        TypeError
            If 'jobs' is not a list, 'func' or 'callback' is not callable,
            'ncores' is not an integer or None, or 'gain' is not a real
            number.

        ValueError
            If 'ncores' is smaller than 1 or 'gain' is not greater than 1.

    -----
    Notes
    -----
        1. SELF-TUNING
            - With 'ncores' = None the concurrency is tuned on the fly: it
              starts with 2 workers and doubles (2, 4, 8, ...) while the
              measured throughput (jobs / s) improves by at least 'gain', up
              to the cores allocated to the process ('_allocated_cores').
              When a step brings no gain, the previous level is restored and
              kept, so CPU-bound runs grow to the allocated cores and
              I/O-bound runs (e.g. network file systems) stop where the disk
              saturates.

            - After each level change one generation of jobs is not measured
              (transition), and the level is then measured over 2 x workers
              (>= 4) completed jobs.

        2. ERRORS
            - Unexpected exceptions of 'func' are re-raised in the main
              thread. On an interruption (e.g. Ctrl+C) the queue is dropped
              and the running jobs finish cleanly.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
    import time


    # ##### PANORAMA #####
    from ..units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(jobs, list):
        raise TypeError(f"Jobs 'jobs' must be a list. Got {type(jobs).__name__}.")
    #
    if not callable(func):
        raise TypeError(f"Job function 'func' must be callable. Got {type(func).__name__}.")
    #
    ncores = _check_ncores(ncores)
    #
    if callback is not None and not callable(callback):
        raise TypeError(f"Progress function 'callback' must be callable or None. Got {type(callback).__name__}.")
    #
    gain = _check_real(gain, "gain", "Throughput gain", unit="ratio", constraint="positive")
    if gain <= 1:
        raise ValueError(f"Throughput gain 'gain' must be greater than 1. Got {gain}.")


    # ####################  PARAMETERS  ####################
    tuned = ncores is None
    max_w = min(ncores or _allocated_cores(), len(jobs))  # Never more workers than jobs
    limit = min(2, max_w) if tuned else max_w  # Jobs kept in flight


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: SERIAL RUN  ##########
    if max_w <= 1:
        for job in jobs:
            func(job)
            if callback is not None:
                callback()
        return min(1, len(jobs)), tuned


    # ##########  OPTION 2: PARALLEL RUN  ##########
    queue      = list(reversed(jobs))  # pop() from the end keeps the original order
    running    = set()
    tuning     = tuned and limit < max_w
    prev_rate  = None  # Throughput of the previous level (jobs / s)
    prev_limit = None  # Number of workers of the previous level
    warmup     = 0  # Completions ignored after a level change
    win_start  = time.perf_counter()
    win_done   = 0
    pool       = ThreadPoolExecutor(max_workers=max_w)
    try:
        while queue or running:
            while queue and len(running) < limit:  # Keep 'limit' jobs in flight
                running.add(pool.submit(func, queue.pop()))
            done, running = wait(running, return_when=FIRST_COMPLETED)
            for future in done:
                future.result()  # Re-raise unexpected errors
                if callback is not None:
                    callback()
            if warmup > 0:  # Transition to the new level: not measured
                warmup -= len(done)
                if warmup <= 0:
                    win_start = time.perf_counter()
                    win_done  = 0
                continue
            win_done += len(done)


            # ##### SELF-TUNING #####
            if tuning and win_done >= max(4, 2 * limit):
                rate = win_done / (time.perf_counter() - win_start)
                if prev_rate is None or rate >= gain * prev_rate:  # Worth it: keep doubling
                    prev_rate  = rate
                    prev_limit = limit
                    limit      = min(2 * limit, max_w)
                    tuning     = limit > prev_limit
                else:  # No gain: go back (less contention) and stay
                    limit  = prev_limit
                    tuning = False
                warmup = limit  # Let one generation of jobs at the new level start
    except BaseException:  # E.g. Ctrl+C: drop the queue; running jobs finish cleanly
        queue.clear()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)


    # ####################  OUTPUT  ####################
    return limit, tuned




# ==============================
def _gzip_backend(
    ) -> tuple[ModuleType, tuple[type, ...]]:
    """
    Return the ISA-L gzip module used to read and write the PRISMAS archives and the exceptions it raises on corrupted data.

    -------
    Returns
    -------
        gz : ModuleType
            The 'isal.igzip' module (same interface as the standard 'gzip'
            module, including 'open').

        gz_errors : tuple[type, ...]
            Exceptions raised on corrupted data: zlib.error and the ISA-L
            error ('isal.isal_zlib.error').

    -----
    Notes
    -----
        - ISA-L ('isal', the Python bindings of the Intel ISA-L library) is a
          required dependency of PANORAMA: it compresses much faster than
          zlib and releases the GIL, so that the archives can be rewritten in
          parallel threads. The result is a standard '.tar.gz', readable by
          any gzip implementation.

        - zlib.error is kept among the exceptions as a safeguard, so that any
          corrupted-data error is reported as an 'error' status by
          '_rename_archive' instead of stopping the run.
    """

    # ####################  MODULES & PACKAGES  ####################
    from isal import igzip
    from isal.isal_zlib import error as isal_error
    import zlib


    # ####################  COMPUTATION  ####################
    gz        = igzip
    gz_errors = (zlib.error, isal_error)


    # ####################  OUTPUT  ####################
    return gz, gz_errors




# ==============================
def _rename_member(
                   name : str,
                   old  : str,
                   new  : str
    ) -> str:
    """
    Rename one archive member path ('old/old_wave.npy' -> 'new/new_wave.npy').

    ----------
    Parameters
    ----------
        name : str
            Member path inside the archive.

        old : str
            Old model name (folder name and file prefix).

        new : str
            New model name.

    -------
    Returns
    -------
        renamed : str
            Member path with every component equal to 'old', or starting with
            'old_', renamed.

    ------
    Raises
    ------
        TypeError
            If 'name', 'old' or 'new' is not a string.

        ValueError
            If 'old' or 'new' is empty.
    """

    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(name, str):
        raise TypeError(f"Archive member path 'name' must be a string. Got {type(name).__name__}.")
    #
    for arg_name, arg_value in (("old", old), ("new", new)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Model name '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value:
            raise ValueError(f"Model name '{arg_name}' must be a non-empty string. Got ''.")


    # ####################  COMPUTATION  ####################
    parts = name.split("/")  # Tar member paths always use '/'
    for i, part in enumerate(parts):
        if part == old:  # Folder
            parts[i] = new
        elif part.startswith(old + "_"):  # Files ('old_wave.npy')
            parts[i] = new + part[len(old):]


    # ####################  OUTPUT  ####################
    return "/".join(parts)




# ==============================
def _rewrite_tar(
                 src : str | Path,
                 dst : str | Path,
                 old : str,
                 new : str
    ) -> None:
    """
    Rewrite a '.tar.gz' archive with its folder and files renamed ('old' -> 'new').

    ----------
    Parameters
    ----------
        src : str | Path
            Archive to read.

        dst : str | Path
            Archive to write (replaced only when the rewrite is complete).

        old : str
            Old model name (folder name and file prefix).

        new : str
            New model name.

    -------
    Returns
    -------
        None
            The renamed archive is written to 'dst'.

    ------
    Raises
    ------
        TypeError
            If 'src' or 'dst' is not a path, or 'old' or 'new' is not a
            string.

        ValueError
            If 'old' or 'new' is empty.

        tarfile.TarError
            If the source archive is not a valid tar stream.

        EOFError
            If the source archive is truncated.

        OSError
            If the source archive is not a valid gzip file or the output
            cannot be written.

        zlib.error
            If the compressed data are corrupted (also the ISA-L error
            'isal.isal_zlib.error'; see '_gzip_backend').

    -----
    Notes
    -----
        1. STREAMING
            - The rewrite is streamed member by member (nothing is extracted
              to disk, compresslevel=1 as in 'prismas.py'), written to a
              temporary '.part' file and moved into place only when
              complete. On any error the '.part' file is removed and the
              error re-raised.

            - The whole gzip stream is read, so its CRC is checked: a
              corrupted or truncated archive (e.g. a killed HTCondor job)
              raises an error and never replaces anything.

        2. GZIP BACKEND
            - The gzip layer is handled outside 'tarfile' (which only sees an
              uncompressed stream) with '_gzip_backend' (ISA-L): compressing
              is most of the cost of a rewrite.
    """

    # ####################  MODULES & PACKAGES  ####################
    import os
    from pathlib import Path
    import tarfile


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("src", src), ("dst", dst)):
        if not isinstance(arg_value, (str, os.PathLike)):
            raise TypeError(f"Archive path '{arg_name}' must be a path (str or pathlib.Path). Got {type(arg_value).__name__}.")
    src = Path(src)
    dst = Path(dst)
    #
    for arg_name, arg_value in (("old", old), ("new", new)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Model name '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value:
            raise ValueError(f"Model name '{arg_name}' must be a non-empty string. Got ''.")


    # ####################  PARAMETERS  ####################
    gz, _  = _gzip_backend()
    buffer = _TAR_BUFFER
    tmp    = dst.with_name(dst.name + ".part")  # Moved into place only when complete


    # ####################  COMPUTATION  ####################
    try:
        with gz.open(src, "rb") as fin, gz.open(tmp, "wb", compresslevel=1) as fout:
            with tarfile.open(fileobj=fin, mode="r|", bufsize=buffer, copybufsize=buffer) as tin, tarfile.open(fileobj=fout, mode="w|", bufsize=buffer, copybufsize=buffer) as tout:
                for member in tin:  # Streaming: one sequential pass
                    fileobj     = tin.extractfile(member) if member.isfile() else None
                    member.name = _rename_member(member.name, old, new)
                    tout.addfile(member, fileobj)
            # 'tarfile' stops at the end-of-archive blocks, before the gzip trailer: read it so that the CRC is checked
            while fin.read(buffer):
                pass
        os.replace(tmp, dst)  # Atomic on the same file system
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise




# ==============================
def _rename_archive(
                    archive        : dict,
                    rename_content : bool,
                    keep_original  : bool
    ) -> None:
    """
    Rename one PRISMAS archive (run in a worker thread of 'rename_tars'; updates only its own entry).

    ----------
    Parameters
    ----------
        archive : dict
            Entry of 'rename_tars', updated in place. A dictionary containing
            at least the following key-value pairs:
                - 'source' : Path
                    Original archive.

                - 'target' : Path
                    Renamed archive.

                - 'status' : str
                    Set to 'renamed', or to 'error' if the archive cannot be
                    read or written.

                - 'reason' : str
                    Set to the error message when 'status' is 'error'.

        rename_content : bool
            If True, rewrite the archive with its folder and files renamed
            too ('_rewrite_tar'); if False, only rename (or copy) the file.

        keep_original : bool
            If True, keep the original archive.

    -------
    Returns
    -------
        None
            The result is stored in 'archive' ('status' and 'reason').

    ------
    Raises
    ------
        TypeError
            If 'archive' is not a dictionary with Path values for 'source'
            and 'target', or 'rename_content' or 'keep_original' is not a
            boolean.

    -----
    Notes
    -----
        - Expected failures (truncated or corrupted archives, unwritable
          output) are stored as 'error' instead of being raised, so that one
          bad archive does not stop the others.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from pathlib import Path
    import shutil
    import tarfile


    # ##### PANORAMA #####
    from ..obs.obs_support import _check_bool


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(archive, dict) or not all(isinstance(archive.get(k), Path) for k in ('source', 'target')):
        raise TypeError(f"Archive entry 'archive' must be a dictionary with pathlib.Path values for 'source' and 'target'. Got {type(archive).__name__}.")
    #
    rename_content = _check_bool(rename_content, "rename_content", "Content renaming flag")
    #
    keep_original = _check_bool(keep_original, "keep_original", "Original archive flag")


    # ####################  PARAMETERS  ####################
    _, gz_errors = _gzip_backend()
    src          = archive['source']
    dst          = archive['target']


    # ####################  COMPUTATION  ####################
    try:
        if rename_content:  # Rewrite with the folder and files renamed
            _rewrite_tar(src, dst, src.name[:-len(".tar.gz")], dst.name[:-len(".tar.gz")])
            if not keep_original:
                src.unlink()
        elif keep_original:  # Copy under the new name
            shutil.copy2(src, dst)
        else:  # Move under the new name
            shutil.move(src, dst)
        archive['status'] = "renamed"
    except (tarfile.TarError, EOFError, OSError, *gz_errors) as err:  # Truncated / corrupted / unwritable
        archive.update({
                        'status' : "error",
                        'reason' : f"{type(err).__name__}: {err}",
                       })




# ==================================================================
# ============================== TEXT ==============================
# ==================================================================

# ==============================
def _format_value(
                  value : float,
                  spec  : str
    ) -> str:
    """
    Format a stellar parameter for a summary table, writing '-' for the parameters not defined for the grid (NaN).

    ----------
    Parameters
    ----------
        value : float
            Parameter value (NaN if not defined for the grid).

        spec : str
            Format specification (e.g. '6.2f').

    -------
    Returns
    -------
        text : str
            The formatted value, or '-' if 'value' is NaN.

    ------
    Raises
    ------
        TypeError
            If 'value' is not a real number or 'spec' is not a string.

        ValueError
            If 'spec' is not a valid format specification for a float.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(value, numbers.Real) or isinstance(value, bool):
        raise TypeError(f"Parameter value 'value' must be a real number (NaN if not defined). Got {type(value).__name__}.")
    #
    if not isinstance(spec, str):
        raise TypeError(f"Format specification 'spec' must be a string. Got {type(spec).__name__}.")
    try:
        format(0.0, spec)
    except ValueError as e:
        raise ValueError(f"Format specification 'spec' must be valid for a float (e.g. '6.2f'). Got '{spec}'.") from e


    # ####################  OUTPUT  ####################
    return "-" if np.isnan(value) else format(value, spec)




# ==============================
def _short_list(
                names     : list[str],
                max_names : int = 5
    ) -> str:
    """
    Format a list of names for a summary line, showing at most 'max_names' of them.

    ----------
    Parameters
    ----------
        names : list[str]
            Names to show.

        max_names : int, optional
            Maximum number of names shown (>= 1). Default is 5.

    -------
    Returns
    -------
        text : str
            Comma-separated names, followed by ', ...' if there are more than
            'max_names' ('' for an empty list).

    ------
    Raises
    ------
        TypeError
            If 'names' is not a list of strings or 'max_names' is not an
            integer.

        ValueError
            If 'max_names' is smaller than 1.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(names, list) or not all(isinstance(n, str) for n in names):
        raise TypeError(f"Names 'names' must be a list of strings. Got {type(names).__name__}.")
    #
    if not isinstance(max_names, numbers.Integral) or isinstance(max_names, bool):
        raise TypeError(f"Maximum number of names 'max_names' must be an integer. Got {type(max_names).__name__}.")
    elif max_names < 1:
        raise ValueError(f"Maximum number of names 'max_names' must be >= 1. Got {max_names}.")


    # ####################  OUTPUT  ####################
    return ", ".join(names[:max_names]) + (", ..." if len(names) > max_names else '')




# ==============================
def _message_block(
                   messages : dict[str, list[str]],
                   title    : str,
                   width    : int = 100
    ) -> str:
    """
    Format the messages of several input files as a titled block, with the file names aligned and long messages wrapped.

    ----------
    Parameters
    ----------
        messages : dict[str, list[str]]
            A dictionary keyed by input file name; each value is the
            non-empty list of messages of that file (e.g. SPAMMS warnings or
            the reason of a failed run).

        title : str
            Title of the block (e.g. 'SPAMMS warnings').

        width : int, optional
            Maximum length of the lines, in characters (>= 40). Default is
            100.

    -------
    Returns
    -------
        block : str
            The block, without trailing newline, e.g.:
                >>> SPAMMS failures (1 input file):
                >>>   input.txt : SPAMMS exited with status 1
                >>>               (ModuleNotFoundError: No module named 'phoebe')

    ------
    Raises
    ------
        TypeError
            If 'messages' is not a dictionary of lists of strings, 'title'
            is not a string or 'width' is not an integer.

        ValueError
            If 'messages' is empty or has an empty list, 'title' is empty,
            or 'width' is smaller than 40.

    -----
    Notes
    -----
        - Used by 'needed_models' for the SPAMMS warnings and failures, so
          that they are readable whatever the Python warning filters.

        - Each message starts on its own line at the column after the file
          names; its continuation lines are indented two more characters.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import textwrap


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(messages, dict):
        raise TypeError(f"Messages of the input files 'messages' must be a dictionary of lists of strings. Got {type(messages).__name__}.")
    elif len(messages) == 0:
        raise ValueError("Messages of the input files 'messages' must contain at least one input file. Got an empty dictionary.")
    for file_name, file_messages in messages.items():
        if not isinstance(file_name, str) or not isinstance(file_messages, list) or not all(isinstance(m, str) for m in file_messages):
            raise TypeError(f"Messages of the input files 'messages' must map each file name (str) to a list of strings. Got {file_messages!r} for {file_name!r}.")
        elif len(file_messages) == 0:
            raise ValueError(f"Messages of the input files 'messages' must contain at least one message per file. Got none for '{file_name}'.")
    #
    if not isinstance(title, str):
        raise TypeError(f"Block title 'title' must be a string. Got {type(title).__name__}.")
    elif not title.strip():
        raise ValueError(f"Block title 'title' must be a non-empty string. Got '{title}'.")
    #
    if not isinstance(width, numbers.Integral) or isinstance(width, bool):
        raise TypeError(f"Line width 'width' (in characters) must be an integer. Got {type(width).__name__}.")
    elif width < 40:
        raise ValueError(f"Line width 'width' (in characters) must be >= 40. Got {width}.")


    # ####################  PARAMETERS  ####################
    name_width = max(len(file_name) for file_name in messages)
    indent     = " " * (2 + name_width + 3)  # Column of the messages: '  <name> : '
    text_width = max(width - len(indent), 20)  # Room left for the messages


    # ####################  COMPUTATION  ####################
    n_files = len(messages)
    lines   = [f"{title} ({n_files} input file{'s' if n_files > 1 else ''}):"]
    for file_name, file_messages in messages.items():
        for i, message in enumerate(file_messages):
            wrapped = textwrap.wrap(message, width=text_width, subsequent_indent="  ", break_long_words=False, break_on_hyphens=False) or ['']
            prefix  = f"  {file_name:<{name_width}} : " if i == 0 else indent  # File name only on its first message
            lines.append(prefix + wrapped[0])
            lines.extend(indent + line for line in wrapped[1:])


    # ####################  OUTPUT  ####################
    return "\n".join(lines)
