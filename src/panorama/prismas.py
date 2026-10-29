"""
Model atmosphere utilities for preparing and post-processing the runs of
PRISMAS, the code that computes grids of synthetic spectra from model
atmospheres ('prismas.py'), within the spectroscopic analysis of (massive)
stars: selection of the model atmospheres of a grid, grid check of SPAMMS and
renaming of the archives produced in HPC mode.

Conventions shared by every function in this module:
    - Model grids are identified by the same codes as '--mode' in
      'prismas.py' ('KC', 'KM', 'T'; case-insensitive), and the models by the
      names (without file extension) listed in the catalogs distributed with
      PANORAMA ('support/prismas/models_grids/').
    - Stellar parameters follow the PANORAMA units: 'teff' in K, 'logg' as
      log10(g / [cm s^-2]), 'met' as [M/H] in dex (Kurucz) or Z / Zsun
      (TLUSTY), 'alpha' and 'carbon' in dex; parameters not defined for a
      grid are NaN. Microturbulences are in km/s.
    - Directory and file arguments are paths, given either as strings or as
      pathlib.Path objects, and are returned as pathlib.Path objects. Output
      directories are never created.
    - All functions validate their inputs before computing.
"""

from __future__ import annotations

# ====================================================================
# ============================== MODELS ==============================
# ====================================================================

# ==============================
def find_models(
                grid        : str,
                teff        : int | float | list[float] | tuple[float | None, float | None] | None = None,
                logg        : int | float | list[float] | tuple[float | None, float | None] | None = None,
                met         : int | float | list[float] | tuple[float | None, float | None] | None = None,
                alpha       : int | float | list[float] | tuple[float | None, float | None] | None = None,
                carbon      : int | float | list[float] | tuple[float | None, float | None] | None = None,
                modeldir    : str | Path | None = None,
                recursive   : bool = False,
                output_dir  : str | Path | None = None,
                output_file : str | Path | None = None,
                verbose     : bool = True
    ) -> list[dict]:
    """
    Search the PRISMAS model catalog of a grid for the model atmospheres whose stellar parameters fall within the requested ranges.

    ----------
    Parameters
    ----------
        grid : str
            Model atmosphere grid (same code as '--mode' in 'prismas.py';
            case-insensitive): 'KC' (Kurucz-Castelli), 'KM'
            (Kurucz-Mészáros) or 'T' (TLUSTY OSTAR2002/BSTAR2006). See Notes.

        teff : int | float | list[float] | tuple[float | None, float | None] | None, optional
            Effective temperature selection in K (values > 0). See Notes for
            the accepted formats. Default is None (no filter).

        logg : int | float | list[float] | tuple[float | None, float | None] | None, optional
            Surface gravity selection, log10(g / [cm s^-2]). Default is None
            (no filter).

        met : int | float | list[float] | tuple[float | None, float | None] | None, optional
            Metallicity selection: [M/H] in dex for 'KC' / 'KM', and Z / Zsun
            (values >= 0) for 'T'. Default is None (no filter).

        alpha : int | float | list[float] | tuple[float | None, float | None] | None, optional
            Alpha enhancement [alpha/M] selection in dex. Only defined for
            'KC' and 'KM'. Default is None (no filter).

        carbon : int | float | list[float] | tuple[float | None, float | None] | None, optional
            Carbon abundance [C/M] selection in dex. Only defined for 'KM'.
            Default is None (no filter).

        modeldir : str | Path | None, optional
            Existing directory where the model files are stored (same as
            '--modeldir' in 'prismas.py'). If given, the files of the
            selected models are located in it. Default is None (only the
            catalog is searched).

        recursive : bool, optional
            If True, also search the subdirectories of 'modeldir'. Default is
            False.

        output_dir : str | Path | None, optional
            Existing, writable directory, different from 'modeldir', where
            the files of the selected models found complete in 'modeldir'
            (every file of the grid, e.g. '.7', '.5' and '.nst' for 'T') are
            copied; existing files with the same name are overwritten.
            Requires 'modeldir'. Default is None (nothing copied).

        output_file : str | Path | None, optional
            Path and filename of a text file where the names of the selected
            models (without file extension) are written, one per line. Its
            directory must exist and an existing file is overwritten. Default
            is None (no file written).

        verbose : bool, optional
            If True, print a summary table of the selected models (the first
            50 rows when more models are selected), with the file(s) of each
            model found in 'modeldir' when it is given. Default is True.

    -------
    Returns
    -------
        models : list[dict]
            Models that fulfill every criterion, sorted by (teff, logg, met,
            alpha, carbon).
            A list of dictionaries, each containing the following key-value
            pairs:
                - 'model' : str
                    Model name, without file extension (as listed in the
                    catalog).

                - 'grid' : str
                    Model grid code ('KC', 'KM' or 'T').

                - 'subgrid' : str
                    'OSTAR2002' or 'BSTAR2006' for 'T'; '' for Kurucz.

                - 'teff' : float
                    Effective temperature in K.

                - 'logg' : float
                    Surface gravity, log10(g / [cm s^-2]).

                - 'met' : float
                    [M/H] in dex for 'KC' / 'KM'; Z / Zsun for 'T'.

                - 'alpha' : float
                    Alpha enhancement [alpha/M] in dex (NaN for 'T').

                - 'carbon' : float
                    Carbon abundance [C/M] in dex (NaN except for 'KM').

                - 'files' : list[Path]
                    Files of the model found in 'modeldir', in the order of
                    the extensions of the grid ('.dat' for 'KC', '.mod' for
                    'KM', '.7', '.5' and '.nst' for 'T'). They are copied to
                    'output_dir' only if every file of the model is found.
                    [] if 'modeldir' is None or no file of the model is
                    found.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If 'grid' is not supported, a criterion is invalid or not defined
            for the grid, or 'output_dir' is given without 'modeldir', is the
            same directory as 'modeldir', or lies inside it with 'recursive'.

        FileNotFoundError
            If 'modeldir', 'output_dir' or the directory of 'output_file' does
            not exist, or the catalog of the grid is not installed.

        NotADirectoryError
            If 'modeldir', 'output_dir' or the parent of 'output_file' exists
            but is not a directory.

        IsADirectoryError
            If 'output_file' is an existing directory.

        PermissionError
            If 'modeldir' cannot be read, 'output_dir' cannot be written, or
            the catalog cannot be read.

        RuntimeError
            If the catalog of the grid contains no valid model name.

    -----
    Notes
    -----
        1. GRIDS AND CATALOGS
            - The models are searched in the catalogs stored in
              'support/prismas/models_grids/' ('KC_models.dat',
              'KM_models.dat', 'T_OSTAR02+BSTAR06models.dat'), one model name
              per line without file extension, and the stellar parameters are
              decoded from the names:
                >>> 'KC' : 'm0.50t03500g0.0am05ak2odfnew'
                >>> 'KM' : 'amm03cm03om03t10000g20v20'
                >>> 'T'  : 'G30000g400v10' (OSTAR2002), 'BG20000g300v2' (BSTAR2006)

            - Names that do not follow the naming of the grid are skipped
              with a warning. For 'KC', an 'a' before 'k2' marks the
              alpha-enhanced models ([alpha/M] = +0.4). For 'KM', the
              abundance codes ending in 3 or 8 are quarter steps ('m03' =
              -0.25). For 'T', Z / Zsun is read from the composition letter
              with the values of 'prismas.py' (e.g. 'G' = 1, 'V' = 0.032).

            - Each catalog is read once per Python session and kept in
              memory; the first 'KM' request (about 850000 models) takes a
              few seconds.

            - The 'T' catalog lists BSTAR2006 from 15000 to 28000 K and
              OSTAR2002 from 27500 to 55000 K: in the overlap of both subgrids
              (27500-30000 K) the OSTAR2002 models are recommended, and only
              the BSTAR2006 models at 28000 K are kept, as the upper node of
              the interpolation for 27000 K < Teff < 27500 K (e.g. 27400 K is
              interpolated between BSTAR2006 at 27000 and 28000 K, so that
              the two subgrids are never mixed). The full listings of both
              subgrids are only used for the figure of 'needed_models'.

            - New grids are added in '_grid_registry' (support module),
              without changing this function.

        2. SELECTION CRITERIA
            - Every criterion ('teff', 'logg', 'met', 'alpha', 'carbon')
              accepts:
                >>> None              : no filter.
                >>> 25000             : single value.
                >>> [25000, 27500]    : list of discrete values.
                >>> (25000, 35000)    : closed range [min, max].
                >>> (None, 3.5)       : open range (any limit can be None).

            - Comparisons use an absolute tolerance of 1e-6, so fractional
              TLUSTY metallicities (0.032, 0.020, ...) can also be selected
              with a range, e.g. met=(0.03, 0.04).

        3. MODEL FILES
            - The files of a model are its name followed by each extension
              of the grid: '<model>.dat' for 'KC', '<model>.mod' for 'KM', and
              '<model>.7', '<model>.5' and '<model>.nst' for 'T' (all needed
              by 'prismas.py').

            - Only complete models (every file of the grid found in
              'modeldir') are copied to 'output_dir'. Selected models with no
              file in 'modeldir' keep 'files' = [], and TLUSTY models with
              only some of their files keep the files found but are not
              copied; each case is reported in a single warning.

            - With 'recursive', the subdirectories are visited in
              alphabetical order; a file name found in more than one
              directory is taken from the first one and reported in a single
              warning, since all the files are copied to 'output_dir'.

        4. TYPICAL USE
            - Select the models of a PRISMAS run and write the list of names
              for the HTCondor submission file ('output_file'), or copy them
              from the full grid to a working directory ('modeldir' +
              'output_dir'):
                >>> models = find_models("T", teff=(25000, 35000), logg=(3.0, 4.0), met=1, output_file="T_run.dat")

    ----------
    References
    ----------
        - Castelli, F. & Kurucz, R. L. 2003, IAUS, 210, A20.
        - Mészáros, Sz., Allende Prieto, C., Edvardsson, B. et al. 2012,
          AJ, 144, 120.
        - Lanz, T. & Hubeny, I. 2003, ApJS, 146, 417.
        - Lanz, T. & Hubeny, I. 2007, ApJS, 169, 83.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numpy as np
    from pathlib import Path
    import shutil
    import warnings


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_bool
    from .support.prismas.prismas_support import _check_criterion, _check_directory, _check_grid, _format_value, _grid_registry, _locate_files, _PARAM_TOL, _passes, _read_catalog
    from .support.tex.tex_support import _check_output_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]
    #
    for crit_name, crit_value in (("teff", teff), ("logg", logg), ("met", met), ("alpha", alpha), ("carbon", carbon)):
        if crit_value is not None and crit_name not in grid_info['params']:
            raise ValueError(f"Selection criterion '{crit_name}' must be None for the '{grid}' grid ({grid_info['name']}), which defines only {', '.join(repr(p) for p in grid_info['params'])}. Got {crit_value!r}.")
    #
    teff   = _check_criterion(teff, "teff", "Effective temperature", unit="in K", constraint="positive")
    logg   = _check_criterion(logg, "logg", "Surface gravity", unit="log10(g / [cm s^-2])")
    met    = _check_criterion(met, "met", "Metallicity", unit=grid_info['met_unit'], constraint=grid_info['met_constraint'])
    alpha  = _check_criterion(alpha, "alpha", "Alpha enhancement", unit="[alpha/M] in dex")
    carbon = _check_criterion(carbon, "carbon", "Carbon abundance", unit="[C/M] in dex")
    #
    modeldir = _check_directory(modeldir, "modeldir", "Model directory", "model directory", allow_none=True)
    #
    recursive = _check_bool(recursive, "recursive", "Recursive search flag")
    #
    output_dir = _check_directory(output_dir, "output_dir", "Output directory", "output directory", readable=False, writable=True, allow_none=True)
    if output_dir is not None:
        if modeldir is None:
            raise ValueError("Output directory 'output_dir' must be given together with the model directory 'modeldir' (files are copied from it). Got modeldir=None.")
        elif output_dir.resolve() == modeldir.resolve():
            raise ValueError(f"Output directory 'output_dir' must be different from 'modeldir'. Got '{output_dir}'.")
        elif recursive and modeldir.resolve() in output_dir.resolve().parents:
            raise ValueError(f"Output directory 'output_dir' must not lie inside 'modeldir' when 'recursive' is True (copied files would be found again). Got '{output_dir}'.")
    #
    output_file = _check_output_file(output_file, allow_none=True)
    #
    verbose = _check_bool(verbose, "verbose", "Summary table flag")


    # ####################  PARAMETERS  ####################
    # ##########  CRITERIA  ##########
    criteria = {
                'teff'   : teff,
                'logg'   : logg,
                'met'    : met,
                'alpha'  : alpha,
                'carbon' : carbon,
               }


    # ##########  TABLE ROWS  ##########
    max_rows = 50  # Rows of the summary table


    # ####################  COMPUTATION  ####################
    # ##########  READING  ##########
    catalog = _read_catalog(grid)  # Columns of the whole catalog (cached)


    # ##########  SELECTION  ##########
    mask = np.ones(len(catalog['model']), dtype=bool)
    for key, crit in criteria.items():
        mask &= _passes(catalog[key], crit, _PARAM_TOL)  # Absolute tolerance of 1e-6, in the units of each parameter
    selected = np.flatnonzero(mask)


    # ##########  SORTING  ##########
    # np.lexsort sorts by the last key first and is stable (ties keep the catalog order); NaN values go last
    order    = np.lexsort(tuple(catalog[key][selected] for key in ('carbon', 'alpha', 'met', 'logg', 'teff')))
    selected = selected[order]


    # ##########  RESULTS LIST  ##########
    models = []
    for i in selected:
        models.append({
                       'model'   : catalog['model'][i],
                       'grid'    : grid,
                       'subgrid' : catalog['subgrid'][i],
                       'teff'    : float(catalog['teff'][i]),
                       'logg'    : float(catalog['logg'][i]),
                       'met'     : float(catalog['met'][i]),
                       'alpha'   : float(catalog['alpha'][i]),
                       'carbon'  : float(catalog['carbon'][i]),
                       'files'   : [],
                      })


    # ##########  FILES OF EACH MODEL  ##########
    if modeldir is not None and models:
        model_files = _locate_files([m['model'] for m in models], modeldir, recursive, grid_info['extensions'])
        for m in models:
            m['files'] = model_files[m['model']]


        # ##### CHECKS #####
        n_ext      = len(grid_info['extensions'])
        missing    = [m['model'] for m in models if not m['files']]
        incomplete = [m['model'] for m in models if 0 < len(m['files']) < n_ext]
        if missing:
            warnings.warn(f"{len(missing)} of the {len(models)} selected model(s) have no file in 'modeldir' ('{modeldir}') (e.g. '{missing[0]}').", stacklevel=2)
        if incomplete:
            not_copied = " and are not copied to 'output_dir'" if output_dir is not None else ''
            warnings.warn(f"{len(incomplete)} of the {len(models)} selected model(s) lack some of their files ({', '.join(grid_info['extensions'])}) in 'modeldir' ('{modeldir}'){not_copied} (e.g. '{incomplete[0]}').", stacklevel=2)


    # ####################  OUTPUT  ####################
    # ##########  WRITE OUTPUT  ##########
    if output_file is not None:
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(''.join(m['model'] + "\n" for m in models))  # One name per line, without header


    # ##########  COPY  ##########
    copied_models = []
    copied_files  = []
    if output_dir is not None:
        for m in models:
            if len(m['files']) < len(grid_info['extensions']):  # Incomplete model: 'prismas.py' needs every file, so it is not copied
                continue
            copied_models.append(m['model'])
            for src in m['files']:
                copied_files.append(Path(shutil.copy2(src, output_dir / src.name)))  # Metadata (e.g. modification time) preserved


    # ##########  SUMMARY  ##########
    if verbose:
        print(f"\n{len(models)} {grid} model(s) selected from the {grid_info['name']} catalog")


        # ##### TABLE #####
        if models:
            # With 'modeldir', the last column gives the file(s) of each model found there: the path of the first one
            # relative to 'modeldir', followed by the extensions of the others (e.g. 'BSTAR2006/BG20000g300v2.7, .5, .nst')
            width  = max(len("model"), max(len(m['model']) for m in models[:max_rows]))
            found  = "  found in 'modeldir'" if modeldir is not None else ''
            header = f"  {'Teff':>7} {'logg':>5} {grid_info['met_label']:>7} {'alpha':>6} {'C':>6}  {'model':<{width}}{found}".rstrip()
            print(header)
            print("  " + "-" * (len(header) - 2))
            for m in models[:max_rows]:
                if modeldir is None:
                    found_text = ''
                elif m['files']:
                    found_text = "  " + ", ".join([str(m['files'][0].relative_to(modeldir))] + [p.suffix for p in m['files'][1:]])
                else:
                    found_text = "  -"
                print(f"  {m['teff']:>7.0f} {m['logg']:>5.2f} {m['met']:>7.3f} {_format_value(m['alpha'], '6.2f'):>6} {_format_value(m['carbon'], '6.2f'):>6}  {m['model']:<{width}}{found_text}".rstrip())
            if len(models) > max_rows:
                print(f"  ... and {len(models) - max_rows} more model(s) (see the returned list)")


        # ##### FILES #####
        if output_file is not None:
            print(f"\n  Model names written to : {output_file}")
        if output_dir is not None:
            print(f"  {len(copied_models)} complete model(s) ({len(copied_files)} file(s)) copied to : {output_dir}")


    # ##########  RESULTS  ##########
    return models




# ====================================================================
# ============================== SPAMMS ==============================
# ====================================================================

# ==============================
def needed_models(
                  spamms_inputs : str | Path | list[str | Path],
                  spamms_dir    : str | Path,
                  grid          : str,
                  output_name   : str | Path = "spamms",
                  show          : bool = False,
                  verbose       : bool = True
    ) -> dict[str, list | dict | Path | Figure | None]:
    """
    Find the grid entries needed by SPAMMS, match them with the model atmospheres of the grid and plot them on the grid.

    ----------
    Parameters
    ----------
        spamms_inputs : str | Path | list[str | Path]
            SPAMMS input files: a directory (every '*.txt' file in it is
            used, in alphabetical order), a single file, or a list of files
            (a tuple is also accepted). The file names must be unique.

        spamms_dir : str | Path
            Directory where the SPAMMS scripts ('spamms.py',
            'spamms_disk.py', 'spamms_noreflec.py') are located; each input
            file is run with the script that corresponds to its name (see
            Notes). SPAMMS is run from this directory.

        grid : str
            Model atmosphere grid (same code as '--mode' in 'prismas.py';
            case-insensitive): 'KC', 'KM' or 'T'. See Notes.

        output_name : str | Path, optional
            Path and base name of the two output files,
            '<output_name>_<grid>_tocompute.dat' (model atmospheres to
            compute with PRISMAS) and '<output_name>_<grid>_missing.dat'
            (needed grid entries, with their model, '-' if outside the grid
            or 'dummy' for dummy entries), e.g. "results/HD_91316" ->
            'results/HD_91316_T_tocompute.dat'. A file extension is removed
            ("results/HD_91316.dat" gives the same names). Their directory
            must exist and existing files are overwritten. See Notes.
            Default is "spamms" (current directory).

        show : bool, optional
            If True, display the figure with plt.show(). The figure is
            returned in any case. Default is False.

        verbose : bool, optional
            If True, print the progress of the SPAMMS runs, the SPAMMS
            warnings and failures, and a summary; if False, the SPAMMS
            warnings and failures are issued as Python warnings. Default is
            True.

    -------
    Returns
    -------
        results_dict : dict[str, list | dict | Path | Figure | None]
            A dictionary containing the following key-value pairs:
                - 'needed' : list[str]
                    Alphabetically sorted grid entries needed by all the
                    SPAMMS input files, without duplicates (e.g.
                    'T30000_G4.00_Z1.000_V10').

                - 'per_file' : dict[str, list[str]]
                    A dictionary keyed by SPAMMS input file name; each value
                    is the list of grid entries needed by that file ([] if
                    none is missing). Files whose run failed are not
                    included.

                - 'models' : list[dict]
                    Needed grid entries whose model is in the grid (they can
                    be computed with PRISMAS).
                    A list of dictionaries, each containing the following
                    key-value pairs:
                        - 'name' : str
                            SPAMMS grid entry.

                        - 'model' : str
                            Model atmosphere, without file extension (e.g.
                            'G30000g400v10').

                        - 'subgrid' : str
                            'OSTAR2002' or 'BSTAR2006' for 'T'; '' for
                            Kurucz.

                        - 'teff' : float
                            Effective temperature in K.

                        - 'logg' : float
                            Surface gravity, log10(g / [cm s^-2]).

                        - 'met' : float
                            [M/H] in dex for 'KC' / 'KM'; Z / Zsun for 'T'.

                        - 'alpha' : float
                            Alpha enhancement [alpha/M] in dex (NaN for
                            'T').

                        - 'carbon' : float
                            Carbon abundance [C/M] in dex (NaN except for
                            'KM').

                        - 'vmic' : float
                            Microturbulence of the synthesis in km/s.

                - 'not_in_grid' : list[str]
                    Entries without model atmosphere in the catalog of the
                    grid (outside the grid).

                - 'dummy' : list[str]
                    Dummy entries (disk patches of SPAMMS), neither plotted
                    nor matched, but listed in 'missing_file'.

                - 'invalid' : list[str]
                    Entries that do not follow the SPAMMS naming of the
                    grid.

                - 'tocompute_file' : Path
                    File with the model atmospheres to compute with PRISMAS
                    ('<output_name>_<grid>_tocompute.dat').

                - 'missing_file' : Path
                    File with the needed grid entries and their models, '-'
                    for the entries outside the grid and 'dummy' for the
                    dummy entries ('<output_name>_<grid>_missing.dat').

                - 'fig' : Figure | None
                    Figure of the needed entries on the grid (None if no
                    entry is needed), to be saved e.g. with 'savefig_pkl'.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If 'grid' is not valid, 'output_name' is empty, no SPAMMS input
            file is found, or two input files have the same name.

        FileNotFoundError
            If an input file, 'spamms_dir', the SPAMMS script needed by an
            input file, the directory of the output files, or the catalog of
            the grid does not exist.

        NotADirectoryError
            If 'spamms_dir', or the parent of the output files, exists but
            is not a directory.

        IsADirectoryError
            If an element of a 'spamms_inputs' list or an output file is a
            directory.

        PermissionError
            If an input file, a SPAMMS script or 'spamms_dir' cannot be
            read.

        RuntimeError
            If SPAMMS fails for every input file (the message gives the
            reason of each failure).

    -----
    Notes
    -----
        1. WORKFLOW
            - The function runs the two steps of a SPAMMS grid check:
                >>> 1. SPAMMS check : needed grid entries and their models ('needed', 'models', 'tocompute_file', 'missing_file').
                >>> 2. Figure       : needed entries on the grid nodes ('fig').

            - Example:
                >>> res = needed_models("SPAMMS_inputs/", "/home/user/SPAMMS/", "T")
                >>> savefig_pkl(res['fig'], "Gridneeded")

        2. OUTPUT FILES
            - The two files are named after 'output_name' and the grid:
                >>> '<output_name>_<grid>_tocompute.dat' : model atmospheres to compute with PRISMAS.
                >>> '<output_name>_<grid>_missing.dat'   : needed grid entries and their models.

            - '<output_name>_<grid>_tocompute.dat': the model atmospheres
              (without extension) of the needed entries inside the grid, one
              per line, alphabetically sorted and without duplicates or
              header (e.g. the input of 'find_models' or of an HTCondor
              submission).

            - '<output_name>_<grid>_missing.dat': every needed entry inside
              or outside the model atmosphere grid, including the dummy
              entries, sorted by increasing Teff (then log g, composition and
              microturbulence), one per line with two columns:
              SPAMMS grid entry and model atmosphere (without extension), '-'
              for the entries outside the grid or 'dummy' for the dummy
              entries. Comment lines start with '#', so that the first column
              can be read as a list of names.

            - Invalid entries are in neither file (see 'invalid'), and dummy
              entries are not in the '_tocompute.dat' file. Both files are
              always written (empty, or with only the header, if no entry is
              needed).

            - The last file extension of 'output_name' is removed (e.g.
              'HD_91316.dat' -> 'HD_91316', 'HD_91316.txt' -> 'HD_91316'), so
              that the names always end in '_<grid>_tocompute.dat' and
              '_<grid>_missing.dat'.

        3. GRIDS AND ENTRY NAMES
            - The SPAMMS grid entries of each grid are:
                >>> 'T'  : 'T{teff}_G{logg}_Z{Z/Zsun}_V{vmic}'
                >>> 'KC' : 'T{teff}_G{logg}_M[a]{±[M/H]}_V{vmic}' or 'T{teff}_G{logg}_M{±[M/H]}_a{±alpha}_V{vmic}'
                >>> 'KM' : 'T{teff}_G{logg}_M{±[M/H]}_a{±alpha}_c{±carbon}_V{vmic}'

            - Signed values are written 'm' (-) / 'p' (+). For 'KC', the
              SPAMMS naming 'M[a]' means [alpha/M] = +0.4 with 'a' and 0
              without it.

        4. SPAMMS
            - The SPAMMS script is chosen from the name of each input file
              (without extension, case-insensitive): 'spamms_disk.py' if it
              ends in '_disk', 'spamms_noreflec.py' if it ends in 'noreflec',
              and 'spamms.py' otherwise, all in 'spamms_dir'. Only the
              scripts needed by the input files must exist.

            - Each input file is run with the Python interpreter of the
              current session as 'python <script> -i <input file> -c -n
              <n_cores>' from 'spamms_dir', with every core allocated to the
              process (all the cores of the computer, or only those assigned
              to the job by HTCondor / Slurm); with '-c', SPAMMS only checks
              the grid and prints the missing grid points. If the input file
              is not in 'spamms_dir', it is copied there for the run and
              always removed afterwards, also if SPAMMS fails or the run is
              interrupted; an input file already there with the same content
              is used and kept, and a different file with the same name is
              never overwritten (that run fails).

            - The progress is printed as one line per input file (with the
              number of missing grid points and of SPAMMS warnings). After
              all the runs, the warnings printed by SPAMMS and the reason of
              each failed run (exit status and last error line of SPAMMS) are
              printed as aligned blocks with 'verbose', and issued as one
              Python warning per input file otherwise. Failed input files are
              skipped.

            - If SPAMMS fails for every input file, the RuntimeError lists
              the reason of each failure (whatever 'verbose').

            - Entries with Teff = 0 or 1000 K are dummy entries (disk
              patches), neither plotted nor matched, but listed in the
              '_missing.dat' file with model 'dummy'.

        5. MODEL ATMOSPHERES
            - The models are found with 'find_models' (same catalogs and
              decoding), matching Teff, log g and the composition ([M/H],
              [alpha/M] and, for 'KM', [C/M]; Z / Zsun for 'T') of each
              entry.

            - TLUSTY: OSTAR2002 is used for Teff >= 27500 K (recommended in
              the overlap of both subgrids, 27500-30000 K) and BSTAR2006
              otherwise, falling back to the other subgrid when the model
              only exists there. The exception is 28000 K, taken from
              BSTAR2006: a Teff between 27000 and 27500 K is interpolated
              between the BSTAR2006 models at 27000 and 28000 K, so that the
              two subgrids are never mixed. From 27500 K the OSTAR2002 nodes
              are used (e.g. 28500 K between 27500 and 30000 K), so SPAMMS
              never requests 29000 K, which is not in the selection catalog.

        6. FIGURE
            - Grid nodes and needed entries (black crosses in the grid, red
              circled crosses outside it, both larger than the nodes), with
              reversed axes (hot / low gravity to the left / top) and the
              legend above the axes. For 'T', the nodes are read from the full
              OSTAR2002 and BSTAR2006 listings ('T_OSTAR02_models.dat' and
              'T_BSTAR06_models.dat'), which also contain the BSTAR2006
              models at 29000 and 30000 K, the only ones missing from the
              catalog used for the selection; for Kurucz, from the catalog of
              the grid. The figure is returned open, so that it can be saved
              with 'savefig_pkl'.

    ----------
    References
    ----------
        - Abdul-Masih, M., Sana, H., Conroy, K. E. et al. 2020, A&A, 636, A59.
        - Galán-Diéguez, D., Abdul-Masih, M., Allende Prieto, C. et al. 2026,
          arXiv:2609.15331.
        - Lanz, T. & Hubeny, I. 2003, ApJS, 146, 417.
        - Lanz, T. & Hubeny, I. 2007, ApJS, 169, 83.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from pathlib import Path
    import warnings


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_bool
    from .support.prismas.prismas_support import _allocated_cores, _check_directory, _check_grid, _decode_spamms, _DUMMY_TEFF, _grid_registry, _match_key, _message_block, _PARAM_TOL, _plot_needed_models, _run_spamms, _short_list, _sort_key, _SPAMMS_SCRIPTS
    from .support.tex.tex_support import _check_input_file, _check_output_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    # ##### OPTION 1: LIST OF FILES #####
    if isinstance(spamms_inputs, (list, tuple)):
        if len(spamms_inputs) == 0:
            raise ValueError("SPAMMS inputs 'spamms_inputs' must contain at least one input file. Got an empty list.")
        input_files = [_check_input_file(f, "spamms_inputs", "SPAMMS input file", "SPAMMS input file") for f in spamms_inputs]


    # ##### OPTION 2: DIRECTORY #####
    elif isinstance(spamms_inputs, (str, os.PathLike)) and Path(spamms_inputs).is_dir():
        input_files = [_check_input_file(f, "spamms_inputs", "SPAMMS input file", "SPAMMS input file") for f in sorted(Path(spamms_inputs).glob("*.txt"))]
        if len(input_files) == 0:
            raise ValueError(f"SPAMMS inputs 'spamms_inputs' must contain at least one input file ('*.txt'). Got none in '{spamms_inputs}'.")


    # ##### OPTION 3: SINGLE FILE #####
    elif isinstance(spamms_inputs, (str, os.PathLike)):
        input_files = [_check_input_file(spamms_inputs, "spamms_inputs", "SPAMMS input file", "SPAMMS input file or directory")]
    else:
        raise TypeError(f"SPAMMS inputs 'spamms_inputs' must be a path (str or pathlib.Path) to an input file or directory, or a list of paths. Got {type(spamms_inputs).__name__}.")
    #
    file_names = [f.name for f in input_files]
    repeated   = sorted({n for n in file_names if file_names.count(n) > 1})
    if repeated:
        raise ValueError(f"SPAMMS inputs 'spamms_inputs' must have unique file names (the results are keyed by file name). Got repeated names {repeated}.")
    #
    spamms_dir = _check_directory(spamms_dir, "spamms_dir", "SPAMMS directory", "SPAMMS directory")
    #
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]
    #
    if not isinstance(output_name, (str, os.PathLike)):
        raise TypeError(f"Output base name 'output_name' must be a path (str or pathlib.Path) without the '_<grid>_tocompute.dat' ending. Got {type(output_name).__name__}.")
    elif not os.fspath(output_name).strip() or not Path(output_name).name or os.fspath(output_name).endswith(("/", os.sep)):
        raise ValueError(f"Output base name 'output_name' must contain a file name (e.g. 'results/HD_91316'). Got '{output_name}'.")
    output_base    = Path(output_name).with_suffix('')  # 'HD_91316.dat' -> 'HD_91316'
    tocompute_file = _check_output_file(Path(f"{output_base}_{grid}_tocompute.dat"))
    missing_file   = _check_output_file(Path(f"{output_base}_{grid}_missing.dat"))
    #
    show = _check_bool(show, "show", "Display flag")
    #
    verbose = _check_bool(verbose, "verbose", "Progress and summary flag")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'needed'         : [],
                    'per_file'       : {},
                    'models'         : [],
                    'not_in_grid'    : [],
                    'dummy'          : [],
                    'invalid'        : [],
                    'tocompute_file' : None,
                    'missing_file'   : None,
                    'fig'            : None,
                   }


    # ##########  SPAMMS SCRIPTS  ##########
    run_scripts = {}  # Input file -> SPAMMS script
    for input_file in input_files:
        ending                  = next(e for e in _SPAMMS_SCRIPTS if input_file.stem.lower().endswith(e))  # '' matches any name
        run_scripts[input_file] = _check_input_file(spamms_dir / _SPAMMS_SCRIPTS[ending], "spamms_dir", "SPAMMS script", "Python script")


    # ####################  COMPUTATION  ####################
    # ##########  SPAMMS RUNS  ##########
    if verbose:
        print(f"SPAMMS grid check ('-c', {_allocated_cores()} cores) of {len(input_files)} input file(s)")
    name_width      = max(len(f.name) for f in input_files)  # Columns of the progress lines
    script_width    = max(len(p.name) for p in run_scripts.values())
    count_width     = len(str(len(input_files)))
    needed          = set()
    errors          = {}  # Input file name -> reason of the failure
    spamms_warnings = {}  # Input file name -> warnings printed by SPAMMS
    for n, input_file in enumerate(input_files, start=1):
        run = _run_spamms(input_file, run_scripts[input_file], spamms_dir)
        if run['warnings']:
            spamms_warnings[input_file.name] = run['warnings']
        if run['entries'] is None:
            errors[input_file.name] = run['error']
            status                  = "FAILED"
        else:
            results_dict['per_file'][input_file.name] = run['entries']
            needed.update(run['entries'])
            status = f"{len(run['entries'])} missing grid point(s)"
        if verbose:
            print(f"  [{n:>{count_width}}/{len(input_files)}]  {input_file.name:<{name_width}}  {run_scripts[input_file].name:<{script_width}}  {status:<26}  {len(run['warnings'])} SPAMMS warning(s)")


    # ##### SPAMMS WARNINGS AND FAILURES #####
    # Reported after the runs, so that they do not interrupt the progress lines: printed as aligned blocks with 'verbose'
    # (Python warnings may be filtered or shown only once per session), and issued as one Python warning per input file otherwise
    if len(errors) == len(input_files):  # Nothing to analyze: the reasons go in the exception, whatever 'verbose'
        if verbose and spamms_warnings:
            print("\n" + _message_block(spamms_warnings, "SPAMMS warnings"))
        raise RuntimeError(f"SPAMMS grid check failed for every input file ({len(input_files)}).\n" + _message_block({name: [error] for name, error in errors.items()}, "Reason of each failure"))
    elif verbose:
        if spamms_warnings:
            print("\n" + _message_block(spamms_warnings, "SPAMMS warnings"))
        if errors:
            print("\n" + _message_block({name: [error] for name, error in errors.items()}, "SPAMMS failures, input files skipped"))
    else:
        for file_name, messages in spamms_warnings.items():
            extra = f" (and {len(messages) - 1} more)" if len(messages) > 1 else ''
            warnings.warn(f"SPAMMS issued {len(messages)} warning(s) for '{file_name}': {messages[0]}{extra}", stacklevel=2)
        for file_name, error in errors.items():
            warnings.warn(f"SPAMMS grid check failed for '{file_name}' (input file skipped): {error}.", stacklevel=2)
    needed = sorted(needed)


    # ##########  GRID ENTRIES  ##########
    decoded = {}  # Grid entry -> stellar parameters and microturbulence
    for name in needed:
        params = _decode_spamms(name, grid)
        if params is None:
            results_dict['invalid'].append(name)
        elif any(abs(params['teff'] - t) < _PARAM_TOL for t in _DUMMY_TEFF):
            results_dict['dummy'].append(name)
        else:
            decoded[name] = params
    if results_dict['invalid']:
        warnings.warn(f"{len(results_dict['invalid'])} grid entr{'y does' if len(results_dict['invalid']) == 1 else 'ies do'} not follow the SPAMMS naming of the '{grid}' grid (e.g. '{results_dict['invalid'][0]}').", stacklevel=2)


    # ##########  MODEL ATMOSPHERES  ##########
    if decoded:
        criteria   = {p: sorted({params[p] for params in decoded.values()}) for p in grid_info['params']}  # Discrete values of every parameter
        candidates = find_models(grid, **criteria, verbose=False)
        index      = {}  # Matching key -> models of the grid (one per subgrid)
        for m in candidates:
            index.setdefault(_match_key(m, grid_info['params']), []).append(m)


        # ##### MATCHING #####
        for name, params in decoded.items():
            options = index.get(_match_key(params, grid_info['params']), [])
            if grid_info['subgrid_order'] is not None:  # Preferred subgrid first
                options = [m for subgrid in grid_info['subgrid_order'](params['teff']) for m in options if m['subgrid'] == subgrid]
            if not options:
                results_dict['not_in_grid'].append(name)
                continue
            chosen = options[0]
            results_dict['models'].append({
                                           'name'    : name,
                                           'model'   : chosen['model'],
                                           'subgrid' : chosen['subgrid'],
                                           'teff'    : chosen['teff'],
                                           'logg'    : chosen['logg'],
                                           'met'     : chosen['met'],
                                           'alpha'   : chosen['alpha'],
                                           'carbon'  : chosen['carbon'],
                                           'vmic'    : params['vmic'],
                                          })


    # ####################  PLOTTING  ####################
    if needed:
        outside             = set(results_dict['not_in_grid'])
        plot_teff           = np.array([params['teff'] for params in decoded.values()], dtype=float)
        plot_logg           = np.array([params['logg'] for params in decoded.values()], dtype=float)
        plot_in_grid        = np.array([name not in outside for name in decoded], dtype=bool)
        results_dict['fig'] = _plot_needed_models(plot_teff, plot_logg, plot_in_grid, grid)
        if show:
            plt.show()


    # ####################  OUTPUT  ####################
    # ##########  OUTPUT FILES  ##########
    with open(tocompute_file, "w", encoding="utf-8") as f:  # Model atmospheres, without header
        f.write(''.join(model + "\n" for model in sorted({m['model'] for m in results_dict['models']})))
    #
    entry_models = {m['name']: m['model'] for m in results_dict['models']}  # Grid entry -> model atmosphere (inside the grid)
    entry_models.update({name: "-" for name in results_dict['not_in_grid']})  # '-': no model atmosphere in the grid
    entry_models.update({name: "dummy" for name in results_dict['dummy']})  # Dummy entries (disk patches): listed, but not computed
    entry_width  = max([len(name) for name in entry_models] + [len("# Grid entry")])
    entry_keys   = {name: (_sort_key(_decode_spamms(name, grid), ('teff', 'logg', 'met', 'alpha', 'carbon', 'vmic')), name) for name in entry_models}  # NaN values last
    entry_order  = sorted(entry_models, key=entry_keys.get)  # Increasing Teff, then log g, composition and vmic (not alphabetical: 'T3500' < 'T10000')
    with open(missing_file, "w", encoding="utf-8") as f:
        f.write(f"# Needed grid entries of the {grid_info['name']} grid: {len(entry_models)} ({len(results_dict['models'])} inside the grid, {len(results_dict['not_in_grid'])} outside it, model '-'; {len(results_dict['dummy'])} dummy, model 'dummy')\n")
        f.write(f"{'# Grid entry':<{entry_width}}  Model\n")
        f.write(''.join(f"{name:<{entry_width}}  {entry_models[name]}\n" for name in entry_order))


    # ##########  SUMMARY  ##########
    if verbose:
        rows = [
                ("Grid entries needed", len(needed), []),
                ("  In the grid", len(results_dict['models']), []),
                ("  Not in the grid", len(results_dict['not_in_grid']), results_dict['not_in_grid']),
                ("  Dummy", len(results_dict['dummy']), []),
                ("  Invalid", len(results_dict['invalid']), results_dict['invalid']),
               ]
        if errors:
            rows.append(("SPAMMS failures", len(errors), list(errors)))
        label_width = max(len(label) for label, _, _ in rows)
        value_width = max(len(str(value)) for _, value, _ in rows)
        print()
        for label, value, names in rows:
            print(f"{label:<{label_width}} : {value:>{value_width}}" + (f"  ({_short_list(names)})" if names else ''))


        # ##### FILES #####
        files      = (("Models to compute written to", tocompute_file), ("Needed entries written to", missing_file))
        file_width = max(len(label) for label, _ in files)
        print()
        for label, path in files:
            print(f"{label:<{file_width}} : {path}")


    # ##########  RESULTS  ##########
    results_dict.update({
                         'needed'         : needed,
                         'tocompute_file' : tocompute_file,
                         'missing_file'   : missing_file,
                        })
    return results_dict




# ======================================================================
# ============================== ARCHIVES ==============================
# ======================================================================

# ==============================
def rename_tars(
                tardir         : str | Path,
                grid           : str,
                vmic           : int | float,
                recursive      : bool = False,
                output_dir     : str | Path | None = None,
                rename_content : bool = True,
                keep_original  : bool = False,
                overwrite      : bool = False,
                dry_run        : bool = False,
                ncores         : int | None = None
    ) -> list[dict]:
    """
    Rename the '.tar.gz' archives produced by PRISMAS in HPC mode (--opt 2) to the names that 'prismas.py' produces in modes 0 / 1.

    ----------
    Parameters
    ----------
        tardir : str | Path
            Directory where the '.tar.gz' archives of HPC mode are stored
            (the '--numpydir' used in 'prismas.py').

        grid : str
            Model atmosphere grid (same code as '--mode' in 'prismas.py';
            case-insensitive): 'KC', 'KM' or 'T'.

        vmic : int | float
            Microturbulence of the synthesis in km/s (the single value passed
            to '--vmic' in 'prismas.py'), >= 0. Only its integer part is
            used in the names, as in 'prismas.py'.

        recursive : bool, optional
            If True, also search the subdirectories of 'tardir'. Default is
            False.

        output_dir : str | Path | None, optional
            Existing directory, different from 'tardir', where the renamed
            archives are written. Default is None (each archive is renamed
            in its own directory).

        rename_content : bool, optional
            If True, rewrite the archive so that the inner folder and '.npy'
            files are renamed too. If False, only the '.tar.gz' file is
            renamed (fast, but its content keeps the HPC names). Default is
            True.

        keep_original : bool, optional
            If True, keep the original HPC archive after renaming. If False,
            it is removed once the renamed archive has been written. Default
            is False.

        overwrite : bool, optional
            If True, overwrite a renamed archive that already exists. If
            False, the archive is not renamed. Default is False.

        dry_run : bool, optional
            If True, only report the renaming that would be done, without
            touching any file. Default is False.

        ncores : int | None, optional
            Number of archives processed in parallel (threads), >= 1. See
            Notes. Default is None (self-tuned).

    -------
    Returns
    -------
        archives : list[dict]
            Archives of the grid found in 'tardir' (HPC names and names
            already of modes 0 / 1), sorted by (teff, logg, met, alpha,
            carbon). A summary table is also printed (sources relative to
            'tardir').
            A list of dictionaries, each containing the following key-value
            pairs:
                - 'source' : Path
                    Original archive.

                - 'target' : Path
                    Renamed archive (the source itself for 'done').

                - 'teff' : float
                    Effective temperature in K.

                - 'logg' : float
                    Surface gravity, log10(g / [cm s^-2]).

                - 'met' : float
                    [M/H] in dex for 'KC' / 'KM'; Z / Zsun for 'T'.

                - 'alpha' : float
                    Alpha enhancement [alpha/M] in dex (NaN for 'T').

                - 'carbon' : float
                    Carbon abundance [C/M] in dex (NaN except for 'KM').

                - 'status' : str
                    One of:
                        >>> 'renamed' : renamed (and rewritten, with 'rename_content').
                        >>> 'dry_run' : would be renamed ('dry_run' = True).
                        >>> 'done'    : already named as in modes 0 / 1.
                        >>> 'exists'  : target exists and 'overwrite' is False.
                        >>> 'skipped' : same target as another archive of the same subgrid.
                        >>> 'error'   : unreadable archive (e.g. truncated) or unwritable target.

                - 'reason' : str
                    Explanation for 'exists', 'skipped' and 'error', and for
                    the archives renamed with a subgrid tag (see Notes); ''
                    otherwise.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If 'grid' is not supported, 'vmic' is negative or not finite,
            'ncores' is smaller than 1, or 'output_dir' is the same
            directory as 'tardir'.

        FileNotFoundError
            If 'tardir' or 'output_dir' does not exist.

        NotADirectoryError
            If 'tardir' or 'output_dir' exists but is not a directory.

        PermissionError
            If 'tardir' cannot be read, or the directory of the renamed
            archives cannot be written (unless 'dry_run').

    -----
    Notes
    -----
        1. NAMING
            - In HPC mode, 'prismas.py' names each archive after its input
              model without extension, and the microturbulence is not part of
              the name:
                >>> 'KC' : 'm0.50t49000g4.5am05ak2odfnew.tar.gz'
                >>> 'KM' : 'amm03cm03om03t10000g20v20.tar.gz'
                >>> 'T'  : 'G30000g400v10.tar.gz' / 'BG20000g300v2.tar.gz'

            - In modes 0 / 1 the name is built from the stellar parameters,
              with signed values written as 'm' (-) / 'p' (+):
                >>> 'KC' : 'T{teff}_G{logg}_M{±[M/H]}_a{±alpha}_V{vmic}.tar.gz'
                >>> 'KM' : 'T{teff}_G{logg}_M{±[M/H]}_a{±alpha}_c{±carbon}_V{vmic}.tar.gz'
                >>> 'T'  : 'T{teff}_G{logg}_Z{Z/Zsun}_V{vmic}.tar.gz'
              e.g. 'T49000_G4.50_Mm0.500_ap0.400_V2.tar.gz'.

            - The HPC names are decoded with the same decoders as the model
              catalogs ('find_models'), and the names of modes 0 / 1 (the
              same format as the SPAMMS grid entries) with those of
              'needed_models'. Archives that follow neither naming of the
              grid are ignored.

            - TLUSTY archives of different subgrids with the same Teff, log g
              and Z / Zsun (e.g. 'G30000g300v10' and 'BG30000g300v2', in the
              overlap of OSTAR2002 and BSTAR2006) would get the same name of
              modes 0 / 1. Both are kept: a warning is issued and the subgrid
              tag is appended to each name, e.g.
              'T30000_G3.00_Z1.000_V10_OSTAR02.tar.gz' and
              'T30000_G3.00_Z1.000_V10_BSTAR06.tar.gz' ('reason' gives the
              other archive). Archives of the same subgrid with the same
              target (e.g. the same file in two subdirectories with
              'recursive') are 'skipped', keeping the first one.

        2. MICROTURBULENCE
            - The microturbulence cannot be recovered from the HPC names (in
              TLUSTY names, 'v' is the microturbulence of the model
              atmosphere, not of the synthesis), so it is a required input.

            - '--opt 2' admits a single microturbulence per run, so every
              archive in 'tardir' is assumed to share it.

        3. CONTENT
            - Each archive contains a folder and five '.npy' files carrying
              the HPC name ('{name}/{name}_wave.npy', ...). With
              'rename_content' = True they are renamed too, so the result is
              identical to a mode 0 / 1 output.

            - The rewrite is streamed member by member (nothing is extracted
              to disk), written to a temporary '.part' file and moved into
              place only when complete. The whole gzip stream is read, so its
              CRC is checked: a corrupted or truncated archive (e.g. a killed
              HTCondor job) is reported as 'error' and never replaces
              anything.

        4. PERFORMANCE
            - The gzip layer uses ISA-L (package 'isal', a required
              dependency), which compresses much faster than zlib, so disk /
              network I/O usually becomes the limiting factor.

            - The archives are processed in parallel threads (ISA-L and zlib
              release the GIL). With 'ncores' = None the number of threads
              tunes itself: it starts at 2 and doubles while the throughput
              improves by >= 10 %, up to the cores allocated to the process
              (CPU affinity, 'SLURM_CPUS_PER_TASK' / 'OMP_NUM_THREADS'). An
              integer fixes it (e.g. 1 for a serial run).

        5. RE-RUNS
            - Archives already named as in modes 0 / 1 (also with a subgrid
              tag) are reported as 'done', so the function can safely be run
              several times on the same directory. Temporary '.part' files
              are never listed.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from functools import partial
    import os
    from tqdm import tqdm
    import warnings


    # ##### PANORAMA #####
    from .support.obs.obs_support import _check_bool
    from .support.prismas.prismas_support import _archive_name, _check_directory, _check_grid, _check_ncores, _decode_model, _decode_spamms, _format_value, _grid_registry, _rename_archive, _run_parallel, _sort_key
    from .support.units.units_support import _check_real


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    tardir = _check_directory(tardir, "tardir", "Archive directory", "archive directory")
    #
    grid      = _check_grid(grid)
    grid_info = _grid_registry()[grid]
    #
    vmic = _check_real(vmic, "vmic", "Microturbulence", unit="in km/s", constraint="non-negative")
    if vmic != int(vmic):
        warnings.warn(f"Microturbulence 'vmic' = {vmic} km/s is not an integer; as in 'prismas.py', only its integer part is used in the names ('_V{int(vmic)}').", stacklevel=2)
    vmic = int(vmic)  # Integer part, as in 'prismas.py'
    #
    recursive = _check_bool(recursive, "recursive", "Recursive search flag")
    #
    output_dir = _check_directory(output_dir, "output_dir", "Output directory", "output directory", readable=False, allow_none=True)  # Write permission checked before renaming (not needed with 'dry_run')
    if output_dir is not None and output_dir.resolve() == tardir.resolve():
        raise ValueError(f"Output directory 'output_dir' must be different from 'tardir' (use None to rename in place). Got '{output_dir}'.")
    #
    rename_content = _check_bool(rename_content, "rename_content", "Content renaming flag")
    #
    keep_original = _check_bool(keep_original, "keep_original", "Original archive flag")
    #
    overwrite = _check_bool(overwrite, "overwrite", "Overwrite flag")
    #
    dry_run = _check_bool(dry_run, "dry_run", "Dry run flag")
    #
    ncores = _check_ncores(ncores)


    # ####################  COMPUTATION  ####################
    # ##########  FILES  ##########
    paths = sorted(tardir.rglob("*.tar.gz") if recursive else tardir.glob("*.tar.gz"))  # Temporary '.part' files do not match


    # ##########  DECODING  ##########
    archives = []
    subgrids = {}  # Source -> subgrid of its model ('' for Kurucz and for the archives already renamed)
    tags     = tuple(f"_{tag}" for tag in (grid_info['subgrid_tags'] or {}).values())  # E.g. '_OSTAR02'
    for path in paths:
        stem      = path.name[:-len(".tar.gz")]
        base_stem = next((stem[:-len(tag)] for tag in tags if stem.endswith(tag)), stem)  # Name without subgrid tag


        # ##### OPTION 1: MODES 0 / 1 NAME #####
        params = _decode_spamms(base_stem, grid)  # Same format as the SPAMMS grid entries
        if params is not None:
            target = path
            status = "done"


        # ##### OPTION 2: HPC NAME #####
        else:
            params = _decode_model(stem, grid)
            if params is None:  # Not an archive of the grid
                continue
            target         = (output_dir if output_dir is not None else path.parent) / f"{_archive_name(params, grid, vmic)}.tar.gz"
            status         = ''  # Still to be processed
            subgrids[path] = params['subgrid']
        archives.append({
                         'source' : path,
                         'target' : target,
                         'teff'   : float(params['teff']),
                         'logg'   : float(params['logg']),
                         'met'    : float(params['met']),
                         'alpha'  : float(params['alpha']),
                         'carbon' : float(params['carbon']),
                         'status' : status,
                         'reason' : '',
                        })


    # ##########  DUPLICATED TARGETS  ##########
    # ##### DIFFERENT SUBGRIDS #####
    # Models of different subgrids with the same parameters (e.g. 'G30000g300v10' and 'BG30000g300v2') would get the same
    # name: both are kept, with the subgrid tag appended to their names, so that no model replaces the other
    tagged = []
    if grid_info['subgrid_tags'] is not None:
        groups = {}  # Target -> archives still to be processed with that target
        for a in archives:
            if a['status'] == '':
                groups.setdefault(a['target'], []).append(a)
        for target, group in groups.items():
            if len({subgrids[a['source']] for a in group}) < 2:  # One subgrid: handled below
                continue
            for a in group:
                tag    = grid_info['subgrid_tags'][subgrids[a['source']]]
                others = ", ".join(str(b['source'].relative_to(tardir)) for b in group if b is not a)
                a.update({
                          'target' : target.with_name(f"{target.name[:-len('.tar.gz')]}_{tag}.tar.gz"),
                          'reason' : f"tag '_{tag}' added: same name as {others}",
                         })
                tagged.append(a)
    if tagged:
        warnings.warn(f"{len(tagged)} archive(s) of different subgrids of the {grid_info['name']} grid share the same stellar parameters and would get the same name; the subgrid tag is appended to their names (e.g. '{tagged[0]['target'].name}').", stacklevel=2)


    # ##### SAME SUBGRID #####
    seen = {}  # Target -> source of the first archive with that target
    for a in archives:
        if a['status'] == '':
            if a['target'] in seen:
                a.update({
                          'status' : "skipped",
                          'reason' : f"same target as {seen[a['target']].relative_to(tardir)}",
                         })
            else:
                seen[a['target']] = a['source']


    # ##########  RENAMING  ##########
    jobs = []
    for a in archives:
        if a['status'] != '':
            continue


        # ##### OPTION 1: EXISTING TARGET #####
        if a['target'].exists() and not overwrite:
            a.update({
                      'status' : "exists",
                      'reason' : "target already exists (use overwrite=True)",
                     })


        # ##### OPTION 2: DRY RUN #####
        elif dry_run:
            a['status'] = "dry_run"


        # ##### OPTION 3: TO RENAME #####
        else:
            if not os.access(a['target'].parent, os.W_OK):
                raise PermissionError(f"Permission denied to write to the directory of the renamed archives: '{a['target'].parent}'.")
            jobs.append(a)


    # ##### PARALLEL RUN #####
    n_workers = 0
    tuned     = False
    pbar      = tqdm(total=len(jobs), desc="Renaming archives", unit="archive", disable=not jobs)
    try:
        if jobs:
            n_workers, tuned = _run_parallel(jobs, partial(_rename_archive, rename_content=rename_content, keep_original=keep_original), ncores, callback=pbar.update)  # Expected errors are stored as 'error' status
    finally:
        pbar.close()


    # ##########  SORTING  ##########
    archives.sort(key=_sort_key)  # (teff, logg, met, alpha, carbon); NaN values last


    # ####################  OUTPUT  ####################
    # ##########  SUMMARY  ##########
    counts = {status: sum(a['status'] == status for a in archives) for status in ("renamed", "dry_run", "done", "exists", "skipped", "error")}
    print(f"\n{len(archives)} {grid} archive(s) found in {tardir} (vmic = {vmic} km/s)")


    # ##### TABLE #####
    if archives:
        header = f"  {'Teff':>7} {'logg':>5} {grid_info['met_label']:>7} {'alpha':>6} {'C':>6}  {'status':<8}  source -> target"
        print(header)
        print("  " + "-" * (len(header) + 30))
        for a in archives:
            reason = f"  [{a['reason']}]" if a['reason'] else ''
            print(f"  {a['teff']:>7.0f} {a['logg']:>5.2f} {a['met']:>7.3f} {_format_value(a['alpha'], '6.2f'):>6} {_format_value(a['carbon'], '6.2f'):>6}  {a['status']:<8}  {a['source'].relative_to(tardir)} -> {a['target'].name}{reason}")


    # ##### COUNTS #####
    print("\n  " + ", ".join(f"{n} {status}" for status, n in counts.items() if n))
    if jobs:
        print(f"  Processed with {n_workers} parallel worker(s)" + (" (self-tuned)." if tuned else "."))
    if dry_run:
        print("  Dry run: no file was modified.")
    elif output_dir is not None:
        print(f"  Renamed archives written to : {output_dir}")


    # ##########  RESULTS  ##########
    return archives
