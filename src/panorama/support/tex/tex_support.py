"""
Internal helpers of the 'tex' module: validation of the file arguments shared
by the public functions (also reused by 'spec' and 'obs'), reading and
normalization of the lists of object names used by the catalog queries
('queryviz', 'querysimbad'), and reading of the SIMBAD result fields.

Conventions shared by every function in this module:
    - Functions are private (leading underscore); they are used by the
      public functions of 'panorama.tex' (and by those of 'panorama.spec'
      and 'panorama.obs' for the file arguments) and are not exposed at the
      package level.
    - File arguments are paths (path + filename), given either as strings
      or as pathlib.Path objects, and are returned as pathlib.Path objects.
      Output directories are never created.
    - All functions validate their inputs before computing.
"""

from __future__ import annotations

# ===================================================================
# ============================== FILES ==============================
# ===================================================================

# ==============================
def _check_input_file(
                      file  : str | Path,
                      name  : str,
                      label : str,
                      kind  : str
    ) -> Path:
    """
    Validate an input file argument (path, existence, type and read permission) and convert it to a Path.

    ----------
    Parameters
    ----------
        file : str | Path
            Path and filename of the input file.

        name : str
            Name of the argument in the error messages (e.g. 'aux_file').

        label : str
            Role of the argument in the error messages (e.g. 'LaTeX auxiliary
            file').

        kind : str
            Kind of file in the error messages (e.g. '.aux file', 'pickle
            file').

    -------
    Returns
    -------
        file : Path
            The input file as a Path.

    ------
    Raises
    ------
        TypeError
            If 'file' is not a path, or 'name', 'label' or 'kind' are not
            strings.

        ValueError
            If 'name', 'label' or 'kind' are empty.

        FileNotFoundError
            If the file does not exist.

        IsADirectoryError
            If the path points to a directory.

        PermissionError
            If the file cannot be read.

    -----
    Notes
    -----
        - Shared by 'nobib', 'missingbib', 'dupbib', 'read_pkl' and
          'panorama.spec.readspec', so that the same argument gives the same
          error messages.
    """

    # ####################  MODULES & PACKAGES  ####################
    import os
    from pathlib import Path


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label), ("kind", kind)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")
    #
    if not isinstance(file, (str, os.PathLike)):
        raise TypeError(f"{label} '{name}' must be a path (str or pathlib.Path) to the {kind}. Got {type(file).__name__}.")
    file = Path(file)
    if not file.exists():
        raise FileNotFoundError(f"{label} not found: '{file}'.")
    elif not file.is_file():
        raise IsADirectoryError(f"Expected a {kind} but found a directory: '{file}'.")
    elif not os.access(file, os.R_OK):
        raise PermissionError(f"Permission denied to read {kind}: '{file}'.")


    # ####################  OUTPUT  ####################
    return file




# ==============================
def _check_output_file(
                       output_file : str | Path | None,
                       allow_none  : bool = False
    ) -> Path | None:
    """
    Validate an output file argument (path, existing parent directory) and convert it to a Path.

    ----------
    Parameters
    ----------
        output_file : str | Path | None
            Path and filename of the output file, or None if 'allow_none'.

        allow_none : bool, optional
            If True, None is accepted and returned unchanged (nothing is
            written). Default is False.

    -------
    Returns
    -------
        output_file : Path | None
            The output file as a Path, or None.

    ------
    Raises
    ------
        TypeError
            If 'output_file' is not a path (nor None when 'allow_none'), or
            'allow_none' is not a boolean.

        IsADirectoryError
            If 'output_file' is an existing directory.

        FileNotFoundError
            If the directory of 'output_file' does not exist.

        NotADirectoryError
            If the parent of 'output_file' exists but is not a directory.

    -----
    Notes
    -----
        - Shared by 'queryviz', 'querysimbad', 'panorama.spec.writespec',
          'rmvspec', 'normspec' and 'panorama.obs.fits_catalog', so that the
          same argument gives the same error messages. Output directories are
          never created.
    """

    # ####################  MODULES & PACKAGES  ####################
    import os
    from pathlib import Path


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(allow_none, bool):
        raise TypeError(f"None flag 'allow_none' must be a boolean. Got {type(allow_none).__name__}.")
    #
    if output_file is None and allow_none:
        return None
    elif not isinstance(output_file, (str, os.PathLike)):
        expected = "a path (str or pathlib.Path) or None" if allow_none else "a path (str or pathlib.Path)"
        raise TypeError(f"Output file 'output_file' must be {expected}. Got {type(output_file).__name__}.")
    output_file = Path(output_file)
    if output_file.is_dir():
        raise IsADirectoryError(f"Output file 'output_file' is an existing directory, not a file: '{output_file}'.")
    elif not output_file.parent.exists():
        raise FileNotFoundError(f"Output directory does not exist: '{output_file.parent}'.")
    elif not output_file.parent.is_dir():
        raise NotADirectoryError(f"Output directory is not a directory: '{output_file.parent}'.")


    # ####################  OUTPUT  ####################
    return output_file




# =====================================================================
# ============================== OBJECTS ==============================
# =====================================================================

# ==============================
def _read_object_names(
                       input_file : str | Path | list[str]
    ) -> list[str]:
    """
    Validate a list of object names, or read it from a text file, following the PANORAMA convention for object lists.

    ----------
    Parameters
    ----------
        input_file : str | Path | list[str]
            Either a non-empty list of object names or the path and filename
            of a text file (txt, dat, ascii...) with the object names in the
            first column.

    -------
    Returns
    -------
        object_names : list[str]
            Object names, in the input order and as given by the user (not
            normalized; see '_normalize_name').

    ------
    Raises
    ------
        TypeError
            If 'input_file' is neither a list nor a path, or the list
            contains non-string elements.

        FileNotFoundError
            If 'input_file' is a path that does not exist.

        IsADirectoryError
            If 'input_file' is an existing directory.

        PermissionError
            If 'input_file' is a file that cannot be read.

        ValueError
            If the list is empty or contains empty names, or the file
            contains no object names.

    -----
    Notes
    -----
        - In a file, lines starting with '#' and empty lines are ignored,
          and only the first whitespace-separated column is read (names with
          spaces must use underscores, e.g. 'HD_91316').

        - A list given by the user is returned as a new list (the input is
          not modified).
    """

    # ####################  MODULES & PACKAGES  ####################
    import os
    from pathlib import Path


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if isinstance(input_file, (str, os.PathLike)):
        input_file = Path(input_file)
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: '{input_file}'.")
        elif not input_file.is_file():
            raise IsADirectoryError(f"Expected an input file but found a directory: '{input_file}'.")
        elif not os.access(input_file, os.R_OK):
            raise PermissionError(f"Permission denied to read input file: '{input_file}'.")
    elif isinstance(input_file, list):
        if len(input_file) == 0:
            raise ValueError("Object list 'input_file' must contain at least one object name. Got an empty list.")
        for name in input_file:
            if not isinstance(name, str):
                raise TypeError(f"Object list 'input_file' must contain only strings (object names). Got {type(name).__name__}.")
            elif not name.strip():
                raise ValueError(f"Object list 'input_file' must not contain empty object names. Got '{name}'.")
    else:
        raise TypeError(f"Object list 'input_file' must be either a list of object names or a path (str or pathlib.Path) to a file containing object names. Got {type(input_file).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  OPTION 1: LIST OF OBJECT NAMES  ##########
    if isinstance(input_file, list):
        object_names = list(input_file)


    # ##########  OPTION 2: FILE PATH  ##########
    else:
        object_names = []
        with open(input_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):  # Skip empty lines and comments
                    continue
                object_names.append(line.split()[0])  # First column only
        #
        if len(object_names) == 0:
            raise ValueError(f"Object list 'input_file' must contain at least one object name. Got none in '{input_file}'.")


    # ####################  OUTPUT  ####################
    return object_names




# ==============================
def _normalize_name(
                    name : str
    ) -> str:
    """
    Normalize an object name before querying a catalog (SIMBAD name resolver).

    ----------
    Parameters
    ----------
        name : str
            Object name as given by the user (e.g. 'HD_91316', 'BD+60_2522').

    -------
    Returns
    -------
        norm_name : str
            Name with underscores replaced by spaces, repeated whitespace
            collapsed and no surrounding whitespace (e.g. 'HD 91316',
            'BD+60 2522').

    ------
    Raises
    ------
        TypeError
            If 'name' is not a string.

        ValueError
            If 'name' is empty.

    -----
    Notes
    -----
        - Underscores stand for spaces in the PANORAMA object lists, so they
          are replaced by spaces, not removed: removing them would merge the
          words of multi-word names ('Cyg_OB2_12' -> 'CygOB212', which SIMBAD
          does not resolve).
    """

    # ####################  MODULES & PACKAGES  ####################
    import re


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(name, str):
        raise TypeError(f"Object name 'name' must be a string. Got {type(name).__name__}.")
    elif not name.strip():
        raise ValueError(f"Object name 'name' must be a non-empty string. Got '{name}'.")


    # ####################  COMPUTATION  ####################
    norm_name = name.replace("_", " ")  # Underscores stand for spaces ('HD_91316' → 'HD 91316')
    norm_name = re.sub(r"\s+", " ", norm_name)  # Collapse repeated whitespace


    # ####################  OUTPUT  ####################
    return norm_name.strip()




# ===================================================================
# ============================== QUERY ==============================
# ===================================================================

# ==============================
def _get_simbad_field(
                      row        : Row,
                      field      : str,
                      candidates : dict[str, tuple[str, ...]]
    ) -> object:
    """
    Read a field of a SIMBAD result row, trying the column names of the TAP-based and of the legacy astroquery APIs.

    ----------
    Parameters
    ----------
        row : Row
            astropy.table.Row returned by Simbad.query_object.

        field : str
            Key of 'candidates' (e.g. 'ra', 'sp_type').

        candidates : dict[str, tuple[str, ...]]
            A dictionary keyed by field; each value is the tuple of column
            names to try, in order (TAP name first, legacy name second).

    -------
    Returns
    -------
        value : object
            Value of the first column present in the row (bytes decoded to
            str), or None if no column is present or the value is masked.

    ------
    Raises
    ------
        TypeError
            If 'row' is not an astropy.table.Row, 'field' is not a string or
            'candidates' is not a dictionary of tuples of strings.

        ValueError
            If 'field' is not a key of 'candidates'.

    -----
    Notes
    -----
        - Used by 'panorama.tex.querysimbad' to support both astroquery
          >= 0.4.8 (TAP service, lowercase field names) and older versions
          (uppercase legacy names).
    """

    # ####################  MODULES & PACKAGES  ####################
    from astropy.table import Row
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(row, Row):
        raise TypeError(f"SIMBAD result row 'row' must be an astropy.table.Row. Got {type(row).__name__}.")
    #
    if not isinstance(candidates, dict):
        raise TypeError(f"SIMBAD column lookup 'candidates' must be a dictionary of tuples of column names. Got {type(candidates).__name__}.")
    for key, names in candidates.items():
        if not isinstance(names, tuple) or not all(isinstance(n, str) for n in names):
            raise TypeError(f"SIMBAD column lookup 'candidates' must map each field to a tuple of column names (strings). Got {names!r} for '{key}'.")
    #
    if not isinstance(field, str):
        raise TypeError(f"SIMBAD field 'field' must be a string. Got {type(field).__name__}.")
    elif field not in candidates:
        raise ValueError(f"SIMBAD field 'field' must be a key of 'candidates' ({', '.join(candidates)}). Got '{field}'.")


    # ####################  COMPUTATION  ####################
    for col in candidates[field]:  # TAP name first, legacy name second
        if col in row.colnames:
            value = row[col]
            if np.ma.is_masked(value):  # Masked cell: SIMBAD has no value for this object
                return None
            return value.decode() if isinstance(value, bytes) else value  # Old astropy versions return bytes for string columns


    # ####################  OUTPUT  ####################
    return None  # None of the candidate columns is in the table




# ==============================
def _to_float(
              x : object
    ) -> float:
    """
    Convert a value read from a catalog to float, returning NaN when it cannot be converted.

    ----------
    Parameters
    ----------
        x : object
            Value to convert (number, numeric string or None).

    -------
    Returns
    -------
        value : float
            float(x), or NaN if 'x' is None or cannot be converted.

    -----
    Notes
    -----
        - Any object is a valid input (missing catalog values arrive as None),
          so the argument is not validated: conversion failures are the
          expected way of signaling missing data.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numpy as np


    # ####################  COMPUTATION  ####################
    try:
        value = float(x)
    except (TypeError, ValueError):  # None, empty or non-numeric strings
        value = np.nan


    # ####################  OUTPUT  ####################
    return value
