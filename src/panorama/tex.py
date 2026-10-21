"""
Manuscript, catalog and figure utilities for the final stages of the
spectroscopic analysis of (massive) stars.

Conventions shared by every function in this module:
    - File arguments are paths (path + filename), given either as strings
      or as pathlib.Path objects.
    - Bibliography checks work on BibTeX and biblatex documents.
    - Figures are saved following the PANORAMA conventions and stored
      as pickle files so that they can be edited later.
    - Catalog queries (VizieR, SIMBAD) require internet access.
    - All functions validate their inputs before computing.
"""

from __future__ import annotations

# =================================================================
# ============================== BIB ==============================
# =================================================================

# ==============================
def nobib(
          aux_file : str | Path,
          bib_file : str | Path
    ) -> list[str]:
    """
    Identify the bibliography entries of a .bib file that are not cited in a LaTeX .aux file.

    ----------
    Parameters
    ----------
        aux_file : str | Path
            Path and filename of the LaTeX .aux file (produced during
            compilation).

        bib_file : str | Path
            Path and filename of the bibliography .bib file.

    -------
    Returns
    -------
        unused_keys : list[str]
            Alphabetically sorted list of the keys defined in the .bib file
            that are never cited in the .aux file.

    ------
    Raises
    ------
        TypeError
            If any input is not a string or a pathlib.Path.

        FileNotFoundError
            If any file does not exist.

        IsADirectoryError
            If any path points to a directory.

        PermissionError
            If any file cannot be read.

    -----
    Notes
    -----
        - Supports both BibTeX (\\citation{...}) and biblatex
          (\\abx@aux@cite{...}) citation formats.

        - @string, @comment and @preamble blocks of the .bib file are not
          bibliography entries and are ignored.

        - If the document uses \\nocite{*}, every entry is formally cited;
          a warning is issued and the comparison ignores the wildcard.

        - Run it after a full compilation (LaTeX + BibTeX/Biber + LaTeX), so
          that the .aux file contains every citation. Typical use: cleaning an
          ADS-exported .bib file before submitting a paper.

        - See 'missingbib' for the opposite check (cited but undefined keys)
          and 'dupbib' for duplicated entries.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import re
    import warnings


    # ##### PANORAMA #####
    from .support.tex.tex_support import _check_input_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    aux_file = _check_input_file(aux_file, "aux_file", "LaTeX auxiliary file", ".aux file")
    #
    bib_file = _check_input_file(bib_file, "bib_file", "Bibliography file", ".bib file")


    # ####################  COMPUTATION  ####################
    # ##########  AUX FILE  ##########
    # ##### READ #####
    with open(aux_file, "r", encoding="utf-8") as f:
        aux_content = f.read()  # Whole file as a single string


    # ##### EXTRACT CITATIONS #####
    cited_keys = set()  # A set: a key cited several times is counted once

    for match in re.findall(r"\\citation{([^}]+)}", aux_content):  # Standard BibTeX citations
        cited_keys.update(k.strip() for k in match.split(","))  # \cite{a,b} produces a single \citation{a,b}

    for match in re.findall(r"\\abx@aux@cite(?:{[^}]*})?{([^}]+)}", aux_content):  # Citations of biblatex (with or without the refsection argument)
        cited_keys.update(k.strip() for k in match.split(","))


    # ##### CHECKS #####
    if "*" in cited_keys:  # \nocite{*} cites the whole .bib file
        warnings.warn(
                      "The document uses \\nocite{*}: every .bib entry is formally cited. "
                      "The wildcard is ignored and only explicit citations are compared."
                     )
        cited_keys.discard("*")  # The wildcard is not a key


    # ##########  BIB FILE  ##########
    # ##### READ #####
    with open(bib_file, "r", encoding="utf-8") as f:
        bib_content = f.read()


    # ##### EXTRACT KEYS #####
    bib_keys = set(re.findall(r"@(?!string\b|comment\b|preamble\b)\w+\s*{\s*([^,\s]+)\s*,", bib_content, flags=re.IGNORECASE))  # Entry keys ('@type{key,'), skipping non-entry blocks


    # ####################  OUTPUT  ####################
    unused_keys = sorted(bib_keys - cited_keys)  # Defined in the .bib file but never cited
    return unused_keys




# ==============================
def missingbib(
               aux_file : str | Path,
               bib_file : str | Path
    ) -> list[str]:
    """
    Identify the citation keys of a LaTeX .aux file that are not defined in a .bib file.

    ----------
    Parameters
    ----------
        aux_file : str | Path
            Path and filename of the LaTeX .aux file (produced during
            compilation).

        bib_file : str | Path
            Path and filename of the bibliography .bib file.

    -------
    Returns
    -------
        missing_keys : list[str]
            Alphabetically sorted list of the keys cited in the document but
            not defined in the .bib file.

    ------
    Raises
    ------
        TypeError
            If any input is not a string or a pathlib.Path.

        FileNotFoundError
            If any file does not exist.

        IsADirectoryError
            If any path points to a directory.

        PermissionError
            If any file cannot be read.

    -----
    Notes
    -----
        - Complementary to 'nobib'. Missing keys produce "Citation undefined"
          warnings and "?" marks in the compiled PDF.

        - Supports both BibTeX (\\citation{...}) and biblatex
          (\\abx@aux@cite{...}) citation formats. The \\nocite{*} wildcard
          is ignored.

        - BibTeX keys are compared case-sensitively, as LaTeX does when
          resolving citations.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import re


    # ##### PANORAMA #####
    from .support.tex.tex_support import _check_input_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    aux_file = _check_input_file(aux_file, "aux_file", "LaTeX auxiliary file", ".aux file")
    #
    bib_file = _check_input_file(bib_file, "bib_file", "Bibliography file", ".bib file")


    # ####################  COMPUTATION  ####################
    # ##########  AUX FILE  ##########
    # ##### READ #####
    with open(aux_file, "r", encoding="utf-8") as f:
        aux_content = f.read()


    # ##### EXTRACT CITATIONS #####
    cited_keys = set()

    for match in re.findall(r"\\citation{([^}]+)}", aux_content):  # Standard BibTeX citations
        cited_keys.update(k.strip() for k in match.split(","))

    for match in re.findall(r"\\abx@aux@cite(?:{[^}]*})?{([^}]+)}", aux_content):  # Citations of biblatex
        cited_keys.update(k.strip() for k in match.split(","))

    cited_keys.discard("*")  # \nocite{*} wildcard is not a key


    # ##########  BIB FILE  ##########
    # ##### READ #####
    with open(bib_file, "r", encoding="utf-8") as f:
        bib_content = f.read()


    # ##### EXTRACT KEYS #####
    bib_keys = set(re.findall(r"@(?!string\b|comment\b|preamble\b)\w+\s*{\s*([^,\s]+)\s*,", bib_content, flags=re.IGNORECASE))  # Entry keys, skipping non-entry blocks


    # ####################  OUTPUT  ####################
    missing_keys = sorted(cited_keys - bib_keys)  # Cited in the document but not defined in the .bib file
    return missing_keys




# ==============================
def dupbib(
           bib_file : str | Path
    ) -> dict[str, dict]:
    """
    Find duplicated entries in a .bib file, either by citation key or by DOI.

    ----------
    Parameters
    ----------
        bib_file : str | Path
            Path and filename of the bibliography .bib file.

    -------
    Returns
    -------
        results_dict : dict[str, dict]
            A dictionary containing the following key-value pairs:
                - 'keys' : dict[str, int]
                    A dictionary keyed by citation key (defined more than
                    once); each value is its number of occurrences. Keys
                    differing only in case are grouped (BibTeX treats them as
                    the same key), and the first spelling found is reported.
                    Empty dictionary if there are no duplicated keys.

                - 'dois' : dict[str, list[str]]
                    A dictionary keyed by normalized DOI (shared by more than
                    one entry); each value is the list of citation keys
                    pointing to it (same paper cited under different keys).
                    Empty dictionary if there are no duplicated DOIs.

    ------
    Raises
    ------
        TypeError
            If 'bib_file' is not a string or a pathlib.Path.

        FileNotFoundError
            If the file does not exist.

        IsADirectoryError
            If the path points to a directory.

        PermissionError
            If the file cannot be read.

    -----
    Notes
    -----
        - Duplicated keys make BibTeX stop with "Repeated entry" errors and
          biblatex silently keep only the first entry.

        - Duplicated DOIs are common when merging .bib files from different
          projects (e.g. the same paper exported from ADS as
          "2014A&A...564A..52L" and stored by hand as "Langer2014"). They do
          not break the compilation but produce duplicated references.

        - DOIs are compared case-insensitively, after removing any
          "https://doi.org/" prefix.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from collections import Counter
    import re


    # ##### PANORAMA #####
    from .support.tex.tex_support import _check_input_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    bib_file = _check_input_file(bib_file, "bib_file", "Bibliography file", ".bib file")


    # ##########  DICTIONARY  ##########
    results_dict = {
                    'keys' : {},
                    'dois' : {},
                   }


    # ####################  COMPUTATION  ####################
    # ##########  BIB FILE  ##########
    # ##### READ #####
    with open(bib_file, "r", encoding="utf-8") as f:
        bib_content = f.read()


    # ##### SPLIT ENTRIES #####
    entries = re.split(r"(?=@\w+\s*{)", bib_content)  # One chunk per '@type{' block


    # ##### EXTRACT KEYS & DOIS #####
    keys    = []
    doi_map = {}

    for entry in entries:
        key_match = re.match(r"@(?!string\b|comment\b|preamble\b)\w+\s*{\s*([^,\s]+)\s*,", entry, flags=re.IGNORECASE)
        if key_match is None:  # Not a bibliography entry (@string, @comment, text between entries...)
            continue
        key = key_match.group(1)
        keys.append(key)  # Every occurrence is kept, so that duplicates can be counted

        doi_match = re.search(r"\bdoi\s*=\s*[{\"]\s*([^}\"]+?)\s*[}\"]", entry, flags=re.IGNORECASE)  # doi = {...} or doi = "..."
        if doi_match is not None:
            doi = re.sub(r"^https?://(dx\.)?doi\.org/", '', doi_match.group(1).strip(), flags=re.IGNORECASE).lower()  # Normalized DOI
            doi_map.setdefault(doi, []).append(key)  # Normalized DOI -> keys of the entries that share it


    # ##########  DUPLICATES  ##########
    # ##### KEYS #####
    key_counts     = Counter(k.lower() for k in keys)  # BibTeX keys are case-insensitive for duplicates
    first_spelling = {}
    for k in keys:
        first_spelling.setdefault(k.lower(), k)  # Report the first spelling found
    dup_keys = {first_spelling[k]: n for k, n in key_counts.items() if n > 1}


    # ##### DOIS #####
    dup_dois = {doi: ks for doi, ks in doi_map.items() if len(ks) > 1}  # DOIs shared by more than one entry


    # ####################  OUTPUT  ####################
    results_dict.update({
                         'keys' : dup_keys,
                         'dois' : dup_dois,
                        })
    return results_dict




# ===================================================================
# ============================== QUERY ==============================
# ===================================================================

# ==============================
def queryviz(
             input_file  : str | Path | list[str],
             catalog     : str,
             output_file : str | Path = "vizier_output.txt",
             col_width   : int = 15
    ) -> None:
    """
    Query a VizieR catalog for a list of objects and save the results as a fixed-width text table.

    ----------
    Parameters
    ----------
        input_file : str | Path | list[str]
            Either a non-empty list of object names or the path and filename
            of a text file (txt, dat, ascii...) with the object names in the
            first column.

        catalog : str
            VizieR catalog identifier (e.g. 'I/239/hip_main' for Hipparcos,
            'I/355/gaiadr3' for Gaia DR3, 'B/mk/mktypes' for MK spectral
            types).

        output_file : str | Path, optional
            Path and filename of the output file. Default is
            "vizier_output.txt".

        col_width : int, optional
            Fixed width of each column of the output file, in characters.
            Must be positive. Default is 15.

    -------
    Returns
    -------
        None
            The results are written to 'output_file'.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        FileNotFoundError
            If 'input_file' is a path that does not exist, or the directory
            of 'output_file' does not exist.

        IsADirectoryError
            If 'input_file' or 'output_file' is an existing directory.

        NotADirectoryError
            If the parent of 'output_file' exists but is not a directory.

        PermissionError
            If 'input_file' is a file that cannot be read.

        ValueError
            If the list of object names is empty or contains empty names, the
            file contains no object names, 'catalog' is empty, or
            'col_width' is not positive.

        RuntimeError
            If every query fails or returns no data.

    -----
    Notes
    -----
        1. INPUT
            - In a file, lines starting with '#' and empty lines are ignored,
              and only the first whitespace-separated column is read (names
              with spaces must use underscores, e.g. 'HD_91316').

            - Names are normalized before querying: underscores are replaced
              by spaces and repeated spaces collapsed ('HD_91316' →
              'HD 91316', 'BD+60_2522' → 'BD+60 2522'), as expected by SIMBAD
              (used by VizieR to resolve names).

            - The list or file of names is read with '_read_object_names' and
              the names are normalized with '_normalize_name' (see
              'support/tex/tex_support.py'), shared with 'querysimbad'.

        2. QUERY
            - All catalog columns are requested (columns=["**"]) and only the
              first row returned for each object is kept (row_limit=1), i.e.
              the first match within the default VizieR search radius. In
              crowded fields check that the match is the intended star.

            - Failed queries issue a warning and are kept as empty rows.

        3. OUTPUT
            - The order of the input objects is preserved.

            - The header contains 'Object' followed by every catalog column.

            - Columns are left-aligned and padded with spaces to 'col_width'
              characters. Values longer than 'col_width' are not truncated
              and will merge with the next column: increase 'col_width' for
              catalogs with long fields (e.g. bibcodes, notes).

            - Objects without data have 'nan' in every catalog column.

    ----------
    References
    ----------
        - Ochsenbein, F., Bauer, P. & Marcout, J. 2000, A&AS, 143, 23.
        - Ginsburg, A., Sipőcz, B. M., Brasseur, C. E. et al. 2019,
          AJ, 157, 98.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astroquery.vizier import Vizier
    import numbers
    from tqdm import tqdm
    import warnings


    # ##### PANORAMA #####
    from .support.tex.tex_support import _check_output_file, _normalize_name, _read_object_names


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    object_names = _read_object_names(input_file)  # Validates the list or the file and reads the names
    #
    if not isinstance(catalog, str):
        raise TypeError(f"VizieR catalog identifier 'catalog' must be a string. Got {type(catalog).__name__}.")
    elif not catalog.strip():
        raise ValueError(f"VizieR catalog identifier 'catalog' must be a non-empty string (e.g. 'I/239/hip_main'). Got '{catalog}'.")
    #
    output_file = _check_output_file(output_file)
    #
    if not isinstance(col_width, numbers.Integral) or isinstance(col_width, bool):
        raise TypeError(f"Column width 'col_width' (in characters) must be an integer. Got {type(col_width).__name__}.")
    elif col_width <= 0:
        raise ValueError(f"Column width 'col_width' (in characters) must be positive. Got {col_width}.")


    # ####################  COMPUTATION  ####################
    # ##########  NORMALIZE NAMES FOR SEARCHING  ##########
    norm_names = [_normalize_name(n) for n in object_names]  # Names used for the queries; the original ones label the output rows


    # ##########  CONFIG  ##########
    viz     = Vizier(columns=["**"], row_limit=1)  # All catalog columns, first matching row only
    results = {}  # Original name -> first catalog row (None if not found or failed)
    columns = None  # Catalog column names, unknown until the first successful query


    # ##########  QUERY  ##########
    # ##### MAKE QUERIES #####
    for orig, norm in tqdm(
                           zip(object_names, norm_names),
                           total=len(object_names),
                           desc="Querying VizieR..."
                          ):
        try:
            res = viz.query_object(norm, catalog=catalog)  # TableList with one table per matching catalog
            if res:  # Empty TableList if the object has no entry in the catalog
                table = res[0]
                if columns is None:  # Column names taken from the first successful query
                    columns = table.colnames
                results[orig] = table[0]  # First matching row (row_limit=1)
            else:
                results[orig] = None
        except Exception as e:
            warnings.warn(f"VizieR query for '{orig}' failed ({e}).")  # Non-fatal: the loop goes on with the next object
            results[orig] = None


    # ##### CHECKS #####
    if columns is None:  # No query returned data, so there is no table to write
        raise RuntimeError("No successful VizieR queries.")


    # ####################  OUTPUT  ####################
    with open(output_file, "w", encoding="utf-8") as f:
        header = ''.join(f"{c:<{col_width}}" for c in ["Object"] + columns)  # Left-aligned, fixed-width columns
        f.write(header + "\n")

        for obj in object_names:  # Input order, not query completion order
            row = results[obj]
            if row is None:
                values = ["nan"] * len(columns)  # Keep the row so that the input order is preserved
            else:
                values = [str(row[c]) for c in columns]  # Masked values are written as '--' by astropy

            line = ''.join(f"{v:<{col_width}}" for v in [obj] + values)
            f.write(line + "\n")




# ==============================
def querysimbad(
                input_file  : str | Path | list[str],
                output_file : str | Path | None = None
    ) -> dict[str, dict]:
    """
    Query SIMBAD for the coordinates, spectral type, B and V magnitudes and parallax of a list of objects.

    ----------
    Parameters
    ----------
        input_file : str | Path | list[str]
            Either a non-empty list of object names or the path and filename
            of a text file with the object names in the first column (same
            format as 'queryviz').

        output_file : str | Path | None, optional
            Path and filename of a tab-separated output file. Default is None
            (results are only returned).

    -------
    Returns
    -------
        results_dict : dict[str, dict]
            A dictionary keyed by the input object names; each value is a
            dictionary containing the following key-value pairs:
                - 'main_id' : str
                    SIMBAD main identifier ('' if not found).

                - 'ra' : float
                    Right ascension (ICRS, epoch J2000) in degrees (NaN if
                    not available).

                - 'dec' : float
                    Declination (ICRS, epoch J2000) in degrees (NaN if not
                    available).

                - 'ra_hms' : str
                    Right ascension as 'hh:mm:ss.ss' ('' if not available).

                - 'dec_dms' : str
                    Declination as '±dd:mm:ss.s' ('' if not available).

                - 'sp_type' : str
                    Spectral type ('' if not available).

                - 'B' : float
                    Johnson B magnitude (mag; NaN if not available).

                - 'V' : float
                    Johnson V magnitude (mag; NaN if not available).

                - 'B-V' : float
                    Observed B - V color (mag; NaN if any magnitude is not
                    available).

                - 'plx' : float
                    Parallax in mas (NaN if not available).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        FileNotFoundError
            If 'input_file' is a path that does not exist, or the directory
            of 'output_file' does not exist.

        IsADirectoryError
            If 'input_file' or 'output_file' is an existing directory.

        NotADirectoryError
            If the parent of 'output_file' exists but is not a directory.

        PermissionError
            If 'input_file' is a file that cannot be read.

        ValueError
            If the list of object names is empty or contains empty names, or
            the file contains no object names.

        RuntimeError
            If the SIMBAD service cannot be configured (e.g. no internet
            connection).

    -----
    Notes
    -----
        1. INPUT
            - Names are read and normalized as in 'queryviz' (underscores
              replaced by spaces, repeated spaces collapsed).

        2. QUERY
            - Supports both the TAP-based astroquery API (>= 0.4.8; fields
              'sp_type', 'B', 'V', 'plx_value') and the legacy one ('sp',
              'flux(B)', 'flux(V)', 'plx').

            - Objects not found or failed queries issue a warning and keep
              the empty values of the template, so that the input order is
              preserved.

        3. PHYSICAL INTERPRETATION
            - The intrinsic color of OB stars is (B - V)_0 ≈ -0.26 to -0.32
              mag, so 'B-V' gives a first estimate of the color excess
              E(B - V) = (B - V) - (B - V)_0 and of the extinction
              A_V ≈ 3.1 * E(B - V).

            - Spectral types in SIMBAD are compilations of heterogeneous
              quality; for OB stars prefer homogeneous classifications (e.g.
              GOSSS, IACOB) when available.

        4. OUTPUT
            - If 'output_file' is given, a tab-separated table with the
              column 'Object' followed by the keys above is written.

    ----------
    References
    ----------
        - Wenger, M., Ochsenbein, F., Egret, D. et al. 2000, A&AS, 143, 9.
        - Ginsburg, A., Sipőcz, B. M., Brasseur, C. E. et al. 2019,
          AJ, 157, 98.
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.simbad import Simbad
    import numpy as np
    from tqdm import tqdm
    import warnings


    # ##### PANORAMA #####
    from .support.tex.tex_support import _check_output_file, _get_simbad_field, _normalize_name, _read_object_names, _to_float


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    object_names = _read_object_names(input_file)  # Validates the list or the file and reads the names
    #
    output_file = _check_output_file(output_file, allow_none=True)  # None: results are only returned


    # ##########  DICTIONARY  ##########
    results_dict = {}  # Filled with one entry per input object (keys are the input names)


    # ####################  PARAMETERS  ####################
    # ##########  COLUMN LOOKUP  ##########
    # Result key -> column names to try in the SIMBAD table (TAP name of astroquery >= 0.4.8 first, legacy name second)
    candidates = {
                  'main_id' : ("main_id", "MAIN_ID"),
                  'ra'      : ("ra", "RA"),
                  'dec'     : ("dec", "DEC"),
                  'sp_type' : ("sp_type", "SP_TYPE"),
                  'B'       : ("B", "FLUX_B"),
                  'V'       : ("V", "FLUX_V"),
                  'plx'     : ("plx_value", "PLX_VALUE"),
                 }


    # ##########  OUTPUT COLUMNS  ##########
    # Order of the columns of 'output_file' (after 'Object')
    output_keys = ['main_id', 'ra', 'dec', 'ra_hms', 'dec_dms', 'sp_type', 'B', 'V', 'B-V', 'plx']


    # ####################  COMPUTATION  ####################
    # ##########  CONFIG  ##########
    # The astroquery SIMBAD interface changed in version 0.4.8 (TAP service, new field names): the new names are tried first and the legacy ones as a fallback
    simbad = Simbad()  # Independent instance, so the global Simbad configuration is not modified
    try:
        simbad.add_votable_fields("sp_type", "B", "V", "plx_value")  # Field names of astroquery >= 0.4.8
    except Exception:
        try:
            simbad = Simbad()  # Fresh instance: the failed call may have left the previous one half-configured
            simbad.add_votable_fields("sp", "flux(B)", "flux(V)", "plx")  # Field names of astroquery < 0.4.8
        except Exception as e:
            raise RuntimeError(f"SIMBAD service could not be configured (internet connection required): {e}") from e


    # ##########  QUERY  ##########
    for orig in tqdm(object_names, desc="Querying SIMBAD..."):
        # Template with every key set to its "not available" value, so that missing data never leave a key undefined
        entry = {
                 'main_id' : '',
                 'ra'      : np.nan,
                 'dec'     : np.nan,
                 'ra_hms'  : '',
                 'dec_dms' : '',
                 'sp_type' : '',
                 'B'       : np.nan,
                 'V'       : np.nan,
                 'B-V'     : np.nan,
                 'plx'     : np.nan,
                }
        try:
            table = simbad.query_object(_normalize_name(orig))  # SIMBAD resolves names with spaces, not underscores
            if table is not None and len(table) > 0:
                row = table[0]  # Only the first (best) match is used


                # ##### IDENTIFIER & SPECTRAL TYPE #####
                entry['main_id'] = str(_get_simbad_field(row, 'main_id', candidates) or '').strip()  # '' if missing or masked
                entry['sp_type'] = str(_get_simbad_field(row, 'sp_type', candidates) or '').strip()


                # ##### COORDINATES #####
                ra  = _get_simbad_field(row, 'ra', candidates)
                dec = _get_simbad_field(row, 'dec', candidates)
                if ra is not None and dec is not None:
                    if isinstance(ra, str):  # Legacy API: sexagesimal strings ('hh mm ss.s', '±dd mm ss')
                        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
                    else:  # TAP API: decimal degrees
                        coord = SkyCoord(float(ra) * u.deg, float(dec) * u.deg, frame="icrs")
                    entry['ra']      = float(coord.ra.deg)
                    entry['dec']     = float(coord.dec.deg)
                    entry['ra_hms']  = coord.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)  # 'hh:mm:ss.ss'
                    entry['dec_dms'] = coord.dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True)  # '±dd:mm:ss.s'


                # ##### PHOTOMETRY & PARALLAX #####
                entry['B']   = _to_float(_get_simbad_field(row, 'B', candidates))  # NaN if missing or masked
                entry['V']   = _to_float(_get_simbad_field(row, 'V', candidates))
                entry['B-V'] = entry['B'] - entry['V']  # NaN if any magnitude is missing
                entry['plx'] = _to_float(_get_simbad_field(row, 'plx', candidates))  # mas
            else:
                warnings.warn(f"Object '{orig}' not found in SIMBAD.")
        except Exception as e:
            warnings.warn(f"SIMBAD query for '{orig}' failed ({e}).")  # Non-fatal: the empty template is kept

        results_dict[orig] = entry  # Insertion order = input order


    # ####################  OUTPUT  ####################
    # ##########  WRITE OUTPUT  ##########
    if output_file is not None:
        with open(output_file, "w", encoding="utf-8") as f:
            f.write("\t".join(["Object"] + output_keys) + "\n")  # Header row
            for obj, entry in results_dict.items():
                f.write("\t".join([obj] + [str(entry[k]) for k in output_keys]) + "\n")  # One tab-separated row per object


    # ##########  RESULTS  ##########
    return results_dict




# ========================================================================
# ============================== MATPLOTLIB ==============================
# ========================================================================

# ==============================
def read_pkl(
             filename : str | Path,
             show     : bool = True
    ) -> Figure:
    """
    Open a matplotlib figure stored in a pickle file, so that it can be displayed and edited again.

    ----------
    Parameters
    ----------
        filename : str | Path
            Path and filename of the pickle file containing the matplotlib
            figure.

        show : bool, optional
            If True, display the figure with plt.show(). Default is True.

    -------
    Returns
    -------
        fig : Figure
            The loaded matplotlib.figure.Figure, which can be modified
            (labels, limits, annotations...) and saved again.

    ------
    Raises
    ------
        TypeError
            If 'filename' is not a string or a pathlib.Path, or 'show' is not
            a boolean.

        FileNotFoundError
            If the file does not exist.

        IsADirectoryError
            If the path points to a directory.

        PermissionError
            If the file cannot be read.

        RuntimeError
            If the file cannot be unpickled or does not contain a figure.

    -----
    Notes
    -----
        - Only open pickle files from trusted sources: unpickling can execute
          arbitrary code.

        - Pickled figures are tied to the matplotlib version that created
          them; loading them with a different version may fail or render
          differently.

        - Figures made with LaTeX text keep text.usetex=True, so LaTeX and the
          packages of the preamble used to create them must be available to
          draw them again.

        - With a non-interactive backend (e.g. 'Agg'), show=True does not
          open any window.

        - See 'savefig_pkl' to save a figure together with its pickle file.

    ----------
    References
    ----------
        - Hunter, J. D. 2007, CSE, 9, 90.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    import pickle


    # ##### PANORAMA #####
    from .support.tex.tex_support import _check_input_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    filename = _check_input_file(filename, "filename", "Pickle file", "pickle file")
    #
    if not isinstance(show, bool):
        raise TypeError(f"Display flag 'show' must be a boolean. Got {type(show).__name__}.")


    # ####################  COMPUTATION  ####################
    try:
        with open(filename, "rb") as f:
            fig = pickle.load(f)  # Unpickled figures are registered with pyplot
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the figure from '{filename}': {e}") from e

    if not isinstance(fig, Figure):  # A valid pickle file may contain any Python object
        raise RuntimeError(f"File '{filename}' does not contain a matplotlib figure. Got {type(fig).__name__}.")


    # ####################  OUTPUT  ####################
    if show:
        plt.show()  # Inline backend: displays the figure; interactive backends: opens a window
    return fig




# ==============================
def savefig_pkl(
                fig      : Figure,
                filename : str | Path,
                formats  : list[str] | tuple[str, ...] = ("pdf",),
                dpi      : int | float = 500,
                close    : bool = True
    ) -> list[Path]:
    """
    Save a figure following the PANORAMA conventions, together with a pickle copy for later editing.

    ----------
    Parameters
    ----------
        fig : Figure
            matplotlib.figure.Figure to save.

        filename : str | Path
            Path and base filename, without extension (a known extension, if
            given, is removed). E.g. 'figs/HD91316_SiIII' produces
            'figs/HD91316_SiIII.pdf' and 'figs/HD91316_SiIII.pkl'.

        formats : list[str] | tuple[str, ...], optional
            Image formats to write (any format supported by Figure.savefig,
            e.g. 'pdf', 'png', 'eps'). A single format given as a string is
            also accepted. Default is ('pdf',).

        dpi : int | float, optional
            Resolution in dots per inch for raster elements and formats.
            Must be finite and positive. Default is 500.

        close : bool, optional
            If True, close the figure after saving to free memory. Default
            is True.

    -------
    Returns
    -------
        written_files : list[Path]
            Paths of the files written (image formats first, pickle last).

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If 'formats' is empty, contains an empty format or 'pkl', or
            'dpi' is not finite or not positive.

        FileNotFoundError
            If the output directory does not exist.

        IsADirectoryError
            If 'filename' is an existing directory.

        NotADirectoryError
            If the parent of 'filename' exists but is not a directory.

    -----
    Notes
    -----
        - Saving convention:
                fig.align_xlabels()
                fig.align_ylabels()
                fig.savefig(filename, format='pdf', dpi=dpi, bbox_inches='tight', pad_inches=0.01)
                plt.close(fig)

        - The pickle copy (written before closing) can be reopened with
          'read_pkl' to change labels, limits or annotations after referee
          comments without rerunning the analysis (e.g. MCMC fits).

    ----------
    References
    ----------
        - Hunter, J. D. 2007, CSE, 9, 90.
    """

    # ####################  MODULES & PACKAGES  ####################
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    import numbers
    import numpy as np
    import os
    from pathlib import Path
    import pickle


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(fig, Figure):
        raise TypeError(f"Figure 'fig' must be a matplotlib Figure. Got {type(fig).__name__}.")
    #
    if not isinstance(filename, (str, os.PathLike)):
        raise TypeError(f"Base filename 'filename' must be a path (str or pathlib.Path). Got {type(filename).__name__}.")
    filename = Path(filename)
    if filename.is_dir():
        raise IsADirectoryError(f"Base filename 'filename' is an existing directory, not a base filename: '{filename}'.")
    elif not filename.parent.exists():
        raise FileNotFoundError(f"Output directory does not exist: '{filename.parent}'.")
    elif not filename.parent.is_dir():
        raise NotADirectoryError(f"Output directory is not a directory: '{filename.parent}'.")
    #
    if isinstance(formats, str):  # Allow a single format given as a string
        formats = (formats,)
    if not isinstance(formats, (list, tuple)):
        raise TypeError(f"Image formats 'formats' must be a list or tuple of strings. Got {type(formats).__name__}.")
    elif len(formats) == 0:
        raise ValueError("Image formats 'formats' must contain at least one format. Got an empty sequence.")
    for fmt in formats:
        if not isinstance(fmt, str):
            raise TypeError(f"Image formats 'formats' must contain only strings. Got {type(fmt).__name__}.")
        elif not fmt.strip().lstrip("."):
            raise ValueError(f"Image formats 'formats' must not contain empty formats. Got '{fmt}'.")
        elif fmt.strip().lower().lstrip(".") == "pkl":
            raise ValueError(f"Image formats 'formats' must not include 'pkl' (the pickle copy is always written). Got {formats}.")
    #
    if not isinstance(dpi, numbers.Real) or isinstance(dpi, bool):
        raise TypeError(f"Figure resolution 'dpi' must be a real number (dots per inch). Got {type(dpi).__name__}.")
    elif not np.isfinite(dpi):
        raise ValueError(f"Figure resolution 'dpi' must be finite. Got {dpi}.")
    elif dpi <= 0:
        raise ValueError(f"Figure resolution 'dpi' must be positive. Got {dpi}.")
    #
    if not isinstance(close, bool):
        raise TypeError(f"Close flag 'close' must be a boolean. Got {type(close).__name__}.")


    # ##########  BASE FILENAME  ##########
    formats = [fmt.strip().lower().lstrip(".") for fmt in formats]  # Normalized formats (e.g. '.PDF' → 'pdf')
    if filename.suffix.lower().lstrip(".") in formats + ["pkl"]:  # Remove a known extension
        filename = filename.with_suffix('')


    # ####################  OUTPUT  ####################
    # ##########  SAVE  ##########
    # ##### AXIS ALIGNMENT #####
    fig.align_xlabels()  # Align the x labels of the subplots of the same row
    fig.align_ylabels()  # Align the y labels of the subplots of the same column

    written_files = []  # Paths in the order they are written


    # ##### IMAGE FORMATS #####
    for fmt in formats:
        path = filename.with_name(f"{filename.name}.{fmt}")  # with_name keeps dots of the base name (with_suffix would replace them)
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches="tight", pad_inches=0.01)  # PANORAMA saving convention
        written_files.append(path)


    # ##### PICKLE #####
    path = filename.with_name(f"{filename.name}.pkl")
    with open(path, "wb") as f:
        pickle.dump(fig, f)  # Written before closing, so that the figure can be reopened with 'read_pkl'
    written_files.append(path)


    # ##### CLOSE #####
    if close:
        plt.close(fig)  # Free memory


    # ##########  RESULTS  ##########
    return written_files
