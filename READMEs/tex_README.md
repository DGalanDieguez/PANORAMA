# `panorama.tex` — Manuscript, Catalog and Figure Utilities

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](../LICENSE)

> Bibliography checks for LaTeX manuscripts, VizieR and SIMBAD queries for lists of targets, and pickled figures that can be reopened and edited.

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

`tex` collects small tools for the final stages of a spectroscopic project on (massive) stars:

* **Bibliography:** entries of a `.bib` file never cited in the compiled manuscript, citations without a `.bib` entry, and duplicated entries (same key, or same paper stored under different keys).
* **Catalog queries:** all the columns of a VizieR catalog, or the coordinates, spectral type, B and V magnitudes and parallax from SIMBAD, for a list of targets.
* **Figures:** saving a figure together with a pickle copy, and reopening it later to change labels, limits or annotations without rerunning the analysis.

<br>

**Conventions shared by the module:**
* File arguments (inputs and outputs) accept a `str` or a `pathlib.Path`; functions that return paths return `pathlib.Path` objects. Output directories are never created: they must exist.
* Bibliography checks read the `.aux` file of a full compilation (LaTeX + BibTeX/Biber + LaTeX) and support both BibTeX (`\citation`) and biblatex (`\abx@aux@cite`, with or without the refsection argument). `@string`, `@comment` and `@preamble` blocks of the `.bib` file are not entries and are ignored.
* Object names can be given as a non-empty Python list or as a text file with the names in the first column (`#` comments, empty lines and any further columns are ignored; names with spaces are written with underscores, e.g. `HD_91316`). Names are normalized before querying: underscores are replaced by spaces and repeated spaces are collapsed (`HD_91316` → `HD 91316`, `BD+60_2522` → `BD+60 2522`).
* Catalog queries require internet access. Objects not found and failed queries are kept with empty values (`nan` / `''`), so the output preserves the input order; failed queries (and, in SIMBAD, objects not found) issue a warning.
* Figures are saved following the PANORAMA A&A conventions: labels aligned, `bbox_inches='tight'`, `pad_inches=0.01`, 500 dpi by default, plus a pickle copy.
* All functions validate their inputs before computing (`TypeError` for wrong types, `FileNotFoundError` / `IsADirectoryError` / `NotADirectoryError` / `PermissionError` for paths, `ValueError` for invalid values, including empty or non-finite ones). Error messages name the meaning of the argument, the argument and, when relevant, its unit (e.g. `Column width 'col_width' (in characters) must be positive.`).

---
<br>

## 2. Function Summary

### Bibliography (`BIB`)

| Function | Description |
|---|---|
| [`nobib`](#nobib) | `.bib` entries never cited in the `.aux` file. |
| [`missingbib`](#missingbib) | Citations in the `.aux` file without a `.bib` entry. |
| [`dupbib`](#dupbib) | Duplicated `.bib` entries, by key and by DOI. |

<br>

### Catalog queries (`QUERY`)

| Function | Description |
|---|---|
| [`queryviz`](#queryviz) | Query a VizieR catalog for a list of objects and save a fixed-width table. |
| [`querysimbad`](#querysimbad) | Coordinates, spectral type, B, V, B−V and parallax from SIMBAD. |

<br>

### Figures (`MATPLOTLIB`)

| Function | Description |
|---|---|
| [`read_pkl`](#read_pkl) | Reopen a pickled figure to display or edit it. |
| [`savefig_pkl`](#savefig_pkl) | Save a figure in one or more formats plus a pickle copy. |

---
<br>

## 3. Function Reference

Each function lists its arguments in an **Inputs** table and its outputs in a **Returns** table. For functions that return a dictionary, the first row is the dictionary itself and the following rows are its keys (nested keys are written as `[name]['key']` or `'outer'[name]['key']`).

### `nobib`
```python
nobib(
      aux_file : str | Path,
      bib_file : str | Path
    ) -> list[str]
```
Compares the citation keys of the compiled document (`.aux` file) with the keys defined in the `.bib` file and returns those never cited. Keys are compared case-sensitively, as LaTeX does. If the document uses `\nocite{*}`, a warning is issued and only explicit citations are compared. Typical use: cleaning an ADS-exported `.bib` file before submission. See [`missingbib`](#missingbib) for the opposite check and [`dupbib`](#dupbib) for duplicated entries.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `aux_file` | `str \| Path` | — | LaTeX `.aux` file produced during compilation. |
| `bib_file` | `str \| Path` | — | Bibliography `.bib` file. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `unused_keys` | `list[str]` | Alphabetically sorted keys present in the `.bib` file but never cited. |

**Raises:** `TypeError`, `FileNotFoundError`, `IsADirectoryError` or `PermissionError` for invalid paths.

<br>

### `missingbib`
```python
missingbib(
           aux_file : str | Path,
           bib_file : str | Path
    ) -> list[str]
```
Opposite check of `nobib`: keys cited in the document but not defined in the `.bib` file, which produce "Citation undefined" warnings and `?` marks in the PDF. Keys are compared case-sensitively, as LaTeX does (`sana2012` does not match `Sana2012`). The `\nocite{*}` wildcard is ignored without warning.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `aux_file` | `str \| Path` | — | LaTeX `.aux` file produced during compilation. |
| `bib_file` | `str \| Path` | — | Bibliography `.bib` file. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `missing_keys` | `list[str]` | Alphabetically sorted keys cited in the document but not defined in the `.bib` file. |

**Raises:** `TypeError`, `FileNotFoundError`, `IsADirectoryError` or `PermissionError` for invalid paths.

<br>

### `dupbib`
```python
dupbib(
       bib_file : str | Path
    ) -> dict[str, dict]
```
Finds keys defined more than once (case-insensitive, as BibTeX treats them; BibTeX stops with "Repeated entry" and biblatex silently keeps the first one) and DOIs shared by several keys, a common situation when merging `.bib` files from different projects. DOIs, given in braces or in quotes, are compared case-insensitively after removing any `https://doi.org/` (or `dx.doi.org`) prefix.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bib_file` | `str \| Path` | — | Bibliography `.bib` file. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, dict]` | Dictionary with the keys below. |
| `'keys'` | `dict[str, int]` | `{key: number of occurrences}` for keys defined more than once (case-insensitive; first spelling found). Empty dictionary if there are no duplicated keys. |
| `'dois'` | `dict[str, list[str]]` | `{normalized DOI: [keys]}` for DOIs shared by more than one entry. Empty dictionary if there are no duplicated DOIs. |

**Raises:** `TypeError`, `FileNotFoundError`, `IsADirectoryError` or `PermissionError` for invalid paths.

<br>

### `queryviz`
```python
queryviz(
         input_file  : str | Path | list[str],
         catalog     : str,
         output_file : str | Path = "vizier_output.txt",
         col_width   : int = 15
    ) -> None
```
Queries a VizieR catalog for each object, requesting all its columns and keeping the first row returned (the first match within the default search radius; check crowded fields). Writes a table with `Object` followed by every catalog column (names taken from the first successful query), left-aligned and padded to `col_width` characters, in the order of the input. Values longer than `col_width` are not truncated and merge with the next column, so use wider columns for catalogs with long fields. Objects without data get `nan` in every column; queries that raise an error also issue a warning.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_file` | `str \| Path \| list[str]` | — | Non-empty list of object names, or text file with the names in the first column. |
| `catalog` | `str` | — | VizieR catalog identifier (e.g. `I/239/hip_main`, `I/355/gaiadr3`, `B/mk/mktypes`). |
| `output_file` | `str \| Path` | `"vizier_output.txt"` | Output file (its directory must exist). |
| `col_width` | `int` | `15` | Fixed column width [characters]; must be positive. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `None` | `None` | The results are written to `output_file`. |

**Raises:** `TypeError` for wrong types; `FileNotFoundError`, `IsADirectoryError`, `NotADirectoryError` or `PermissionError` for invalid paths; `ValueError` if the list is empty or contains empty names, the file contains no names, `catalog` is empty or `col_width` is not positive; `RuntimeError` if every query fails or returns no data.

<br>

### `querysimbad`
```python
querysimbad(
            input_file  : str | Path | list[str],
            output_file : str | Path | None = None
    ) -> dict[str, dict]
```
Retrieves from SIMBAD the main identifier, ICRS coordinates (epoch J2000; degrees and sexagesimal), spectral type, Johnson B and V magnitudes, the observed B − V color and the parallax of each object. Works with both the TAP-based astroquery interface (≥ 0.4.8) and the legacy one. Objects not found and failed queries issue a warning and keep empty values. Since (B − V)₀ ≈ −0.26 to −0.32 mag for OB stars, `B-V` gives a first estimate of E(B − V) and A_V ≈ 3.1 E(B − V). SIMBAD spectral types are heterogeneous compilations: for OB stars prefer homogeneous classifications (e.g. GOSSS, IACOB) when available.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_file` | `str \| Path \| list[str]` | — | Non-empty list of object names, or text file with the names in the first column. |
| `output_file` | `str \| Path \| None` | `None` | Tab-separated output file (`Object` followed by the keys below); if `None`, results are only returned. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, dict]` | Dictionary keyed by the input object names (`name`); each value is a dictionary with the keys below. |
| `[name]['main_id']` | `str` | SIMBAD main identifier (`''` if not found). |
| `[name]['ra']` | `float` | Right ascension (ICRS, epoch J2000) [deg] (`NaN` if not available). |
| `[name]['dec']` | `float` | Declination (ICRS, epoch J2000) [deg] (`NaN` if not available). |
| `[name]['ra_hms']` | `str` | Right ascension as `hh:mm:ss.ss` (`''` if not available). |
| `[name]['dec_dms']` | `str` | Declination as `±dd:mm:ss.s` (`''` if not available). |
| `[name]['sp_type']` | `str` | Spectral type (`''` if not available). |
| `[name]['B']` | `float` | Johnson B magnitude [mag] (`NaN` if not available). |
| `[name]['V']` | `float` | Johnson V magnitude [mag] (`NaN` if not available). |
| `[name]['B-V']` | `float` | Observed B − V color [mag] (`NaN` if any magnitude is not available). |
| `[name]['plx']` | `float` | Parallax [mas] (`NaN` if not available). |

**Raises:** `TypeError` for wrong types; `FileNotFoundError`, `IsADirectoryError`, `NotADirectoryError` or `PermissionError` for invalid paths; `ValueError` if the list is empty or contains empty names, or the file contains no names; `RuntimeError` if SIMBAD cannot be configured (e.g. no internet access).

<br>

### `read_pkl`
```python
read_pkl(
         filename : str | Path,
         show     : bool = True
    ) -> Figure
```
Loads a pickled `matplotlib` figure (registered again with `pyplot`), optionally displays it, and returns it so that it can be edited and saved again (e.g. with `savefig_pkl`). Only open pickle files from trusted sources, since unpickling can execute arbitrary code. Pickled figures are tied to the `matplotlib` version that created them, and figures made with LaTeX text need LaTeX and the same preamble packages to be drawn again. With a non-interactive backend (e.g. `Agg`), `show=True` does not open any window; with the inline backend of Jupyter, the figure is shown in the notebook.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `filename` | `str \| Path` | — | Pickle file containing the figure. |
| `show` | `bool` | `True` | Display the figure with `plt.show()`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `fig` | `Figure` | The loaded `matplotlib.figure.Figure`, which can be modified (labels, limits, annotations...) and saved again. |

**Raises:** `TypeError` for wrong types; `FileNotFoundError`, `IsADirectoryError` or `PermissionError` for invalid paths; `RuntimeError` if the file cannot be unpickled or does not contain a figure.

<br>

### `savefig_pkl`
```python
savefig_pkl(
            fig      : Figure,
            filename : str | Path,
            formats  : list[str] | tuple[str, ...] = ("pdf",),
            dpi      : int | float = 500,
            close    : bool = True
    ) -> list[Path]
```
Aligns the axis labels, saves the figure in each requested format with the PANORAMA conventions (`bbox_inches='tight'`, `pad_inches=0.01`) and writes a pickle copy before closing it. A known extension in `filename` (`.pkl` or one of the requested formats) is removed, so `"figs/HD91316_SiIII.pdf"` and `"figs/HD91316_SiIII"` are equivalent; any other suffix is kept as part of the base name (e.g. `"figs/SiIII_4552.62"` → `SiIII_4552.62.pdf`). A single format can be given as a string (e.g. `formats="png"`), and formats are case-insensitive (`".PDF"` → `pdf`).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fig` | `Figure` | — | `matplotlib.figure.Figure` to save. |
| `filename` | `str \| Path` | — | Path and base filename (extension optional; its directory must exist). |
| `formats` | `list[str] \| tuple[str, ...]` | `("pdf",)` | Image formats (`pdf`, `png`, `eps`, `svg`...); `pkl` is always written and must not be listed. |
| `dpi` | `int \| float` | `500` | Resolution [dots per inch]; must be finite and positive. |
| `close` | `bool` | `True` | Close the figure after saving to free memory. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `written_files` | `list[Path]` | Paths of the files written (image formats first, pickle last). |

**Raises:** `TypeError` for wrong types; `FileNotFoundError`, `IsADirectoryError` or `NotADirectoryError` for invalid paths; `ValueError` if `formats` is empty, contains an empty format or `pkl`, or `dpi` is not finite or not positive.

---
<br>

## 4. Dependencies

| Package | Used for |
|---|---|
| `astropy` | Coordinate conversions (`querysimbad`). |
| `astroquery` | VizieR and SIMBAD queries (`queryviz`, `querysimbad`; require internet access). |
| `matplotlib` | Saving and reopening figures (`savefig_pkl`, `read_pkl`). |
| `numpy` | Numerical checks and missing values. |
| `tqdm` | Progress bars of the catalog queries. |

---
<br>

## 5. Internal Structure

| Path | Content |
|---|---|
| `src/panorama/tex.py` | Public functions, grouped in three sections: `BIB`, `QUERY`, `MATPLOTLIB`. |
| `src/panorama/support/tex/tex_support.py` | Shared helpers, in three groups. `FILES`: validation of input files (`_check_input_file`, used by `nobib`, `missingbib`, `dupbib` and `read_pkl`, and reused by `panorama.spec.readspec`) and of output files (`_check_output_file`, used by `queryviz` and `querysimbad`, and reused by `writespec`, `rmvspec` and `normspec` of `spec` and by `fits_catalog` of `obs`), so that file arguments give the same error messages in the whole package. `OBJECTS`: reading and validation of the lists (or files) of object names (`_read_object_names`) and normalization of the names before querying SIMBAD / VizieR (`_normalize_name`), shared by `queryviz` and `querysimbad`. `QUERY`: reading of a field of a SIMBAD result row with the column names of both the TAP-based (astroquery >= 0.4.8) and the legacy APIs (`_get_simbad_field`), and tolerant conversion of catalog values to float, with NaN for missing data (`_to_float`), used by `querysimbad`. |

---
<br>

## 6. Demo

The Jupyter notebook [`Demos/tex_demo.ipynb`](../Demos/tex_demo.ipynb) runs every function of the module in 9 examples, one section per example:

* **Examples 1-3:** uncited, missing and duplicated entries of a mock manuscript, with BibTeX and biblatex citations (with and without the refsection argument), `@string`, `@comment` and `@preamble` blocks, a key cited with the wrong case, a key repeated in another case and two papers stored under two keys each (DOIs with prefix, in other case and in quotes), including the `\nocite{*}` warning of `nobib` and the wildcard ignored by `missingbib`.
* **Examples 4-6:** VizieR queries of OB stars from a list (Hipparcos) and from a file (Gaia DR3; comments, empty lines and an extra column ignored, and an unknown object kept as a row of `nan`), and SIMBAD queries with and without output file (an unknown object kept with empty values and its warning), with an estimate of the color excess.
* **Examples 7-9:** saving, reopening and editing a pickled line-profile figure (the original, reopened and edited figures are shown inline), and the special cases of `savefig_pkl` (single format string, known extension removed, other suffixes kept in the base name, `pkl` rejected).

Input files are written to `Demos/data_and_results/tex/data/` and outputs to `Demos/data_and_results/tex/results/`. Examples 4-6 require internet access (and `astroquery`); without it they print `Skipped: ...` and the rest of the notebook still runs. Figures are shown inline (`%matplotlib inline`) and also saved with `savefig_pkl`.

```bash
pip install -e ".[demos]"     # PANORAMA + JupyterLab and PyQt6 (demo extras)
cd Demos
jupyter lab tex_demo.ipynb    # Or: jupyter notebook tex_demo.ipynb
```

---
<br>

## 7. References

* Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022, *ApJ*, 935, 167 — `astropy` (coordinate conversions).
* Ginsburg, A., Sipőcz, B. M., Brasseur, C. E. et al. 2019, *AJ*, 157, 98 — `astroquery`.
* Hunter, J. D. 2007, *CSE*, 9, 90 — `matplotlib`.
* Ochsenbein, F., Bauer, P. & Marcout, J. 2000, *A&AS*, 143, 23 — the VizieR database.
* Wenger, M., Ochsenbein, F., Egret, D. et al. 2000, *A&AS*, 143, 9 — the SIMBAD database.
