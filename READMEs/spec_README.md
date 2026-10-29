# `panorama.spec` — Reading, Writing, Plotting, Cleaning and Normalization of Spectra

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](../LICENSE)

> Reading and writing of 1D spectra in ASCII and FITS formats, interactive comparison with the diagnostic lines of OB stars, interactive removal of spectral regions, and IRAF-like continuum normalization, automatic or interactive.

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

`spec` covers the first steps of the work with the observed spectra of (massive) stars, before any line measurement or model fitting:

* **Reading:** 1D spectra stored as ASCII tables, FITS binary tables or FITS images with WCS keywords (linear, log10 or ln λ scales, e.g. HERMES and IACOB spectra), optionally gzip-compressed, returned as `(wave, flux)` arrays with increasing wavelengths.
* **Writing:** `(wave, flux)` arrays or spectrum files written as ASCII tables or FITS binary tables (format chosen from the extension, optionally gzip-compressed), with optional header lines; every file written can be read back with `readspec`, so the same function also converts between formats.
* **Plotting:** several spectra in one interactive window, with check buttons to show or hide each one and the diagnostic lines of O and B stars (H, He I, He II, C, N, O, Mg and Si) drawn and labeled.
* **Cleaning:** removal of wavelength spans (cosmic rays, telluric or interstellar features, bad pixels) selected with the mouse, with undo / redo and export of the removed spans.
* **Normalization:** continuum fitting in the style of the IRAF `continuum` task (cubic or linear splines, Legendre or Chebyshev polynomials, sample regions, pixel averaging and iterative sigma rejection), in batch mode or in an interactive window with automatic and manual continua.

<br>

**Conventions shared by the module:**
* A **spectrum** is given either as a `(wave, flux)` tuple of 1D array-like objects with the same length, or as the path to a spectrum file read with [`readspec`](#readspec) (`str` or `pathlib.Path`).
* **Wavelengths** keep the units of the input spectrum (usually Å) and are never converted or corrected for the barycentric velocity (see [`bary_corr`](obs_README.md#bary_corr)); plot axes are labeled in Å. The diagnostic lines of OB stars are rest-frame **air** wavelengths in Å, so spectra not corrected for the radial velocity of the star show their lines shifted by v_rad / c · λ from the marks.
* **Non-finite values:** `readspec` rejects non-finite wavelengths and keeps non-finite fluxes with a warning, and `writespec` does the same. `rmvspec` discards points with a non-finite wavelength or flux, and `normspec` discards points with a non-finite wavelength and keeps non-finite fluxes, which are never used as sample points; both warn and sort the spectrum by wavelength if needed.
* **Output spectra** are written with [`writespec`](#writespec) by every function of the module (`rmvspec` and `normspec` included): ASCII (`.txt`, `.dat`, `.asc`, `.ascii`; two tab-separated columns after a `# wave flux` header line, also tab-separated) or FITS (`.fits`, `.fit`, `.fts`; binary table with `WAVE` and `FLUX` columns), chosen from the extension of the output file, with an optional final `.gz`. `normspec` adds header lines describing the fit. Output directories are never created: they must exist.
* **Interactive functions** (`plotspec`, `rmvspec` and the interactive mode of `normspec`) need an interactive Matplotlib backend (e.g. `QtAgg`, `TkAgg`; `%matplotlib qt` in Jupyter). With a non-interactive backend (e.g. `Agg`, the inline backend) `plotspec` draws the figure without opening a window (with the inline backend of Jupyter it is shown in the notebook, without working check buttons), and `rmvspec` and `normspec` issue a `RuntimeWarning` and return at once (see their descriptions). With an interactive backend, `rmvspec` and the interactive mode of `normspec` return when their window is closed, also in Jupyter.
* All functions validate their inputs before computing (`TypeError` for wrong types, `FileNotFoundError` / `IsADirectoryError` / `NotADirectoryError` / `PermissionError` for paths, `ValueError` for invalid values, including empty or non-finite ones, `FileExistsError` when `writespec` would replace a file without `overwrite=True`). Error messages name the meaning of the argument, the argument and, when relevant, its unit (e.g. `Wavelength range 'wave_lim' must satisfy min < max. Got (4600.0, 4440.0).`).

---
<br>

## 2. Function Summary

### Spectra (`SPECTRA`)

| Function | Description |
|---|---|
| [`readspec`](#readspec) | Read a 1D spectrum file (ASCII or FITS, optionally gzip-compressed). |
| [`writespec`](#writespec) | Write a 1D spectrum to an ASCII or FITS file (optionally gzip-compressed). |
| [`plotspec`](#plotspec) | Plot one or several spectra with toggles and the diagnostic lines of OB stars, and return the figure. |
| [`rmvspec`](#rmvspec) | Remove regions of a spectrum interactively and return the cleaned spectrum. |
| [`normspec`](#normspec) | Fit the continuum of a spectrum and normalize it, automatically or interactively. |

---
<br>

## 3. Function Reference

Each function lists its arguments in an **Inputs** table and its outputs in a **Returns** table.

### `readspec`
```python
readspec(
         file : str | Path
    ) -> tuple[np.ndarray[float], np.ndarray[float]]
```
Reads a 1D spectrum file and returns its wavelength and flux arrays. The format is chosen from the extension (a final `.gz` is accepted):

* **ASCII** (`.txt`, `.dat`, `.asc`, `.ascii`): two or more columns (wavelength, flux and optional extra columns, which may contain text and trigger a warning), separated by whitespace, commas or semicolons. Lines starting with `#` and empty lines are ignored, non-numeric lines before the first data line (column names, e.g. `wavelength flux`) are skipped with a warning, and every data line must have the same number of columns.
* **FITS** (`.fits`, `.fit`, `.fts`): first, a binary table with wavelength and flux columns (`WAVE` / `WAVELENGTH` / `LAMBDA` / `LAM` and `FLUX` / `INTENSITY` / `FLX` / `SPEC` / `FLAMBDA` / `COUNTS` / `CTS`, case-insensitive), searched across every extension (tables without these columns are skipped); otherwise, an image with WCS keywords (`CRVAL1`, `CDELT1` or `CD1_1`, `CRPIX1`) in the primary HDU or, if it is empty, in the first image extension with data. For multi-dimensional images (e.g. flux + error, IRAF multispec) the first row is used, and a missing `CRPIX1` means pixel 1 (FEROS / IRAF convention).
* **Logarithmic scales:** `DC-FLAG = 1` (IRAF, log10 λ), `CTYPE1` ending in `-LOG` (FITS standard; Greisen et al. 2006), or ln λ values in files with `log` in their name (e.g. HERMES `*_log_merged_cf.fits` spectra; the name check is needed because their headers do not flag the logarithmic scale).

Descending wavelength scales (e.g. negative `CDELT1`) are reversed. Every file written by [`writespec`](#writespec) can be read back, since both functions use the same extensions. Each format is read by a private function (`_read_ascii`, `_read_fits`; see [Section 5](#5-internal-structure)) that only returns the raw arrays, and `readspec` applies the generic checks to all of them.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `file` | `str \| Path` | — | Spectrum file (ASCII or FITS, optionally `.gz`). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `wave` | `np.ndarray[float]` | Wavelengths, strictly increasing [units of the file, usually Å]. |
| `flux` | `np.ndarray[float]` | Fluxes (same length as `wave`; non-finite values kept, with a warning). |

**Raises:** `TypeError` if `file` is not a path; `FileNotFoundError`, `IsADirectoryError` or `PermissionError` for invalid paths; `ValueError` if the extension of the file is not supported or the file does not contain a valid spectrum (no data, missing WCS keywords, arrays of different lengths, non-finite or non-monotonic wavelengths, e.g. overlapping échelle orders); `RuntimeError` if the file cannot be read for any other reason (e.g. a corrupt FITS file).

<br>

### `writespec`
```python
writespec(
          spec        : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
          output_file : str | Path,
          overwrite   : bool = False,
          header      : str | None = None
    ) -> Path
```
Writes a 1D spectrum and returns the path of the file. The format is chosen from the extension of `output_file` (case-insensitive), and a final `.gz` writes the file gzip-compressed (e.g. `HD886.fits.gz`):

* **ASCII** (`.txt`, `.dat`, `.asc`, `.ascii`): two tab-separated columns (wavelength, flux) after the `#` header lines, i.e. the lines of `header` (if given) and the column names `wave` and `flux`. Wavelengths are written with `%.12g` and fluxes with `%.10g`, so absolute fluxes (~1e-15) keep their precision; non-finite fluxes are written as `nan` / `inf`.
* **FITS** (`.fits`, `.fit`, `.fts`): an empty primary HDU plus a binary table extension `SPECTRUM` with the `float64` columns `WAVE` and `FLUX`, and the lines of `header` (if given, ASCII characters only) as `COMMENT` cards of the extension. No units (`TUNITn`) are written, because wavelengths keep the units of the input.

The extensions are those accepted by [`readspec`](#readspec), so every file written can be read back; giving a spectrum file as `spec` converts it to the format of `output_file` (e.g. `writespec("HD886.fits", "HD886.txt")`). Descending wavelength arrays are reversed, repeated or unsorted wavelengths are rejected, and non-finite fluxes are kept with a warning. An existing file is only replaced with `overwrite=True`. [`rmvspec`](#rmvspec) and [`normspec`](#normspec) write their output spectra with this function (replacing the previous file, since they save repeatedly from their windows).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `spec` | `tuple[np.ndarray[float], np.ndarray[float]] \| list[np.ndarray[float]] \| str \| Path` | — | Spectrum to write: `(wave, flux)` tuple, `[wave, flux]` list or spectrum file (format conversion). |
| `output_file` | `str \| Path` | — | Output spectrum file: ASCII (`.txt`, `.dat`, `.asc`, `.ascii`) or FITS (`.fits`, `.fit`, `.fts`), optionally followed by `.gz`. Its directory must exist. |
| `overwrite` | `bool` | `False` | Replace `output_file` if it already exists. |
| `header` | `str \| None` | `None` | Text stored with the spectrum (e.g. how it was processed): one `#` line per line of the text (ASCII, given without the `#`) or one `COMMENT` card per line (FITS, ASCII characters only). `None` = no extra lines. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `output_file` | `Path` | Path of the written spectrum file. |

**Raises:** `TypeError` for wrong types; `ValueError` if the spectrum is not valid (arrays not 1D, of different lengths, empty, with non-finite or non-monotonic wavelengths), the extension of `output_file` is not supported, or `header` contains non-ASCII characters for a FITS file; `IsADirectoryError`, `FileNotFoundError` or `NotADirectoryError` for an invalid `output_file`; `FileExistsError` if `output_file` exists and `overwrite` is `False`; `RuntimeError` if the file cannot be written for any other reason (e.g. no write permission). Spectrum files raise the exceptions of [`readspec`](#readspec).

<br>

### `plotspec`
```python
plotspec(
         spec       : list[tuple[np.ndarray[float], np.ndarray[float]] | str | Path] | tuple[np.ndarray[float], np.ndarray[float]] | str | Path,
         wave_lim   : tuple[float, float] | None = None,
         show_lines : bool = True,
         show       : bool = True
    ) -> Figure
```
Plots one or several spectra in an interactive window, with a check button per spectrum to show or hide it. Spectrum files are labeled with their file name (without extension) and `(wave, flux)` tuples as `Spectrum N`; duplicated labels get a ` (k)` suffix. If `show_lines` is `True`, the 40 diagnostic lines of O and B stars inside the wavelength range of the spectra (Balmer lines, He I, He II, and C, N, O, Mg and Si lines used in spectral classification; Walborn & Fitzpatrick 1990; Gray & Corbally 2009; Sota et al. 2011) are drawn as dotted vertical lines, colored by species (H I red, He I blue, He II purple; metals by element: C brown, N orange, O cyan, Mg olive, Si green) and labeled at the bottom of the axes; the labels follow their lines when zooming or panning in wavelength and stay at the bottom when zooming in flux. An extra check button, `Lines (OB-type stars)`, shows or hides them all. The lines are drawn at rest, and blends are labeled by their main component (e.g. `He I+II 4026`, `C III 4650` for C III 4647-4651). The figure is returned without saving it (it can be saved with [`savefig_pkl`](tex_README.md#savefig_pkl)), and it keeps a reference to the check buttons so that they keep working while it is open. With a non-interactive backend the figure is drawn but no window is opened; with the inline backend of Jupyter it is shown in the notebook (without working check buttons).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `spec` | `list[tuple[np.ndarray[float], np.ndarray[float]] \| str \| Path] \| tuple[np.ndarray[float], np.ndarray[float]] \| str \| Path` | — | A spectrum (`(wave, flux)` tuple or spectrum file) or a non-empty list of them. |
| `wave_lim` | `tuple[float, float] \| None` | `None` | Wavelength range `(min, max)` shown [units of the spectra]; the flux axis is scaled to the points inside it. `None` = automatic limits. |
| `show_lines` | `bool` | `True` | Draw and label the diagnostic lines of OB stars. |
| `show` | `bool` | `True` | Display the figure with `plt.show()`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `fig` | `Figure` | The `matplotlib` figure with the spectra and the check buttons. |

**Raises:** `TypeError` if `spec` or one of its items has the wrong type, `wave_lim` is not a pair of real numbers, or `show_lines` or `show` is not a boolean; `ValueError` if `spec` is an empty list, a spectrum is not valid or `wave_lim` is not a finite range with min < max. Spectrum files raise the exceptions of [`readspec`](#readspec).

<br>

### `rmvspec`
```python
rmvspec(
        spec        : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
        output_file : str | Path | None = None,
        wave_lim    : tuple[float, float] | None = None,
        autosave    : bool = True,
        show_lines  : bool = True
    ) -> tuple[np.ndarray[float], np.ndarray[float]]
```
Opens an interactive window to remove regions of a spectrum: spans are selected by dragging the left mouse button (shown in red) and removed with `r`; removed points stay visible in light gray, and every edit can be undone or redone. The cleaned spectrum is returned when the window is closed and, if `autosave` is `True` and `output_file` is given, also written then with [`writespec`](#writespec) (unless it has not changed since the last save). The zoom is kept after every removal, and the diagnostic lines of OB stars are drawn as in [`plotspec`](#plotspec). The keys of the tool take priority over Matplotlib's shortcuts, which keep working otherwise (`o` zoom, `p` pan, `q` close, arrows back / forward); no global Matplotlib setting is modified. Points with a non-finite wavelength or flux are discarded with a warning, and the spectrum is sorted by wavelength if needed. With a non-interactive backend a `RuntimeWarning` is issued, the figure is closed without saving and the spectrum is returned unmodified (inside `wave_lim`).

In the table below, `<stem>` is `output_file` without its extension (and `.gz`), e.g. `results/HD36822_clean` for `results/HD36822_clean.fits.gz`.

| Key | Action |
|---|---|
| mouse drag | Select a span (red). |
| `r` | Remove all selected spans. |
| `x` / `c` | Discard the last selected span / clear all selected spans. |
| `a` | Restore all removed points. |
| `u` / `U` | Undo the last edit / redo. |
| `g` | Toggle the grid. |
| `l` | Show / hide the OB-type star lines (replaces Matplotlib's log-scale toggle in this window). |
| `h` | Print the keyboard help. |
| `s` / `S` | Save the cleaned spectrum to `output_file` (with `writespec`) / and the figure as `<stem>.png` plus its `<stem>.pkl` copy (see [`savefig_pkl`](tex_README.md#savefig_pkl)). |
| `d` / `D` | Export the removed spans (merged) to `<stem>_rmvspans.txt` / the pending selected spans to `<stem>_selectedspans.txt` (columns `wave_min`, `wave_max`). |

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `spec` | `tuple[np.ndarray[float], np.ndarray[float]] \| list[np.ndarray[float]] \| str \| Path` | — | Spectrum to clean: `(wave, flux)` tuple, `[wave, flux]` list or spectrum file. |
| `output_file` | `str \| Path \| None` | `None` | File where the cleaned spectrum is written with [`writespec`](#writespec) (ASCII or FITS, chosen from its extension, optionally `.gz`); exported spans and figures are written next to it. Its directory must exist. `None` = nothing written. |
| `wave_lim` | `tuple[float, float] \| None` | `None` | Wavelength range `(min, max)` displayed and cleaned [units of the spectrum]; points outside it are excluded from the tool and from the returned spectrum. `None` = whole spectrum. |
| `autosave` | `bool` | `True` | Write the cleaned spectrum to `output_file` when the window is closed (if it changed since the last save). |
| `show_lines` | `bool` | `True` | Draw and label the diagnostic lines of OB stars (toggled with `l`). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `wave_clean` | `np.ndarray[float]` | Wavelengths of the kept points [units of the spectrum]. |
| `flux_clean` | `np.ndarray[float]` | Fluxes of the kept points. |

**Raises:** `TypeError` for wrong types; `IsADirectoryError`, `FileNotFoundError` or `NotADirectoryError` for an invalid `output_file`; `ValueError` if the spectrum or `wave_lim` are not valid, fewer than two finite points remain inside `wave_lim`, or the extension of `output_file` is not supported by `writespec` (checked before the window opens). Spectrum files raise the exceptions of [`readspec`](#readspec). Errors while saving inside the window (e.g. permissions, or repeated wavelengths, which `writespec` rejects) are printed, not raised.

<br>

### `normspec`
```python
normspec(
         spec        : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
         output_file : str | Path | None = None,
         wave_lim    : tuple[float, float] | None = None,
         function    : str = "spline3",
         order       : int = 1,
         sample      : str | list[tuple[float, float]] | tuple[float, float] | None = "*",
         naverage    : int = 1,
         low_reject  : float = 2.,
         high_reject : float = 0.,
         niterate    : int = 10,
         grow        : float = 1.,
         output_type : str = "ratio",
         replace     : bool = False,
         markrej     : bool = True,
         interactive : bool = True
    ) -> tuple[np.ndarray[float], np.ndarray[float] | None]
```
Fits the continuum of a spectrum as the IRAF `continuum` / `icfit` tasks do (Tody 1986, 1993): sample points are taken from the wavelength ranges in `sample`, combining `naverage` consecutive pixels into each point (groups never cross region borders); the function is least-squares fitted to them; points with residuals below −`low_reject` · rms or above `high_reject` · rms are rejected, together with the points within `grow` sample points of them, and the function is refitted, up to `niterate` times or until no new points are rejected (rejected points stay rejected; the rms uses N − 1 degrees of freedom). The output spectrum is then computed from the final continuum; ratios where the continuum is ≤ 0 are set to NaN, with a warning. The defaults (`low_reject = 2`, `high_reject = 0`) suit absorption-line spectra. If some spline pieces contain no sample points (e.g. gaps between sample regions), the knots are moved to quantiles of the sample points, with a warning, so that the Schoenberg-Whitney conditions hold (Schoenberg & Whitney 1953).

The output spectrum is written with [`writespec`](#writespec) (ASCII or FITS, chosen from the extension of `output_file`), replacing the previous file; the fit is described in `#` lines before the `wave` / `flux` column names (ASCII) or in `COMMENT` cards (FITS), and non-finite output fluxes are written as they are, with a warning. With `interactive=False` the continuum is fitted, saved (if `output_file` is given) and returned, which is useful for batch processing. With `interactive=True` a window with two panels opens: the spectrum, sample regions and points, rejected and deleted points and continuum (top), and the output spectrum (bottom); its status line and the header of `output_file` always describe the continuum shown, even if a later fit failed. Matplotlib's keyboard shortcuts are disabled in this window (the toolbar buttons still work). With a non-interactive backend a `RuntimeWarning` is issued, the initial automatic fit is returned (shown in the notebook with the inline backend of Jupyter) and nothing is written. With an interactive backend the function returns when the window is closed, also in Jupyter / IPython.

| Key | Action |
|---|---|
| `f` | Fit the continuum with the current parameters. |
| `+` / `-` | Increase / decrease the order. |
| `1` / `3` | Linear / cubic spline. |
| `l` / `c` | Legendre / Chebyshev polynomial. |
| `s` + `s` | Mark a sample region (press `s` at both ends); the first region replaces the whole-spectrum sample. |
| `t` | Reset the sample to the whole spectrum. |
| `d` / `u` | Delete / undelete the sample point nearest to the cursor. |
| left-click / right-click | Add a manual continuum point / remove it. |
| `enter` | Fit a cubic spline through the manual points (manual continuum, extrapolated outside them). |
| `r` | Reset to the initial parameters. |
| `w` | Save the output spectrum to `output_file` (with `writespec`). |
| `q` / `?` | Close the window / print the controls. |

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `spec` | `tuple[np.ndarray[float], np.ndarray[float]] \| list[np.ndarray[float]] \| str \| Path` | — | Spectrum to normalize: `(wave, flux)` tuple, `[wave, flux]` list or spectrum file. |
| `output_file` | `str \| Path \| None` | `None` | Output file, written with [`writespec`](#writespec) (ASCII or FITS, chosen from its extension, optionally `.gz`) with a header describing the fit; written with `w` in interactive mode and automatically otherwise. Its directory must exist. `None` = nothing written. |
| `wave_lim` | `tuple[float, float] \| None` | `None` | Wavelength range `(min, max)` normalized [units of the spectrum]; points outside it are excluded from the fit and from the returned spectrum. `None` = whole spectrum. |
| `function` | `str` | `"spline3"` | Continuum function: `"spline3"` (cubic spline), `"spline1"` (linear spline), `"legendre"` or `"chebyshev"`. |
| `order` | `int` | `1` | Number of spline pieces of equal length, or number of polynomial terms (degree + 1); ≥ 1. |
| `sample` | `str \| list[tuple[float, float]] \| tuple[float, float] \| None` | `"*"` | Wavelength ranges used in the fit [units of the spectrum]: `"*"` or `None` (whole spectrum), an IRAF range string (e.g. `"4000:4090,4150:4300"`), a list of `(min, max)` tuples or a single `(min, max)` tuple. |
| `naverage` | `int` | `1` | Consecutive pixels combined into each sample point: their average if positive, their median if negative; non-zero. |
| `low_reject` | `float` | `2.` | Rejection limit below the fit [rms units]; ≥ 0, 0 disables it. |
| `high_reject` | `float` | `0.` | Rejection limit above the fit [rms units]; ≥ 0, 0 disables it. |
| `niterate` | `int` | `10` | Maximum number of rejection iterations; ≥ 1. |
| `grow` | `float` | `1.` | Rejection growing radius [sample points; pixels if `naverage = 1`]; ≥ 0, only its integer part is used. |
| `output_type` | `str` | `"ratio"` | Output spectrum: `"ratio"` (flux / continuum), `"difference"` (flux − continuum), `"fit"` (continuum) or `"data"` (flux); unique abbreviations are accepted (e.g. `"r"`, `"di"`). |
| `replace` | `bool` | `False` | Replace the rejected and deleted pixels by the fit before computing the `"ratio"`, `"difference"` and `"data"` outputs. |
| `markrej` | `bool` | `True` | Mark the rejected points in the interactive window. |
| `interactive` | `bool` | `True` | Open the interactive window (`True`) or just fit, save and return (`False`). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `wave` | `np.ndarray[float]` | Wavelengths (inside `wave_lim`, if given) [units of the spectrum]. |
| `result` | `np.ndarray[float] \| None` | Output spectrum (`output_type`) of the last continuum fitted (automatic or manual); `None` if no continuum could be fitted in interactive mode. |

**Raises:** `TypeError` for wrong types; `IsADirectoryError`, `FileNotFoundError` or `NotADirectoryError` for an invalid `output_file`; `ValueError` for invalid values (unknown `function` or `output_type`, ambiguous abbreviation, `order` < 1, `naverage` = 0, negative rejection limits or `grow`, malformed `sample` or sample regions outside the spectrum, extension of `output_file` not supported by `writespec`), fewer than two points inside `wave_lim`, or, in non-interactive mode, a continuum that cannot be fitted (e.g. fewer sample points than free coefficients) or written (e.g. repeated wavelengths); `RuntimeError` if, in non-interactive mode, `output_file` cannot be written for any other reason (e.g. no write permission). Spectrum files raise the exceptions of [`readspec`](#readspec). Errors while saving inside the window are printed, not raised.

---
<br>

## 4. Dependencies

| Package | Used for |
|---|---|
| `astropy` | Reading and writing FITS files (`readspec`, `writespec`). |
| `matplotlib` | Interactive windows and figures (`plotspec`, `rmvspec`, `normspec`). |
| `numpy` | Spectra, ASCII output (`writespec`), sample points and numerical operations. |
| `scipy` | Least-squares splines and the cubic spline of the manual continuum (`normspec`). |
| `PyQt6` (optional) | Qt backend for the interactive windows (`rmvspec`, `normspec`, `plotspec`); installed with the `demos` extra (`pip install -e ".[demos]"`). Any other interactive Matplotlib backend also works. |

The figures saved by `rmvspec` (`S` key) use [`savefig_pkl`](tex_README.md#savefig_pkl) from the `tex` module.

---
<br>

## 5. Internal Structure

| Path | Content |
|---|---|
| `src/panorama/spec.py` | Public functions, in one section: `SPECTRA`. |
| `src/panorama/support/spec/spec_support.py` | Input validation shared by the public functions (`_check_array`, `_check_wave_lim`, `_prepare_spectrum`, `_load_spectrum`, `_check_fit_parameters`) and merging of overlapping wavelength ranges (`_merge_ranges`); keyboard handling of the interactive windows (`_take_keyboard`, `_release_keyboard`, `_check_interactive_backend`) and drawing of the OB-type star lines (`_draw_spectral_lines`); the supported extensions (`_ASCII_EXTENSIONS`, `_FITS_EXTENSIONS`) and the format of a spectrum file chosen from them (`_split_spectrum_name`, `_spectrum_file_format`, used by `readspec`, `writespec`, `rmvspec` and `normspec`); the readers of `readspec` (`_read_ascii`, `_read_fits`, `_read_fits_table`, `_read_fits_image`); the callback of the check buttons of `plotspec` (`_SpectrumToggle`, connected through its `_on_clicked` method); the window of `rmvspec` (`_InteractiveSpectrumCleaner`, run with `_run`); and the continuum fitting of `normspec` (`_parse_sample`, `_sample_points`, `_fit_function`, `_fit_continuum`, `_pixel_mask`, `_output_spectrum`; `_ContinuumFit`, whose `_fit`, `_fit_manual` and `_save` methods are also called by `normspec` and by the window; `_InteractiveNormalization`, run with `_run`). Every method of these classes is private (leading underscore). |
| `src/panorama/support/tex/tex_support.py` | Validation of the file arguments, shared with the `tex` and `obs` modules: input files (`_check_input_file`, used by `readspec`) and output files (`_check_output_file`, used by `writespec`, `rmvspec` and `normspec`). |
| `src/panorama/support/units/units_support.py` | Validation of scalar real numbers (`_check_real`, shared with the `units` and `obs` modules), used for wavelength limits and sample regions. |
| `src/panorama/support/spec/spec_OBlines.py` | Diagnostic lines of OB stars (`_OB_LINES`: rest-frame air wavelength [Å] and species of each line) and their colors (`_LINE_COLORS`), used by `plotspec` and `rmvspec`. |

**Adding a spectrum format.** Each format is read by a private function of `spec_support.py` that takes the path of the file and returns the raw `(wave, flux)` arrays; the generic checks (1D arrays, equal lengths, finite wavelengths, ordering) are done by `readspec` for every format. A new format needs:

1. its extensions, as a constant in the `INIT` section of `spec_support.py` (like `_ASCII_EXTENSIONS` and `_FITS_EXTENSIONS`);
2. a reading function (like `_read_ascii`), which raises `ValueError` for invalid content (any other exception is chained to a `RuntimeError` with the file name);
3. a branch for the new extensions in `_spectrum_file_format` and in `readspec`, and, if the files must also be written, in `writespec`.

Instrument-specific variants of a FITS layout (e.g. a pipeline that stores the spectrum in another way) are added to `_read_fits`, as another reading function tried before `_read_fits_table` and `_read_fits_image`.

---
<br>

## 6. Demo

The Jupyter notebook [`Demos/spec_demo.ipynb`](../Demos/spec_demo.ipynb) runs every function of the module in 5 examples, with three spectra read from `Demos/data_and_results/spec/data/`: a HERMES (Mercator) merged spectrum, not normalized, on a ln λ scale (`01179148_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fits`), and two normalized spectra from the IACOB database, one as a FITS image (`HD886_20160709_050356_M_V85000.fits`) and one as an ASCII file (`HD36822_20111109_015842_M_V85000_RV.ascii`). These spectra are distributed with the repository for demonstration purposes only; cite the HERMES spectrograph (Raskin et al. 2011) and the IACOB spectroscopic database if you use them.

* **Example 1:** reading the three formats with `readspec` (the ln λ scale of the HERMES file is detected from `log` in its name), with the wavelength coverage, sampling and median flux of each spectrum, and the warnings issued while reading (e.g. non-finite fluxes).
* **Example 2:** removal of bad points from the IACOB ASCII spectrum with `rmvspec`, between 4440 and 4600 Å (He I 4471, Mg II 4481, He II 4542 and the Si III 4553-4575 triplet); in interactive mode the cleaned spectrum is written to `HD36822_clean.txt` with `writespec` when the window is closed (`autosave`).
* **Example 3:** continuum normalization of the HERMES spectrum between 4000 and 5000 Å with `normspec`, with a set of recommended parameters (`NORM_PARAMS` in the setup block: cubic spline of 13 pieces, median of 10 pixels per sample point, rejection −2σ / +3σ, 10 iterations, growing radius of 1 sample point), each one printed with the reason for its value; the normalized spectrum is written to `HERMES_norm_4000-5000.txt` with `writespec` (in interactive mode, with `w`), and its header lines describing the fit are printed.
* **Example 4:** comparison with `plotspec` of the normalized HERMES spectrum (example 3), the IACOB FITS spectrum and the cleaned IACOB ASCII spectrum (example 2) in the same range, with the diagnostic lines of OB stars drawn and labeled; the returned figure is saved with `savefig_pkl` (PDF + pickle).
* **Example 5:** writing spectra with `writespec`: the normalized HERMES spectrum (example 3) as a gzip-compressed FITS table whose `COMMENT` cards repeat the header lines of `HERMES_norm_4000-5000.txt` (the fit written by `normspec`), the IACOB FITS image converted to ASCII (format conversion from a file), and a round trip with `readspec` that checks that the data are kept; then the particular cases of its Notes, with a mock 5-point spectrum: descending wavelengths (reversed), non-finite fluxes (kept, with a warning), an existing file without `overwrite=True` (`FileExistsError`), an unsupported extension and repeated wavelengths (`ValueError`).

**Interactive mode:** the flag `INTERACTIVE` in the first code cell selects how the notebook runs (restart the kernel after changing it). With `INTERACTIVE = False` (default) the inline backend is used and figures are shown in the notebook: `rmvspec` returns the spectrum unmodified with a `RuntimeWarning` and writes nothing, `normspec` runs in batch mode (`interactive=False`) and writes the automatic fit, and `plotspec` shows the comparison inline. With `INTERACTIVE = True` the Qt backend is selected (`%matplotlib qt`; a Qt binding such as PyQt5 / PyQt6 and a display are needed) and examples 2-4 open Matplotlib windows; `rmvspec` and `normspec` wait until their window is closed. Outputs are written to `Demos/data_and_results/spec/results/`. No internet access is needed.

```bash
pip install -e ".[demos]"     # PANORAMA + JupyterLab and PyQt6 (demo extras)
cd Demos
jupyter lab spec_demo.ipynb   # Or: jupyter notebook spec_demo.ipynb
```

---
<br>

## 7. References

* Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022, *ApJ*, 935, 167 — `astropy` (FITS input/output of `readspec` and `writespec`).
* Gray, R. O. & Corbally, C. J. 2009, *Stellar Spectral Classification* (Princeton University Press) — diagnostic lines of OB stars (`plotspec`, `rmvspec`).
* Greisen, E. W., Calabretta, M. R., Valdes, F. G. et al. 2006, *A&A*, 446, 747 — FITS spectral coordinates, logarithmic wavelength scales (`readspec`).
* Hunter, J. D. 2007, *CSE*, 9, 90 — `matplotlib`.
* Pence, W. D., Chiappetti, L., Page, C. G. et al. 2010, *A&A*, 524, A42 — the FITS standard (`readspec`, `writespec`).
* Raskin, G., van Winckel, H., Hensberge, H. et al. 2011, *A&A*, 526, A69 — the HERMES spectrograph (demo spectra).
* Schoenberg, I. J. & Whitney, A. 1953, *Trans. Amer. Math. Soc.*, 74, 246 — conditions for least-squares splines (`normspec`).
* Sota, A., Maíz Apellániz, J., Walborn, N. R. et al. 2011, *ApJS*, 193, 24 — spectral classification of O stars (`plotspec`, `rmvspec`).
* Tody, D. 1986, *Proc. SPIE*, 627, 733 — the IRAF system (`normspec`).
* Tody, D. 1993, *ASP Conf. Ser.*, 52, 173 — IRAF in the nineties (`normspec`).
* Virtanen, P., Gommers, R., Oliphant, T. E. et al. 2020, *Nature Methods*, 17, 261 — `scipy` (splines of `normspec`).
* Walborn, N. R. & Fitzpatrick, E. L. 1990, *PASP*, 102, 379 — spectral atlas of OB stars (`plotspec`, `rmvspec`).
