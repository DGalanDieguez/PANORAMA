# `panorama.prismas` — Model Atmospheres, SPAMMS Grid Checks and PRISMAS Archives

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](../LICENSE)

> Selection of the Kurucz and TLUSTY model atmospheres of a PRISMAS run, grid check of SPAMMS (grid entries needed and their model atmospheres), and renaming of the archives produced by PRISMAS in HPC mode.

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

`prismas` prepares and post-processes the runs of PRISMAS (`prismas.py`), the code that computes grids of synthetic spectra from model atmospheres, used to model (massive) stars with SPAMMS:

* **Model atmospheres:** selection of the models of a grid by effective temperature, surface gravity and composition, from the model catalogs distributed with PANORAMA, with the list of names written to a file (e.g. for an HTCondor submission) and the model files located and copied from a local copy of the grid.
* **SPAMMS:** grid check of SPAMMS for one or several input files, the grid entries needed by them, the model atmospheres that can be computed with PRISMAS and those outside the grid, and a figure of the needed entries on the grid.
* **Archives:** renaming of the `.tar.gz` archives produced by PRISMAS in HPC mode (`--opt 2`), named after their input model, to the names of modes 0 / 1, built from the stellar parameters.

<br>

**Conventions shared by the module:**
* **Model grids** are given by the same codes as `--mode` in `prismas.py` (case-insensitive): `'KC'` (Kurucz-Castelli ATLAS9 ODFNEW), `'KM'` (Kurucz-Mészáros ATLAS9 APOGEE) and `'T'` (TLUSTY OSTAR2002 / BSTAR2006). Models are identified by their names without file extension, as listed in the catalogs of `src/panorama/support/prismas/models_grids/`:
  * `'KC'`: `m0.50t03500g0.0am05ak2odfnew` (file `.dat`); an `a` before `k2` marks the alpha-enhanced models ([α/M] = +0.4).
  * `'KM'`: `amm03cm03om03t10000g20v20` (file `.mod`); abundance codes ending in 3 or 8 are quarter steps (`m03` = −0.25 dex).
  * `'T'`: `G30000g400v10` (OSTAR2002) and `BG20000g300v2` (BSTAR2006) (files `.7`, `.5` and `.nst`); the composition letter gives Z / Z☉ with the values of `prismas.py` (`C` = 2, `G` = 1, `L` = 1/2, `S` = 1/5, `T` = 1/10, `V` = 0.032, `W` = 1/50, `X` = 1/100, `Y` = 1/1000, `Z` = 0).
* **Stellar parameters** use the PANORAMA units: `teff` [K], `logg` = log₁₀(g / [cm s⁻²]), `met` = [M/H] [dex] for Kurucz or Z / Z☉ for TLUSTY, `alpha` = [α/M] [dex] and `carbon` = [C/M] [dex]. Parameters not defined for a grid are `NaN` (`alpha` and `carbon` for `'T'`, `carbon` for `'KC'`). Microturbulences are in km/s.
* **SPAMMS grid entries** and the **archives of modes 0 / 1** share the same names, with signed values written `m` (−) / `p` (+): `T{teff}_G{logg}_Z{Z/Zsun}_V{vmic}` for `'T'`, `T{teff}_G{logg}_M{±[M/H]}_a{±alpha}_V{vmic}` for `'KC'` (or the SPAMMS naming `T{teff}_G{logg}_M[a]{±[M/H]}_V{vmic}`) and `T{teff}_G{logg}_M{±[M/H]}_a{±alpha}_c{±carbon}_V{vmic}` for `'KM'` (e.g. `T30000_G4.00_Z1.000_V10`, `T49000_G4.50_Mm0.500_ap0.400_V2`).
* The **TLUSTY catalog** used for the selection (`T_OSTAR02+BSTAR06models.dat`) merges both subgrids: BSTAR2006 from 15000 to 28000 K and OSTAR2002 from 27500 to 55000 K (see [TLUSTY subgrids and their overlap](#tlusty-subgrids-and-their-overlap)). The full listings of both subgrids are only used for the figure of `needed_models`.
* Each catalog is read once per Python session and kept in memory (the first `'KM'` request, about 850000 models, takes a few seconds).
* Directory and file arguments accept a `str` or a `pathlib.Path`; returned paths are `pathlib.Path` objects. Output directories are never created: they must exist.
* All functions validate their inputs before computing (`TypeError` for wrong types, `FileNotFoundError` / `IsADirectoryError` / `NotADirectoryError` / `PermissionError` for paths, `ValueError` for invalid values, including empty or non-finite ones). Error messages name the meaning of the argument, the argument and, when relevant, its unit (e.g. `Effective temperature range 'teff' must satisfy min <= max. Got (35000.0, 25000.0).`).


<br>

### TLUSTY subgrids and their overlap

The TLUSTY grid (`'T'`) combines two subgrids that overlap in effective temperature:

| Subgrid | Teff range [K] | Teff step [K] | Model names |
|---|---|---|---|
| OSTAR2002 (Lanz & Hubeny 2003) | 27500-55000 | 2500 | `G30000g400v10` |
| BSTAR2006 (Lanz & Hubeny 2007) | 15000-30000 | 1000 | `BG20000g300v2` |

* **Overlap (27500-30000 K):** the OSTAR2002 models are recommended, so they are used from 27500 K, whenever a model exists in that subgrid.
* **Interpolation:** when a required Teff does not fall exactly on a model in SPAMMS, the spectra are weighted between the two closest models, Ti (lower) and Ts (upper). Both nodes must belong to the same subgrid, so that models of different grids are never mixed in the weighting.
* **Exception at 28000 K:** only for 27000 K < Teff < 27500 K (e.g. 27400 K), the nodes are the BSTAR2006 models at Ti = 27000 K and Ts = 28000 K. This is the only case in which a BSTAR2006 model is used inside the overlap. From 27500 K, the OSTAR2002 nodes are used (e.g. 27700 K or 28500 K between OSTAR2002 at 27500 and 30000 K).
* **Selection catalog:** for this reason, `T_OSTAR02+BSTAR06models.dat` lists BSTAR2006 from 15000 to 28000 K (including all the BSTAR2006 models at 28000 K) and OSTAR2002 from 27500 to 55000 K. The other BSTAR2006 models of the overlap (29000 and 30000 K) are never needed and are left out: from 27500 K the OSTAR2002 nodes are used, so SPAMMS never requests 29000 K.
* **Subgrid of each SPAMMS entry** (`needed_models`): OSTAR2002 for Teff ≥ 27500 K except 28000 K, and BSTAR2006 otherwise (also at 28000 K); if the model only exists in the other subgrid, that one is used.

| Teff [K] | Subgrid used |
|---|---|
| 15000-27000 | BSTAR2006 |
| 27500 | OSTAR2002 |
| 28000 | BSTAR2006 (upper node only for 27000 K < Teff < 27500 K) |
| 29000 | — (never requested by SPAMMS; not in the selection catalog) |
| 30000-55000 | OSTAR2002 |

---
<br>

## 2. Function Summary

### Model atmospheres (`MODELS`)

| Function | Description |
|---|---|
| [`find_models`](#find_models) | Models of a grid within the requested ranges of the stellar parameters, optionally copied from a local copy of the grid. |

<br>

### SPAMMS (`SPAMMS`)

| Function | Description |
|---|---|
| [`needed_models`](#needed_models) | Grid entries needed by SPAMMS, their model atmospheres and a figure of the entries on the grid. |

<br>

### Archives (`ARCHIVES`)

| Function | Description |
|---|---|
| [`rename_tars`](#rename_tars) | Rename the archives of PRISMAS HPC mode to the names of modes 0 / 1. |

---
<br>

## 3. Function Reference

Each function lists its arguments in an **Inputs** table and its outputs in a **Returns** table. For functions that return a dictionary (or a list of dictionaries), the first row is the container itself and the following rows are its keys (keys of a nested container are written as `'outer'[i]['key']`).

### `find_models`
```python
find_models(
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
    ) -> list[dict]
```
Searches the model catalog of a grid for the models that fulfill every criterion and returns them sorted by (`teff`, `logg`, `met`, `alpha`, `carbon`). Each criterion accepts `None` (no filter), a single value (`25000`), a list of discrete values (`[25000, 27500]`), a closed range (`(25000, 35000)`) or an open range with any limit set to `None` (`(None, 3.5)`). Comparisons use an absolute tolerance of 10⁻⁶, so the fractional TLUSTY metallicities can also be selected with a range (e.g. `met=(0.03, 0.04)` for Z / Z☉ = 0.032). A criterion not defined for the grid (e.g. `alpha` for `'T'`) must be `None`. With `modeldir`, the files of each model (its name followed by each extension of the grid: `.dat` for `'KC'`, `.mod` for `'KM'`, and `.7`, `.5` and `.nst` for `'T'`) are located in that directory (and in its subdirectories with `recursive`, visited in alphabetical order) and, with `output_dir`, the complete models (every file of the grid found, as `prismas.py` needs them all) are copied there. Models without files, TLUSTY models with only some of their files (not copied) and file names repeated in several subdirectories (the first one is used) are each reported in a single warning. Typical use: writing the list of models of a PRISMAS run for an HTCondor submission (`output_file`), or copying them from the full grid to a working directory (`modeldir` + `output_dir`).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `grid` | `str` | — | Model grid: `'KC'`, `'KM'` or `'T'` (case-insensitive). |
| `teff` | `int \| float \| list[float] \| tuple[float \| None, float \| None] \| None` | `None` | Effective temperature selection [K]; positive values. |
| `logg` | `int \| float \| list[float] \| tuple[float \| None, float \| None] \| None` | `None` | Surface gravity selection [log₁₀(g / [cm s⁻²])]. |
| `met` | `int \| float \| list[float] \| tuple[float \| None, float \| None] \| None` | `None` | Metallicity selection: [M/H] [dex] for `'KC'` / `'KM'`, Z / Z☉ (non-negative) for `'T'`. |
| `alpha` | `int \| float \| list[float] \| tuple[float \| None, float \| None] \| None` | `None` | Alpha enhancement [α/M] selection [dex]; only for `'KC'` and `'KM'`. |
| `carbon` | `int \| float \| list[float] \| tuple[float \| None, float \| None] \| None` | `None` | Carbon abundance [C/M] selection [dex]; only for `'KM'`. |
| `modeldir` | `str \| Path \| None` | `None` | Existing directory with the model files (`--modeldir` of `prismas.py`); `None` = only the catalog is searched. |
| `recursive` | `bool` | `False` | Also search the subdirectories of `modeldir`. |
| `output_dir` | `str \| Path \| None` | `None` | Existing, writable directory, different from `modeldir` (and not inside it with `recursive`), where the files of the complete models are copied (existing files are overwritten); requires `modeldir`. |
| `output_file` | `str \| Path \| None` | `None` | Text file where the names of the selected models are written, one per line (its directory must exist). |
| `verbose` | `bool` | `True` | Print a summary table (first 50 models), with the files found in `modeldir`. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `models` | `list[dict]` | Selected models, sorted by (`teff`, `logg`, `met`, `alpha`, `carbon`); each element is a dictionary with the keys below. |
| `[i]['model']` | `str` | Model name, without file extension. |
| `[i]['grid']` | `str` | Model grid code (`'KC'`, `'KM'` or `'T'`). |
| `[i]['subgrid']` | `str` | `'OSTAR2002'` or `'BSTAR2006'` for `'T'`; `''` for Kurucz. |
| `[i]['teff']` | `float` | Effective temperature [K]. |
| `[i]['logg']` | `float` | Surface gravity [log₁₀(g / [cm s⁻²])]. |
| `[i]['met']` | `float` | [M/H] [dex] for `'KC'` / `'KM'`; Z / Z☉ for `'T'`. |
| `[i]['alpha']` | `float` | Alpha enhancement [α/M] [dex] (`NaN` for `'T'`). |
| `[i]['carbon']` | `float` | Carbon abundance [C/M] [dex] (`NaN` except for `'KM'`). |
| `[i]['files']` | `list[Path]` | Files of the model found in `modeldir`, in the order of the extensions of the grid (copied to `output_dir` only if every file of the model is found); `[]` if `modeldir` is `None` or no file is found. |

**Raises:** `TypeError` for wrong types; `ValueError` for an unsupported grid, invalid criteria (non-finite values, empty lists, ranges without limits or with min > max, criteria not defined for the grid), or `output_dir` given without `modeldir`, equal to it or inside it with `recursive`; `FileNotFoundError`, `NotADirectoryError`, `IsADirectoryError` or `PermissionError` for invalid paths or a catalog not installed; `RuntimeError` if the catalog contains no valid model name.

<br>

### `needed_models`
```python
needed_models(
              spamms_inputs : str | Path | list[str | Path],
              spamms_dir    : str | Path,
              grid          : str,
              output_name   : str | Path = "spamms",
              show          : bool = False,
              verbose       : bool = True
    ) -> dict[str, list | dict | Path | Figure | None]
```
Runs the grid check of SPAMMS (`python <script> -i <input file> -c -n <n_cores>`, with the Python interpreter of the current session, from `spamms_dir`, using every core allocated to the process: all the cores of the computer, or only those assigned to the job by HTCondor / Slurm) for each input file, collects the grid entries that SPAMMS reports as missing, and matches them with the model atmospheres of the grid (same catalogs and decoding as [`find_models`](#find_models)). The SPAMMS script is chosen from the name of each input file (without extension, case-insensitive): `spamms_disk.py` if it ends in `_disk`, `spamms_noreflec.py` if it ends in `noreflec`, and `spamms.py` otherwise; only the scripts needed by the input files must exist in `spamms_dir`. Input files outside `spamms_dir` are copied there for the run and always removed afterwards, also if SPAMMS fails, the copy is partial or the run is interrupted (an input file already there with the same content is used and kept; a different file with the same name is never overwritten: that run fails). The progress is printed as one line per input file; after all the runs, the warnings printed by SPAMMS and the reason of each failed run (exit status and last error line of SPAMMS; the input file is skipped) are printed as aligned blocks with `verbose`, and issued as one Python warning per input file otherwise. If SPAMMS fails for every input file, the `RuntimeError` lists the reason of each failure. Entries with Teff = 0 or 1000 K are dummy entries (disk patches of SPAMMS): they are neither matched nor plotted, but are listed in the `_missing.dat` file with model `dummy`. For TLUSTY, OSTAR2002 is used for Teff ≥ 27500 K and BSTAR2006 otherwise, falling back to the other subgrid, except at 28000 K, taken from BSTAR2006 as the upper node of the interpolation for 27000 K < Teff < 27500 K (see [TLUSTY subgrids and their overlap](#tlusty-subgrids-and-their-overlap)). The figure shows the grid nodes (for `'T'`, the full OSTAR2002 and BSTAR2006 listings) and the needed entries on top of them, with larger markers (black crosses in the grid, red crosses surrounded by a red circle outside it), with reversed axes (hot / low gravity to the left / top) and the legend above the axes; it is returned open, to be saved e.g. with [`savefig_pkl`](tex_README.md#savefig_pkl).

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `spamms_inputs` | `str \| Path \| list[str \| Path]` | — | SPAMMS input files: a directory (every `*.txt` file in it, in alphabetical order), a single file, or a list (or tuple) of files with unique names. |
| `spamms_dir` | `str \| Path` | — | Directory where SPAMMS is installed (`spamms.py`, `spamms_disk.py`, `spamms_noreflec.py`); SPAMMS is run from it. |
| `grid` | `str` | — | Model grid of the input files: `'KC'`, `'KM'` or `'T'` (case-insensitive). |
| `output_name` | `str \| Path` | `"spamms"` | Path and base name of the two output files, `<output_name>_<grid>_tocompute.dat` and `<output_name>_<grid>_missing.dat` (a file extension is removed; their directory must exist and existing files are overwritten). |
| `show` | `bool` | `False` | Display the figure with `plt.show()` (the figure is returned in any case). |
| `verbose` | `bool` | `True` | Print the progress of the SPAMMS runs, the SPAMMS warnings and failures, and a summary (if `False`, the warnings and failures are issued as Python warnings). |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `results_dict` | `dict[str, list \| dict \| Path \| Figure \| None]` | Dictionary with the keys below. |
| `'needed'` | `list[str]` | Alphabetically sorted grid entries needed by all the input files, without duplicates (e.g. `T30000_G4.00_Z1.000_V10`). |
| `'per_file'` | `dict[str, list[str]]` | `{input file name: [grid entries]}` (`[]` if none is missing); files whose run failed are not included. |
| `'models'` | `list[dict]` | Needed entries whose model is in the grid (they can be computed with PRISMAS); each element is a dictionary with the keys below. |
| `'models'[i]['name']` | `str` | SPAMMS grid entry. |
| `'models'[i]['model']` | `str` | Model atmosphere, without file extension (e.g. `G30000g400v10`). |
| `'models'[i]['subgrid']` | `str` | `'OSTAR2002'` or `'BSTAR2006'` for `'T'`; `''` for Kurucz. |
| `'models'[i]['teff']` | `float` | Effective temperature [K]. |
| `'models'[i]['logg']` | `float` | Surface gravity [log₁₀(g / [cm s⁻²])]. |
| `'models'[i]['met']` | `float` | [M/H] [dex] for `'KC'` / `'KM'`; Z / Z☉ for `'T'`. |
| `'models'[i]['alpha']` | `float` | Alpha enhancement [α/M] [dex] (`NaN` for `'T'`). |
| `'models'[i]['carbon']` | `float` | Carbon abundance [C/M] [dex] (`NaN` except for `'KM'`). |
| `'models'[i]['vmic']` | `float` | Microturbulence of the synthesis [km/s]. |
| `'not_in_grid'` | `list[str]` | Entries without model atmosphere in the catalog of the grid. |
| `'dummy'` | `list[str]` | Dummy entries (disk patches, Teff = 0 or 1000 K), neither matched nor plotted, but listed in `'missing_file'` (model `dummy`). |
| `'invalid'` | `list[str]` | Entries that do not follow the SPAMMS naming of the grid (also reported in a warning). |
| `'tocompute_file'` | `Path` | `<output_name>_<grid>_tocompute.dat`: model atmospheres to compute with PRISMAS (entries inside the grid), one per line, sorted, without duplicates or header (e.g. the input of an HTCondor submission). |
| `'missing_file'` | `Path` | `<output_name>_<grid>_missing.dat`: every needed entry inside or outside the grid and the dummy entries, sorted by increasing Teff (then log g, composition and microturbulence; not alphabetically, so that `T3500` comes before `T10000`), two columns (SPAMMS grid entry and model atmosphere, `-` outside the grid, `dummy` for dummy entries), with `#` comment lines. |
| `'fig'` | `Figure \| None` | Figure of the needed entries on the grid (`None` if no entry is needed). |

**Raises:** `TypeError` for wrong types; `ValueError` for an unsupported grid, an empty `output_name`, no input file or repeated input file names; `FileNotFoundError`, `NotADirectoryError`, `IsADirectoryError` or `PermissionError` for invalid paths, a missing SPAMMS script or a catalog not installed; `RuntimeError` if SPAMMS fails for every input file (with the reason of each failure).

<br>

### `rename_tars`
```python
rename_tars(
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
    ) -> list[dict]
```
In HPC mode (`--opt 2`), `prismas.py` names each archive after its input model (e.g. `G30000g400v10.tar.gz`), without the microturbulence of the synthesis; in modes 0 / 1 the name is built from the stellar parameters (e.g. `T30000_G4.00_Z1.000_V10.tar.gz`). This function renames the former to the latter, decoding the HPC names with the same decoders as the model catalogs. The microturbulence cannot be recovered from the HPC names (in TLUSTY names, `v` is the microturbulence of the model atmosphere), so it is a required input. With `rename_content`, the inner folder and the five `.npy` files (`{name}/{name}_wave.npy`, ...) are renamed too, so the result is identical to a mode 0 / 1 output of `prismas.py`: the archive is streamed member by member (nothing is extracted to disk), written to a temporary `.part` file and moved into place only when complete; the whole gzip stream is read, so a corrupted or truncated archive (e.g. a killed HTCondor job) is reported as `'error'` and never replaces anything. The archives are processed in parallel threads; with `ncores=None` the number of threads tunes itself (it starts at 2 and doubles while the throughput improves by ≥ 10 %, up to the cores allocated to the process: CPU affinity, `SLURM_CPUS_PER_TASK` / `OMP_NUM_THREADS`). The gzip layer uses ISA-L (package `isal`, a required dependency), which compresses much faster than `zlib`. TLUSTY archives of different subgrids with the same Teff, log g and Z / Z☉ (e.g. `G30000g300v10` and `BG30000g300v2`, in the overlap of OSTAR2002 and BSTAR2006) would get the same name: both are kept, a warning is issued and the subgrid tag is appended to each name (`T30000_G3.00_Z1.000_V10_OSTAR02.tar.gz` and `T30000_G3.00_Z1.000_V10_BSTAR06.tar.gz`), so that neither model replaces the other; archives of the same subgrid with the same target (e.g. the same file in two subdirectories with `recursive`) are `'skipped'`, keeping the first one. Archives already named as in modes 0 / 1 (also with a subgrid tag) are reported as `'done'`, so the function can be run several times on the same directory. A summary table is printed.

**Inputs:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tardir` | `str \| Path` | — | Directory with the `.tar.gz` archives of HPC mode (`--numpydir` of `prismas.py`). |
| `grid` | `str` | — | Model grid: `'KC'`, `'KM'` or `'T'` (case-insensitive). |
| `vmic` | `int \| float` | — | Microturbulence of the synthesis (`--vmic` of `prismas.py`) [km/s]; non-negative, only its integer part is used. |
| `recursive` | `bool` | `False` | Also search the subdirectories of `tardir`. |
| `output_dir` | `str \| Path \| None` | `None` | Existing directory, different from `tardir`, where the renamed archives are written; `None` = each archive is renamed in its own directory. |
| `rename_content` | `bool` | `True` | Also rename the inner folder and `.npy` files (if `False`, only the `.tar.gz` file is renamed: fast, but its content keeps the HPC names). |
| `keep_original` | `bool` | `False` | Keep the original archive (if `False`, it is removed once the renamed archive has been written). |
| `overwrite` | `bool` | `False` | Overwrite renamed archives that already exist. |
| `dry_run` | `bool` | `False` | Only report the renaming, without touching any file. |
| `ncores` | `int \| None` | `None` | Number of archives processed in parallel (threads), ≥ 1; `None` = self-tuned. |

**Returns:**

| Return | Type | Description |
|---|---|---|
| `archives` | `list[dict]` | Archives of the grid found in `tardir` (HPC names and names already of modes 0 / 1), sorted by (`teff`, `logg`, `met`, `alpha`, `carbon`); each element is a dictionary with the keys below. Other files are ignored. |
| `[i]['source']` | `Path` | Original archive. |
| `[i]['target']` | `Path` | Renamed archive (the source itself for `'done'`). |
| `[i]['teff']` | `float` | Effective temperature [K]. |
| `[i]['logg']` | `float` | Surface gravity [log₁₀(g / [cm s⁻²])]. |
| `[i]['met']` | `float` | [M/H] [dex] for `'KC'` / `'KM'`; Z / Z☉ for `'T'`. |
| `[i]['alpha']` | `float` | Alpha enhancement [α/M] [dex] (`NaN` for `'T'`). |
| `[i]['carbon']` | `float` | Carbon abundance [C/M] [dex] (`NaN` except for `'KM'`). |
| `[i]['status']` | `str` | `'renamed'`, `'dry_run'` (would be renamed), `'done'` (already named as in modes 0 / 1), `'exists'` (target exists and `overwrite` is `False`), `'skipped'` (same target as another archive of the same subgrid) or `'error'` (unreadable archive or unwritable target). |
| `[i]['reason']` | `str` | Explanation for `'exists'`, `'skipped'` and `'error'`, and for the archives renamed with a subgrid tag (the other archive with the same name); `''` otherwise. |

**Raises:** `TypeError` for wrong types; `ValueError` for an unsupported grid, a negative or non-finite `vmic`, `ncores` < 1 or `output_dir` equal to `tardir`; `FileNotFoundError`, `NotADirectoryError` or `PermissionError` for invalid directories, or if the directory of the renamed archives cannot be written (unless `dry_run`).

---
<br>

## 4. Dependencies

| Package | Used for |
|---|---|
| `matplotlib` | Figure of the needed grid entries (`needed_models`). |
| `numpy` | Catalog columns and selection masks. |
| `tqdm` | Progress bar of the renaming (`rename_tars`). |
| `isal` | Fast gzip compression and decompression of the archives with ISA-L (`rename_tars`). |

`needed_models` also needs a local installation of SPAMMS (not a Python dependency of PANORAMA), run with the Python interpreter of the current session. All the functions work offline.

---
<br>

## 5. Internal Structure

| Path | Content |
|---|---|
| `src/panorama/prismas.py` | Public functions, grouped in three sections: `MODELS`, `SPAMMS`, `ARCHIVES`. |
| `src/panorama/support/prismas/prismas_support.py` | Constants (`_CATALOG_CACHE`, `_PARAM_TOL`, `_SPAMMS_SCRIPTS`, `_SPAMMS_PHRASE`, `_DUMMY_TEFF`, `_TAR_BUFFER`) and helpers, in nine groups. `GRIDS`: registry of the supported grids, with their catalogs, naming patterns, decoders, criteria, file extensions and archive names and subgrid tags (`_grid_registry`). `VALIDATION`: grid codes (`_check_grid`), directories (`_check_directory`, shared by the three public functions, `_locate_files` and `_run_spamms`) and number of parallel workers (`_check_ncores`). `CATALOGS`: reading, decoding and caching of the model catalogs (`_read_catalog`). `DECODING`: model names (`_decode_model`, `_decode_kc`, `_decode_km`, `_decode_tlusty`, with `_signed`, `_km_value` and `_tlusty_metallicity`), SPAMMS grid entries and names of modes 0 / 1 (`_decode_spamms`, `_spamms_kc`, `_spamms_km`, `_spamms_tlusty`) and preferred TLUSTY subgrid (`_tlusty_subgrid_order`). `SELECTION`: validation and application of the selection criteria (`_check_criterion`, `_passes`), and matching and sorting keys (`_match_key`, `_sort_key`). `SPAMMS`: SPAMMS grid check of one input file (`_run_spamms`). `FIGURES`: needed entries on the grid nodes (`_plot_needed_models`). `FILES`: location of the model files (`_locate_files`). `ARCHIVES`: names of modes 0 / 1 (`_archive_name`), allocated cores (`_allocated_cores`, also used for the `-n` option of SPAMMS in `_run_spamms`) and parallel runs (`_run_parallel`), gzip backend (`_gzip_backend`) and streaming renaming of an archive (`_rename_member`, `_rewrite_tar`, `_rename_archive`). `TEXT`: formatting of the summary tables (`_format_value`, `_short_list`) and of the SPAMMS warnings and failures (`_message_block`). |
| `src/panorama/support/prismas/models_grids/` | Model catalogs (one model name per line, without extension): `KC_models.dat`, `KM_models.dat`, `T_OSTAR02+BSTAR06models.dat` (selection catalog of TLUSTY: BSTAR2006 from 15000 to 28000 K and OSTAR2002 from 27500 to 55000 K) and the full TLUSTY listings `T_OSTAR02_models.dat` and `T_BSTAR06_models.dat` (figure of `needed_models`). |
| `src/panorama/support/tex/tex_support.py` | Validation of input and output files (`_check_input_file`, `_check_output_file`) and reading of the catalog files (`_read_object_names`), reused by the three public functions. |
| `src/panorama/support/units/units_support.py` | Validation of real numbers (`_check_real`), reused for the selection criteria, `vmic` and the TLUSTY temperatures. |
| `src/panorama/support/obs/obs_support.py` | Validation of boolean flags (`_check_bool`). |
| `src/panorama/support/spec/spec_support.py` | Validation of NumPy arrays (`_check_array`), reused by `_passes` and `_plot_needed_models`. |

---
<br>

## 6. Demo

The Jupyter notebook [`Demos/prismas_demo.ipynb`](../Demos/prismas_demo.ipynb) runs every function of the module in 6 examples, one section per example:

* **Examples 1-3:** TLUSTY models of a PRISMAS run selected by a range of temperatures and a list of gravities (BSTAR2006 up to 28000 K and OSTAR2002 from 27500 K, names written to a file) and a fractional Z / Z☉ selected with a range; alpha-enhanced Kurucz-Castelli models and Kurucz-Mészáros models with a quarter-step metallicity and an open range in [C/M]; and the checks of the selection criteria (a criterion not defined for the grid, an inverted range, a negative Z / Z☉, a range without limits, an unknown grid and a copy without `modeldir`).
* **Example 4:** copy of the selected models from a *mock* local grid (empty files in subdirectories): only the complete OSTAR2002 model is copied, while a BSTAR2006 model with only one of its three files (repeated in another subdirectory) and a model without files are reported in their warnings and not copied.
* **Example 5:** SPAMMS grid check of `input_SPAMMS.txt` (run with `spamms.py`), with the model atmospheres to compute, the content of the two output files (`_tocompute.dat`, and `_missing.dat` with every needed entry sorted by increasing Teff, its model, `-` outside the grid or `dummy`) and the figure of the needed entries on the grid (shown inline and saved with `savefig_pkl`). It uses the SPAMMS installation given by the environment variable `SPAMMS_DIR` (e.g. `export SPAMMS_DIR=/home/user/SPAMMS` before launching Jupyter; `SPAMMS_GRID` in the setup cell sets the grid of the input file) or, if it is not set, a *mock* SPAMMS script written to `results/mock_spamms/`, which issues a warning and prints a fixed list of missing grid points (entries inside the BSTAR2006 grid, two outside it and a dummy entry).
* **Example 6:** renaming of *mock* TLUSTY archives of HPC mode (folder and five `.npy` files each): a dry run and the renaming with the content renamed too and 2 parallel workers (a truncated archive reported as `'error'`, two archives of different subgrids with the same parameters, renamed with their subgrid tags and reported in a warning, an archive already renamed and a file that is not a model archive).

The SPAMMS input file is read from `Demos/data_and_results/prismas/data/`, and the outputs (list of selected models, mock and copied model files, mock and renamed archives, the mock SPAMMS script, and the output files and the figure of the grid check, saved as PDF + pickle) are written to `Demos/data_and_results/prismas/results/`; the mock model files, archives and SPAMMS script are rebuilt at every run. No internet access is needed. The notebook takes a few seconds (mostly the first reading of the `'KM'` catalog).

```bash
pip install -e ".[demos]"         # PANORAMA + JupyterLab and PyQt6 (demo extras)
cd Demos
jupyter lab prismas_demo.ipynb    # Or: jupyter notebook prismas_demo.ipynb
```

---
<br>

## 7. References

* Abdul-Masih, M., Sana, H., Conroy, K. E. et al. 2020, *A&A*, 636, A59 — SPAMMS (Spectroscopic PAtch Model for Massive Stars).
* Castelli, F. & Kurucz, R. L. 2003, *IAUS*, 210, A20 — Kurucz-Castelli ATLAS9 ODFNEW model atmospheres (`'KC'`).
* Galán-Diéguez, D., Abdul-Masih, M., Allende Prieto, C. et al. 2026, arXiv:2609.15331 — PRISMAS and its use with SPAMMS.
* Hunter, J. D. 2007, *CSE*, 9, 90 — `matplotlib`.
* Lanz, T. & Hubeny, I. 2003, *ApJS*, 146, 417 — TLUSTY OSTAR2002 grid (`'T'`).
* Lanz, T. & Hubeny, I. 2007, *ApJS*, 169, 83 — TLUSTY BSTAR2006 grid (`'T'`).
* Mészáros, Sz., Allende Prieto, C., Edvardsson, B. et al. 2012, *AJ*, 144, 120 — Kurucz-Mészáros ATLAS9 APOGEE model atmospheres (`'KM'`).
