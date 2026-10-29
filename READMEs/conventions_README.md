# PANORAMA — Contributing and Code Conventions

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](../LICENSE)
> How to add a module to PANORAMA and the code, documentation, demo and check conventions shared by every module.

[← Back to the main README](../README.md)

---
<br>

## Table of Contents
- [1. Adding a Module](#1-adding-a-module)
- [2. Code Conventions](#2-code-conventions)
  - [2.1. Language](#21-language)
  - [2.2. Public vs. Internal Code](#22-public-vs-internal-code)
  - [2.3. Module Layout](#23-module-layout)
  - [2.4. Function Signatures](#24-function-signatures)
  - [2.5. Function Body](#25-function-body)
  - [2.6. Imports](#26-imports)
  - [2.7. Code Style](#27-code-style)
  - [2.8. Files and Directories](#28-files-and-directories)
  - [2.9. Input Validation](#29-input-validation)
  - [2.10. Docstrings](#210-docstrings)
  - [2.11. Returned Results](#211-returned-results)
  - [2.12. Figures](#212-figures)
  - [2.13. Scientific Consistency](#213-scientific-consistency)
  - [2.14. Module READMEs](#214-module-readmes)
  - [2.15. Demos](#215-demos)
  - [2.16. Checks](#216-checks)

---
<br>

## 1. Adding a Module
1. Create `src/panorama/<module>.py` with its public functions. It is detected automatically by `__init__.py` and its functions are exposed as `panorama.<function>`: no other code needs to be changed.
2. Put its internal code in `src/panorama/support/<module>/<module>_support.py` (with an `__init__.py` in `support/<module>/`) and its data files in `src/panorama/support/<module>/`. Data files are installed with the package if their extension is listed under `[tool.setuptools.package-data]` in `pyproject.toml`.
3. Add any new third-party dependency (only with a technical justification) to `dependencies` in `pyproject.toml`, and to the *Dependencies* table of the module README. Dependencies used only by the demos (e.g. JupyterLab, PyQt6 for the interactive windows) go in the optional `demos` extra (`pip install -e ".[demos]"`) instead.
4. Write `READMEs/<module>_README.md`, following the structure of the existing module READMEs (see [Section 2.14](#214-module-readmes)).
5. Add a row to the table of [Section 4 of the main README](../README.md#4-modules).
6. Add the Jupyter notebook `Demos/<module>_demo.ipynb`, with its input data in `Demos/data_and_results/<module>/data/` and its outputs in `Demos/data_and_results/<module>/results/` (see [Section 2.15](#215-demos)).

Function names must be unique across modules, since all of them share the `panorama` namespace (a warning is issued on import otherwise). Modules whose name starts with an underscore are not exposed.

---
<br>

## 2. Code Conventions
The following conventions are mandatory for every function: public functions, internal helpers under `support/`, class methods and nested functions. When existing code is modified, they apply to all of the delivered code, not only to the new lines. The `tex` module (`tex.py`, `tex_README.md` and `tex_demo.ipynb`) is the reference implementation: cases not covered below follow what it does.

<br>

### 2.1. Language
* All code, comments, docstrings, variable names, error messages, READMEs and demos must be written in English, using American spelling consistently (*modeling*, *color*, *normalization*, *center*).

<br>

### 2.2. Public vs. Internal Code
* Only the functions that belong to the package API are public: those exported by `panorama` (`from panorama import ...`) and documented in the *Function Reference* of the module README. **Everything else always starts with an underscore**, even inside something that is already private: the helpers and classes of `support/` (`_check_input_file`, `_ContinuumFit`), **every method of the private classes** (also when it is called from outside the class, e.g. `_ContinuumFit._save`), the methods connected as Matplotlib or other callbacks (`_on_key`, `_on_click`), nested functions (`_normalize`), and class and static methods.
* The only names without an underscore outside the API are Python special methods (`__init__`, `__call__`, ...), methods that override a method of an external base class and are called by name by the library (e.g. `do_GET` and `log_message` of `http.server.BaseHTTPRequestHandler`, marked with a comment such as `# Overrides BaseHTTPRequestHandler.do_GET; name fixed by the library`), and names required by external tools (e.g. pytest `test_*` functions). A double leading underscore (`__name`) is never used outside special methods, since it triggers Python's name mangling.
* Renaming a private function or method to follow this rule does not change the public API, but every call in the package (including other `support/` files) and the *Internal Structure* table of the README must be updated. Demos and usage READMEs never call private functions, classes or methods: they only appear in the *Internal Structure* table.
* Everything under `support/` is internal and may change without notice.
* **Helpers always live in `support/`.** The files of `src/panorama/` (the public modules, e.g. `tex.py`, `spec.py`, and the package `__init__.py`) define only public functions: no private (`_name`) function or class, no nested function and no class inside a function. Every helper is a module-level function (or a private class, when it needs shared state or callbacks) of `src/panorama/support/<module>/<module>_support.py`, imported in the `PANORAMA` sub-subsection of `MODULES & PACKAGES` of the function that uses it; what a former nested function captured from its enclosing function is passed as an explicit argument. The private helpers of the package `__init__.py` live in `src/panorama/support/package/package_support.py`. Nested functions are only allowed inside the functions and methods of `support/`, for logic that makes no sense outside them. In the public modules, `lambda` is only allowed as the argument of a trivial one-expression call (e.g. `key=lambda row: row['ra']`), never assigned to a name.
* Logic shared by several functions (e.g. the normalization of object names, the validation of an argument or an astronomical computation) must be moved to a helper in `src/panorama/support/<module>/<module>_support.py` instead of being duplicated, and the *Internal Structure* table of the module README must be updated. Existing helpers of other modules are reused rather than rewritten: `_check_input_file` / `_check_output_file` (`support/tex/tex_support.py`) for file arguments and `_check_real` (`support/units/units_support.py`) for scalar real numbers.
* The existing API (function names, signatures, dictionary keys and output formats) must be preserved unless a change is explicitly required. PANORAMA supports Python ≥ 3.10.

<br>

### 2.3. Module Layout
* Every module starts with a module docstring describing its purpose, followed, after a blank line, by the list `Conventions shared by every function in this module:` (`... every object ...` in data-only modules such as `spec_OBlines.py`).
* The docstring is followed by `from __future__ import annotations`. It is required because imports live inside the functions, so annotations such as `Path`, `Figure` or `np.ndarray` must not be evaluated at definition time.
* Exactly one blank line separates `from __future__ import annotations` from the `INIT` section or the first group title.
* Any other module-level import goes in an `INIT` section (a group title), under a `MODULES & PACKAGES` header, with an explanatory string stating why it cannot live inside the function. Module-level constants also go in `INIT` (also in data-only modules such as `spec_OBlines.py`), each block under its own level-1 header.
* Functions are grouped by topic (e.g. `BIB`, `QUERY`, `MATPLOTLIB` in `tex.py`). Each group starts with a three-line title in uppercase, with exactly 30 `=` on each side of the name, and top and bottom lines of the same length as the central one:

```python
# =================================================================
# ============================== BIB ==============================
# =================================================================
```

* Every module-level `def` or `class` is immediately preceded by the separator `# ==============================` (30 `=`). Consecutive functions are separated by four blank lines; a group title is preceded by four blank lines (except the first one of the module, see above) and followed by one.
* Nested functions carry no separator; consecutive nested functions are separated by a line containing only `#`.

<br>

### 2.4. Function Signatures
Every function, without exception, uses this layout:

```python
# ==============================
def queryviz(
             input_file  : str | Path | list[str],
             catalog     : str,
             output_file : str | Path = "vizier_output.txt",
             col_width   : int = 15
    ) -> None:
```

* `def function_name(` on its own line, and one argument per line, even for single-argument functions (including `self` and `cls` in methods).
* Arguments start in the column right after the opening parenthesis, and names are padded so that the colons are vertically aligned.
* Exactly one space on each side of the colon (`name : type`) and of the `=` of default values (`col_width : int = 15`).
* No trailing comma after the last argument.
* The closing parenthesis goes on its own line, indented four spaces more than the `def`, followed by ` -> return_type:`. The return type is always annotated, including `-> None`, and full types are kept (`list[str] | tuple[str, ...]`, `dict[str, dict]`).

<br>

### 2.5. Function Body
After the docstring, the body is divided into the following level-1 sections, in this order and omitting only those that are not needed: `MODULES & PACKAGES`, `INPUTS VALIDATION`, `PARAMETERS`, `COMPUTATION`, `PLOTTING` and `OUTPUT`.

* `INPUTS VALIDATION` contains, when needed and in this order, the subsections `VALUES` (validation of the arguments, one argument per block separated by a line containing only `#`), `FUNCTIONS` (all nested helper functions are defined here; only in the functions and methods of `support/`), `DICTIONARY` (initialization of `results_dict`) and any further preparation of the validated inputs (e.g. `BASE FILENAME`).
* `COMPUTATION`, `PLOTTING` and `OUTPUT` are divided into descriptive subsections (e.g. `READING`, `CHECKS`, `QUERY`, `SAVE`, `RESULTS`); alternative branches are named `OPTION 1: ...`, `OPTION 2: ...`. Loose comments (e.g. `# Title & legend`) are never used as block titles: a header of the right level is used instead.
* Section headers have three levels, in uppercase and with the exact number of `#` and spaces:

```python
# ####################  SECTION  ####################
# ##########  SUBSECTION  ##########
# ##### SUB-SUBSECTION #####
```

* Header names are in uppercase, with words separated by spaces (`MODULES & PACKAGES`, not `MODULES&PACKAGES`). Each header is preceded by two blank lines, except when it directly follows another header (no blank line) or the docstring (one blank line).
* Nested functions carry no section headers: their body follows the docstring directly.
* **Classes:** the class line fits on one line (`class _Name:` or `class _Name(Base):`), and the class docstring documents the constructor *Parameters* (and *Notes* when needed). Methods are grouped by level-1 headers (`CONFIGURATION` for class attributes, `REGISTRATION`, `INIT`, then topics such as `READING`, `FIGURE`, `EDITING`, `SAVING`, `EVENTS`, `HELPERS`, `PUBLIC API`), so method bodies use level-2 sections and level-3 subsections in the order above (e.g. `# ##### VALUES #####` inside `# ##########  INPUTS VALIDATION  ##########`). Consecutive methods of a group are separated by a line containing only `#`, every method starts with an underscore (see [Section 2.2](#22-public-vs-internal-code)), and every method documents its *Parameters* (if it takes arguments besides `self` / `cls`) and *Returns* (including `None`).

<br>

### 2.6. Imports
* Imports are placed inside each function, in the `MODULES & PACKAGES` section. Nested functions use the imports of the enclosing function.
* They form a single block sorted alphabetically by full module name, case-insensitively and regardless of whether the statement is `import x` or `from x import y`; the standard library is not separated from third-party packages. Names within each `from ... import ...` are also sorted alphabetically.
* If the function uses internal PANORAMA imports (including imports between `support/` files), the section is split into the sub-subsections `PYTHON` (standard library and third-party) and `PANORAMA`.
* Only what the function actually uses is imported (including types used in `isinstance`, such as `Figure`). New dependencies require a technical justification.

<br>

### 2.7. Code Style
* Ordinary assignments have exactly one space on each side of `=`; keyword arguments in calls have none (`row_limit=1`), also in multi-line calls.
* The `=` of consecutive single-line assignments forming a logical block (same indentation, no blank line or comment between them) is vertically aligned, unless the alignment hurts readability. Multi-line assignments are never aligned with their neighbors (exactly one space before their `=`).
* Multi-line calls, lists and dictionaries start their elements right after the opening bracket and close it on its own line, aligned under the opening bracket, with no trailing comma after the last argument of a call. The 4-space hanging indentation of PEP 8 is not used:

```python
    results_dict.update({
                         'keys' : dup_keys,
                         'dois' : dup_dois,
                        })
```

* Literal dictionaries with several entries use one key per line, the `'key' : value` format with aligned colons and a trailing comma after every entry. Dictionary comprehensions keep the conventional syntax.
* Dictionary keys use single quotes, also when accessed (`row['ra']`, `data.get('ra')`) and in lists of key names; all other strings (messages, paths, file modes, formats, values) use double quotes. The empty string is always written `''` (as in the `results_dict` templates). Messages with variables use f-strings.
* Comments start with an uppercase letter (unless they start with a symbol, a variable or a code identifier, e.g. `# sec(z) = 1 / sin(alt)`), do not end with a period, are separated from inline code by two spaces, and explain the why, the units or the expected format rather than what the code already says. Inside functions, bare strings are never used as comments: a multi-line explanation is a block of `#` lines right below the header it refers to.
* At module level, a block of constants (or of `INIT` imports) is preceded, right below its header, by an explanatory string made of `-` items (wrapped as described in [Section 2.10](#210-docstrings)), with no blank line before the code it describes:

```python
# ####################  LINE COLORS  ####################
"""
- Matplotlib colors of the lines: H and He by species, metals by element
"""
_LINE_COLORS = {
                'H I'  : "tab:red",
                'He I' : "tab:blue",
               }
```

* Indentation is four spaces, with no trailing whitespace. Every text file ends with exactly one newline character after its last line, with no blank lines at the end. Names are descriptive and consistent with astrophysical terminology (`object_names`, `cited_keys`, `written_files`).

<br>

### 2.8. Files and Directories
* Every file or directory argument, input or output, accepts both a `str` and a `pathlib.Path`. It is validated with `isinstance(x, (str, os.PathLike))` and converted to `Path` immediately afterward.
* Paths are manipulated with `Path`, and functions that return paths return `Path` objects (or `list[Path]`), never strings.
* Text files are always opened with `encoding="utf-8"`.
* Output directories are never created automatically: if the parent directory does not exist, an error is raised.
* Functions accepting a list of object names or a file of names follow the package convention: in the file, empty lines and lines starting with `#` are ignored, only the first column is read, and names containing spaces are written with underscores (e.g. `HD_91316`). The existing helpers `_read_object_names` and `_normalize_name` of `support/tex/tex_support.py` are reused instead of being rewritten.

<br>

### 2.9. Input Validation
All arguments must be validated, in the order of the signature, before any computation. This applies to public functions, to private module-level helpers (`support/`), and to class constructors and the methods called from outside their class. Nested functions (which receive already validated values) and Matplotlib / HTTP callbacks (which must never raise, and report errors instead) are exempt. Validation repeated in several functions (e.g. `wave_lim`, `output_file`) is moved to a `support/` helper.

| Situation | Exception |
|---|---|
| Wrong type | `TypeError` |
| Correct type, invalid value (out of range, empty, non-finite) | `ValueError` |
| Required file or directory does not exist | `FileNotFoundError` |
| A file was expected but a directory was found | `IsADirectoryError` |
| The parent of an output exists but is not a directory | `NotADirectoryError` |
| File without read permission | `PermissionError` |
| External failure not attributable to the arguments (service down, corrupted file, all queries failed) | `RuntimeError` (chained with `from e`) |

* Non-fatal failures of a single element inside a loop (e.g. one catalog query) issue a `warnings.warn(...)` and the loop continues, preserving the input order.
* Error messages never use the generic `Variable '<name>'`: they identify the argument by its physical quantity, scientific meaning or role, include its actual name in single quotes and the expected type, value or constraint (with units when relevant), and end with what was received: `Got {type(x).__name__}.` for types, `Got {x}.` for values, and the path between single quotes for files (e.g. `Wavelength range 'wave_lim' must satisfy min < max. Got (5000, 4000).`). The same argument gives the same message across the package.
* Real numbers are checked with `isinstance(x, numbers.Real) and not isinstance(x, bool)` and integers with `isinstance(x, numbers.Integral) and not isinstance(x, bool)` (which also accepts NumPy integers). Non-finite values raise `ValueError`, and physical and mathematical constraints (positivity, intervals, inclusive or exclusive limits) are checked separately.
* Mandatory strings must also be non-empty; NumPy arrays are checked with `isinstance(x, np.ndarray)` and, when relevant, their dimension and numeric type; lists and tuples are validated as containers and then element by element; optional arguments (`... | None = None`) are validated only when not `None`; functions are checked with `callable(x)`; dictionaries are checked for their mandatory keys, value types and constraints.

<br>

### 2.10. Docstrings
* Every function and method has a NumPy-style docstring. Nested functions may be limited to the description, *Parameters* and *Returns*.
* The first line is a single sentence starting with an imperative verb (`Identify ...`, `Query ...`, `Save ...`) and ending with a period, never split over several lines.
* Sections appear in the order *Parameters*, *Returns*, *Raises*, *Notes* and *References*, only when relevant. Each title is placed between two lines of dashes of its own length, and its content is indented four spaces:

```text
    ----------
    Parameters
    ----------
        show : bool, optional
            If True, display the figure with plt.show(). Default is True.
```

* *Parameters* use the `name : type` format with exactly the same type as the signature; arguments with a default value add `, optional` and end with `Default is <value>.`. Units and examples (`e.g.`) are given when they help.
* In dash lists of docstrings and of module-level explanatory strings, continuation lines are aligned with the start of the item text, not with the dash, and lines are wrapped to a length similar to the rest of the docstring (about 79 characters). This does not apply to the descriptions of dictionary keys, indented four spaces below the key, nor to `>>>` examples, formulas or URLs that cannot be split.
* *Returns* uses the name of the variable actually returned; for `-> None`, it states where the results are written.
* Every returned or nested dictionary is documented key by key (never grouping keys), with the type of each value, its physical or computational meaning, its units and sign convention when relevant, the order of the elements of tuples, and the value taken when the data is missing.
* *Raises* documents every exception the function raises explicitly. *Notes* documents assumptions, limitations, conventions, algorithms, particular cases, dependency-version compatibility (e.g. the `astroquery` 0.4.8 API change) and relations with other functions, split into numbered uppercase topics when there is more than one.
* *References* use the A&A format (three authors at most, followed by `et al.`) with ADS journal abbreviations (e.g. `Hunter, J. D. 2007, CSE, 9, 90`) and cite the equations, models, algorithms, databases and essential external packages used (e.g. VizieR, SIMBAD, `astroquery`). A reference longer than about 79 characters is split after the year, with the continuation aligned with the text of the item. References are never invented: unverifiable details are reported to the maintainer instead.

<br>

### 2.11. Returned Results
* A single simple value (a list, a figure, a path) is returned directly; results written to a file with nothing to return use `-> None`.
* Several related results are returned in a single `results_dict`. It is initialized in the `DICTIONARY` subsection with every key that will be returned (`np.nan` for numbers, `''` for strings, `None` for tuples or objects, `{}`/`[]` for containers), updated at the end with `results_dict.update({...})` and returned.
* If its keys are identifiers only known during the computation (e.g. object names), it is initialized empty, with a comment explaining why, and a template with every key initialized is created inside the loop before filling it.
* Keys, types and meanings are kept consistent with the docstring and the module README.

<br>

### 2.12. Figures
* Figures are saved with `savefig_pkl` (label alignment, `bbox_inches='tight'`, `pad_inches=0.01`, 500 dpi and a `.pkl` copy) instead of calling `fig.savefig` directly.
* They follow the A&A conventions of PANORAMA (e.g. a single-column width of 3.46 in), with axis labels giving units between square brackets (`Wavelength [Å]`).
* Functions that create figures accept a `show : bool` argument and return the figure without saving it; they work both with the inline backend of Jupyter (figure in the notebook) and with interactive backends (window). Interactive tools that wait for their window to be closed use `plt.show(block=True)`, which also blocks in IPython / Jupyter.

<br>

### 2.13. Scientific Consistency
* Physical definitions, units, coordinate systems (e.g. ICRS) and sign conventions of PANORAMA are preserved, and units are made explicit in names, docstrings or comments whenever they avoid ambiguity. Documented ranges and units match the ones validated in the code.
* Functions depending on external services (SIMBAD, VizieR, ...) or on specific dependency versions support every version accepted by the package and document it in *Notes*.

<br>

### 2.14. Module READMEs
Each public module has a `READMEs/<module>_README.md` (model: `tex_README.md`) with: title, badges, a one-sentence summary and a link back to the main README; *Table of Contents*; *Overview* (including the conventions shared by the module); *Function Summary* (one table per group of functions); *Function Reference* (signature, description, and the labels **Inputs:**, **Returns:** and **Raises:**, each on its own line: the first two followed by a blank line and their table, the parameter table with units and defaults, and the last one by the exceptions); *Dependencies*; *Internal Structure*; *Demo*; and *References*. Signatures are copied from the code (without `def` and the final colon, closing parenthesis indented four spaces), and parameter types, return types, defaults, keys and exceptions must match the code exactly.

<br>

### 2.15. Demos
* Each module has a demo notebook `Demos/<module>_demo.ipynb` (model: `tex_demo.ipynb`), delivered executed (with its outputs and figures), so that it can be read directly on GitHub.
* The first cell is a Markdown cell with the title (`` # `panorama.<module>` — Demo ``), the purpose of the notebook, what each example covers, where data and results are written and which examples require internet access.
* The notebook is organized in sections, each one a Markdown cell with a level-2 heading (`## Imports`, `## Setup`, `## Example 1: <title>`, ...) and a one-line sentence with its goal, followed by its code cell(s). Code cells follow the code conventions of this document. Functions are imported from `panorama` in the imports cell, in alphabetical order with a comment indicating their group.
* Figures are shown inline (`%matplotlib inline` in the imports cell): the plotting functions are called with `show=True`, or the figure is displayed with `IPython.display.display(fig)` before saving it with `savefig_pkl` (which closes it). Figures that are not shown are closed with `plt.close(fig)`, so that they are not displayed at the end of the cell.
* Input data are read from `Demos/data_and_results/<module>/data/` and results written to `Demos/data_and_results/<module>/results/` (paths relative to `Demos/`, the folder of the notebook), created with `Path.mkdir(parents=True, exist_ok=True)`.
* Examples depending on internet or external services are wrapped in `try/except` and print `Skipped: ...` on failure without stopping the notebook.
* Notebooks cover every public function of the module and the particular cases documented in *Notes*, and run from top to bottom (*Run All*) without user interaction. Interactive tools (e.g. `rmvspec` and `normspec` of `spec`) are run non-interactively by default; a flag in the first code cell (e.g. `INTERACTIVE = False`) selects the Qt backend (`%matplotlib qt`, with PyQt6 from the `demos` extra) to use them.
* Each example prints its results with tabs. Expected warnings are captured with `warnings.catch_warnings(record=True)` and only their message is printed, so that no internal paths of the library appear. Mock (invented) test data are described as *mock* in the notebook and in the README. Demos never call private functions, classes or methods.
* Before committing a notebook: restart the kernel and *Run All* with the package installed from the repository (`pip install -e ".[demos]"`), so that the execution counts start at 1; check that the outputs contain no local paths or personal information, that the kernel metadata is generic (`"Python 3"`) and that no old files in `results/` alter the outputs.

<br>

### 2.16. Checks
* No code formatter (e.g. `black`) is used, since it would undo the alignment of signatures, dictionaries and assignments.
* Every modified Python file must pass `python -m py_compile <file>.py`, and the corresponding demo notebook must be run from top to bottom (e.g. `jupyter nbconvert --to notebook --execute --inplace Demos/<module>_demo.ipynb`) before delivering changes. When changes affect the API, the documented behavior or the internal structure, the module README and demo are updated in the same change.
* Before delivering, the whole file is checked against these conventions (layout, signatures, headers, naming, imports, quotes, alignment, validation, docstrings, READMEs and demos), not only the modified lines. Syntactic correctness, scientific validity and compatibility with PANORAMA take priority over purely aesthetic preferences.
