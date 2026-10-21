# PANORAMA — Package for the ANalysis Of RApidly rotating MAssive stars

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT + citation](https://img.shields.io/badge/License-MIT%20%2B%20citation-yellow.svg)](LICENSE)

**PANORAMA** (*Package for the ANalysis Of RApidly rotating MAssive stars*) is a Python framework for the physical and spectroscopic analysis of rotating massive stars.

It brings together, in a single and consistent environment, the tools needed along the whole analysis chain: from planning the observations, through handling and measuring the spectra, to the physical modeling of the rotating star and the preparation of the results for publication.

---
<br>

## Table of Contents
- [1. Scientific Context](#1-scientific-context)
- [2. Installation](#2-installation)
- [3. Quick Start](#3-quick-start)
- [4. Modules](#4-modules)
- [5. Repository Structure](#5-repository-structure)
- [6. Demos](#6-demos)
- [7. Contributing and Code Conventions](#7-contributing-and-code-conventions)
- [8. How to Cite](#8-how-to-cite)
- [9. License](#9-license)
  - [9.1. Summary](#91-summary)
  - [9.2. Citation Requirement](#92-citation-requirement)
- [10. Authors and Contact](#10-authors-and-contact)
- [11. Project Links](#11-project-links)

---
<br>

## 1. Scientific Context

Massive stars are frequently characterized by very rapid rotation and a high multiplicity fraction. Rotation distorts the star away from spherical symmetry (Roche geometry), makes the surface temperature and gravity depend on latitude (gravity darkening), and broadens and reshapes the spectral lines. As a consequence, standard 1D, non-rotating approaches are often insufficient to interpret their observables.

**PANORAMA** addresses these challenges with a modular design: each module covers one stage or area of the analysis (observations, spectra, line measurements, stellar physics, publication) and all of them share the same conventions, so that they can be combined freely. It is equally suited to the analysis of a single object and to large observing campaigns.

The modules currently available are listed in [Section 4](#4-modules).

---
<br>

## 2. Installation

### Requirements
PANORAMA requires **Python 3.10 or later** and has been developed and tested on Linux. Its third-party dependencies and their minimum versions are declared in [`pyproject.toml`](pyproject.toml) and are installed automatically by `pip`. Module-specific dependencies are described in the *Dependencies* section of each module README.

### Install
```bash
git clone https://github.com/DGalanDieguez/PANORAMA
cd PANORAMA
pip install -e .            # Editable installation (recommended): code changes take effect without reinstalling
```

Use `pip install .` instead for a regular, non-editable installation.

> **Note:** PANORAMA is installed from this repository. The name `panorama` on PyPI belongs to an unrelated project, so `pip install panorama` does **not** install this package.

### Verification
```bash
python -c "import panorama; print(panorama.__version__); panorama.functions()"
```

The version and the list of public functions, grouped by module, are printed if the installation succeeded.

---
<br>

## 3. Quick Start

Every public function of every module is available directly at package level:

```python
import panorama

panorama.functions()              # All public functions, grouped by module
panorama.functions("<module>")    # Public functions of a single module (e.g. "obs")
help(panorama.<function>)         # Full documentation of a function (e.g. panorama.visibility)

from panorama import <function>   # Import a function directly
```

A short startup message is printed when PANORAMA is imported. To silence it (e.g. in scripts or batch jobs), set the environment variable `PANORAMA_QUIET=1`.

Worked examples for each module are given in [Demos](#6-demos).

---
<br>

## 4. Modules

Each module has its own README in [`READMEs/`](READMEs), all with the same structure: overview, function summary and reference, dependencies, internal structure, demo and references.

| Module | Description | Status | Documentation |
|---|---|:---:|---|
| `fitline` | Spectral-line profile fitting and measurement of line parameters | 🚧 | — |
| `obs` | Target visibility, observation scheduling, timing (JD/HJD and ephemeris phases) and observing logs | ✅ | [obs_README.md](READMEs/obs_README.md) |
| `prismas` | Model-atmosphere grid selection and output handling for the PRISMAS pipeline | 🚧 | — |
| `pyiacob_broad` | Python implementation of the IACOB-BROAD method: projected rotational velocity (*v* sin *i*) and macroturbulence from combined Fourier-transform and goodness-of-fit analyses | 🚧 | — |
| `rot` | Physics of rapidly rotating stars: Roche-model geometry and gravity darkening | 🚧 | — |
| `sii` | Joint fitting of stellar intensity interferometric and spectroscopic data | 🚧 | — |
| `spec` | Reading, visualization, cleaning and continuum normalization of observed spectra | ✅ | [spec_README.md](READMEs/spec_README.md) |
| `specblur` | Signal-to-noise (SNR) estimation and degradation of spectra to a target resolution, SNR, rotation and macroturbulence | 🚧 | — |
| `tex` | LaTeX manuscript and bibliography checks, VizieR/SIMBAD catalog queries and A&A-ready figures | ✅ | [tex_README.md](READMEs/tex_README.md) |
| `units` | Physical and astronomical unit conversions and derived stellar parameters | ✅ | [units_README.md](READMEs/units_README.md) |

**Development Status:** ✅ *Stable*: available and documented; no major changes are expected&nbsp;&nbsp;|&nbsp;&nbsp;🚧 *Under Development*: not publicly available; access can be requested from the authors

---
<br>

## 5. Repository Structure

```text
PANORAMA/
├── pyproject.toml                  # Build configuration, metadata and dependencies
├── README.md                       # This file
├── LICENSE                         # MIT License with a citation requirement
├── CITATION.cff                    # Citation metadata (GitHub "Cite this repository")
├── .gitignore
├── READMEs/
│   ├── <module>_README.md          # Documentation of each module
│   ├── BinarAlt_UserGuide.pdf      # User guide for BinarAlt
│   └── conventions_README.md       # How to add a module and code conventions
├── Demos/
│   ├── <module>_demo.ipynb         # Jupyter notebook with worked examples of each module
│   └── data_and_results/
│       └── <module>/
│           ├── data/               # Example input data used by each demo
│           └── results/            # Outputs written by each demo (not tracked by git)
└── src/
    └── panorama/
        ├── __init__.py             # Detects every module and exposes its public functions
        ├── <module>.py             # User-facing functions of each module
        └── support/
            ├── <module>/           # Internal code and data files of each module
            └── package/            # Internal helpers of __init__.py
```

Every module follows the same pattern, so its files can be located from its name alone:

| Path | Content |
|---|---|
| `src/panorama/<module>.py` | Public (user-facing) functions. |
| `src/panorama/support/<module>/` | Internal helpers and data files (every private function and class of the module lives here). |
| `READMEs/<module>_README.md` | Documentation. |
| `Demos/<module>_demo.ipynb` | Worked examples, as a Jupyter notebook (when available). |
| `Demos/data_and_results/<module>/data/` | Example input data of the demo. |
| `Demos/data_and_results/<module>/results/` | Outputs of the demo. |

---
<br>

## 6. Demos

The Jupyter notebooks `Demos/<module>_demo.ipynb` exercise every public function of each module end to end, one section per example, with annotated output and the figures shown inline. They are distributed already executed, so they can be read directly on GitHub. JupyterLab and PyQt6 (for the interactive Matplotlib windows) are optional dependencies, installed with the `demos` extra; open the notebooks from `Demos/`, since their paths are relative to that folder:

```bash
pip install -e ".[demos]"     # PANORAMA + JupyterLab and PyQt6 (demo extras)
cd Demos
jupyter lab                   # Or: jupyter notebook
```

Some notebooks need extra resources, such as internet access; the examples that need them print `Skipped: ...` if they are not available, and the rest of the notebook still runs. This is stated in the first cell of each notebook and in the *Demo* section of the corresponding module README. Input data and outputs are written to `Demos/data_and_results/<module>/`.

---
<br>

## 7. Contributing and Code Conventions

New modules and changes to existing code must follow the PANORAMA development conventions, described in detail in [`READMEs/conventions_README.md`](READMEs/conventions_README.md): how to add a module, module and function layout, imports, code style, file handling, input validation, docstrings, returned results, figures, module READMEs, demos and checks.

In short:
* **Public vs. internal code:** only the functions of `src/panorama/<module>.py` without a leading underscore are public (they are exposed as `panorama.<function>`); everything else (helpers and classes under `support/`, their methods, callbacks and nested functions) starts with an underscore, is internal and may change without notice.
* **Unique names:** function names must be unique across modules, since all of them share the `panorama` namespace.
* **Reference implementation:** the `tex` module (`tex.py`, `tex_README.md` and `tex_demo.ipynb`) is the model for new code, documentation and demos.

---
<br>

## 8. How to Cite

If you use PANORAMA in your research, please cite the reference publication (placeholder, to be replaced by the ADS entry once the paper is published):

```bibtex
@ARTICLE{PANORAMA_Reference,
       author = {{Gal{\'a}n-Di{\'e}guez}, D. and others},
        title = "{PANORAMA: A Package for the Analysis of Rapidly Rotating Massive Stars}",
      journal = {xxx},
         year = {202X},
       volume = {xxx},
        pages = {xxx},
          doi = {xxx/xxx}
}
```

The same information is available in machine-readable form in [`CITATION.cff`](CITATION.cff) (GitHub shows it as *Cite this repository*).

Some modules implement methods published by other authors (e.g. IACOB-BROAD) or rely on external services and packages (e.g. VizieR, SIMBAD, Astropy). When you use them, please also cite the original works listed in the *References* section of the corresponding module README.

---
<br>

## 9. License

### 9.1. Summary
PANORAMA is distributed under the terms of the **MIT License with an additional citation requirement**. The full text is given in the [LICENSE](LICENSE) file; in short:

* **You may** use, copy, modify, merge, publish, distribute, sublicense and sell copies of the software, free of charge.
* **You must**:
  * keep the copyright notice and the license text in all copies or substantial portions of the software;
  * cite the reference publication given in [Section 8](#8-how-to-cite) in any academic, scientific or technical publication, presentation or public dissemination of results obtained in whole or in part with PANORAMA.
* **No warranty:** the software is provided "as is", without warranty of any kind, and the authors are not liable for any claim or damage arising from its use.

This summary is for convenience only; the [LICENSE](LICENSE) file is the legally binding text. Because of the citation clause, the license is not the standard MIT License; it is declared in `pyproject.toml` as `LicenseRef-MIT-Citation`.

<br>

### 9.2. Citation Requirement
The reference publication that the [LICENSE](LICENSE) file requires to be cited is given, with its BibTeX entry, in [Section 8](#8-how-to-cite). Methods implemented from other authors' work should also be cited as indicated in the corresponding module README.

---
<br>

## 10. Authors and Contact

* **D. Galán-Diéguez** – *Lead developer* – [dgalandieguez.astro@gmail.com](mailto:dgalandieguez.astro@gmail.com)
* **M. Arévalo-Morales** – *Co-developer of `BinarAlt`* – [marevalo@ull.edu.es](mailto:marevalo@ull.edu.es)

---
<br>

## 11. Project Links

* **Personal website:** [dgalandieguez.com](https://dgalandieguez.com)
* **PANORAMA website:** [dgalandieguez.com/PANORAMA](https://dgalandieguez.com/PANORAMA)
