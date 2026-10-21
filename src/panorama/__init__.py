"""
PANORAMA: Package for the ANalysis Of RApidly rotating MAssive stars.

Every public function defined in the top-level submodules of the package is
exposed directly at the package level (e.g. 'panorama.visibility' for
'panorama.obs.visibility'). Run 'panorama.functions()' to list them, grouped
by submodule.

Conventions shared by every function in this module:
    - Functions whose names start with an underscore are private and are not
      exposed at the package level; the private helpers used to build the
      namespace are defined in 'support/package/package_support.py'.
    - Submodules are detected automatically: adding a new public module to
      'src/panorama/' is enough to expose its functions.
    - The startup message can be silenced with the environment variable
      PANORAMA_QUIET=1.
"""

from __future__ import annotations

# ==================================================================
# ============================== INIT ==============================
# ==================================================================

# ####################  MODULES & PACKAGES  ####################
"""
- The package namespace (version, submodules and startup message) is built
  when 'panorama' is imported, outside any function, so these imports are
  needed at module level
- The private helpers that detect and expose the public functions are defined
  in 'support/package/package_support.py', since no private function is
  defined in the files of 'src/panorama/'
"""

# ##### PYTHON #####
from importlib.metadata import PackageNotFoundError, version
import os
import pkgutil


# ##### PANORAMA #####
from .support.package.package_support import _public_functions


# ####################  PACKAGE METADATA  ####################
"""
- The version is defined only once, in 'pyproject.toml', and read from the
  installed package metadata; if PANORAMA is imported without being installed
  (e.g. directly from 'src/'), a placeholder version is used instead
"""
try:
    __version__ = version("panorama")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__author__  = "D. Galán-Diéguez"
__license__ = "LicenseRef-MIT-Citation"  # MIT License with an additional citation clause (see LICENSE)


# ####################  CONFIGURATION  ####################
"""
- Submodules never imported automatically, environment variable that silences
  the startup message, and package-level names that must not be overwritten
"""
_EXCLUDED_SUBMODULES = frozenset()
_QUIET_ENV_VARIABLE  = "PANORAMA_QUIET"  # Set to "1" to silence the startup message
_RESERVED_NAMES      = frozenset({'functions', 'submodules'})




# ===================================================================
# ============================== UTILS ==============================
# ===================================================================

# ==============================
def functions(
              module : str | None = None
    ) -> None:
    """
    Print the public functions available in PANORAMA, grouped by submodule.

    ----------
    Parameters
    ----------
        module : str | None, optional
            Restrict the listing to a single submodule (e.g. "units").
            Default is None (all submodules are listed).

    -------
    Returns
    -------
        None
            The function names and the first line of their docstrings are
            printed to the standard output.

    ------
    Raises
    ------
        TypeError
            If 'module' is not a string.

        ValueError
            If 'module' is not one of the available submodules.

    -----
    Notes
    -----
        - Example:
            >>> import panorama
            >>> panorama.functions("obs")
    """

    # ####################  MODULES & PACKAGES  ####################
    import inspect


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if module is not None:
        if not isinstance(module, str):
            raise TypeError(f"Submodule name 'module' must be a string. Got {type(module).__name__}.")
        elif module not in _functions_by_module:
            raise ValueError(f"Submodule name 'module' must be one of {submodules}. Got '{module}'.")


    # ####################  PARAMETERS  ####################
    selected = [module] if module is not None else submodules  # Requested submodule or all available ones


    # ####################  OUTPUT  ####################
    for i, submodule in enumerate(selected):
        leading_blank_lines = "\n\n" if i != 0 else ''
        print(f"{leading_blank_lines}   {'-' * (len(submodule) + 10)}")
        print(f"   {submodule} SUBMODULE")
        print(f"   {'-' * (len(submodule) + 10)}")

        for function_name, function in _functions_by_module[submodule].items():
            docstring  = inspect.getdoc(function)
            first_line = docstring.strip().splitlines()[0] if docstring else ''  # First line of the docstring
            print(f"\t{function_name:<25} {first_line}")
    print()




# =======================================================================
# ============================== NAMESPACE ==============================
# =======================================================================

# ####################  SUBMODULES  ####################
"""
- Public, non-package submodules detected automatically
"""
submodules = sorted(name for _, name, is_pkg in pkgutil.iter_modules(__path__) if not is_pkg and not name.startswith("_") and name not in _EXCLUDED_SUBMODULES)


# ####################  PUBLIC API  ####################
"""
- Public functions exposed at the package level, and names available through
  'from panorama import *'
"""
_functions_by_module = _public_functions(__name__, submodules, globals(), _RESERVED_NAMES)  # Also adds every public function to this namespace
public_functions     = sorted({name for funcs in _functions_by_module.values() for name in funcs} - _RESERVED_NAMES)
__all__              = [*public_functions, "functions", "__version__"]


# ####################  STARTUP MESSAGE  ####################
"""
- Printed once when the package is imported, unless PANORAMA_QUIET=1
"""
if os.environ.get(_QUIET_ENV_VARIABLE, "0") != "1":
    print(
          f"\tPANORAMA v{__version__} loaded...\n"
          "\t\t- To see the available functions, run:\n"
          "\t\t\tpanorama.functions()          # Lists all functions, grouped by submodule\n"
          "\t\t\tpanorama.functions('obs')     # Lists the functions of a single submodule\n"
          "\t\t- To see the docstring of a specific function (e.g. 'visibility'), run:\n"
          "\t\t\thelp(panorama.visibility)\n"
          f"\t\t- To silence this message, set the environment variable {_QUIET_ENV_VARIABLE}=1\n"
         )
