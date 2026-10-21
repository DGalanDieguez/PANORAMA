"""
Internal helpers of the package initialization file ('panorama/__init__.py'):
detection of the public functions of each submodule and their exposure at
the package level.

Conventions shared by every function in this module:
    - Functions are private (leading underscore); they are used only by
      'panorama/__init__.py' when the package is imported, and are not
      exposed at the package level.
    - Public functions are the functions defined in a submodule whose names
      do not start with an underscore.
    - All functions validate their inputs before computing.
"""

from __future__ import annotations

# =======================================================================
# ============================== DISCOVERY ==============================
# =======================================================================

# ==============================
def _get_public_functions(
                          module : ModuleType
    ) -> dict[str, FunctionType]:
    """
    Return the public functions defined in a module.

    ----------
    Parameters
    ----------
        module : ModuleType
            Module to inspect.

    -------
    Returns
    -------
        module_functions : dict[str, FunctionType]
            A dictionary keyed by function name, in alphabetical order; each
            value is the function object.

    ------
    Raises
    ------
        TypeError
            If 'module' is not a module.

    -----
    Notes
    -----
        - Functions whose names start with an underscore are private and are
          not included.

        - Functions imported into the namespace of the module from elsewhere
          are not included.
    """

    # ####################  MODULES & PACKAGES  ####################
    import inspect
    from types import ModuleType


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(module, ModuleType):
        raise TypeError(f"Submodule 'module' must be a Python module. Got {type(module).__name__}.")


    # ####################  COMPUTATION  ####################
    # inspect.getmembers returns the members sorted by name; only the functions defined in the module itself are kept
    module_functions = {name: obj for name, obj in inspect.getmembers(module, inspect.isfunction) if not name.startswith("_") and obj.__module__ == module.__name__}


    # ####################  OUTPUT  ####################
    return module_functions




# ==============================
def _public_functions(
                      package_name    : str,
                      submodule_names : list[str],
                      namespace       : dict[str, object],
                      reserved_names  : frozenset[str]
    ) -> dict[str, dict[str, FunctionType]]:
    """
    Import every detected submodule and expose its public functions directly at the package level.

    ----------
    Parameters
    ----------
        package_name : str
            Name of the package (e.g. 'panorama'), used to import its
            submodules.

        submodule_names : list[str]
            Names of the submodules to import (e.g. ['obs', 'spec']).

        namespace : dict[str, object]
            Namespace of the package (the globals() of its '__init__.py'),
            where the public functions are added.

        reserved_names : frozenset[str]
            Package-level names that must not be overwritten (e.g.
            'functions', 'submodules').

    -------
    Returns
    -------
        functions_by_module : dict[str, dict[str, FunctionType]]
            A dictionary keyed by submodule name; each value is a dictionary
            keyed by function name whose values are the function objects.

    ------
    Raises
    ------
        TypeError
            If 'package_name' is not a string, 'submodule_names' is not a
            list of strings, 'namespace' is not a dictionary or
            'reserved_names' is not a frozenset of strings.

        ValueError
            If 'package_name' is empty.

    -----
    Notes
    -----
        - Kept in a function (rather than a loop at module level) so that its
          loop variables stay local and never leak into the 'panorama'
          namespace as spurious public attributes.

        - If two submodules define a function with the same name, the one
          from the submodule that comes later alphabetically is kept and a
          warning is issued.

        - Functions whose names are reserved are not exposed at the package
          level (a warning is issued); they remain available as
          'panorama.<submodule>.<function>'.
    """

    # ####################  MODULES & PACKAGES  ####################
    import importlib
    import warnings


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(package_name, str):
        raise TypeError(f"Package name 'package_name' must be a string. Got {type(package_name).__name__}.")
    elif not package_name.strip():
        raise ValueError(f"Package name 'package_name' must be a non-empty string. Got '{package_name}'.")
    #
    if not isinstance(submodule_names, list) or not all(isinstance(name, str) for name in submodule_names):
        raise TypeError(f"Submodule names 'submodule_names' must be a list of strings. Got {type(submodule_names).__name__}.")
    #
    if not isinstance(namespace, dict):
        raise TypeError(f"Package namespace 'namespace' must be a dictionary (the globals() of the package). Got {type(namespace).__name__}.")
    #
    if not isinstance(reserved_names, frozenset) or not all(isinstance(name, str) for name in reserved_names):
        raise TypeError(f"Reserved names 'reserved_names' must be a frozenset of strings. Got {type(reserved_names).__name__}.")


    # ##########  DICTIONARY  ##########
    functions_by_module = {}  # Filled with one entry per submodule (keys are the submodule names)


    # ####################  COMPUTATION  ####################
    origin = {}  # Submodule that each exposed name comes from (to detect duplicated names)

    for submodule in submodule_names:  # Alphabetical order: on duplicated names, the last submodule wins
        module           = importlib.import_module(f"{package_name}.{submodule}")  # E.g. 'panorama.obs'
        module_functions = _get_public_functions(module)

        for function_name, function in module_functions.items():


            # ##### RESERVED NAMES #####
            if function_name in reserved_names:
                warnings.warn(
                              f"'{package_name}.{submodule}.{function_name}' is not exposed at package level because "
                              f"'{function_name}' is reserved. Use '{package_name}.{submodule}.{function_name}' instead.",
                              stacklevel=2
                             )
                continue


            # ##### DUPLICATED NAMES #####
            if function_name in origin:
                warnings.warn(
                              f"Function name '{function_name}' is defined in both '{package_name}.{origin[function_name]}' and "
                              f"'{package_name}.{submodule}'; '{package_name}.{function_name}' now refers to the one in '{package_name}.{submodule}'.",
                              stacklevel=2
                             )


            # ##### EXPOSE #####
            namespace[function_name] = function  # Available as 'panorama.<function>'
            origin[function_name]    = submodule

        functions_by_module[submodule] = {name: func for name, func in module_functions.items() if name not in reserved_names}  # Listed by 'panorama.functions'


    # ####################  OUTPUT  ####################
    return functions_by_module
