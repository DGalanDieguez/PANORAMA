"""
Internal helpers of the 'units' module: validation of the scalar physical
quantities shared by the public functions.

Conventions shared by every function in this module:
    - Functions are private (leading underscore); they are used by the
      public functions of 'panorama.units' (and 'panorama.obs',
      'panorama.spec' and 'panorama.prismas') and are not exposed at the
      package level.
    - Validated values are returned as Python floats, so that every
      function works with the same numerical type.
    - All functions validate their inputs before computing.
"""

from __future__ import annotations

# =====================================================================
# ============================== NUMBERS ==============================
# =====================================================================

# ==============================
def _check_real(
                value      : object,
                name       : str,
                label      : str,
                unit       : str | None = None,
                constraint : str | tuple[float, float, str] | None = None,
                allow_none : bool = False,
                allow_inf  : bool = False
    ) -> float | None:
    """
    Validate a scalar real argument (type, finiteness and physical constraint) and return it as a float.

    ----------
    Parameters
    ----------
        value : object
            Argument to check.

        name : str
            Name of the argument in the error messages (e.g. 'teff').

        label : str
            Physical quantity of the argument in the error messages (e.g.
            'Effective temperature').

        unit : str | None, optional
            Unit or meaning shown between parentheses in the type error
            (e.g. 'in K', 'dimensionless', 'log10(g / [cm s^-2])'). Default
            is None (nothing shown).

        constraint : str | tuple[float, float, str] | None, optional
            Physical constraint: 'positive' (> 0), 'non-negative' (>= 0), or
            an interval (low, high, bounds) where 'bounds' is '[]', '[)',
            '(]' or '()'. Default is None (any finite value).

        allow_none : bool, optional
            If True, None is accepted and returned unchanged (optional
            arguments). Default is False.

        allow_inf : bool, optional
            If True, -inf and +inf are accepted (NaN never is); used for open
            wavelength limits in 'spec'. Default is False.

    -------
    Returns
    -------
        value : float | None
            The argument as a Python float (possibly infinite if
            'allow_inf'), or None if it was None and 'allow_none' is True.

    ------
    Raises
    ------
        TypeError
            If 'value' is not a real number (booleans excluded; None only
            when 'allow_none'), or 'name', 'label', 'unit', 'constraint',
            'allow_none' or 'allow_inf' have the wrong type.

        ValueError
            If 'value' is NaN, infinite (unless 'allow_inf') or does not
            satisfy 'constraint', or 'name', 'label' or 'constraint' are not
            valid.

    -----
    Notes
    -----
        - Shared by every public function of 'units' (and by 'obs', 'spec'
          and 'prismas'), so that the same argument gives the same error
          messages in the whole package, e.g.:
                >>> Effective temperature 'teff' must be a real number (in K). Got str.
                >>> Effective temperature 'teff' must be positive. Got -5.0.
                >>> Hydrogen mass fraction 'H_frac' must lie in [0, 1]. Got 1.2.
    """

    # ####################  MODULES & PACKAGES  ####################
    import numbers
    import numpy as np


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    for arg_name, arg_value in (("name", name), ("label", label)):
        if not isinstance(arg_value, str):
            raise TypeError(f"Error-message text '{arg_name}' must be a string. Got {type(arg_value).__name__}.")
        elif not arg_value.strip():
            raise ValueError(f"Error-message text '{arg_name}' must be a non-empty string. Got '{arg_value}'.")
    #
    if unit is not None and not isinstance(unit, str):
        raise TypeError(f"Unit text 'unit' must be a string or None. Got {type(unit).__name__}.")
    #
    if constraint is not None:
        if isinstance(constraint, str):
            if constraint not in ("positive", "non-negative"):
                raise ValueError(f"Constraint 'constraint' must be 'positive', 'non-negative' or a (low, high, bounds) tuple. Got '{constraint}'.")
        elif isinstance(constraint, tuple):
            if len(constraint) != 3 or constraint[2] not in ("[]", "[)", "(]", "()"):
                raise ValueError(f"Constraint 'constraint' must be a (low, high, bounds) tuple with bounds '[]', '[)', '(]' or '()'. Got {constraint!r}.")
            elif not all(isinstance(x, numbers.Real) and not isinstance(x, bool) and np.isfinite(x) for x in constraint[:2]) or constraint[0] >= constraint[1]:
                raise ValueError(f"Constraint 'constraint' must have finite real limits with low < high. Got {constraint!r}.")
        else:
            raise TypeError(f"Constraint 'constraint' must be a string, a (low, high, bounds) tuple or None. Got {type(constraint).__name__}.")
    #
    if not isinstance(allow_none, bool):
        raise TypeError(f"None flag 'allow_none' must be a boolean. Got {type(allow_none).__name__}.")
    #
    if not isinstance(allow_inf, bool):
        raise TypeError(f"Infinity flag 'allow_inf' must be a boolean. Got {type(allow_inf).__name__}.")


    # ####################  COMPUTATION  ####################
    # ##########  TYPE  ##########
    if value is None and allow_none:  # Optional argument not given
        return None
    if not isinstance(value, numbers.Real) or isinstance(value, bool):  # numbers.Real also accepts NumPy scalars; bool is excluded explicitly
        unit_text = f" ({unit})" if unit else ''  # E.g. ' (in K)'
        none_text = " or None" if allow_none else ''
        raise TypeError(f"{label} '{name}' must be a real number{unit_text}{none_text}. Got {type(value).__name__}.")


    # ##########  FINITENESS  ##########
    if np.isnan(value):  # NaN is never accepted
        expected = "must not be NaN" if allow_inf else "must be finite"
        raise ValueError(f"{label} '{name}' {expected}. Got {value}.")
    elif not allow_inf and np.isinf(value):  # ±inf only accepted for open limits
        raise ValueError(f"{label} '{name}' must be finite. Got {value}.")


    # ##########  CONSTRAINT  ##########
    if constraint == "positive" and value <= 0:
        raise ValueError(f"{label} '{name}' must be positive. Got {value}.")
    elif constraint == "non-negative" and value < 0:
        raise ValueError(f"{label} '{name}' must be non-negative. Got {value}.")
    elif isinstance(constraint, tuple):
        low, high, bounds = constraint
        above             = value >= low if bounds[0] == "[" else value > low  # '[': closed lower bound, '(': open
        below             = value <= high if bounds[1] == "]" else value < high  # ']': closed upper bound, ')': open
        if not (above and below):
            raise ValueError(f"{label} '{name}' must lie in {bounds[0]}{low:g}, {high:g}{bounds[1]}. Got {value}.")


    # ####################  OUTPUT  ####################
    return float(value)  # Python float (also for NumPy scalars and integers)
