"""
Reference data of the 'spec' module: rest-frame diagnostic lines of O and B
stars and the colors used to draw them.

Conventions shared by every object in this module:
    - Constants are private (leading underscore); they are used by the
      functions of 'panorama.spec' and are not exposed at the package level.
    - Wavelengths are rest-frame air wavelengths in Å (convert them with
      'air2vac' for vacuum spectra).
    - Lines are keyed by their plot label: the Greek letter for the Balmer
      lines and '<species> <rounded wavelength>' otherwise (e.g.
      'Si III 4553').
    - Species are written in spectroscopic notation ('He I' = neutral,
      'He II' = singly ionized, ...).
"""

from __future__ import annotations

# ==================================================================
# ============================== INIT ==============================
# ==================================================================

# ####################  OB LINES  ####################
"""
- Diagnostic lines of O and B stars; each value is:
    >> rest-frame air wavelength [Å],
    >> species
"""
_OB_LINES = {
             # Hydrogen: Balmer series
             'Hε'           : (3970.07, "H I"),
             'Hδ'           : (4101.74, "H I"),
             'Hγ'           : (4340.47, "H I"),
             'Hβ'           : (4861.33, "H I"),
             'Hα'           : (6562.80, "H I"),
             # Neutral helium: Teff and log g diagnostics of B stars
             'He I 4009'    : (4009.26, "He I"),
             'He I+II 4026' : (4026.19, "He I"),  # He I 4026.19 blended with He II 4025.6 in O stars
             'He I 4144'    : (4143.76, "He I"),
             'He I 4388'    : (4387.93, "He I"),
             'He I 4471'    : (4471.48, "He I"),
             'He I 4713'    : (4713.15, "He I"),
             'He I 4922'    : (4921.93, "He I"),
             'He I 5016'    : (5015.68, "He I"),
             'He I 5876'    : (5875.62, "He I"),
             'He I 6678'    : (6678.15, "He I"),
             # Ionized helium: Teff diagnostics of O stars (He I / He II ionization balance)
             'He II 4200'   : (4199.83, "He II"),
             'He II 4542'   : (4541.59, "He II"),
             'He II 4686'   : (4685.71, "He II"),
             'He II 5412'   : (5411.52, "He II"),
             # Carbon
             'C II 4267'    : (4267.26, "C II"),
             'C III 4650'   : (4650.25, "C III"),
             'C IV 5801'    : (5801.33, "C IV"),
             'C IV 5812'    : (5811.98, "C IV"),
             # Nitrogen
             'N II 3995'    : (3994.99, "N II"),
             'N III 4634'   : (4634.14, "N III"),
             'N III 4641'   : (4640.64, "N III"),
             'N IV 4058'    : (4057.76, "N IV"),
             'N V 4604'     : (4603.73, "N V"),
             'N V 4620'     : (4619.97, "N V"),
             # Oxygen
             'O II 4415'    : (4414.90, "O II"),
             'O II 4662'    : (4661.63, "O II"),
             'O III 5592'   : (5592.25, "O III"),
             # Magnesium: v sin i and spectral type of B stars
             'Mg II 4481'   : (4481.13, "Mg II"),
             # Silicon: Teff diagnostics of B stars (Si II / Si III / Si IV ionization balance)
             'Si II 4128'   : (4128.07, "Si II"),
             'Si II 4131'   : (4130.89, "Si II"),
             'Si III 4553'  : (4552.62, "Si III"),
             'Si III 4568'  : (4567.84, "Si III"),
             'Si III 4575'  : (4574.76, "Si III"),
             'Si IV 4089'   : (4088.86, "Si IV"),
             'Si IV 4116'   : (4116.10, "Si IV"),
            }


# ####################  LINE COLORS  ####################
"""
- Matplotlib colors of the lines: H and He by species, metals by element (first
  word of the species, e.g. 'Si III' -> 'Si')
"""
_LINE_COLORS = {
                'H I'   : "tab:red",
                'He I'  : "tab:blue",
                'He II' : "tab:purple",
                'C'     : "tab:brown",
                'N'     : "tab:orange",
                'O'     : "tab:cyan",
                'Mg'    : "tab:olive",
                'Si'    : "tab:green",
               }
