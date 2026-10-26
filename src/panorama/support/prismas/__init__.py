"""
Internal (private) helpers and data files of the 'panorama.prismas' module;
their contents may change without notice.

Conventions shared by every object in this module:
    - Nothing is imported here: the public functions of 'panorama.prismas'
      import the helpers they need explicitly from 'prismas_support.py'.
    - Names are private (leading underscore) and are not exposed at the
      package level.
    - The model catalogs ('.dat' files named in '_grid_registry') are stored
      in 'models_grids/', list one model name per line without file
      extension, and are read by '_read_catalog'.
"""
