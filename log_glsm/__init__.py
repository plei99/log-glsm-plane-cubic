"""Exact log-GLSM calculations, with ordinary Python module namespaces.

Sage implementations keep their .sage source and exact arithmetic. Importing
the Fraction-based curve modules does not require Sage to be installed.
"""

from ._sage import install_importer

install_importer()
del install_importer
