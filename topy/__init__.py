"""
# ==============================================================================
# ToPy -- Topology optimization with Python.
# Copyright (C) 2012, 2015, 2016, 2017 William Hunter.
# Copyright (C) 2020 Tarcísio L. de Oliveira
# ==============================================================================
"""

from .topology import *
from .visualisation import *
from .elements import *
from .optimisation import *

__version__ = "1.0.0"
__author__  = "Tarcisio L. de Oliveira"

__all__ = (
	topology.__all__ +
	visualisation.__all__ +
	elements.__all__ +
	optimisation.__all__
)
