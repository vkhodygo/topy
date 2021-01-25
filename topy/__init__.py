"""
# ==============================================================================
# ToPy -- Topology optimization with Python.
# Copyright (C) 2012, 2015, 2016, 2017 William Hunter.
# Copyright (C) 2020, 2021, Tarcísio L. de Oliveira
# ==============================================================================
"""

from .topology_trad import *
from .topology_gen import *
from .visualisation import *
from .elements import *
from .optimisation import *
from .pathfinding import * 

__version__ = "1.0.0"
__author__  = "Tarcisio L. de Oliveira"

__all__ = (
	topology_trad.__all__ +
	topology_gen.__all__ +
	visualisation.__all__ +
	elements.__all__ +
	optimisation.__all__ +
        pathfinding.__all__
)
