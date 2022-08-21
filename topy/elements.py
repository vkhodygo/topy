# -*- coding: utf-8 -*-
"""
# =============================================================================
# Finite element stiffness matrices.
#
# To define your own finite elements, see Python scripts in 'data' directory.
#
# Author: William Hunter, Tarcísio L. de Oliveira
# Copyright (C) 2008, 2015, William Hunter.
# Copyright (C) 2020, 2021, Tarcísio L. de Oliveira
# =============================================================================
"""
from os import path, system

from numpy import array, linspace, unique, sqrt, round, load
from numpy.linalg import eigvalsh

from .utils import get_logger
from .data import *

logger = get_logger(__name__)

__all__ = ['create_element']

# ===================================================
# === Messages used for errors, information, etc. ===
# ===================================================
MSG0 = "finite element stiffness matrix."

MSG1 = "Element stiffness matrices does not exist.\n Created... Please re-run \
your last attempt."


# ==================================================
# === Dictionary of elements - for easy creation ===
# ==================================================
create_element = {
    # 2D elements
    "Q4bar": Q4bar_K.create_K, # KBar of Q4, see De Klerk and Groenwold.
    "Q4"   : Q4_K.create_K,    # Stiffness matrix of a square 4 node plane stress bi-linear element.
    "Q5B"  : Q5B_K.create_K,   # Stiffness matrix of a square 4 node plane stress '5-beta' element.
    "Q4T"  : Q4T_K.create_K,   # Matrix for an element used in 2D thermal problems.
    "Q4a5B": Q4a5B_K.create_K, # Stiffness matrix of a square 4 node 'Q4a5B' element.
    # 3D elements
    "H8"   : H8_K.create_K,    # Stiffness matrix for a hexahedron 8 node tri-linear 3D element.
    "H18B" : H18B_K.create_K,  # Stiffness matrix of a cubic 8 node '18-beta' element.
    "H8T"  : H8T_K.create_K    # Stiffness matrix for a hexahedron 8 node tri-linear 3D element for thermal problems.
}

# EOF elements.py
