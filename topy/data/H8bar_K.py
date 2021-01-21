"""
# =============================================================================
# Creates the stiffness matrix as requested, using the material properties 
# provided in the TPD file (for v2020 files).
#
# Author: William Hunter, Tarcísio L. de Oliveira
# Copyright (C) 2008, 2015, William Hunter.
# Copyright (C) 2020, Tarcísio L. de Oliveira
# =============================================================================
"""
from __future__ import division

import os

from sympy import symbols, Matrix, diff, integrate, zeros
from numpy import abs, array

from ..utils import get_logger

logger = get_logger(__name__)

def create_K(_L, _E, _nu, _k, _t):
    # Initialize variables
    _a, _b, _c = _L, _L, _L  # element dimensions (half-lengths)
    _G = _E / (2 * (1 + _nu))  # modulus of rigidity
    _g = _E /  ((1 + _nu) * (1 - 2 * _nu))

    # SymPy symbols:
    x, y, z = symbols('x y z')
    N1, N2, N3, N4 = symbols('N1 N2 N3 N4')
    N5, N6, N7, N8 = symbols('N5 N6 N7 N8')
    o = symbols('o') #  dummy symbol
    xlist = [x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x]
    ylist = [y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y, y]
    zlist = [z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z]
    yxlist = [y, x, o, y, x, o, y, x, o, y, x, o, y, x, o, y, x, o, y, x, o, y, x, o]
    zylist = [o, z, y, o, z, y, o, z, y, o, z, y, o, z, y, o, z, y, o, z, y, o, z, y]
    zxlist = [z, o, x, z, o, x, z, o, x, z, o, x, z, o, x, z, o, x, z, o, x, z, o, x]

    # Shape functions:
    N1 = (_a - x) * (_b - y) * (_c - z) / (8 * _a * _b * _c)
    N2 = (_a + x) * (_b - y) * (_c - z) / (8 * _a * _b * _c)
    N3 = (_a + x) * (_b + y) * (_c - z) / (8 * _a * _b * _c)
    N4 = (_a - x) * (_b + y) * (_c - z) / (8 * _a * _b * _c)
    N5 = (_a - x) * (_b - y) * (_c + z) / (8 * _a * _b * _c)
    N6 = (_a + x) * (_b - y) * (_c + z) / (8 * _a * _b * _c)
    N7 = (_a + x) * (_b + y) * (_c + z) / (8 * _a * _b * _c)
    N8 = (_a - x) * (_b + y) * (_c + z) / (8 * _a * _b * _c)

    # Create strain-displacement matrix B:
    B0 = tuple(map(diff, [N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0,\
                    N5, 0, 0, N6, 0, 0, N7, 0, 0, N8, 0, 0], xlist))
    B1 = tuple(map(diff, [0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0,\
                    0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8, 0], ylist))
    B2 = tuple(map(diff, [0, 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4,\
                    0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8], zlist))
    B3 = tuple(map(diff, [N1, N1, N1, N2, N2, N2, N3, N3, N3, N4, N4, N4,\
                    N5, N5, N5, N6, N6, N6, N7, N7, N7, N8, N8, N8], yxlist))
    B4 = tuple(map(diff, [N1, N1, N1, N2, N2, N2, N3, N3, N3, N4, N4, N4,\
                    N5, N5, N5, N6, N6, N6, N7, N7, N7, N8, N8, N8], zylist))
    B5 = tuple(map(diff, [N1, N1, N1, N2, N2, N2, N3, N3, N3, N4, N4, N4,\
                    N5, N5, N5, N6, N6, N6, N7, N7, N7, N8, N8, N8], zxlist))
    B = Matrix([B0, B1, B2, B3, B4, B5])

    # Create constitutive (material property) matrix:
    C = Matrix([[(1 - _nu) * _g, _nu * _g, _nu * _g, 0, 0, 0],
                [_nu * _g, (1 - _nu) * _g, _nu * _g, 0, 0, 0],
                [_nu * _g, _nu * _g, (1 - _nu) * _g, 0, 0, 0],
                [0, 0, 0,                           _G, 0, 0],
                [0, 0, 0,                           0, _G, 0],
                [0, 0, 0,                           0, 0, _G]])

    CB = C * B

    # Create delB matrix:
    delCB0x = array(tuple(map(diff, CB[0, :], xlist)))
    delCB0y = array(tuple(map(diff, CB[3, :], ylist)))
    delCB0z = array(tuple(map(diff, CB[5, :], zlist)))
    delCB0 = delCB0x + delCB0y + delCB0z

    delCB1y = array(tuple(map(diff, CB[1, :], ylist)))
    delCB1x = array(tuple(map(diff, CB[3, :], xlist)))
    delCB1z = array(tuple(map(diff, CB[4, :], zlist)))
    delCB1 = delCB1x + delCB1y + delCB1z

    delCB2z = array(tuple(map(diff, CB[2, :], zlist)))
    delCB2y = array(tuple(map(diff, CB[4, :], ylist)))
    delCB2x = array(tuple(map(diff, CB[5, :], xlist)))
    delCB2 = delCB2x + delCB2y + delCB2z

    Bbar = Matrix([delCB0.tolist(), delCB1.tolist(), delCB2.tolist()])

    dKbar = Bbar.T * Bbar #  a matrix of constants, i.e., no x, y or z vals

    # Integration:
    logger.info('SymPy is integrating: K for H8bar...')
    Kbar = dKbar.integrate((x, -_a, _a),(y, -_b, _b),(z, -_c, _c))

    # Convert SymPy Matrix to NumPy array:
    K = array(Kbar, dtype='double')
    C = array(C, dtype='double')

    # Set small (<< 0) values equal to zero:
    K[abs(K) < 1e-6] = 0

    # Return result:
    logger.info('Created stiffness matrix.')
    return K, B, C

# EOF H8bar_K.py
