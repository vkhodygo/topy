# -*- coding: utf-8 -*-
"""
# =============================================================================
# Creates the stiffness matrix as requested, using the material properties 
# provided in the TPD file (for v2020 files).
#
# Author: William Hunter, Tarcísio L. de Oliveira
# Copyright (C) 2008, 2015, William Hunter.
# Copyright (C) 2020, 2021, Tarcísio L. de Oliveira
# =============================================================================
"""
from __future__ import division
from __future__ import print_function

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
    x, y = symbols('x y')
    N1, N2, N3, N4 = symbols('N1 N2 N3 N4')
    xlist = [x, x, x, x, x, x, x, x]
    ylist = [y, y, y, y, y, y, y, y]
    yxlist = [y, x, y, x, y, x, y, x]

    # Shape functions:
    N1 = (_a - x) * (_b - y) / (4 * _a * _b)
    N2 = (_a + x) * (_b - y) / (4 * _a * _b)
    N3 = (_a + x) * (_b + y) / (4 * _a * _b)
    N4 = (_a - x) * (_b + y) / (4 * _a * _b)

    # Create strain-displacement matrix B:

    B0 = my_map(diff, [N1,  0, N2,  0, N3,  0, N4,  0],  xlist)
    B1 = my_map(diff, [0,  N1,  0, N2,  0, N3,  0, N4],  ylist)
    B2 = my_map(diff, [N1, N1, N2, N2, N3, N3, N4, N4], yxlist)


    B = Matrix([B0, B1, B2])

    # Create constitutive (material property) matrix for plane stress:
    C = (_E / (1 - _nu**2)) * Matrix([[1, _nu, 0],
                                      [_nu, 1, 0],
                                      [0,  0, (1 - _nu) / 2]])

    dK = B.T * C * B

    # Integration:
    logger.info('SymPy is integrating: K for Q4...')
    K = dK.integrate((x, -_a, _a),(y, -_b, _b))

    # Convert SymPy Matrix to NumPy array:
    K = _t * array(K, dtype='double')
    C = array(C, dtype='double')

    # Set small (<< 0) values equal to zero:
    K[abs(K) < 1e-6] = 0

    # Return result:
    logger.info('Created stiffness matrix.')
    return K, B, C

    print('Created {0} (stiffness matrix).'.format(fname))


# EOF Q4_K.py
