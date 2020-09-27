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
    a, b, x, y = symbols('a b x y')
    E, nu = symbols('E nu')
    N1, N2, N3, N4 = symbols('N1 N2 N3 N4')
    k = symbols('k')
    xlist = [x, x, x, x]
    ylist = [y, y, y, y]

    # Shape functions:
    N1 = (a - x) * (b - y) / (4 * a * b)
    N2 = (a + x) * (b - y) / (4 * a * b)
    N3 = (a + x) * (b + y) / (4 * a * b)
    N4 = (a - x) * (b + y) / (4 * a * b)

    # Create strain-displacement matrix B:
    B0 = tuple(map(diff, [N1, N2, N3, N4], xlist))
    B1 = tuple(map(diff, [N1, N2, N3, N4], ylist))
    B = Matrix([B0, B1])

    # Create conductivity matrix:
    C = Matrix([[k, 0],
                [0, k]])

    dK = B.T * C * B

    # Integration:
    logger.info('SymPy is integrating: K for Q4T...')
    K = dK.integrate((x, -a, a),(y, -b, b))

    # Convert SymPy Matrix to NumPy array:
    K = _t * array(K.subs({a:_a, b:_b, k:_k})).astype('double')
    B = B.subs({a:_a, b:_b, k:_k})
    C = array(C.subs({a:_a, b:_b, k:_k})).astype('double')

    # Set small (<< 0) values equal to zero:
    K[abs(K) < 1e-6] = 0

    # Return result:
    logger.info('Created stiffness matrix.')
    return K, B, C

# EOF Q4T_K.py
