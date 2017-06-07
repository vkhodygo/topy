# -*- coding: utf-8 -*-
"""
# =============================================================================
# Write the stiffness matrix of finite element to file. The created file name
# is equal to the string between the underscores of *this* file's name, plus a
# 'K' extension, e.g.,
#
#     python ELEM_K.py
#
# gives a file named ELEM.K in the same directory.
#
# Author: William Hunter
# Copyright (C) 2008, 2015, William Hunter.
# =============================================================================
"""

from __future__ import division
from __future__ import print_function

from sympy import symbols, Matrix, diff, integrate, zeros, eye

from numpy import abs, array

from .matlcons import *
from ..helper_functions import my_map

# Get file name:
fname = __file__.split('_')[0] + '.K'
fname = __file__[:-5] + '.K'
print("working on filename: {0}".format(fname))

try:
    f = open(fname)
    print('{0} (stiffness matrix) exists!'.format(fname))
    f.close()
except IOError:
    # SymPy symbols:
    a, b, x, y = symbols('a b x y')
    E, nu = symbols('E nu')
    N1, N2, N3, N4 = symbols('N1 N2 N3 N4')
    xlist = [x, x, x, x, x, x, x, x]
    ylist = [y, y, y, y, y, y, y, y]
    yxlist = [y, x, y, x, y, x, y, x]

    # Shape functions:
    N1 = (a - x) * (b - y) / (4 * a * b)
    N2 = (a + x) * (b - y) / (4 * a * b)
    N3 = (a + x) * (b + y) / (4 * a * b)
    N4 = (a - x) * (b + y) / (4 * a * b)

    # Create strain-displacement matrix B:
    B0 = my_map(diff, [N1, 0, N2, 0, N3, 0, N4, 0], xlist)
    B1 = my_map(diff, [0, N1, 0, N2, 0, N3, 0, N4], ylist)
    B2 = my_map(diff, [N1, N1, N2, N2, N3, N3, N4, N4], yxlist)
    B = Matrix([B0, B1, B2])

    # Create constitutive (material property) matrix for plane stress:
    C = (E / (1 - nu**2)) * Matrix([[1, nu, 0],
                                    [nu, 1, 0],
                                    [0,  0, (1 - nu) / 2]])

    PI = eye(3)
    PH2 = Matrix([[y / b, 0], [0, x / a], [0, 0]])
    P = PI.row_join(PH2)
    tP = P.transpose()

    dJ = tP * B
    dH = tP * C.inv() * P

    # Integration:
    print('SymPy is integrating: K for Q5B...')
    J = dJ.integrate((x, -a, a),(y, -b, b))
    H = dH.integrate((x, -a, a),(y, -b, b))

    # Convert SymPy Matrix to NumPy array:
    K = J.transpose() * H.inv() * J

    # Convert SymPy Matrix to NumPy array:
    K = array(K.subs({a:_a, b:_b, E:_E, nu:_nu})).astype('double')

    # Set small (<< 0) values equal to zero:
    K[abs(K) < 1e-6] = 0

    # Create file:
    K.dump(fname)
    print('Created {0} (stiffness matrix).'.format(fname))

# EOF Q5B_K.py
