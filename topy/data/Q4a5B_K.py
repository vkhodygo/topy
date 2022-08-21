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

from . import Q4_K, Q4bar_K
from ..utils import get_logger

logger = get_logger(__name__)

# ===========================================================
# === Stiffness matrix of a square 4 node 'Q4a5B' element ===
# ===========================================================
# This element is based on the '5-beta' assumed stress element for plane
# stress, but elemental parameters are introduced and selected such that
# spurious zero energy modes are not introduced, for which an investigation
# of characteristic equations of the elemental stiffness matrix is needed.
# Element thickness set = 1. See De Klerk and Groenwold for details.
def create_K(_L, _E, _nu, _k, _t):
    logger.info('Generating K for Q4a5B...')
    logger.info('Q4a5B uses Q4 and Q4bar for its creation, so they will also be created...')

    # Symbolic value of alpha_opt for bending:
    alpha2D = (2 * _L**2 * (1 - _nu) * (2 * _nu**2 - _nu + 1)) \
              / (3 * (_nu + 1) * _E**2)

    K_Q4, B_Q4, C_Q4 = Q4_K.create_K(_L, _E, _nu, _k, _t)
    K_Q4bar, B_Q4bar, C_Q4bar = Q4bar_K.create_K(_L, _E, _nu, _k, _t)

    # Stiffness matrix
    K = K_Q4 - alpha2D * _E * K_Q4bar
    # I'm not sure if this is correct
    B = B_Q4 - alpha2D * _E * B_Q4bar
    C = C_Q4 - alpha2D * _E * C_Q4bar

    logger.info("Created stiffness matrix for Q4a5B.")
    return K, B, C
