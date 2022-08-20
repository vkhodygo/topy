# -*- coding: utf-8 -*-
"""Common utilities."""
import itertools
import logging
import sys
import os

from datashape import dshape
from scipy.sparse import lil_matrix
from scipy.sparse import linalg

def get_logger(name: str) -> logging.Logger:
    """Return a `Logger` instance for `name`."""
    logger = logging.getLogger(name)
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.setLevel(logging.DEBUG)
    return logger



# ===================================
# === Matrix methods and helpers ===
# ===================================


def update_add_mask_sym(
    A: dshape("M... * M * complex"),
    B: dshape("M... * M * complex"),
    ind: dshape("M... * complex"),
    mask: dshape("M... * complex"),
    # symmetric: bool = True,
) -> dshape("M... * M * complex"):
    """
    Assemble a symmetric global finite element matrix. Changes `A` in-place.

    :param A: MxM matrix.
    :param B: MxM matrix.
    :param ind: Mx1 matrix.
    :param mask: Mx1 matrix.
    :returns: A matrix with the same dimensions as A and B.

    source: http://pysparse.sourceforge.net/spmatrix.html#spmatrix.ll_mat.update_add_mask_sym
    """
    for (i, ind_i), (j, ind_j) in itertools.product(enumerate(ind), repeat=2):
        if mask[i]:
            A[ind_i, ind_j] += B[i, j]

    return A


def precondition_sparse_matrix(A: lil_matrix) -> linalg.LinearOperator:
    """Compute an approximate inverse of `A` using incomplete LU decomposition."""
    ilu = linalg.spilu(A)
    Mx = ilu.solve
    return linalg.LinearOperator(A.shape, Mx)




def get_data_file(source_file_name):
    # type: (str) -> (str)
    """Return the data file path to store result."""
    path = list(os.path.split(source_file_name))
    path[-1] = path[-1].split('_')[0] + '.K'
    return os.path.join(*path)
