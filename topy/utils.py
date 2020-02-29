# -*- coding: utf-8 -*-
"""Common utilities."""
import itertools
import logging
import sys

from datashape import dshape


def get_logger(name: str) -> logging.Logger:
    """Return a `Logger` instance for `name`."""
    logger = logging.getLogger(name)
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.setLevel(logging.DEBUG)
    return logger


def update_add_mask_sym(
    A: dshape("M... * M * complex"),
    B: dshape("M... * M * complex"),
    ind: dshape("M... * complex"),
    mask: dshape("M... * complex"),
    # symmetric: bool = True,
) -> dshape("M... * M * complex"):
    """
    Assemble a symmetric global finite element matrix. Changes `A` in-place.

    source: http://pysparse.sourceforge.net/spmatrix.html#spmatrix.ll_mat.update_add_mask_sym
    """
    for (i, ind_i), (j, ind_j) in itertools.product(enumerate(ind), repeat=2):
        if mask[i]:
            A[ind_i, ind_j] += B[i, j]

    return A
