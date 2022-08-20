"""Common utilities."""
import logging
import sys

import numpy as np
import scipy as sp


def get_logger(name: str) -> logging.Logger:
    """Return a `Logger` instance for `name`."""
    logger = logging.getLogger(name)
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.setLevel(logging.DEBUG)
    return logger


def update_add_mask(
    A: sp.sparse.spmatrix,
    B: sp.sparse.spmatrix,
    ind0: np.ndarray,
    ind1: np.ndarray,
    mask0: np.ndarray,
    mask1: np.ndarray,
) -> sp.sparse.csr_matrix:
    """
    Add two sparse matrices with a boolean mask on `B` and a permutation on `A`.

    This function is meant to replace `pysparse.spmatrix.ll_mat.update_add_mask` (http://pysparse.sourceforge.net/spmatrix.html#spmatrix.ll_mat.update_add_mask).

    Parameters
    ----------
    A : sp.sparse.spmatrix
        A sparse matrix whose rows and columns will be permuted by `ind`.
    B : sp.sparse.spmatrix
        A sparse matrix whose rows and columns will be masked by `mask`. Must have the same dimensions as `A`.
    ind0: np.ndarray
        A column vector representing a permutation of the rows a`A`. Must have the same number of rows as `A`.
    ind1: np.ndarray
        A column vector representing a permutation of the columns of `A`. Must have the same number of columns as `A`.
    mask0: np.ndarray
        A column vector of 0's and 1's by which indicates the rows in `B` to be zeroed out. Must have the same number of rows as `B`.
    mask1: np.ndarray
        A column vector of 0's and 1's by which indicates the columns in `B` to be zeroed out. Must have the same number of columns as `B`.

    Returns
    -------
    sp.sparse.csr_matrix
        The matrix result of the operation.
    """
    # Apply mask to B. We convert it to CSR for fast product.
    masked_B = B.tocsr() * mask0.T * mask1
    # Permute the rows of A. We convert to CSR for efficient row slicing and to CSC for efficient column slicing.
    permuted_A = permute(A, ind0, ind1)
    # Add the matrices. The entries of A are still permuted.
    permuted_sum = permuted_A + masked_B
    # Let's revert the permutation and return it.
    return permute(permuted_sum, inverse_permutation(ind0), inverse_permutation(ind1))


def update_add_mask_sym(
    A: sp.sparse.spmatrix, B: sp.sparse.spmatrix, ind: np.ndarray, mask: np.ndarray
) -> sp.sparse.csr_matrix:
    """
    Add two symmetric sparse matrices with a boolean mask on `B` and a permutation on `A`. Does not change the inputs.

    This is the symmetric version of `update_add_mask`. This function is meant to replace `pysparse.spmatrix.ll_mat.update_add_mask_sym` (http://pysparse.sourceforge.net/spmatrix.html#spmatrix.ll_mat.update_add_mask_sym).

    Parameters
    ----------
    A : sp.sparse.spmatrix
        A sparse square matrix whose rows and columns will be permuted by `ind`.
    B : sp.sparse.spmatrix
        A sparse matrix whose rows and columns will be masked by `mask`. Must have the same dimensions as `A`.
    ind : np.ndarray
        A column vector representing a permutation of the indices of rows and columns of `A`. Must have the same number of rows (and columns) as `A`.
    mask : np.ndarray
        A column vector of 0's and 1's by which indicates the rows and columns in `B` to be zeroed out. Must have the same number of rows (and columns) as `A`.

    Returns
    -------
    sp.sparse.csr_matrix
        The matrix result of the operation.
    """
    return update_add_mask(A, B, ind, ind, mask, mask)


def permute(
    A: sp.sparse.spmatrix, ind0: np.ndarray, ind1: np.ndarray,
) -> sp.sparse.csr_matrix:
    """
    Return the matrix A with its rows permuted by `ind0` and columns permuted by `ind1`.

    Converts the matrix to CSR for efficient row slicing and to CSC for efficient column slicing.

    Args:
        A (sp.sparse.spmatrix): The sparse matrix to permute.
        ind0 (np.ndarray): The permutation on the rows of A.
        ind1 (np.ndarray): The permutation on the columns of A.

    Returns:
        sp.sparse.csr_matrix: The permuted sparse matrix.
    """
    return A.tocsc()[:, ind1].tocsr()[ind0, :]


def inverse_permutation(p: np.ndarray) -> np.ndarray:
    """
    Return an array `s`, where `s[i]` gives the index of `i` in `p`.

    The argument `p` is assumed to be some permutation of `0, 1, ..., len(p)-1`.

    Code obtained from https://stackoverflow.com/a/25535723/9954163.
    """
    s = np.empty(p.size, p.dtype)
    s[p] = np.arange(p.size)
    return s
