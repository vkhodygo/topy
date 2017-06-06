# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse as sp_sparse

def my_map(*args):
    return list(map(*args))


################################################################################
## https://stackoverflow.com/questions/13077527/is-there-a-numpy-delete-equivalent-for-sparse-matrices
##
def identity_minus_rows(N, rows):
    if np.isscalar(rows):
        rows = [rows]
    J = sp_sparse.diags(np.ones(N), 0).tocsr()  # make a diag matrix
    #print("J.shape: {0}".format(J.shape))

    indices = []
    for i_ in range(len(rows)):
        if rows[i_] == 0: indices.append(i_)
    J = delete_rows_csr(J, indices)
    ####
    ## delete from back to front to keep ID's intact
    #for r in sorted(rows, reverse=True):
        #J = delete_row_csr(J, r)
    return J

def delete_rows_csr(mat, indices):
    """
    Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
    """
    if not isinstance(mat, sp_sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]

#def delete_row_csr(mat, i):
    ##print("type of matrix: {0}".format(type(mat)))
    #if not isinstance(mat, sp_sparse.csr.csr_matrix):
        #raise ValueError("works only for CSR format -- use .tocsr() first")
    #n = mat.indptr[i+1] - mat.indptr[i]
    #if n > 0:
        #mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        #mat.data = mat.data[:-n]
        #mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        #mat.indices = mat.indices[:-n]
    #mat.indptr[i:-1] = mat.indptr[i+1:]
    #mat.indptr[i:] -= n
    #mat.indptr = mat.indptr[:-1]
    #mat._shape = (mat._shape[0]-1, mat._shape[1])
    #return mat

#def delete_row_lil(mat, i):
    #if not isinstance(mat, scipy.sparse.lil_matrix):
        #raise ValueError("works only for LIL format -- use .tolil() first")
    #mat.rows = np.delete(mat.rows, i)
    #mat.data = np.delete(mat.data, i)
    #mat._shape = (mat._shape[0] - 1, mat._shape[1])
    #return mat


################################################################################
##
##
def update_add_mask_sym(K, updatedKe, e2sdofmap, mask):
    #K.update_add_mask_sym(updatedKe, e2sdofmap, mask)
    #print(e2sdofmap)
    #print(type(K))#[e2sdofmap,e2sdofmap])

    list1 = np.array(e2sdofmap)[:,np.newaxis]
    list2 = np.array(e2sdofmap)
    K[list1,list2] += updatedKe

    #for i in range(len(e2sdofmap)):
        #for j in range(len(e2sdofmap)):
            #if mask[i]:
                #K[e2sdofmap[i],e2sdofmap[j]] += updatedKe[i,j]
    return K