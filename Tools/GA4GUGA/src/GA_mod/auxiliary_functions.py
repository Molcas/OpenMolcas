#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2025, Maru Song                                        *
#***********************************************************************

"""
Auxiliary functions
"""

import numpy as np

def J_mat_from_fcidump(FCIDUMPClass, norb):
    """
    Generate the J matrix from the FCIDUMP class.
    We onlny use exchange integrals.
    """
    J = np.zeros((norb, norb))
    for i in range(norb):
        for j in range(i + 1, norb):
            J[i, j] = -FCIDUMPClass.get_integral(j + 1, i + 1, j + 1, i + 1)
            J[j, i] = -FCIDUMPClass.get_integral(j + 1, i + 1, j + 1, i + 1)

    return J

def X_matrix_openshell_only(d_vec):
    """
    Compute the matrix X_ij for a given step-vector d.
    See equation (A.10) in the appendix (A.2) of Werner Dobrautz's phd thesis.
    The exact form used here will be found in our paper (update this once the
    paper is published).
    **Note that this assumes the step-vector only contains 1's and 2's.**

    Args:
        d_vec (list): step-vector not including 0 and 3

    Returns:
        N x N array where X[i][j] corresponds to the computed value for i < j.
    """

    N = len(d_vec)
    d_vec = np.array(d_vec)
    s_vec = 2 * (d_vec == 1) - 1
    b_vec = np.cumsum(s_vec)

    A = (b_vec - 2 * d_vec + 4) / (b_vec + 2 * d_vec - 2)
    B = (b_vec + 4 * d_vec - 5) / (b_vec + 1)
    X = np.ones((N, N))

    for i, j in zip(*np.triu_indices(N, k=1)):
        X[i, j] = X[i, j - 1] * A[j - 1] * B[j] 
    X = np.sqrt(X)
    np.fill_diagonal(X, 0)
    for i, j in zip(*np.triu_indices(N, k=1)):
        X[i, j] *= s_vec[i] * s_vec[j]
        X[j, i] = X[i, j]

    return X
