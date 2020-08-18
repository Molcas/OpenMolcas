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
# Copyright (C) 2020, Oskar Weser                                      *
#***********************************************************************

from numpy import eye
from scipy.linalg import svd, det

def get_kabsch_transformation(Q, P, force_rotation=False):
    """Calculate the optimal orthogonal transformation from ``P`` onto ``Q``.

    If force_rotation is true, the transformation
    is from the special orthogonal group.
    """
    # Naming of variables follows the wikipedia article:
    # http://en.wikipedia.org/wiki/Kabsch_algorithm
    A = P.T @ Q
    # One can't initialize an array over its transposed. :-(
    U, S, V = svd(A)
    V = V.T
    if force_rotation:
        sign = eye(W)
        # NOTE: W and V are unitary. Up to float noise we know that:
        #   det(W @ V.T) in {-1, +1}
        sign[-1, -1] = 1.0 if det(V @ U.T) > 0.0 else -1.0
        return V @ sign @ U.T
    else:
        return V @ U.T
