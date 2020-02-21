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
