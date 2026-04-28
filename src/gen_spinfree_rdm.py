"""
This module evaluates the diagonal Hamiltonian matrix elements in the GUGA framework.
The equations can be found in Appendix A.2.1 of the PhD thesis of Werner Dobrautz 
(https://dx.doi.org/10.18419/opus-10593).

Copyright (c) 2025 Maru Song
"""

import numpy as np

def stepvec_to_occ(CSF_stepvec: np.array) -> np.array:
    return CSF_stepvec - (CSF_stepvec // 2)

def spin_free_1rdm(CSF_stepvec: np.array) -> np.array:
    CSF_stepvec = np.asarray(CSF_stepvec)
    if CSF_stepvec.ndim != 1:
        raise ValueError(f"CSF_stepvec must be a 1D array, got shape {CSF_stepvec.shape}")
    CSF_occ = stepvec_to_occ(CSF_stepvec)
    # Currently, just return the diagonal elements
    return CSF_occ
    
def spin_free_2rdm(CSF_stepvec: np.array) -> np.array:
    CSF_stepvec = np.asarray(CSF_stepvec)
    if CSF_stepvec.ndim != 1:
        raise ValueError(f"CSF_stepvec must be a 1D array, got shape {CSF_stepvec.shape}")
    CSF_occ = stepvec_to_occ(CSF_stepvec)
    Nmat = N_matrix(CSF_occ)

    # i = j != k = l
    d_iijj = Nmat.copy()
    np.fill_diagonal(d_iijj, d_iijj.diagonal() - CSF_occ)

    # i = l != j = k
    d_ijji = -d_iijj.copy() - X_matrix(CSF_stepvec)
    np.fill_diagonal(d_ijji, d_ijji.diagonal() - CSF_occ)

    return d_iijj / 2, d_ijji / 2

def X_matrix(d_vec):
    N = len(d_vec)

    db = lambda dv: dv - 3 * (dv // 2)
    A2 = lambda b, x, y: (b + x) / (b + y)

    d_vec = np.array(d_vec)
    b_vec = np.cumsum(db(d_vec))
    X = np.zeros((N, N))

    def f(b, d):
        if d == 1:
            return A2(b, 2, 0) * A2(b, -1, 1)
        elif d == 2:
            return A2(b, 0, 2) * A2(b, 3, 1)
        else:
            return 1

    for i, j in zip(*np.triu_indices(N, k=1)):
        di, dj = d_vec[i], d_vec[j]
        bi, bj = b_vec[i], b_vec[j]
        if di == 1 and dj == 1:
            Xij = A2(bi, 2, 0) * A2(bj, -1, 1)
        elif di == 1 and dj == 2:
            Xij = A2(bi, 2, 0) * A2(bj, 3, 1)
        elif di == 2 and dj == 1:
            Xij = A2(bi, 0, 2) * A2(bj, -1, 1)
        elif di == 2 and dj == 2:
            Xij = A2(bi, 0, 2) * A2(bj, 3, 1)
        else:
            Xij = 0

        if Xij != 0:
            for k in range(i + 1, j):
                Xij *= f(b_vec[k], d_vec[k])
            Xij = np.sqrt(Xij)

        if di * dj == 2:
            Xij *= -1.0

        X[i, j] = Xij
        X[j, i] = Xij

    return X

def N_matrix(CSF_occ: np.array) -> np.array:
    return np.outer(CSF_occ, CSF_occ)

def write_spinfree_1rdm(CSF_stepvec: np.array, filename='DMAT.1', thr: float = 1e-12) -> None:
    """
    Write a Spin-free 1-RDM to a file with triangular indexing.

    Parameters:
    -----------
    CSF_stepvec : np.array
        Input CSF step vector
    filename : str
        Output filename
    thr : float
        Threshold for writing matrix elements (only write if |value| > thr)
        
    The output file has two columns:
    - First column: n(n+1)/2 where n is the 1-based index (1, 3, 6, 10, ...)
    - Second column: corresponding element from the input array (if above threshold)
    """
    CSF_stepvec = np.asarray(CSF_stepvec)
    if CSF_stepvec.ndim != 1:
        raise ValueError(f"CSF_stepvec must be a 1D array, got shape {CSF_stepvec.shape}")
    one_rdm = spin_free_1rdm(CSF_stepvec)

    with open(filename, 'w') as f:
        for i, value in enumerate(one_rdm):
            n = i + 1  # 1-based index
            triangular_index = n * (n + 1) // 2
            value = float(value)
            
            # Only write if magnitude is greater than threshold
            if abs(value) > thr:
                f.write(f"{triangular_index:6d}{value:25.17f}\n")

def write_spinfree_2rdm(CSF_stepvec: np.array, thr: float = 1e-12) -> None:
    CSF_stepvec = np.asarray(CSF_stepvec)
    if CSF_stepvec.ndim != 1:
        raise ValueError(f"CSF_stepvec must be a 1D array, got shape {CSF_stepvec.shape}")
    norb = len(CSF_stepvec)
    d_iijj, d_ijji = spin_free_2rdm(CSF_stepvec)
    f = lambda x: x * (x + 1) // 2

    with open('PSMAT.1', 'w') as file:
        counter = 1
        for i in range(1, norb + 1):
            # ijji
            for j in range(1, i + 1):
                # iiii
                if i == j:
                    if CSF_stepvec[i - 1] == 3:
                        value = 1.0
                    else:
                        value = 0.0
                else:
                    value = d_ijji[i - 1, j - 1] / 2
                
                if abs(value) > thr:
                    file.write(f"{f(counter):6d}{value:25.17f}\n")
                counter += 1

            # iijj
            for j in range(1, i):
                value = d_iijj[i - 1, j - 1]
                index = f(f(i) - 1) + f(j)
                
                if abs(value) > thr:
                    file.write(f"{index:6d}{value:25.17f}\n")

    with open('PAMAT.1', 'w') as file:
        counter = 1
        for i in range(1, norb + 1):
            # ijji
            for j in range(1, i + 1):
                if i == j:
                    value = 0.0
                else:
                    value = -d_ijji[i - 1, j - 1] / 2
                
                if abs(value) > thr:
                    file.write(f"{f(counter):6d}{value:25.17f}\n")
                counter += 1


if __name__ == "__main__":
    import sys

    # Read CSF_stepvec from terminal arguments
    if len(sys.argv) < 2:
        print("Usage: python GUGA_diag.py <CSF_stepvec elements separated by space>")
        sys.exit(1)
    CSF_stepvec = np.array([int(x) for x in sys.argv[1:]])
    write_spinfree_1rdm(CSF_stepvec, filename='DMAT.1', thr=1e-12)
    write_spinfree_2rdm(CSF_stepvec, thr=1e-12)
