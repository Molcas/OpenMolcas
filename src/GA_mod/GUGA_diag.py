"""
This module evaluates the diagonal Hamiltonian matrix elements in the GUGA framework.
The equations can be found in Appendix A.2.1 of the PhD thesis of Werner Dobrautz 
(https://dx.doi.org/10.18419/opus-10593).

Copyright (c) 2025 Maru Song
"""

import numpy as np

class DiagElement:
    def __init__(self, norb: int, IntegralClass):
        self.norb = norb
        self.E_core = IntegralClass.core_energy
        self.t_ii = np.array([IntegralClass.get_integral(i, i, 0, 0) for i in range(1, self.norb + 1)])
        self.v_iiii = np.array([IntegralClass.get_integral(i, i, i, i) for i in range(1, self.norb + 1)])
        self.v_iijj = np.array([[IntegralClass.get_integral(i, i, j, j) for j in range(1, self.norb + 1)] for i in range(1, self.norb + 1)])
        self.v_ijji = np.array([[IntegralClass.get_integral(i, j, j, i) for j in range(1, self.norb + 1)] for i in range(1, self.norb + 1)])

    def stepvec_to_occ(self, CSF_stepvec: np.array) -> np.array:
        return CSF_stepvec - (CSF_stepvec // 2)

    def calc_diag_elem(self, CSF_stepvec: np.array, add_core=True) -> float:
        CSF_stepvec = np.asarray(CSF_stepvec)
        if CSF_stepvec.ndim != 1:
            raise ValueError(f"CSF_stepvec must be a 1D array, got shape {CSF_stepvec.shape}")
        CSF_occ = self.stepvec_to_occ(CSF_stepvec)
        one_body = self.one_body_contrib(CSF_occ)
        two_body = self.two_body_contrib(CSF_stepvec, CSF_occ)

        if add_core:
            return one_body + two_body + self.E_core
        else:
            return one_body + two_body

    def one_body_contrib(self, CSF_occ: np.array) -> float:
        contrib = np.dot(CSF_occ, self.t_ii)

        return contrib

    def two_body_contrib(self, CSF_stepvec: np.array, CSF_occ: np.array) -> float:
        Nmat = self.N_matrix(CSF_occ)
        Xmat = self.X_matrix(CSF_stepvec)

        # Case1: i=j=k=l
        case1_contrib = np.sum(self.v_iiii[CSF_stepvec == 3])

        # Case2: i=j!=k=l
        case2_contrib = np.sum(np.triu(self.v_iijj * Nmat, k=1))

        # Case3: i=l!=j=k
        case3_contrib = -np.sum(np.triu(self.v_ijji * (Xmat + Nmat) / 2, k=1))

        return case1_contrib + case2_contrib + case3_contrib

    def X_matrix(self, d_vec):
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

    def N_matrix(self, CSF_occ: np.array) -> np.array:
        return np.outer(CSF_occ, CSF_occ)
