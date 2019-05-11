import numpy as np
from numpy import dot
from numpy.linalg import svd
import scipy as sp
from scipy.linalg import det, norm

from attr import attrs, attrib

def procrust(A, B):
    """Calculate the optimal orthogonal transformation from ``A`` to ``B``.

    The algorithm is described very well in
    `https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem>`_.

    Args:
        A (~numpy.array):
        B (~numpy.array):

    Returns:
        :class:`~numpy.array`: Rotation matrix
    """
    # Naming of variables follows the wikipedia article:
    # https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
    M = B @ A.T
    # One can't initialize an array over its transposed
    U, S, V = svd(M)
    V = V.T
    return U @ V.T

@attrs
class RasOrb:
    orbs = attrib()
    coeff = attrib()
    occ = attrib()
    energy = attrib()
    idx = attrib()

    @classmethod
    def read_orbfile(cls, path):
        coeff, occ, energy, idx = read_orbfile(path)
        orbs = [len(x) for x in coeff]
        return cls(orbs, coeff, occ, energy, idx)


def _forward(f, n):
    for i in range(n):
        line = f.readline()
    return line

def _read_coeffs(f, orbs, cols):
    MO_coeff = [np.empty((n_orbs, n_orbs)) for n_orbs in orbs]
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        for orb in range(n_orbs):
            values = []
            for row in range(rows):
                line = f.readline()
                values.extend([float(x) for x in line.split()])
            MO_coeff[irrep][orb, :] = values
            next(f)
    return MO_coeff

def _read_orb_property(f, orbs, cols):
    MO_v = [np.empty(n_orbs) for n_orbs in orbs]
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        values = []
        for row in range(rows):
            line = f.readline()
            values.extend([float(x) for x in line.split()])
        MO_v[irrep][:] = values
    return MO_v

def _read_CAS_idx(f, orbs, cols):
    idx = []
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        values = []
        next(f)
        for row in range(rows):
            line = f.readline()
            values.extend(line.strip()[2:])
        idx.append(values)
    return idx

def read_orbfile(path):
    with open(path, 'r') as f:
        _forward(f, 5)
        line = f.readline()
        orbs = [int(irrep) for irrep in line.split()]

        line = f.readline()
        while 'ORBITAL' not in line:
            line = f.readline()
        MO_coeff = _read_coeffs(f, orbs, 5)

        line = f.readline()
        while 'OCCUPATION NUMBERS' not in line:
            line = f.readline()
        cols = 5
        MO_occ = _read_orb_property(f, orbs, 5)

        line = f.readline()
        while 'ONE ELECTRON ENERGIES' not in line:
            line = f.readline()
        MO_E = _read_orb_property(f, orbs, 10)

        line = f.readline()
        while 'INDEX' not in line:
            line = f.readline()
        idx = _read_CAS_idx(f, orbs, 10)
    return MO_coeff, MO_occ, MO_E, idx

def _reshape_output(v, cols, formatter):
    rows = (len(v) + cols - 1) // cols
    lines = []
    for i in range(rows - 1):
        current = v[(i * cols):((i + 1) * cols)]
        lines.append(''.join(f'{x:{formatter}}' for x in current))
    current = v[(rows - 1) * cols:]
    lines.append(''.join(f'{x:{formatter}}' for x in current))
    return lines

def _get_header(orbs):
    return f"""#INPORB 2.2
#INFO
* PyOrb written
{0:8d}{len(orbs):8d}{0:8d}
{''.join(f"{x:8d}" for x in orbs)}
{''.join(f"{x:8d}" for x in orbs)}
* PyOrb written
#ORB"""


def _write_MO_coeff(orbs, MO_coeff, file):
    res_out = lambda *args : '\n'.join(_reshape_output(*args))
    for irrep, n_orbs in enumerate(orbs):
        for orb in range(n_orbs):
            print(f'* ORBITAL{irrep + 1:5d}{orb + 1:5d}', file=file)
            print(res_out(MO_coeff[irrep][orb, :] , 5, '22.14E'), file=file)

def _write_MO_occ(orbs, MO_occ, file):
#     res_out = _reshape_output
    res_out = lambda *args : '\n'.join(_reshape_output(*args))
    print('#OCC', file=file)
    print('* OCCUPATION NUMBERS', file=file)
    for irrep, n_orbs in enumerate(orbs):
        print(res_out(MO_occ[irrep][:], 5, '22.14E'), file=file)
    print('#OCHR', file=file)
    print('* OCCUPATION NUMBERS (HUMAN-READABLE)', file=file)
    for irrep, n_orbs in enumerate(orbs):
        print(res_out(MO_occ[irrep][:], 10, '8.4f'), file=file)

def _write_MO_E(orbs, MO_E, file):
#     res_out = _reshape_output
    res_out = lambda *args : '\n'.join(_reshape_output(*args))
    print('#ONE', file=file)
    print('* ONE ELECTRON ENERGIES', file=file)
    for irrep, n_orbs in enumerate(orbs):
        print(res_out(MO_E[irrep][:], 10, '12.4E'), file=file)

def _write_CAS_idx(orbs, idx, file):
    def res_out(v):
        lines = _reshape_output(v, 10, '1')
        return '\n'.join(f'{i % 10} {a}'for i, a in enumerate(lines))

    print('#INDEX', file=file)
    for irrep, n_orbs in enumerate(orbs):
        print('* 1234567890', file=file)
        print(res_out(idx[irrep]), file=file)

def write_orbfile(MO_coeff, MO_occ, MO_E, idx, path):
    orbs = [len(x) for x in MO_coeff]
    with open(path, 'w') as f:
        print(_get_header(orbs), file=f)
        _write_MO_coeff(orbs, MO_coeff, file=f)
        _write_MO_occ(orbs, MO_occ, file=f)
        _write_MO_E(orbs, MO_E, file=f)
        _write_CAS_idx(orbs, idx, file=f)
