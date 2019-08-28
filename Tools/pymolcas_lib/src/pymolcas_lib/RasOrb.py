import numpy as np
from numpy.linalg import svd

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

    def copy(self):
        return self.__class__(
            orbs=self.orbs.copy(),
            coeff=self.coeff.copy(),
            occ=self.occ.copy(),
            energy=self.energy.copy(),
            idx=self.idx.copy())

    @classmethod
    def read_orbfile(cls, path):
        coeff, occ, energy, idx = _read_orbfile(path)
        orbs = [len(x) for x in coeff]
        return cls(orbs, coeff, occ, energy, idx)

    def write_orbfile(self, path):
        with open(path, 'w') as f:
            print(self._get_header(), file=f)
            for line in self._write_MO_coeff():
                print(line, file=f)
            for line in self._write_MO_occ():
                print(line, file=f)
            for line in self._write_MO_E():
                print(line, file=f)
            for line in self._write_CAS_idx():
                print(line, file=f)

    def reindex(self, new_idx, inplace=False):
        if inplace:
            self.coeff = [
                coeff[:, idx] for idx, coeff in zip(new_idx, self.coeff)]
            self.energy = [
                energy[idx] for idx, energy in zip(new_idx, self.energy)]
            self.idx = [kind[idx] for idx, kind in zip(new_idx, self.idx)]
            self.occ = [occ[idx] for idx, occ in zip(new_idx, self.occ)]
            self.orbs = [len(idx) for idx in new_idx]
        else:
            new = self.copy()
            new.reindex(new_idx, inplace=True)
            return new

    def get_index(self):
        return [range(len(energy)) for energy in self.energy]

    def _get_header(self):
        return f"""#INPORB 2.2
#INFO
* PyOrb written
{0:8d}{len(self.orbs):8d}{0:8d}
{''.join(f"{x:8d}" for x in self.orbs)}
{''.join(f"{x:8d}" for x in self.orbs)}
* PyOrb written
#ORB"""

    def _write_MO_coeff(self):
        def res_out(v, cols, formatter):
            return '\n'.join(_reshape_output(v, cols, formatter))

        lines = []
        for irrep, n_orbs in enumerate(self.orbs):
            for orb in range(n_orbs):
                lines.append(f'* ORBITAL{irrep + 1:5d}{orb + 1:5d}')
                lines.append(res_out(self.coeff[irrep][:, orb], 5, '22.14E'))
        return lines

    def _write_MO_occ(self):
        def res_out(v, cols, formatter):
            return '\n'.join(_reshape_output(v, cols, formatter))

        lines = ['#OCC', '* OCCUPATION NUMBERS']
        for irrep, n_orbs in enumerate(self.orbs):
            lines.append(res_out(self.occ[irrep][:], 5, '22.14E'))
        lines.extend(['#OCHR', '* OCCUPATION NUMBERS (HUMAN-READABLE)'])
        for irrep, n_orbs in enumerate(self.orbs):
            lines.append(res_out(self.occ[irrep][:], 10, '8.4f'))
        return lines

    def _write_MO_E(self):
        def res_out(v, cols, formatter):
            return '\n'.join(_reshape_output(v, cols, formatter))

        lines = ['#ONE', '* ONE ELECTRON ENERGIES']
        for irrep, n_orbs in enumerate(self.orbs):
            lines.append(res_out(self.energy[irrep][:], 10, '12.4E'))
        return lines

    def _write_CAS_idx(self):
        def res_out(v):
            lines = _reshape_output(v, 10, '1')
            return '\n'.join(f'{i % 10} {a}'for i, a in enumerate(lines))

        lines = ['#INDEX']
        for irrep, n_orbs in enumerate(self.orbs):
            lines.append('* 1234567890')
            lines.append(res_out(self.idx[irrep]))
        return lines


def _forward(f, n):
    for _ in range(n):
        line = f.readline()
    return line


def _read_coeffs(f, orbs, cols):
    MO_coeff = [np.empty((n_orbs, n_orbs)) for n_orbs in orbs]
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        for orb in range(n_orbs):
            values = []
            for _ in range(rows):
                line = f.readline()
                values.extend([float(x) for x in line.split()])
            MO_coeff[irrep][:, orb] = values
            next(f)
    return MO_coeff


def _read_orb_property(f, orbs, cols):
    MO_v = [np.empty(n_orbs) for n_orbs in orbs]
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        values = []
        for _ in range(rows):
            line = f.readline()
            values.extend([float(x) for x in line.split()])
        MO_v[irrep][:] = values
    return MO_v


def _read_CAS_idx(f, orbs, cols):
    idx = []
    for n_orbs in orbs:
        rows = (n_orbs + cols - 1) // cols
        values = []
        next(f)
        for _ in range(rows):
            line = f.readline()
            values.extend(line.strip()[2:])
        idx.append(np.array(values))
    return idx


def _read_orbfile(path):
    with open(path, 'r') as f:
        _forward(f, 5)
        line = f.readline()
        orbs = [int(irrep) for irrep in line.split()]

        while line:
            line = f.readline()
            if 'ORBITAL' in line:
                MO_coeff = _read_coeffs(f, orbs, 5)
            if 'OCCUPATION NUMBERS' in line and 'HUMAN-READABLE' not in line:
                MO_occ = _read_orb_property(f, orbs, 5)
            if 'ONE ELECTRON ENERGIES' in line:
                MO_E = _read_orb_property(f, orbs, 10)
            if 'INDEX' in line:
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
