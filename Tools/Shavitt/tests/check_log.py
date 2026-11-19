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
# Copyright (C) 2026, Roland Lindh                                     *
#***********************************************************************
"""Check the CSF decoding of shavittlib against the CI printout of RASSCF.

usage: check_log.py rasscf.log rasscf.h5

Every CSF printed by RASSCF (PRWF) is decoded from its index with the graph
rebuilt from the HDF5 file; the occupation string and the coefficient
must agree. Exits with status 1 on any mismatch.
"""

import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from shavittlib.graph import ShavittGraph  # noqa: E402

ROOT_RE = re.compile(r'printout of CI-coefficients larger than\s+\S+\s+for root\s+(\d+)')
ROW_RE = re.compile(r'^\s+(\d+)\s+([0ud2 ]+?)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s*$')


def parse_log(fname):
    """{root (1-based): [(csf index 1-based, occupation, coefficient)]}"""
    out = {}
    root = None
    with open(fname) as f:
        for line in f:
            m = ROOT_RE.search(line)
            if m:
                root = int(m.group(1))
                out.setdefault(root, [])
                continue
            if root is None:
                continue
            if '****' in line:
                root = None
                continue
            m = ROW_RE.match(line)
            if m:
                out[root].append((int(m.group(1)), m.group(2).replace(' ', ''), float(m.group(3))))
    return out


def main(log, h5):
    g = ShavittGraph.from_h5(h5)
    printed = parse_log(log)
    if not printed:
        print('no CI printout found in', log)
        return 1
    nbad = ntot = 0
    for root, rows in printed.items():
        for idx, occ, coef in rows:
            ntot += 1
            steps = g.step_vector(idx - 1)
            mine = g.csf_string(steps)  # RASSCF prints the levels (RAS spaces, then irreps)
            c = g.ci_vectors[root - 1, idx - 1]
            ok = (mine == occ) and abs(c - coef) < 1.0e-4
            g.walk(steps)  # must be a valid walk
            if not ok:
                nbad += 1
                print(f'MISMATCH root {root} CSF {idx}: printed {occ} {coef:+.5f}, decoded {mine} {c:+.5f}')
    if g.count_walks()[g.lsym] != g.ncsf[g.lsym]:
        nbad += 1
        print('MISMATCH walk count', g.count_walks(), g.ncsf)
    if g.ci_vectors is not None and g.ci_vectors.shape[1] != g.ncsf[g.lsym]:
        nbad += 1
        print('MISMATCH CI vector length', g.ci_vectors.shape, g.ncsf)
    print(f'{h5}: {ntot - nbad}/{ntot} printed CSFs reproduced, NCSF {np.asarray(g.ncsf).tolist()}')
    return 1 if nbad else 0


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
