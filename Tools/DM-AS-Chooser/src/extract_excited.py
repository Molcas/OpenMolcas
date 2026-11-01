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
# Copyright (C) 2026, Benjamin W. Kaufold                              *
#               2026, Meredith E. Wanta                                *
#               2026, Sijia S. Dong                                    *
#***********************************************************************

"""
extract_excdata:
    Reads excited-state energies (RASSI's spin-free "SF State" excitation energies, CASPT2 root energies,
    and MC-PDFT total energies) from a single OpenMolcas .log file, and writes a TSV of excitation energies
    (relative to the ground state), one row per excited state:
        CAS   PT2   PDFT   PT2(sorted)   PDFT(sorted)

extract_excerr:
    Compares results from excdata against user-provided high quality references by absolute difference.
"""

import re
from pathlib import Path
from typing import List, Optional

import numpy as np

from src.utils import logger, CONSTANTS


def extract_excdata(log_path: str, out_path: Path = None) -> str:
    """
    Reads excitation-energy data from one OpenMolcas log file and writes the combined TSV.

    :returns out_path : the path to the .exc.tsv output file.
    """
    log_path = Path(log_path)
    if out_path is None:
        out_path = log_path.with_suffix('').with_suffix('.exc.tsv')
    else:
        out_path = Path(out_path)

    # Get excited-state energy values
    casexc = np.array(
        _grep_block_field(log_path, r'SF State', after=7, tail=5,
                         field_1indexed=4),
        dtype=float)
    pt2nrg = np.array(
        _grep_field(log_path, r'CASPT2 Root', field_1indexed=7), dtype=float)
    pftnrg = np.array(
        _grep_field(log_path, r'Total MC-PDFT', field_1indexed=8), dtype=float)

    pt2exc = CONSTANTS['HARTREE_TO_EV'] * (pt2nrg[1:] - pt2nrg[0])
    pftexc = CONSTANTS['HARTREE_TO_EV'] * (pftnrg[1:] - pftnrg[0])
    st2exc = np.sort(pt2exc)
    sftexc = np.sort(pftexc)

    combined = np.column_stack([casexc, pt2exc, pftexc, st2exc, sftexc])
    np.savetxt(out_path, combined, delimiter='\t')

    return str(out_path)

def extract_excerr(exc_path: str, ref_path: Path) -> Optional[str]:
    """
    Computes the absolute error against a user-provided ref-excs.list (high quality calculation references) for a
    given .exc.tsv file.

    :param exc_path: path to a .exc.tsv file (as written by extract_excdata.py).
    :param ref_path: path to the .list file to use for calculation references.

    :returns: the path to the written .excerr.tsv, or None if no ref-excs.list was found alongside exc_path.
    """

    exc_path = Path(exc_path)

    if not ref_path.exists():
        logger.error(f"Calculation references not found at {ref_path}")
        return None

    calc_excs = np.loadtxt(exc_path)
    ref_excs = np.loadtxt(ref_path, ndmin=2)

    exc_err = np.abs(calc_excs - ref_excs)

    err_path = exc_path.parent / exc_path.name.replace('.exc.tsv', '.excerr.tsv', 1)
    np.savetxt(err_path, exc_err)
    return str(err_path)

def _bash_field(line: str, field_1indexed: int) -> str:
    """
    Splits a line into fields separated by runs of one or more spaces, and returns the field at the
    requested position (counting from 1). A line with leading whitespace gets an empty first field,
    shifting every real field's position by one.
    """
    squeezed = re.sub(r' +', ' ', line.rstrip('\n'))
    fields = squeezed.split(' ')
    return fields[field_1indexed - 1]


def _grep_field(path: Path, pattern: str, field_1indexed: int) -> List[str]:
    """
    Extracts one field from every line in the file that contains `pattern`, in file
    order. Meant for data printed as a series of independent single-line facts.
    """
    regex = re.compile(pattern)
    values = []
    with open(path, 'r') as f:
        for line in f:
            if regex.search(line):
                values.append(_bash_field(line, field_1indexed))
    return values


def _grep_block_field(path: Path, pattern: str, after: int, tail: int,
                      field_1indexed: int) -> List[str]:
    """
    Finds the first line matching `pattern`, takes it plus the next `after` lines, keeps only the last `tail` of those,
    and extracts one field from each. Meant for a single contiguous block of data following a matching header line.
    """
    regex = re.compile(pattern)
    with open(path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if regex.search(line):
            block = lines[i:i + after + 1][-tail:]
            return [_bash_field(l, field_1indexed) for l in block]

    return []