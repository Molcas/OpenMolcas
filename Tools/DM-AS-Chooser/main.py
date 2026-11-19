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

import argparse
import json
from typing import List, Optional

from src.workflow_orchestrator import PROTOCOL_SPECS, run_workflow
from src.utils import CONSTANTS


def get_parser() -> argparse.ArgumentParser:
    """
    The only rule for file naming is that the method must be included. Other information is sourced from the file.

    Excluding ground/excited states for specific protocols is managed by the software.

    --methods : which protocols to run.
                Expected format: "--methods method1 method2 ...", where method is a part of PROTOCOL_SPECS.
    --data-dir : the directory to search for candidate and reference files.
    --output-dir : the directory to output calculation matrices and final active space selections.
    -m/--mr_keywords : keyword(s) identifying candidate files by name.
                    Expected format: "-m casci casscf" finds every *.log file in --data-dir whose name
                    contains "casci" or "casscf".
    -r/--ref_keywords : keyword(s) identifying reference sources by name.
                        Expected format: "-r pbe blyp" builds two reference sources (PBE and BLYP)
                        from whichever files in --data-dir contain each keyword.
    --extract-excdata : run excitation-energy extraction ('CASSCF', 'CASPT2', 'CAS-PDFT' only).
    --generate-input : generate OpenMolcas .input files for the selected active spaces.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--methods', nargs='+', required=True, choices=list(PROTOCOL_SPECS),
        help='Which protocol(s) to run, e.g. --methods gdm-as edm-as')
    parser.add_argument(
        '--data-dir', type=str, required=True,
        help='Directory to search for candidate and reference files.')
    parser.add_argument(
        '--output-dir', type=str, required=False,
        help='Directory to output calculation matrices and final active space selections.')
    parser.add_argument(
        '-m', '--mr-keywords', type=str, required=True, nargs='+',
        metavar='KEYWORD',
        help='Keyword(s) identifying candidate files by name, e.g. '
             '-m casci casscf finds every *.log file in --data-dir '
             'whose name contains "casci" or "casscf".')
    parser.add_argument(
        '-r', '--ref-keywords', dest='ref_keywords', required=True, nargs='+',
        metavar='KEYWORD',
        help='Keyword(s) identifying reference sources by name, e.g. '
             '-r b3lyp camb3lyp builds two reference sources (B3LYP '
             'and CAMB3LYP) from whichever files in --data-dir '
             'contain each keyword.')
    parser.add_argument('--debye-per-au', type=float, default=None,
                        help=f"Override DEBYE_PER_AU (default: {CONSTANTS['DEBYE_PER_AU']})")
    parser.add_argument('--hartree-to-ev', type=float, default=None,
                        help=f"Override HARTREE_TO_EV (default: {CONSTANTS['HARTREE_TO_EV']})")
    parser.add_argument(
        '--extract-excdata', action='store_true',
        help='Run excitation-energy extraction.')

    file_gen_group = parser.add_argument_group('--generate-input options')
    file_gen_group.add_argument(
        '--generate-comparison-input', action='store_true',
        help='After running the selected protocols, generate an OpenMolcas timing-comparison input file '
             '(CASCI vs. CASSCF) for each unique winning active space.')
    file_gen_group.add_argument('--coord', default='starting.xyz')
    file_gen_group.add_argument('--basis', default='ANO-RCC-VTZP')
    file_gen_group.add_argument('--integrals', default='$CurrDir/2e-integrals/2e-integrals.OrdInt')
    file_gen_group.add_argument('--fileorb', default='$CurrDir/uhf/uhf.UnaOrb')
    file_gen_group.add_argument('--nroots', type=int, default=6)
    file_gen_group.add_argument('--charge', type=int, default=0)
    file_gen_group.add_argument('--spin', type=int, default=1)
    file_gen_group.add_argument('--cimx', type=int, default=200)
    file_gen_group.add_argument('--caspt2-maxiter', type=int, default=10000)
    file_gen_group.add_argument('--ksdft', default='T:PBE')

    return parser

def main(args: Optional[List[str]] = None):
    parser = get_parser()
    opts = parser.parse_args(args)

    if opts.debye_per_au is not None:
        CONSTANTS['DEBYE_PER_AU'] = opts.debye_per_au
    if opts.hartree_to_ev is not None:
        CONSTANTS['HARTREE_TO_EV'] = opts.hartree_to_ev

    results = run_workflow(opts, parser)

    print(json.dumps(results))


if __name__ == '__main__':
    main()