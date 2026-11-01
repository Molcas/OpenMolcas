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
This script generates multiple OpenMolcas input files by scanning over a range of active spaces
using a template that controls other features of the calculation.

The user can supply command-line arguments to change the active spaces scanned and other calculation details.

In its current form it uses CASCI at the active space selection phase in CASCI-DM-AS,
but this can be changed by editing the templates.
"""

import argparse
from pathlib import Path


def write_file(filename, content):
    """ Write content to filename """
    with open(filename, "w") as f:
        f.write(content)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate OpenMolcas input files for a CASCI active-space scan, "
            "plus the two-electron integral and ground-state orbital inputs "
            "needed to support it."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="The full path to the output directory",
    )
    parser.add_argument(
        "--coord", default="starting.xyz",
        help="Name of Cartesian coordinates file for input system",
    )
    parser.add_argument(
        "--basis", default="ANO-RCC-VTZP",
        help="Basis set in OpenMolcas format",
    )
    parser.add_argument(
        "--nroots", type=int, default=6,
        help="How many roots to include in SA-CASCI",
    )
    parser.add_argument(
        "--charge", type=int, default=0,
        help="Charge",
    )
    parser.add_argument(
        "--spin", type=int, default=1,
        help="Spin as (2S+1)",
    )
    parser.add_argument(
        "--nael-min", type=int, default=6,
        help="Minimum number of active space electrons (n_e) to scan",
    )
    parser.add_argument(
        "--nael-max", type=int, default=14,
        help="Maximum number of active space electrons (n_e) to scan",
    )
    parser.add_argument(
        "--nael-step", type=int, default=2,
        help="Step size between successive n_e values",
    )
    parser.add_argument(
        "--norb-max", type=int, default=14,
        help="Maximum number of active space orbitals (n_o) to scan",
    )
    parser.add_argument(
        "--norb-offset", type=int, default=3,
        help="Offset added to n_e/2 to get the minimum n_o for a given n_e",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Template for generating two-electron integrals
    two_electron_input = f"""&GATEWAY
    Title = DM-AS Two-Electron Integrals
    Coord = {args.coord}
    Basis = {args.basis}
    Group = NoSym

&SEWARD
    grid input
    grid=ultrafine
    end of grid input

>>> Copy $Project.OrdInt $CurrDir/2e-integrals.OrdInt
"""
    write_file(output_dir / "2e-integrals.input", two_electron_input)


    # Template for calculating ground-state orbitals
    ground_orbs_input = f"""&GATEWAY
    Title = DM-AS Ground-State Orbitals
    Coord = {args.coord}
    Basis = {args.basis}
    Group = NoSym

&SEWARD
    OneOnly
    grid input
    grid=ultrafine
    end of grid input

>>> Shell ln -s $CurrDir/2e-integrals.OrdInt $Project.OrdInt

&SCF

&MBPT2
    Prpt

>>> Copy $Project.MP2Orb $CurrDir/ground-orbs.InpOrb
"""
    write_file(output_dir / "ground-orbs.input", ground_orbs_input)

    # Active space scan: nael = number of active electrons,
    #                    norb = number of active orbitals
    for nael in range(args.nael_min, args.nael_max + 1, args.nael_step):
        pael = f"{nael:02d}"
        morb = nael // 2 + args.norb_offset
        for norb in range(morb, args.norb_max + 1):
            porb = f"{norb:02d}"
            inam = f"casci-{pael}-{porb}.input"

            casci_input = f"""&GATEWAY
    Title = DM-AS CASCI Active Space Scan
    Coord = {args.coord}
    Basis = {args.basis}
    Group = NoSym

&SEWARD
    OneOnly
    grid input
    grid=ultrafine
    end of grid input

>>> Shell ln -s $CurrDir/2e-integrals.OrdInt $Project.OrdInt

&RASSCF
    FileOrb = $CurrDir/ground-orbs/ground-orbs.MP2Orb
    nActEl = {nael}
    Ras2 = {norb}
    CIRoot = {args.nroots} {args.nroots} 1
    Charge = {args.charge}
    Spin = {args.spin}
    CIOnly
    CIMX = 200
"""
            write_file(output_dir / inam, casci_input)


if __name__ == "__main__":
    main()