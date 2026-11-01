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
Generates OpenMolcas timing-comparison input files (CASCI vs. CASSCF) for specific active spaces.
"""


def build_comparison_input(nael: int, norb: int, args) -> str:
    """
    Builds the OpenMolcas timing-comparison input file content for one specific active space.

    :param nael: number of active electrons.
    :param norb: number of active orbitals.
    :param args: an object with the same attributes as the standalone script's parsed CLI arguments
        (coord, basis, integrals, fileorb, nroots, charge, spin, cimx, caspt2_maxiter, ksdft).

    :returns: the input file content as a string.
    """
    return f"""&GATEWAY
    Title = CASCI vs. CASSCF Timing
    Coord = {args.coord}
    Basis = {args.basis}
    Group = NoSym

&SEWARD
    OneOnly
    grid input
    grid=ultrafine
    end of grid input

>>> SHELL ln -s {args.integrals} $Project.OrdInt

&RASSCF
    FileOrb = {args.fileorb}
    nActEl = {nael}
    Ras2 = {norb}
    Charge = {args.charge}
    Spin = {args.spin}
    CIRoot = {args.nroots} {args.nroots} 1
    CIMX = {args.cimx}

&GRID_IT
    FileOrb = $Project.RasOrb
    Name = RasOrb
    All

>>> FOREACH A IN (1 .. {args.nroots})
  &GRID_IT
          FileOrb = $Project.RasOrb.$A
          Name = RasOrb.$A
          All
>>> END FOREACH

&CASPT2
    MaxIter = {args.caspt2_maxiter}
    Density
    Prop

>>> FOREACH B IN (1 .. {args.nroots})
  &GRID_IT
          FileOrb = $Project.Pt2Orb.$B
          Name = Pt2Orb.$B
          All
>>> END FOREACH

&MCPDFT
    KSDFT = {args.ksdft}

&RASSI
"""


def unique_winners(method_results: dict) -> list:
    """
    Extracts every unique winning active-space label across every protocol and reference source for ONE
    keyword group -- so the same active space winning multiple times within that group doesn't produce
    duplicate timing input files.

    :param method_results: {method: {source: label_or_None}}, one keyword's slice of run()'s results.

    :returns: sorted list of unique winning labels (e.g. ['06-06', '08-08']).
    """
    winners = set()
    for source_results in method_results.values():
        for label in source_results.values():
            if label is not None:
                winners.add(label)
    return sorted(winners)


def generate_comparison_inputs(method_results: dict, opts, output_dir) -> list:
    """
    Generates one timing-comparison input file per unique winning active space, for one keyword group.

    :param method_results: {method: {source: label_or_None}}, one keyword's slice of run()'s results.
    :param opts: parsed CLI options, providing the timing-input-specific arguments.
    :param output_dir: that keyword's own output directory.

    :returns: list of paths written.
    """
    written = []
    for label in unique_winners(method_results):
        nael, norb = (int(x) for x in label.split('-'))
        content = build_comparison_input(nael, norb, opts)
        path = output_dir / f'comparison-{label}.input'
        with open(path, 'w') as f:
            f.write(content)
        written.append(path)
    return written