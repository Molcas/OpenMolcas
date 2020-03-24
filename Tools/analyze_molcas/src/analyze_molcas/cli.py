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
# Copyright (C) 2020, Oskar Weser                                      *
#***********************************************************************

"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later,
  but that will cause
  problems: the code will get executed twice:

  - When you run `python -manalyze_molcas` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``analyze_molcas.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``analyze_molcas.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
from pathlib import Path
import click
from analyze_molcas import SpatialOrbs


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('input_path', type=click.Path(exists=True))
@click.option('--out', '-o', default=None, type=click.Path(exists=False),
              help='Output file name. '
                   'If omitted, the suffix UhfOrb is appended.')
def spat_to_spin(input_path, out):
    """Transform spatial orbital files to spin orbitals."""
    input_path = Path(input_path)
    if out is None:
        out = input_path.parent / f'{input_path.stem}.UhfOrb'
    SpatialOrbs.read_orbfile(input_path).to_SpinOrbs().write_orbfile(out)


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('reference_path', type=click.File())
@click.argument('target_path', type=click.File())
def assimilate(reference_path, target_path):
    """Orthogonally transform orbitals in target active space
    as close as possible to reference active space."""
    reference = SpatialOrbs.read_orbfile_stream(reference_path).canonicalize()
    target = SpatialOrbs.read_orbfile_stream(target_path).canonicalize()

    for line in (reference
                 .assimilate(target, [target.get_idx_active_space()])
                 .write_orbfile()):
        print(line)


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('orb_files', type=click.File(), nargs=-1)
def combine_orbs(orb_files):
    "Combine several RasOrbs files into one. The symmetry group must match."
    import operator
    from functools import reduce

    orbs = (SpatialOrbs.read_orbfile_stream(path) for path in orb_files)

    for line in reduce(operator.add, orbs).write_orbfile():
        print(line)
