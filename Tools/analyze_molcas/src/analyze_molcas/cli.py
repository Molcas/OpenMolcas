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
