#!/usr/bin/env python3
"""Script to automatically generate grid files from Orbfiles.
"""

import os
from os import remove
from os.path import join, basename, splitext, relpath

import subprocess
from subprocess import run
import click


def create_module_input(module_name, **kwargs):
    indentstep = 2
    indentation = indentstep * ' '

    out = f'&{module_name.upper()}\n'
    for keyword, value in kwargs.items():
        out += f'{1 * indentation}{keyword}\n'
        if value:
            out += f'{2 * indentation}{value}\n'
    out += 'End of Input\n'
    return out


def emil_cmd(cmd, *args):
    return f'> {cmd} ' + ' '.join(args) + '\n'


def get_currdir_path(path):
    return join('$CurrDir', relpath(path, os.getcwd()))


def create_input(coord, basis, inporb, density='sparse', group=None,
                 select=None, total=False):
    inp = create_module_input(
        'GATEWAY',
        coord=get_currdir_path(coord),
        basis=basis, group='full' if group is None else group)

    inp += emil_cmd('COPY', get_currdir_path(inporb), 'INPORB')

    gridit = {'NOLUSCUS': '', 'ASCII': '', density: ''}
    if select:
        gridit['select'] = select
    if total:
        gridit['total'] = ''

    inp += create_module_input('GRIDIT', **gridit)
    return inp


def start_calc(input_str, inp_path='make_grid.inp', molcas_exe='pymolcas'):

    with open(inp_path, 'w') as f:
        f.write(input_str)

    run([molcas_exe, '-f', basename(inp_path)], stdout=subprocess.PIPE,
        check=True)


def get_project_name(inporb):
    return splitext(basename(inporb))[0]


def cleanup(project):
    remove(f'{project}.inp')
    remove(f'{project}.log')
    remove(f'{project}.status')
    remove(f'{project}.err')
    remove('xmldump')


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--coord', '-c', required=True, type=str,
              help='Filepath to xyz-file.')
@click.option('--basis', '-b', required=True, type=str,
              help='Basisset to use in Molcas nomenclature.')
@click.option('--inporb', '-i', required=True, type=str,
              help='Filepath to Molcas orbital file.')
@click.option('--density', '-d', default='sparse', type=str, show_default=True,
              help='Gridit keyword. One of {sparse, normal, dense}.')
@click.option('--group', '-g', default='full', type=str, show_default=True,
              help='Symmetry group in Molcas nomenclature.')
@click.option('--select', '-s', default=None, type=str, show_default=True,
              help=('Gridit keyword to select certain orbitals. '
                    'The default depends on the Molcas implementation. '
                    'Usually the active space is used.'))
@click.option('--total/--no-total', '-t',default=False,
              help='Total density computed from contributions of all orbitals.')
@click.option('--project', '-p', default='make_grid', type=str,
              show_default=True,
              help=('Prefix for the grid calculations. '
                    'The gridfile will be named ${project}.grid.'))
@click.option('--molcas_exe', '-m', default='pymolcas', type=str,
              show_default=True,
              help='Filepath of the Molcas driver script.')
@click.option('--clean/--no-clean', default=True,
              help='Clean temporary files after calculation.')
def main(coord, basis, inporb, density, group, select, total, project, molcas_exe, clean):
    start_calc(
        create_input(coord, basis, inporb, density, group, select, total),
        f'{project}.inp', molcas_exe)
    if clean:
        cleanup(project)



if __name__ == '__main__':
    main()
