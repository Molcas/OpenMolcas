from io import StringIO
import chemcoord as cc
from os.path import join, relpath, normpath, dirname, basename


def create_module_input(module_name, **kwargs):
    indentstep = 2
    indentation = indentstep * ' '

    out = f'&{module_name.upper()}\n'
    for keyword, value in kwargs.items():
        out += f'{1 * indentation}{keyword}\n'
        out += f'{2 * indentation}{value}\n'
    out += 'End of Input'
    return out

def emil_cmd(cmd, *args):
    return f'> {cmd}' + ' '.join(args)


def rel_currdir_path(path):
    return normpath(join(
        '$CurrDir',
        relpath(dirname(path), dirname(el_calc_input)),
        basename(path)))

def _generate_calc_input(
        molecule, hamiltonian, basis, el_calc_input, charge=0,
        forces=True, title=' ', start_orb=None, sym_group=None,
        multiplicity=1):
    """Generate a molcas input file.

    Args:
        molecule (chemcoord.Cartesian or chemcoord.Zmat or str):
            If it is a string, it has to be a valid xyz-file.
        hamiltonian (str): {hamiltonian}
        basis (str): {basis}
        charge (int): {charge}
        forces (bool): {forces}
        title (str): {title}
        multiplicity (int): {multiplicity}
        wfn_symmetry (int): {wfn_symmetry}


    Returns:
        str : molcas input.
    """
    if isinstance(molecule, str):
        molecule = molecule.read_xyz(StringIO(molecule))
    elif isinstance(molecule, cc.Zmat):
        molecule = molecule.get_cartesian()

    inp = create_module_input(
          'GATEWAY',
          title=title, coord=molecule.to_xyz(sorted_index=False),
          basis=basis, group='full' if sym_group is None else sym_group)

    inp += create_module_input('SEWARD')

    if start_orb is not None:
        start_orb_path = normpath(join(
            '$CurrDir',
            relpath(dirname(start_orb), dirname(el_calc_input)),
            basename(start_orb)))

        inp += emil_cmd('COPY', start_orb_path, 'INPORB')

    get_output = """\
&GATEWAY
 Coord
 {geometry}
 Title = {title}
 Basis = {basis}
 {sym_group}

&SEWARD

{hamiltonian_str}

{forces}
""".format

    out = get_output(
        title=title, basis=basis, geometry=molecule.to_xyz(sort_index=False),
        sym_group='' if sym_group is None else 'group = {}'.format(sym_group),
        hamiltonian_str=_get_hamiltonian_str(hamiltonian, charge, multiplicity, start_orb, el_calc_input),
        forces='&ALASKA' if forces else '')
    return out

