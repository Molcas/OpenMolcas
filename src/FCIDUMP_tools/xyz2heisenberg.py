"""
Convert an xyz file to a Heisenberg FCIDUMP or a python array of J values.
The Heisenberg Js are inversely proportional to the distance between the atoms.

Copyright (c) 2025 Maru Song
"""

import numpy as np
from itertools import combinations
from typing import Optional

def xyz_parser(xyzfile: str) -> tuple[list[str], np.ndarray]:
    """
    Parse an XYZ file and extract atomic symbols and coordinates.
    
    Args:
        xyzfile (str): Path to the XYZ file to parse.
        
    Returns:
        tuple[list[str], np.ndarray]: A tuple containing:
            - List of atomic symbols (e.g., ['H', 'O', 'H'])
            - NumPy array of 3D coordinates with shape (natom, 3)
            
    Raises:
        ValueError: If the coordinates are not in 3D format.
    """
    with open(xyzfile, 'r') as f:
        lines = f.readlines()
    natom = int(lines[0])
    atoms: list[str] = []
    coords_list: list[list[float]] = []
    for i in range(2, natom+2):
        atoms.append(lines[i].split()[0])
        coords_list.append(list(map(float, lines[i].split()[1:])))
    coords: np.ndarray = np.array(coords_list)
    if coords.shape[1] != 3:
        raise ValueError("The coordinates should be in 3D.")
    return atoms, coords

def measure_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    return np.linalg.norm(coord1 - coord2)

def gen_combinations(natom: int) -> list[tuple[int, int]]:
    return list(combinations(range(1, natom + 1), 2))

def gen_J_matrix(xyzfile: str, distance_range: Optional[tuple[float, float]] = None) -> np.ndarray:
    """
    Generate a symmetric matrix of J coupling constants from atomic coordinates.
    
    The J values are calculated as the inverse of the distance between atoms,
    representing antiferromagnetic interactions in the Heisenberg model.
    
    Args:
        xyzfile (str): Path to the XYZ file containing atomic coordinates.
        distance_range (tuple[float, float], optional): Tuple of (min_dist, max_dist)
                                                       to filter interactions by distance.
                                                       If None, all pairs are included.
                                                       Defaults to None.
    
    Returns:
        np.ndarray: Symmetric matrix of J coupling constants with shape (natom, natom).
                   J_matrix[i][j] = 1.0 / distance(i, j) for interacting pairs,
                   0.0 otherwise.
                   
    Note:
        The resulting matrix has J[i][j] = J[j][i] and J[i][i] = 0.
    """
    atoms, coords = xyz_parser(xyzfile)
    natom = len(atoms)
    combinations = gen_combinations(natom)
    J_matrix = np.zeros((natom, natom))
    for i, j in combinations:
        distance = measure_distance(coords[i - 1], coords[j - 1])

        if distance_range is not None:
            min_dist, max_dist = distance_range
            if distance < min_dist or distance > max_dist:
                continue

        J = 1.0 / distance # antiferromagnetic interaction
        J_matrix[i - 1][j - 1] = J
        J_matrix[j - 1][i - 1] = J
    return J_matrix

def gen_heisenberg_fcidump(xyzfile: str, fcidumpname: str, nel: int, norb: int, ms: int,
                           distance_range: Optional[tuple[float, float]] = None) -> None:
    """
    Generate a Heisenberg model FCIDUMP file from atomic coordinates.
    
    This function creates an FCIDUMP file containing two-electron integrals
    for the Heisenberg model, where J values are inversely proportional to the
    interaction distances.
    The interactions are antiferromagnetic (negative J values).
    
    Args:
        xyzfile (str): Path to the XYZ file containing atomic coordinates.
        fcidumpname (str): Output filename for the FCIDUMP file.
        nel (int): Number of electrons in the system.
        norb (int): Number of orbitals in the system.
        ms (int): Spin multiplicity (2*S) of the system.
        distance_range (tuple[float, float], optional): Tuple of (min_dist, max_dist)
                                                       to filter interactions by distance.
                                                       If None, all pairs are included.
                                                       Defaults to None.
    
    Returns:
        None: The function writes the FCIDUMP file to disk.
        
    Note:
        The negative sign indicates antiferromagnetic coupling.
    """
    atoms, coords = xyz_parser(xyzfile)
    natom = len(atoms)
    combinations = gen_combinations(natom)
    integral_list = []
    for i, j in combinations:
        distance = measure_distance(coords[i - 1], coords[j - 1])

        if distance_range is not None:
            min_dist, max_dist = distance_range
            if distance < min_dist or distance > max_dist:
                continue

        J = -1.0 / distance # antiferromagnetic interaction
        diag = J / 2.0
        integral_list.append([diag, j, j, i, i])
        integral_list.append([J, j, i, i, j])
    dump_integrals(fcidumpname, nel, norb, ms, integral_list)

#------------------------------------------------------------------------------#
"""
Collection of functions for writing an FCIDUMP file.
This feature is redundant with a feature in IntegralClass.
They have to be combined.
"""

def dump_integrals(fcidumpname, nel, norb, ms, integral_list):
    "integral list has to be [[x1 i1 j1 k1 l1], [x2 i2 j2 k2 l2], ...]"

    header = gen_header(nel, norb, ms)

    with open(fcidumpname, 'w') as f:
        f.write(header)
        for ints in integral_list:
            f.write(gen_int_manual(*ints))

def gen_header(nel, norb, ms):
    """
    Generates the header for a manual fcidump file
    Currently, it is assumed that there is no point group symmetry
    """

    header = f""" &FCI NORB=  {norb}, NELEC=  {nel}, MS2=  {ms},
  ORBSYM= {', '.join(['1'] * norb)}
  ISYM=0
 &END\n"""

    return header
    
def gen_int_manual(x, i, j, k, l, digits=12):
    integral = f"""    {x:.{digits}f}    {i}    {j}    {k}    {l}\n"""
    return integral

#------------------------------------------------------------------------------#
