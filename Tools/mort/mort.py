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
# Copyright (C) 2021, Ignacio Fdez. Galván                             *
#***********************************************************************

import sys
import re
import argparse
import numpy as np
import h5py

version = '1.0'

################################################################################
# FUNCTIONS
################################################################################

# Copy dataset name_in of file_in into dataset name_out of file_out,
# but attributes (description) are taken from name_out.
# If size is specified, no data is copied and a dataset is created with the new size
def copy_dataset(file_in, file_out, name_in, name_out=None, size=None):
  if name_out is None:
    name_out = name_in
  no_attrs = name_in != name_out
  if size is None:
    file_in.copy(file_in[name_in], file_out, name_out, without_attrs=no_attrs)
  else:
    file_out.create_dataset(name_out, size, dtype=file_in[name_out].dtype)
  if no_attrs:
    for attr in file_in[name_out].attrs.keys():
      file_out[name_out].attrs.create(attr, file_in[name_out].attrs[attr])

# Rearrange a symmetry-blocked array into a square array
def sblock_to_square(data, nbas, triangular=False):
  if len(nbas) == 1:
    if triangular:
      mtrx = np.zeros((nbas[0], nbas[0]))
      ii, jj = np.tril_indices(nbas[0])
      mtrx[ii,jj] = data
      mtrx[jj,ii] = data
    else:
      mtrx = np.reshape(data, (nbas[0], nbas[0]))
  else:
    nb = int(np.sum(nbas))
    mtrx = np.zeros((nb, nb))
    j = 0
    for i,n in enumerate(nbas):
      ni = int(np.sum(nbas[:i]))
      if triangular:
        tmp = np.zeros((n, n))
        ii, jj = np.tril_indices(n)
        tmp[ii,jj] = data[j:j+n*(n+1)//2]
        tmp[jj,ii] = data[j:j+n*(n+1)//2]
        mtrx[np.ix_(range(ni,ni+n), range(ni,ni+n))] = tmp[:,:]
        j += n*(n+1)//2
      else:
        mtrx[np.ix_(range(ni,ni+n), range(ni,ni+n))] = np.reshape(data[j:j+n**2], (n, n))
        j += n**2
  return mtrx

# Normalize a rotation matrix,
# i.e., find the closest unitary matrix
def normalize_rot(A, test=True):
  U, s, Vh = np.linalg.svd(A)
  if test:
    if not np.allclose(s, 1.0, atol=1e-3):
      sys.exit('Not a unitary rotation')
  return U @ Vh

# Returns unique permutations (permutations with repeated elements).
# From https://stackoverflow.com/questions/6284396/permutations-with-unique-values
def permute_unique(nums):
  perms = [[]]
  for n in nums:
    new_perm = []
    for perm in perms:
      for i in range(len(perm) + 1):
        new_perm.append(perm[:i]+[n]+perm[i:])
        # handle duplication
        if i < len(perm) and perm[i] == n:
          break
    perms = new_perm
  return perms

# Returns a list of Cartesian factors (0=x, 1=y, 2=z), given l and m
# (note that m goes from -l to +l and beyond)
def factors(l, m):
  ly = int(np.floor((np.sqrt(8*(m+l)+1)-1)/2))
  lz = m+l-ly*(ly+1)//2
  lx = l-ly
  ly -= lz
  f = []
  for i in range(lx):
    f.append(0)
  for i in range(ly):
    f.append(1)
  for i in range(lz):
    f.append(2)
  return f

# Compute the rotation matrix for the functions in a Cartesian shell
# of angular momentum l, given the real-space rotation matrix R.
# (The cases l=0,1 are trivial, but we don't treat them separately)
def cartesian_shell_rotation(l, R):
  assert l >= 0
  num = (l+1)*(l+2)//2
  Rot = np.zeros((num, num))
  # Pre-compute the Cartesian factors for each m
  f = []
  for m in range(-l, num-l):
    f.append(factors(l, m))
  for m1 in range(num):
    perms = permute_unique(f[m1])
    for m2 in range(num):
      # Each element [m2,m1] is of the form sum(prod(R[i,j])),
      # where the sum extends over all distinct permutations of the m1 factors
      # and the prod is over the factors of m1 & m2
      for p in perms:
        K = 1.0
        for i,j in zip(p, f[m2]):
          K *= R[i,j]
        Rot[m2,m1] += K
  return Rot

# Compute the rotation matrix for the functions in a spherical harmonics shell
# of angular momentum l, given the real-space rotation matrix R (implicit).
# Uses the recursive algorithm of Ivanic and Ruedenberg (doi:10.1021/jp953350u, doi:10.1021/jp9833350),
# which is valid for l > 1.
# Note that this builds the inverse/transpose, so indexing may look reversed.
def spherical_shell_rotation(l):
  assert l > 1
  def rP(i,mu,m2):
    if abs(mu) > l-1:
      return 0
    # Instead of referring to the original rotation matrix R, we use RotM[1],
    # which is already reordered according to m (but indices are 0,1,2, not -1,0,1).
    # NB: The first two cases are wrong in the original (m' -> l-1)
    elif m2 == l:
      return RotM[1][2,i+1]*RotM[l-1][2*l-2,l-1+mu]-RotM[1][0,i+1]*RotM[l-1][0,l-1+mu]
    elif m2 == -l:
      return RotM[1][2,i+1]*RotM[l-1][0,l-1+mu]+RotM[1][0,i+1]*RotM[l-1][2*l-2,l-1+mu]
    else:
      return RotM[1][1,i+1]*RotM[l-1][l-1+m2,l-1+mu]
  num = 2*l+1
  Rot = np.zeros((num, num))
  for m2 in range(-l, l+1):
    if abs(m2) == l:
      den = 2*l*(2*l-1)
    else:
      den = (l+m2)*(l-m2)
    for m1 in range(-l, l+1):
      u = np.sqrt((l+m1)*(l-m1)/den)
      rU = rP(0,m1,m2)
      if m1 == 0:
        v = -0.5*np.sqrt(2*(l-1)*l/den)
        rV = rP(1,1,m2)+rP(-1,-1,m2)
        w = 0
        rW = 0
      else:
        s = 1 if m1 > 0 else -1
        absm = abs(m1)
        v = 0.5*np.sqrt((l+absm-1)*(l+absm)/den)
        # NB: The case m < 0 is wrong in the correction (sqrt(1-delta) -> sqrt(1+delta))
        if absm == 1:
          rV = np.sqrt(2)*rP(s,absm-1,m2)
        else:
          rV = rP(1,m1-s,m2)-s*rP(-1,-m1+s,m2)
        w = -0.5*np.sqrt((l-absm-1)*(l-absm)/den)
        rW = rP(1,m1+s,m2)+s*rP(-1,-m1-s,m2)
      Rot[l+m2,l+m1] = u*rU+v*rV+w*rW
  return Rot

################################################################################

# Command-line arguments
parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=24,width=120))
parser.version = f'MORT: molecular orbital rotation and translation (version {version})'
parser.add_argument('-v', '--version', help='print version number', action='version')
parser.add_argument('-r', '--rotate', help='rotation matrix to apply (9 real values)', type=float, nargs=9, metavar=('Rxx','Rxy','Rxz','Ryx','Ryy','Ryz','Rzx','Rzy','Rzz'))
parser.add_argument('-t', '--translate', help='translation vector to apply after rotation (3 real values, in bohr)', type=float, nargs=3, metavar=('Tx','Ty','Tz'))
parser.add_argument('-x', '--exchange', help='reorder atoms (comma-separated list of integers giving order in the output file)', metavar='n1,n2,...')
parser.add_argument('-d', '--desymmetrize', help='desymmetrize data', action='store_true')
parser.add_argument('-a', '--all', help='copy all unhandled datasets and attributes to output file (some items may be wrong!)', action='store_true')
parser.add_argument('infile', help='input file (HDF5 format)', metavar='input_file')
parser.add_argument('outfile', help='output file (HDF5 format, will be overwritten)', metavar='output_file')
if len(sys.argv) < 2 and sys.stdin.isatty():
  parser.print_help(sys.stderr)
  sys.exit(1)
else:
  args = vars(parser.parse_args())

# Read input file

with h5py.File(args['infile'], 'r') as f:
  nsym = f.attrs['NSYM']
  if nsym != 1:
    if not args['desymmetrize']:
      sys.exit('Error: Symmetry is only supported with --desymmetrize')
  try:
    if nsym == 1:
      b_id_ref = f['BASIS_FUNCTION_IDS'][:]
    else:
      b_id_ref = f['DESYM_BASIS_FUNCTION_IDS'][:]
    nb = len(b_id_ref)
  except KeyError:
    sys.exit('Error: Input file does not have basis set information')
  try:
    if nsym == 1:
      nat = f.attrs['NATOMS_UNIQUE']
    else:
      nat = f.attrs['NATOMS_ALL']
  except KeyError:
    sys.exit('Error: Input file does not have number of atoms')

# Compute transformation matrices

if args['rotate']:
  do_rot = True
  R = np.reshape(args['rotate'], (3,3))
  # Normalize to clean up numerical errors
  R = normalize_rot(R)

  # Maximum l in spherical shells
  # (Cartesian shells have negative l and don't need recursion)
  max_l = np.max(b_id_ref[:,2])

  # Group the basis functions by shells
  # (same center, l, and radial part)
  l_idx = {}
  b_num = []
  for i,b in enumerate(b_id_ref):
    l = b[2] 
    for j,bb in enumerate(b_id_ref[:i,:]):
      if tuple(b[0:3]) == tuple(bb[0:3]):
        num = b_num[j] 
        l_idx[l][num].append(i)
        b_num.append(num)
        break
    else:
      if l in l_idx:
        l_idx[l].append([i])
      else:
        l_idx[l] = [[i]]
      b_num.append(len(l_idx[l])-1)
    
  # Make sure each shell is ordered by m
  for l in l_idx:
    for i,shell in enumerate(l_idx[l]):
      l_idx[l][i] = sorted(shell, key=lambda x: b_id_ref[x,3])

  # Compute the rotation matrices for the different shells
  RotM = {}
  # First the Cartesian shells
  for l in l_idx.keys():
    # s and p shells (l=0,1) are always considered Cartesian
    if l < 2:
      RotM[l] = cartesian_shell_rotation(abs(l), R)
  # Reorder the p rotation to match the m ordering (1=x, -1=y, 0=z)
  RotM[1] = RotM[1][np.ix_([1,2,0],[1,2,0])]
  # Then the spherical harmonics shells
  for l in range(2,max_l+1):
    RotM[l] = spherical_shell_rotation(l)
  
  # Build the big (all functions), sparse transformation matrix and its inverse,
  # which is needed because the Cartesian rotations are not unitary
  U = np.zeros((nb, nb))
  V = np.zeros((nb, nb))
  for l in l_idx:
    for shell in l_idx[l]:
      U[np.ix_(shell,shell)] = RotM[l]
      if l < 0:
        V[np.ix_(shell,shell)] = np.linalg.inv(RotM[l])
      else:
        V[np.ix_(shell,shell)] = RotM[l].T
else:
  do_rot = False
  RotM = {}
  R = np.eye(3)
  U = np.eye(nb)
  V = np.eye(nb)

if args['translate']:
  do_trans = True
  T = np.array(args['translate'])
else:
  do_trans = False
  T = np.zeros(3)

if args['exchange']:
  do_exch = True
  # List that specifies how the centers are reordered, its inverse, and the corresponding matrix
  reindex = [int(i)-1 for i in args['exchange'].split(',')]
  assert len(reindex) == nat, 'Error: length of --exchange list does not match number of atoms'
  reindex_inv = []
  for i in range(len(reindex)):
    reindex_inv.append(reindex.index(i))
  reorder = np.eye(nat)[reindex]
else:
  do_exch = False
  reindex = list(range(nat))
  reindex_inv = reindex
  reorder = np.eye(nat)

if args['desymmetrize']:
  do_desym = nsym != 1
else:
  do_desym = False

if do_desym or do_exch:
  # List that specifies how basis functions are reordered and its matrix
  bf_reindex = []
  # Create an auxiliary array to use lexsort on:
  aux_bfid = b_id_ref.copy()
  for b in aux_bfid:
    if do_exch:
      b[0] = reindex_inv[b[0]-1]+1
    b[2] = abs(b[2])
    # p shells are ordered x,y,z
    if b[2] == 1:
      if b[3] == 1:
        b[3] = -1
      elif b[3] == -1:
        b[3] = 0
      elif b[3] == 0:
        b[3] = 1
  bf_reindex = np.lexsort((aux_bfid[:,1], aux_bfid[:,3], aux_bfid[:,2], aux_bfid[:,0]))
  bf_reorder = np.eye(nb)[bf_reindex]
  # U and V are easily reordered
  # (this affects the matrices even if there was no rotation, so enforce rotation)
  do_rot = True
  U = U @ bf_reorder.T
  V = bf_reorder @ V

do_smth = do_rot or do_trans or do_exch or do_desym

# Process the input data and write the output file

# Known attributes and datasets that are generally safe to copy (need no transformation)
safe_to_copy = [
  # attributes
  'A2LEV',
  'CI_TYPE',
  'L2ACT',
  'MOLCAS_MODULE',
  'MOLCAS_VERSION',
  'NACTEL',
  'NCONF',
  'NELEC3',
  'NHOLE1',
  'NPRIM',
  'NROOTS',
  'NSTATE',
  'NSTATES',
  'ORBITAL_TYPE',
  'POTNUC',
  'RASSCF_ITERATIONS',
  'ROOT2STATE',
  'SPINMULT',
  'STATE_SPINMULT',
  # datasets
  'CI_VECTORS',
  'ENERGY',
  'H_EFF',
  'ORIGINAL_OVERLAPS',
  'ROOT_ENERGIES',
  'SFS_COEFFICIENTS',
  'SFS_ENERGIES',
  'STATE_WEIGHT',
  'STATE_PT2_ENERGIES',
  'STATE_REFWF_ENERGIES',
  'STATE_ROOTID',
]

# Attributes and datasets that cannot be transformed with desymmetrization 
# (e.g. state order is not known)
not_with_desym = [
  # attributes
  'NCONF',
  # datasets
  'CI_VECTORS',
  'STATE_ROOTID',
]

warning_item = []
warning_renum = False
with h5py.File(args['infile'], 'r') as fi, h5py.File(args['outfile'], 'w') as fo:

  fo.attrs.create('MORT_version', np.array(f'MORT v{version}', dtype='S'))

  module = fi.attrs['MOLCAS_MODULE'].decode()

  # Pick up the data needed for desymmetrization, which includes reordering MOs
  if do_desym:
    desym = np.reshape(fi['DESYM_MATRIX'], (nb, nb))
    nbas = fi.attrs['NBAS']
    mo_reindex = np.lexsort((fi['MO_ENERGIES'], -fi['MO_OCCUPATIONS'][:]))
    try:
      # In case of UHF
      mo_alpha_reindex = np.lexsort((fi['MO_ALPHA_ENERGIES'], -fi['MO_ALPHA_OCCUPATIONS'][:]))
      mo_beta_reindex = np.lexsort((fi['MO_BETA_ENERGIES'], -fi['MO_BETA_OCCUPATIONS'][:]))
    except KeyError:
      pass
    try:
      # In case of [CR]ASSCF, some data is on active orbitals only
      tpidx = fi['MO_TYPEINDICES'][:]
    except KeyError:
      pass
    else:
      actidx = np.isin(tpidx, [b'1', b'2', b'3']).nonzero()[0]
      nact = actidx.shape[0]
      act_reindex = np.empty(nact)
      for i in range(nact):
        act_reindex[i] = np.where(mo_reindex == actidx[i])[0][0]
      act_reindex = np.argsort(act_reindex)
      act_reorder = np.eye(nact)[act_reindex]
      # In CASPT2 density matrices are stored for all (non-frozen, non-deleted) orbitals
      if module == 'CASPT2':
        pt2idx = np.logical_not(np.isin(tpidx, [b'F', b'D'])).nonzero()[0]
        bins = [0]
        for i in nbas:
          bins.append(bins[-1]+i)
        pt2nbas = np.histogram(pt2idx, bins=bins)[0]
        # TODO: reindex/reorder MOs for CASPT2 (even if nothing happens by default)
    V = V @ desym.T
  else:
    nbas = [nb]

  for name in fi.attrs:

    if name in safe_to_copy:
      if do_desym and name in not_with_desym:
        continue
      fo.attrs[name] = fi.attrs[name]

    elif name == 'LSYM':
      if do_desym:
        fo.attrs[name] = 1
      else:
        fo.attrs[name] = fi.attrs[name]

    elif name == 'NSYM':
      if do_desym:
        fo.attrs[name] = 1
      else:
        fo.attrs[name] = fi.attrs[name]

    elif name == 'IRREP_LABELS':
      if do_desym:
        fo.attrs[name] = np.array(['a  '], dtype='S')
      else:
        fo.attrs[name] = fi.attrs[name]

    elif name == 'STATE_IRREPS':
      fo.attrs[name] = fi.attrs[name]
      if do_desym:
        fo.attrs[name][:] = 1

    elif name == 'NBAS':
      if do_desym:
        fo.attrs[name] = [np.sum(fi.attrs[name])]
      else:
        fo.attrs[name] = fi.attrs[name]

    elif name in ['NATOMS_UNIQUE', 'NATOMS_ALL']:
      if do_desym:
        name2 = 'NATOMS_UNIQUE'
        if name == name2:
          continue
        fo.attrs[name2] = fi.attrs[name]
      else:
        fo.attrs[name] = fi.attrs[name]

    elif args['all']:
      fo.attrs[name] = fi.attrs[name]
      warning_item.append(name)

  for name in fi:

    if name in safe_to_copy:
      if do_desym and name in not_with_desym:
        continue
      copy_dataset(fi, fo, name)

    elif name in ['DESYM_MATRIX']:
      continue

    # Coordinates are transformed with the input rotation
    elif name in ['CENTER_COORDINATES', 'DESYM_CENTER_COORDINATES']:
      name2 = name
      if name.startswith('DESYM_'):
        name2 = name[6:]
      elif do_desym:
        continue
      copy_dataset(fi, fo, name, name2)
      if do_rot or do_trans or do_exch:
        coor = fi[name]
        fo[name2][:] = reorder @ (coor @ R.T + T)

    # Multipole centers too, but they are not reordered
    elif name == 'MLTPL_ORIG':
      copy_dataset(fi, fo, name)
      if do_rot or do_trans:
        coor = fi[name]
        fo[name][:] = coor @ R.T + T

    # Other center data are reordered only
    elif name in ['CENTER_ATNUMS', 'CENTER_CHARGES', 'DESYM_CENTER_ATNUMS', 'DESYM_CENTER_CHARGES']:
      name2 = name
      if name.startswith('DESYM_'):
        name2 = name[6:]
      elif do_desym:
        continue
      copy_dataset(fi, fo, name, name2)
      if do_exch:
        data = fi[name]
        fo[name2][:] = data[:][reindex]

    # Special care with strings (can't just copy dataset because lengths may be different)
    elif name in ['CENTER_LABELS', 'DESYM_CENTER_LABELS']:
      name2 = name
      if name.startswith('DESYM_'):
        name2 = name[6:]
      elif do_desym:
        continue
      copy_dataset(fi, fo, name, name2, size=(nat,))
      if do_exch or do_desym:
        data = fi[name][:]
        # Find out if labels need renumbering
        any_num = False
        do_renum = True
        for i,lab in enumerate(fi[name2]):
          match = re.match(r'[A-Za-z]{1,2}(\d+)\s*', lab.decode())
          # Need renumbering if all labels have a number and all match the index
          if match:
            any_num = True
            if int(match.group(1)) != i+1:
              do_renum = False
          else:
            do_renum = False
        if do_renum:
          for i,lab in enumerate(data[reindex]):
            data[reindex[i]] = re.sub(r'^([A-Za-z]{1,2})(\d+).{4}', r'\g<1>'+f'{i+1}    ', lab.decode())
        elif any_num:
          warning_renum = True
        fo[name2][:] = data[reindex]
      else:
        fo[name2][:] = fi[name][:]

    # Basis info must be reordered and reindexed
    elif name in ['BASIS_FUNCTION_IDS', 'DESYM_BASIS_FUNCTION_IDS']:
      name2 = name
      if name.startswith('DESYM_'):
        name2 = name[6:]
      elif do_desym:
        continue
      copy_dataset(fi, fo, name, name2)
      if do_exch or do_desym:
        bfid = fi[name][:][bf_reindex]
        for b in bfid:
          b[0] = reindex_inv[b[0]-1]+1
        fo[name2][:] = bfid

    # Primitive info is similar to basis, but it was not read before
    elif name in ['PRIMITIVES', 'PRIMITIVE_IDS']:
      # Just to make sure these are processed in the right order
      if name == 'PRIMITIVES':
        continue
      copy_dataset(fi, fo, name)
      name = 'PRIMITIVES'
      copy_dataset(fi, fo, name)
      if do_exch:
        name = 'PRIMITIVE_IDS'
        prid = fi[name][:]
        pr_reindex = []
        for i in reindex:
          pr_reindex.extend([j for j,p in enumerate(prid) if p[0] == i+1])
        for p in prid:
          p[0] = reindex_inv[p[0]-1]+1
        fo[name][:] = prid[pr_reindex]
        name = 'PRIMITIVES'
        prim = fi[name][:][pr_reindex]
        fo[name][:] = prim

    # MO coefficients are transformed with the U matrix
    elif name == 'MO_VECTORS':
      copy_dataset(fi, fo, name, size=(nb**2,))
      if do_rot or do_desym:
        mo = sblock_to_square(fi[name], nbas)
        if do_desym:
          mo = (mo @ desym)[mo_reindex]
        fo[name][:] = (mo @ U).flatten()
    elif name == 'MO_ALPHA_VECTORS':
      copy_dataset(fi, fo, name, size=(nb**2,))
      if do_rot or do_desym:
        mo = sblock_to_square(fi[name], nbas)
        if do_desym:
          mo = (mo @ desym)[mo_alpha_reindex]
        fo[name][:] = (mo @ U).flatten()
    elif name == 'MO_BETA_VECTORS':
      copy_dataset(fi, fo, name, size=(nb**2,))
      if do_rot or do_desym:
        mo = sblock_to_square(fi[name], nbas)
        if do_desym:
          mo = (mo @ desym)[mo_beta_reindex]
        fo[name][:] = (mo @ U).flatten()

    # Other MO data are just reordered if desymmetrized
    elif name in ['MO_ENERGIES', 'MO_OCCUPATIONS', 'MO_TYPEINDICES', 'SUPSYM_IRREP_INDICES']:
      copy_dataset(fi, fo, name)
      if do_desym:
        data = fi[name]
        fo[name][:] = data[:][mo_reindex]
    elif name in ['MO_ALPHA_ENERGIES', 'MO_ALPHA_OCCUPATIONS', 'MO_ALPHA_TYPEINDICES']:
      copy_dataset(fi, fo, name)
      if do_desym:
        data = fi[name]
        fo[name][:] = data[:][mo_alpha_reindex]
    elif name in ['MO_BETA_ENERGIES', 'MO_BETA_OCCUPATIONS', 'MO_BETA_TYPEINDICES']:
      copy_dataset(fi, fo, name)
      if do_desym:
        data = fi[name]
        fo[name][:] = data[:][mo_beta_reindex]

    # AO matrices are transformed with the V matrix
    elif name in ['AO_OVERLAP_MATRIX', 'AO_FOCKINT_MATRIX']:
      copy_dataset(fi, fo, name, size=(nb**2,))
      if do_rot or do_desym:
        mtrx = sblock_to_square(fi[name], nbas)
        fo[name][:] = (V @ mtrx @ V.T).flatten()

    # Multipole AO matrices are further transformed among their components
    elif name in ['AO_MLTPL_X', 'AO_MLTPL_Y', 'AO_MLTPL_Z']:
      comp = ['X', 'Y', 'Z']
      # Do all work only when the first component is found
      if name.endswith(('_'+comp[1], '_'+comp[2])):
        continue
      ao_mltpl1 = np.empty((len(comp), nb, nb))
      for i,lab in enumerate(comp):
        name = 'AO_MLTPL_'+lab
        copy_dataset(fi, fo, name)
        ao_mltpl1[i,:,:] = fi[name]
      if do_rot:
        # Transform each component
        for i in range(len(comp)):
          ao_mltpl1[i,:,:] = V @ ao_mltpl1[i,:,:] @ V.T 
        # Combine components
        ao_mltpl1 = np.tensordot(R, ao_mltpl1, axes=1)
        for i,lab in enumerate(comp):
          name = 'AO_MLTPL_'+lab
          fo[name][:] = ao_mltpl1[i,:,:]
    elif name in ['AO_MLTPL_XX', 'AO_MLTPL_XY', 'AO_MLTPL_XZ', 'AO_MLTPL_YY', 'AO_MLTPL_YZ', 'AO_MLTPL_ZZ']:
      comp = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
      # Do all work only when the first component is found
      if name.endswith(('_'+comp[1], '_'+comp[2], '_'+comp[3], '_'+comp[4], '_'+comp[5])):
        continue
      ao_mltpl2 = np.empty((len(comp), nb, nb))
      for i,lab in enumerate(comp):
        name = 'AO_MLTPL_'+lab
        copy_dataset(fi, fo, name)
        ao_mltpl2[i,:,:] = fi[name]
      if do_rot:
        # Transform each component
        for i in range(len(comp)):
          ao_mltpl2[i,:,:] = V @ ao_mltpl2[i,:,:] @ V.T 
        # Quadrupoles transform like a Cartesian d shell (the inverse)
        if -2 in RotM:
          QRot = np.linalg.inv(RotM[-2])
        else:
          QRot = np.linalg.inv(cartesian_shell_rotation(2, R))
        # Combine components
        ao_mltpl2 = np.tensordot(QRot, ao_mltpl2, axes=1)
        for i,lab in enumerate(comp):
          name = 'AO_MLTPL_'+lab
          fo[name][:] = ao_mltpl2[i,:,:]

    # Density matrices may need to be reordered if desymmetrizing
    elif name in ['DENSITY_MATRIX', 'SPINDENSITY_MATRIX', 'TRANSITION_DENSITY_MATRIX', 'TRANSITION_SPIN_DENSITY_MATRIX']:
      copy_dataset(fi, fo, name)
      if  do_desym:
        if module == 'CASPT2':
          # In CASPT2, DM is in triangular storage
          dmsize = list(fo[name].shape)
          del fo[name]
          n = int(np.sum(pt2nbas))
          dmsize[1] = n*(n+1)//2
          copy_dataset(fi, fo, name, size=dmsize)
          for i in range(dmsize[0]):
            dm = sblock_to_square(fi[name][i,:], pt2nbas, triangular=True)
            fo[name][i,:] = dm[np.tril_indices(n)]
        else:
          dm = fi[name][:]
          for i in range(dm.shape[0]):
            dm[i,:,:] = act_reorder @ dm[i,:,:] @ act_reorder.T
            fo[name][:] = dm

    elif args['all']:
      copy_dataset(fi, fo, name)
      warning_item.append(name)

# In case something unsafe was done
if warning_item and do_smth:
  print('WARNING: The following attributes and/or datasets are not recognized, and were copied verbatim.')
  print('         If they would need transformation, they are probably WRONG!')
  print()
  for name in warning_item:
    print(name)
if warning_renum:
  print('WARNING: Atoms were reordered, but labels are not renumbered')
