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
# Copyright (C) 2021,2024, Ignacio Fdez. Galván                        *
#***********************************************************************

import sys
import re
import argparse
import math
import numpy as np
from fractions import Fraction
import h5py

version = '2.1.1'

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

# Double factorial of 2n-1 (for n >= 0)
def xfact(n):
  if n <= 0:
    return 1 if n == 0 else None
  else:
    return (2*n-1)*xfact(n-1)

# Binomial coefficient as a fraction
# Easy overflow for large arguments, but we are interested in relatively small arguments
def binom(n, k):
  mk = max(k,n-k)
  try:
    binom = Fraction(math.factorial(n), math.factorial(mk))
    binom *= Fraction(1, math.factorial(n-mk))
    assert binom.denominator == 1
  except ValueError:
    binom = Fraction(0, 1)
  return binom

# Compute the coefficient for x^lx * y^ly * z^lz in the expansion of
# the real solid harmonic S(l,±m) = C * r^l*(Y(l,m)±Y(l,-m))
# Since the coefficients are square roots of rational numbers, this
# returns the square of the coefficient as a fraction, with its sign
#
# See:
# Transformation between Cartesian and pure spherical harmonic Gaussians
# doi: 10.1002/qua.560540202
# (note that there appears to be an error in v(4,0), the coefficient 1/4
#  should probably be 3/4*sqrt(3/35) )
def cart_coeff(l, m, lx, ly, lz):
  assert (lx + ly + lz == l) and (lx >= 0) and (ly >= 0) and (lz >= 0)
  am = abs(m)
  assert (am <= l)
  j = lx + ly - am
  if j%2 == 0:
    j = j//2
  else:
    return Fraction(0, 1)
  c = 0
  for i in range((l-am)//2+1):
    c += binom(l, i) * binom(i, j) * Fraction(math.factorial(2*l-2*i), math.factorial(l-am-2*i)) * (-1)**i
  if c == 0:
    return Fraction(0, 1)
  c_sph = c
  c = 0
  for k in range(j+1):
    c += binom(j, k) * binom(am, lx-2*k) * 1j**(am-lx+2*k)
  if m >= 0:
    c = int(np.real(c))
  else:
    c = int(np.imag(c))
  if c == 0:
    return Fraction(0, 1)
  c_sph *= c
  if c_sph < 0:
    c_sph *= -c_sph
  else:
    c_sph *= c_sph
  if (m == 0):
    lm = 1
  else:
    lm = 2
  c = Fraction(math.factorial(l-am), math.factorial(l+am))
  c *= Fraction(lm, math.factorial(l))
  c *= Fraction(1, math.factorial(2*l))
  c_sph *= c
  return c_sph

# Precompute the Cartesian coefficients for spherical harmonics up to l_max
def sph_to_cart(l_max):
  # Get the coefficients for each value of l,m
  sph_c = []
  for l in range(l_max+1):
    s = {}
    for m in range(-l, l+1):
      s[m] = []
      # Go through all possible lx+ly+lz=l
      for lx in range(l+1):
        for ly in range(l-lx+1):
          lz = l-lx-ly
          # Get the coefficient (cart_coeff returns the square as a fraction with the sign)
          c = float(cart_coeff(l, m, lx, ly, lz))
          c = np.sign(c)*np.sqrt(abs(c))
          if (c != 0):
            s[m].append([c, [lx, ly, lz]])
    sph_c.append(s)
  return sph_c

# Scale the primitive coefficients to normalize the radial part
# Note that the coefficients refer to normalized primitives
def normalize_rad(shell):
  e = shell['prims']
  c = shell['coeffs']
  l = shell['l']
  nprim = e.shape[0]
  # Normalized primitives, so just sum the square coefficients
  integral = np.sum(c**2, axis=1)
  # Add the overlap between each pair of normalized primitives
  for i in range(nprim):
    for j in range(i+1, nprim):
      cc = c[:,i]*c[:,j]
      ep = np.sqrt(e[i]*e[j])
      es = (e[i]+e[j])/2
      integral += 2*cc*np.power(ep/es, 1.5+l)
  c /= np.sqrt(integral[:,np.newaxis])

# Compute the Cartesian expansion of r^(2*n), where r^2 = x^2 + y^2 + z^2
def r2_power(n):
  if n == 0:
    return {(0, 0, 0): 1}
  basis = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
  terms = {}
  for i in it.product(basis, repeat=n):
    t = tuple(np.sum(i, axis=0))
    terms[t] = terms.get(t, 0)+1
  return terms

# Find and label contaminant functions
# This is not 100% sure, but it should work in most practical situations
def label_contaminants(b_id, c_labs):
  cont = {}
  for i,b in enumerate(b_id):
    if b[2] < 0:
      continue
    # Only test the first radial function of each shell
    if b[1] != 1:
      continue
    # Has the same "label" already been used?
    for j,c in enumerate(b_id[:i]):
      if all(b == c):
        cont[i] = 0
        break
  # Find out which centers are equivalent
  equiv = []
  for i,lab in enumerate(c_labs):
    this = lab.split(b':')[0]
    for j,lab2 in enumerate(c_labs[:i+1]):
      that = lab2.split(b':')[0]
      if that == this:
        equiv.append(j+1)
        break
  # Now for each duplicated label, find the real shell (loop backwards to find the most likely parent)
  for i in list(cont):
    b = b_id[i]
    for j,c in reversed(list(enumerate(b_id[:i]))):
      if j in cont or c[2] < 0:
        continue
      if (b[0] == c[0]) and (b[1] == c[1]) and (c[2] > b[2]) and ((c[2]-b[2])%2 == 0):
        cont[i] = (c[2]-b[2])//2
        break
    assert cont[i] > 0
    # Label all other consecutive functions in the shell (taking equivalent atoms into account)
    for jj,c in enumerate(b_id[i+1:]):
      if (equiv[c[0]-1] != equiv[b[0]-1]) or (c[1] < b[1]) or (c[2] != b[2]) or (c[3] != b[3]):
        break
      j = jj+i+1
      if j in cont or c[2] < 0:
        continue
      cont[j] = cont[i]
  return cont

# Read shells (functions sharing center, set of primitive exponents and l value)
def read_shells(fileid):
  prefix = 'DESYM_' if fileid.attrs['NSYM'] > 1 else ''
  centers = fileid[prefix+'CENTER_COORDINATES']
  c_labs = fileid[prefix+'CENTER_LABELS']
  b_id = fileid[prefix+'BASIS_FUNCTION_IDS']
  prims = fileid['PRIMITIVES']
  prids = fileid['PRIMITIVE_IDS']
  # Group by shells of functions sharing center and l (different m values)
  shells = []
  shell_idx = {}
  cont = label_contaminants(b_id, c_labs)
  for i,b in enumerate(b_id):
    l = abs(b[2])
    r2 = cont.get(i, 0)
    l = l+2*r2
    label = (b[0], b[1], l, r2)
    if label in shell_idx:
      s = shells[shell_idx[label]]
      assert b[3] not in s['m']
      s['m'].append(b[3])
      s['b_ids'].append(i)
    else:
      shell_idx[label] = len(shells)
      shells.append({'xyz': centers[b[0]-1,:], 'l': l, 'cart': b[2] < 0, 'r2': r2, 'm': [b[3]], 'b_ids': [i], 'prims': [], 'coeffs': [[]]})
  # Assign exponents and coefficients
  for i,j in zip(prids, prims):
    l = i[1]
    # Assign to all shells and contaminants
    for r2 in range(l//2, -1, -1):
      label = (i[0], i[2], l, r2)
      if label in shell_idx:
        idx = shell_idx[label]
        shells[idx]['prims'].append(j[0])
        shells[idx]['coeffs'][0].append(j[1])
  # Merge shells sharing the same set of exponents (and same center, l and m values)
  for i,s in enumerate(shells):
    l = s['l']
    r2 = s['r2']
    num_m = len(s['m'])
    assert s['r2'] == 0 or not s['cart']
    assert (num_m == 2*(l-2*r2)+1 and not s['cart']) or (num_m == (l+1)*(l+2)/2 and s['cart'])
    if not s['b_ids']:
      continue
    for j,t in enumerate(shells[i+1:]):
      if np.any(t['xyz'] != s['xyz']) or any([t[i] != s[i] for i in ['l', 'cart', 'r2', 'm', 'prims']]):
        continue
      s['coeffs'].append(t['coeffs'][0])
      s['b_ids'].extend(t['b_ids'])
      t['b_ids'] = None
  del shell_idx
  # Normalize, assign coefficients for angular part and clean up
  l_max = max([s['l'] for s in shells])
  sph_coeff = sph_to_cart(l_max)
  for s in shells:
    if not s['b_ids']:
      continue
    m_list = s.pop('m')
    # Reorder the basis function ids according to the final arrangement
    s['b_ids'] = np.reshape(s['b_ids'], (len(s['coeffs']), len(m_list))).flatten('F').tolist()
    s['prims'] = np.array(s['prims'])
    s['coeffs'] = np.array(s['coeffs'])
    # Normalize coefficients
    normalize_rad(s)
    # Assign angular part
    l = s['l']
    r2 = s['r2']
    lr = l-2*r2
    m_prims = []
    for lx in range(l,-1,-1):
      for ly in range(l-lx,-1,-1):
        m_prims.append((lx, ly, l-lx-ly))
    s['m_coeffs'] = np.zeros((len(m_list), len(m_prims)))
    if s['cart']:
      # For Cartesian functions, m starts at -l, but the m_prims list starts at 0
      for i,m in enumerate(m_list):
        j = m+l
        s['m_coeffs'][i,j] = np.sqrt(2**l)
    else:
      # Expansion of r^2 power
      r2_terms = r2_power(r2)
      # For spherical harmonics, the coefficients have been precomputed by sph_to_cart above
      for i,m in enumerate(m_list):
        for coeff in sph_coeff[lr][m]:
          for t,c in r2_terms.items():
            lx, ly, lz = (coeff[1][0]+t[0], coeff[1][1]+t[1], coeff[1][2]+t[2])
            j = (ly+lz)*(ly+lz+1)//2+lz
            assert (lx, ly, lz) == m_prims[j]
            s['m_coeffs'][i,j] += coeff[0]*c
      # We include the extra normalization for contaminants here
      if r2 > 0:
        p = 1
        for r in range(lr+1,l+1):
          p *= (2*r+1)/2
        s['m_coeffs'] /= np.sqrt(p)
    s['m_prims'] = np.array(m_prims)
  # Merge contaminant shells with their parents
  for i,s in enumerate(shells):
    l = s['l']
    if (l < 2) or not s['b_ids']:
      continue
    assert s['r2'] == 0
    for j,t in enumerate(shells[i+1:]):
      if not t['b_ids']:
        continue
      if any([not np.array_equal(t[i], s[i]) for i in ['xyz', 'prims', 'coeffs']]) or any([t[i] != s[i] for i in ['l', 'cart']]):
        continue
      assert np.array_equal(t['m_prims'], s['m_prims'])
      s['m_coeffs'] = np.vstack((s['m_coeffs'], t['m_coeffs']))
      s['b_ids'].extend(t['b_ids'])
      t['b_ids'] = None
  # Remove merged shells
  shells = [s for s in shells if s['b_ids']]
  return shells

# Compute overlaps between two shells
def shell_overlap(shell_A, shell_B):
  A = shell_A['xyz']
  B = shell_B['xyz']
  d = A-B
  d = np.dot(d, d)
  # Primitive exponents expanded to 2D shape
  alpha, beta = np.meshgrid(shell_A['prims'], shell_B['prims'], indexing='ij', sparse=True)
  ab1 = alpha*beta
  ab2 = alpha+beta
  # Overlap of s functions
  if d > 0.0:
    E_AB = np.exp(-d*ab1/ab2)
  else:
    E_AB = np.ones_like(ab1)
  E_AB *= np.power(4*ab1/ab2**2, 0.75)
  # Normalization for higher l functions
  l_A = shell_A['l']
  l_B = shell_B['l']
  E_AB *= np.power(2*alpha, l_A/2)*np.power(2*beta, l_B/2)
  # Expand the integrals to all Cartesian power combinations
  nm_A = len(shell_A['m_prims'])
  nm_B = len(shell_B['m_prims'])
  prim_ovlp = np.tile(E_AB, (nm_A, nm_B, 1, 1))
  # Compute the angular contributions with the recurrence relations
  AA = A[:,np.newaxis,np.newaxis]
  BB = B[:,np.newaxis,np.newaxis]
  P = (alpha*AA+beta*BB)/ab2
  max_l = l_A+l_B
  # For each dimension: x, y, z
  for x in [0, 1, 2]:
    # First compute all power factors
    s = {}
    s[(0,0)] = 1.0
    if max_l > 0:
      s[(1,0)] = P[x,:,:]-AA[x,:,:]
      s[(0,1)] = P[x,:,:]-BB[x,:,:]
      for ax in range(2, max_l+1):
        s[(ax,0)] = s[(1,0)]*s[(ax-1,0)]+(ax-1)/(2*ab2)*s[(ax-2,0)]
        for bx in range(1, min(ax, l_B)+1):
          s[(ax-bx,bx)] = s[(ax-bx+1,bx-1)]+(AA[x,:,:]-BB[x,:,:])*s[(ax-bx,bx-1)]
    # Then apply them to the integrals
    for i,m in enumerate(shell_A['m_prims']):
      for j,n in enumerate(shell_B['m_prims']):
        prim_ovlp[i,j,:,:] *= s[(m[x], n[x])]
  del s
  # Finally transform the primitive integrals to the basis functions
  ovlp = np.einsum('ik,ma,jl,nb,abkl->minj', shell_A['coeffs'], shell_A['m_coeffs'], shell_B['coeffs'], shell_B['m_coeffs'], prim_ovlp)
  nf_A = shell_A['coeffs'].shape[0]*shell_A['m_coeffs'].shape[0]
  nf_B = shell_B['coeffs'].shape[0]*shell_B['m_coeffs'].shape[0]
  ovlp = np.reshape(ovlp, (nf_A, nf_B))
  return(ovlp)

def basis_overlap(shells_A, shells_B):
  nbas_A = sum([s['coeffs'].shape[0]*s['m_coeffs'].shape[0] for s in shells_A])
  nbas_B = sum([s['coeffs'].shape[0]*s['m_coeffs'].shape[0] for s in shells_B])
  ovlp = np.empty((nbas_A, nbas_B))
  i = 0
  for s in shells_A:
    j = 0
    for t in shells_B:
      o = shell_overlap(s, t)
      ii, jj = o.shape
      ovlp[i:i+ii,j:j+jj] = o
      j += jj
    i += ii
  # Reorder if needed
  reorder = False
  idx_A = [0 for i in range(nbas_A)]
  ib = 0
  for s in shells_A:
    for b in s['b_ids']:
      if b != ib:
        reorder = True
      idx_A[b] = ib
      ib += 1
  idx_B = [0 for i in range(nbas_B)]
  ib = 0
  for s in shells_B:
    for b in s['b_ids']:
      if b != ib:
        reorder = True
      idx_B[b] = ib
      ib += 1
  if reorder:
    ovlp = ovlp[np.ix_(idx_A, idx_B)]
  return ovlp


################################################################################

# Command-line arguments
parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=24,width=120))
parser.version = f'MORT: molecular orbital rotation and transfer (version {version})'
parser.add_argument('-v', '--version', help='print version number', action='version')
parser.add_argument('-r', '--rotate', help='rotation matrix to apply (9 real values)', type=float, nargs=9, metavar=('Rxx','Rxy','Rxz','Ryx','Ryy','Ryz','Rzx','Rzy','Rzz'))
parser.add_argument('-t', '--translate', help='translation vector to apply after rotation (3 real values, in bohr)', type=float, nargs=3, metavar=('Tx','Ty','Tz'))
parser.add_argument('-x', '--exchange', help='reorder atoms (comma-separated list of integers giving order in the output file)', metavar='n1,n2,...')
parser.add_argument('-d', '--desymmetrize', help='desymmetrize data', action='store_true')
parser.add_argument('-a', '--all', help='copy all unhandled datasets and attributes to output file (some items may be wrong!)', action='store_true')
parser.add_argument('-f', '--transfer', help='transfer (project) orbitals from input file onto output file\n'
                                             'transf_spec: "all", "occupied", "active", or comma-separated list of input orbital indices)\n'
                                             'incompatible with other options', metavar='transf_spec')
parser.add_argument('-s', '--select', help='specify which output orbitals to replace with --transfer ("help" for format)', metavar='select_spec')
parser.add_argument('--threshold', help=argparse.SUPPRESS or 'overlap threshold for orbital projection (default: 0.5)', type=float, default=0.5, metavar='T')

# Partial parse before specifying required arguments
args = vars(parser.parse_known_args()[0])
select_help = '''
select_spec format:
===================

t1:n1[,t2:n2[,...]]

t1, t2 ...: Orbital type labels (any combination of F, I, 1, 2, 3, S, D), case-insensitive
n1, n2 ...: Number of orbitals of the corresponding type

The different specified type labels must not overlap.
The n best-matching orbitals of each specified t will be selected.
Any remaining orbitals will be selected from the remaining types.

Examples:

i:10        Select the 10 best-matching inactive orbitals, the rest from frozen, secondary, etc.
123:4,sd:0  Select the 4 best-matching active orbitals and no virtuals, the rest from frozen and inactive.\
'''
if args['select'] == 'help':
  print(select_help)
  sys.exit(0)

parser.add_argument('infile', help='input file (HDF5 format)', metavar='input_file')
parser.add_argument('outfile', help='output file (HDF5 format, will be overwritten)', metavar='output_file')
if len(sys.argv) < 2 and sys.stdin.isatty():
  parser.print_help(sys.stderr)
  sys.exit(1)
else:
  args = vars(parser.parse_args())

# From here, there are two working modes:
#   A: with transfer (orbital projection)
#   B: without transfer (rotation, reordering, etc.)

# A: Orbital projection
if args['transfer']:

  # TODO: Symmetry support
  #       Alpha/beta orbitals

  # Some input processing

  if (args['rotate'] or args['translate'] or args['exchange'] or args['desymmetrize'] or args['all']):
    sys.exit('Error: The --transfer option is incompatible with --rotate|--translate|--exchange|--desymmetrize|--all')

  if args['transfer'] == 'all':
    numorbs = 'all'
  elif args['transfer'] == 'occupied':
    numorbs = 'occupied'
  elif args['transfer'] == 'active':
    numorbs = 'active'
  else:
    try:
      numorbs = [int(i) for i in args['transfer'].split(',')]
    except ValueError:
      sys.exit('Error: The --transfer option needs "all", "occupied", "active" or a comma-separated list of input orbital indices')

  with h5py.File(args['infile'], 'r') as f:
    nsym = f.attrs['NSYM']
    if nsym != 1:
      sys.exit('Error: Symmetry is not supported with --transfer')
    shells_A = read_shells(f)
    nbas = f.attrs['NBAS']
    ov_A = sblock_to_square(f['AO_OVERLAP_MATRIX'], nbas)
    C_A = sblock_to_square(f['MO_VECTORS'], nbas).T
    tp_A = f['MO_TYPEINDICES'][:]
    en_A = f['MO_ENERGIES'][:]
    oc_A = f['MO_OCCUPATIONS'][:]
    nbas_A = sum(nbas)
    try:
      desym = np.reshape(f['DESYM_MATRIX'], (nbas_A, nbas_A))
      ov_A = desym.T @ ov_B @ desym
      C_A = desym @ C_A
    except KeyError:
      pass
    if not np.allclose(C_A.T @ ov_A @ C_A, np.eye(nbas_A)):
      sys.exit(f'Error: Orbitals in {infile} are not orthonormal!')

  with h5py.File(args['outfile'], 'r') as f:
    nsym = f.attrs['NSYM']
    if nsym != 1:
      sys.exit('Error: Symmetry is not supported with --transfer')
    shells_B = read_shells(f)
    nbas = f.attrs['NBAS']
    ov_B = sblock_to_square(f['AO_OVERLAP_MATRIX'], nbas)
    C_B = sblock_to_square(f['MO_VECTORS'], nbas).T
    tp_B = f['MO_TYPEINDICES'][:]
    nbas_B = sum(nbas)
    try:
      desym = np.reshape(f['DESYM_MATRIX'], (nbas_B, nbas_B))
      ov_B = desym.T @ ov_B @ desym
      C_B = desym @ C_B
    except KeyError:
      pass
    if not np.allclose(C_B.T @ ov_B @ C_B, np.eye(nbas_B)):
      sys.exit(f'Error: Orbitals in {outfile} are not orthonormal!')

  if numorbs == 'all':
    numorbs = [i for i in range(nbas_A)]
  elif numorbs == 'occupied':
    numorbs = [i for i in range(nbas_A) if tp_A[i] in [b'F', b'I', b'1', b'2', b'3']]
  elif numorbs == 'active':
    numorbs = [i for i in range(nbas_A) if tp_A[i] in [b'1', b'2', b'3']]
  else:
    if min(numorbs) < 1 or max(numorbs) > nbas_A:
      sys.exit(f'Error: Orbital indices must be between 1 and {nbas_A}')
    numorbs = [i-1 for i in numorbs]
  if not numorbs:
    sys.exit(f'Error: No orbitals to transfer')

  # Parse the select_spec
  group_spec = [['fi123sd', len(numorbs)]]
  select = args['select'].split(',') if args['select'] else []
  try:
    for s in select:
      t, n = s.split(':')
      t = t.lower()
      n = int(n)
      assert n >= 0
      new = [t, n]
      for i in t:
        if i not in group_spec[-1][0]:
          raise Exception
        group_spec[-1][0] = group_spec[-1][0].replace(i, '')
      group_spec[-1][1] -= n
      if n > 0:
        group_spec.insert(-1, new)
  except:
    print(f'\nError in select_spec: {args["select"]}')
    print(select_help)
    sys.exit(1)

  groups = []
  for g in group_spec:
    new = [g[1], [i for i in range(nbas_B) if tp_B[i].decode().lower() in g[0]]]
    if new[0] > len(new[1]):
      print(f'\nNot enough orbitals to select: {args["select"]}')
      print(select_help)
      sys.exit(1)
    groups.append(new)

  # Sanity check

  S_A = basis_overlap(shells_A, shells_A)
  S_B = basis_overlap(shells_B, shells_B)
  if not np.allclose(S_A, ov_A) or not np.allclose(S_B, ov_B):
    sys.exit(f'Error: Computed overlaps do not match the file(s). This may be a bug')

  # This is our master key

  S_AB = basis_overlap(shells_A, shells_B)

  # Overlap of input orbitals (Cs: "selected" coefficients)
  Cs_A = C_A[:,numorbs]
  M = Cs_A.T @ S_AB @ C_B
  U, s, V = np.linalg.svd(M, full_matrices=False)
  print(f'Overlap between input orbitals and target space: {np.prod(s):8.6f}')
  print()
  num = np.count_nonzero(s > args['threshold'])
  if num < len(numorbs):
    sys.exit(f'Error: only {num} input orbitals have an overlap larger than {args["threshold"]} with the output space')

  # Match the target orbitals as much as possible with the input orbitals
  Cs_B = C_B @ V.T @ U.T

  d = np.diag(Cs_A.T @ S_AB @ Cs_B)
  num = np.count_nonzero(d > args['threshold'])
  if num < len(numorbs):
    print('Overlap between input orbitals and projected orbitals:')
    print(' '.join([f'{i:8.6f}' for i in d]))
    sys.exit(f'Error: only {num} output orbitals have an overlap larger than {args["threshold"]} with the input orbitals')

  # Identify the orbitals in the output space that best match the projected orbitals (approx.)
  # We will be replacing these and changing the rest as little as possible
  sel = []
  rest = []
  for g in groups:
    M = Cs_B.T @ S_B @ C_B[:,g[1]]
    idx = [g[1][i] for i in np.argsort(np.sum(M**2, axis=0))[::-1]]
    sel.extend(idx[:g[0]])
    rest.extend(idx[g[0]:])
  sel.sort()
  rest.sort()
  # project out the selected vectors from the rest
  P = Cs_B.T @ S_B @ C_B[:,rest]
  D = C_B[:,rest] - Cs_B @ P
  # symmetric orthogonalization of the result
  s, U = np.linalg.eigh(D.T @ S_B @ D)
  Dq = D @ U / np.sqrt(s) @ U.T
  if not np.allclose(Cs_A.T @ S_AB @ Dq, 0):
    sys.exit('Error: complementary space is contaminated. This may be a bug')
  # fill up the columns (Cn: new coefficients)
  Cn_B = np.empty((nbas_B, nbas_B))
  Cn_B[:,sel] = Cs_B
  Cn_B[:,rest] = Dq
  if not np.allclose(Cn_B.T @ S_B @ Cn_B, np.eye(nbas_B)):
    sys.exit('Error: output orbitals are not orthonormal. This may be a bug')

  print('Overlap between input orbitals and projected orbitals:')
  print('in    out   overlap')
  print('--------------------')
  for i,v in enumerate(d):
    print(f'{numorbs[i]+1:<5} {sel[i]+1:<5} {d[i]:8.6f}')

  with h5py.File(args['outfile'], 'r+') as f:
    f.attrs.create('MORT_version', np.array(f'MORT v{version}', dtype='S'))
    f['MO_VECTORS'][...] = Cn_B.flatten('F')
    f['MO_TYPEINDICES'][sel] = tp_A[numorbs]
    f['MO_ENERGIES'][sel] = en_A[numorbs]
    f['MO_OCCUPATIONS'][sel] = oc_A[numorbs]

# B: Rotation, reordering, etc.
else:

  if args['select']:
    sys.exit('Error: The --select option is intended for use with --transfer')

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

    translate_mltpl1 = False

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

      # Multipole origins too, but they are not reordered
      elif name == 'MLTPL_ORIG':
        copy_dataset(fi, fo, name)
        if do_rot or do_trans:
          coor = fi[name]
          fo[name][:] = coor @ R.T + T
          # Do not transform origins for order 0 and 1 multipoles if they are zero
          if fi[name].shape[0] > 0 and not np.any(fi[name][0,:]):
            fo[name][0,:] = fi[name][0,:]
          if fi[name].shape[0] > 1 and not np.any(fi[name][1,:]):
            fo[name][1,:] = fi[name][1,:]
            # We will need to translate the integrals instead
            translate_mltpl1 = do_trans

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
        else:
          fo[name][:] = fi[name][:]
      elif name == 'MO_ALPHA_VECTORS':
        copy_dataset(fi, fo, name, size=(nb**2,))
        if do_rot or do_desym:
          mo = sblock_to_square(fi[name], nbas)
          if do_desym:
            mo = (mo @ desym)[mo_alpha_reindex]
          fo[name][:] = (mo @ U).flatten()
        else:
          fo[name][:] = fi[name][:]
      elif name == 'MO_BETA_VECTORS':
        copy_dataset(fi, fo, name, size=(nb**2,))
        if do_rot or do_desym:
          mo = sblock_to_square(fi[name], nbas)
          if do_desym:
            mo = (mo @ desym)[mo_beta_reindex]
          fo[name][:] = (mo @ U).flatten()
        else:
          fo[name][:] = fi[name][:]

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
        else:
          fo[name][:] = fi[name][:]

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

    # If the origin for dipole integrals was not translated, we need to translate the integrals
    if translate_mltpl1:
      if 'AO_OVERLAP_MATRIX' in fo:
        if 'AO_MLTPL_X' in fo:
          ov = sblock_to_square(fo['AO_OVERLAP_MATRIX'], nbas)
          fo['AO_MLTPL_X'][:] += T[0]*ov
          fo['AO_MLTPL_Y'][:] += T[1]*ov
          fo['AO_MLTPL_Z'][:] += T[2]*ov
      else:
        print('WARNING: The AO_MLTPL_X, _Y, _Z integrals need to be translated, but AO_OVERLAP_MATRIX is missing.')
        print('         They are probably WRONG!')
        print()

  # In case something unsafe was done
  if warning_item and do_smth:
    print('WARNING: The following attributes and/or datasets are not recognized, and were copied verbatim.')
    print('         If they would need transformation, they are probably WRONG!')
    print()
    for name in warning_item:
      print(name)
  if warning_renum:
    print('WARNING: Atoms were reordered, but labels are not renumbered')
