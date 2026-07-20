#! /usr/bin/env python3

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
# Copyright (C) 2026, Ignacio Fdez. Galván                             *
#***********************************************************************

import sys
import os
import subprocess

################################################################################
# Sample script to run xTB (https://xtb-docs.readthedocs.io) with              #
# OpenMolcas FALSE                                                             #
#                MAKE SURE TO TEST AND ADAPT TO YOUR NEEDS                     #
################################################################################
# Environment and executable settings
env = os.environ.copy()
command = os.path.expanduser('~/Programs/xtb/build/xtb')

# System and calculation settings
charge = 0       # total charge
multiplicity = 1 # spin multiplicity (1 = singlet, 2 = doublet...)
der = 1          # requested energy derivatives (0 = energy, 1 = gradient)
################################################################################

# Read FALSE input

geometry = []
with open(sys.argv[1], 'r') as f:
  for l in f:
    if l.strip() == '[XYZ]':
      num = int(f.readline())
      f.readline()
      for i in range(num):
        s, x, y, z = f.readline().split()
        geometry.append([s, float(x), float(y), float(z)])

# Write xTB input

tmpfile = sys.argv[2]+'.inp'
with open(tmpfile, 'w') as f:
  f.write(f'$chrg {charge}\n')
  f.write(f'$spin {multiplicity-1}\n')

xyzfile = sys.argv[2]+'.xyz'
with open(xyzfile, 'w') as f:
  f.write(f'{num}\n\n')
  for l in geometry:
    f.write(f'{l[0]} {l[1]} {l[2]} {l[3]}\n')

# Run xTB

args = ['--chrg', f'{charge}']
args += ['--uhf', f'{multiplicity-1}']
if der > 0:
  args += ['--grad']
args += [xyzfile]

with open(sys.argv[2]+'.out', 'w') as f:
  result = subprocess.run([command] + args, env=env, stdout=f)

# Read xTB output

energy = []
dipole = []
gradient = []
with open(sys.argv[2]+'.out', 'r') as f:
  for l in f:
    if l.strip().startswith(':: total energy'):
      energy.append(float(l.split()[3]))
    elif l.strip().startswith('molecular dipole'):
      next(f)
      next(f)
      d = next(f).split()
      dipole.append([float(d[i]) for i in [1, 2, 3]])
if der > 0:
  with open('gradient', 'r') as f:
    for i in range(num+2):
      next(f)
    for i in range(num):
      gradient.extend([float(i) for i in next(f).split()])

# Write FALSE output

with open(sys.argv[2], 'w') as f:
  f.write('[ROOTS]\n')
  f.write('1\n')
  if energy:
    f.write('[ENERGIES]\n')
    for e in energy:
      f.write(f'{e}\n')
  if dipole:
    f.write('[DIPOLES]\n')
    for d in dipole:
      f.write(f'{d[0]} {d[1]} {d[2]}\n')
  if gradient:
    f.write('[GRADIENT]\n')
    f.write('1\n')
    for i in range(num):
      f.write(f'{gradient[3*i]} {gradient[3*i+1]} {gradient[3*i+2]}\n')
