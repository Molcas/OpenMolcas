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
# Sample script to run MOPAC (https://openmopac.net/) with OpenMolcas FALSE    #
#                MAKE SURE TO TEST AND ADAPT TO YOUR NEEDS                     #
################################################################################
# Environment and executable settings
env = os.environ.copy()
env['LD_LIBRARY_PATH'] = os.path.expanduser('~/install/mopac/lib64')
command = os.path.expanduser('~/install/mopac/bin/mopac')

# System and calculation settings
charge = 0       # total charge
multiplicity = 1 # spin multiplicity (1 = singlet, 2 = doublet...)
method = 'PM7'   # method keyword
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

# Write MOPAC input

keyline = ['1SCF'] # Single-point calculation
keyline.append(method)
keyline.append(f'CHARGE={charge}')
keyline.append(f'MS={(multiplicity-1)/2}')
if der > 0:
  keyline.append('GRADIENTS')

tmpfile = sys.argv[2]+'.tmp'
with open(tmpfile, 'w') as f:
  f.write(' '.join(keyline)+'\n\n\n')
  for l in geometry:
    f.write(f'{l[0]} {l[1]} {l[2]} {l[3]}\n')

# Run MOPAC

subprocess.run([command, tmpfile], env=env)

# Read MOPAC output

kjmol = 0.000380879884713        # 1 kJ/mol in Eh
debye = 0.393430269519           # 1 D in e*a0
kcalangstrom = 0.000843297564063 # 1 kcal/mol/Å in Eh/a0

energy = []
dipole = []
gradient = []
with open(tmpfile+'.out', 'r') as f:
  for l in f:
    if l.strip().startswith('FINAL HEAT OF FORMATION'):
      energy.append(float(l.split()[8])*kjmol)
    elif l.strip().startswith('DIPOLE'):
      next(f)
      next(f)
      l = next(f)
      d = l.split()[1:]
      dipole.append([float(d[0])*debye, float(d[1])*debye, float(d[2])*debye])
    elif l.strip().startswith('FINAL  POINT  AND  DERIVATIVES'):
      next(f)
      next(f)
      for i in range(3*num):
        g = next(f).split()[6]
        gradient.append(float(g)*kcalangstrom)

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
