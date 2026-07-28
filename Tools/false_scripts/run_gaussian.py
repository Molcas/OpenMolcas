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
import re

################################################################################
# Sample script to run Gaussian (https://gaussian.com/) with OpenMolcas FALSE  #
#                MAKE SURE TO TEST AND ADAPT TO YOUR NEEDS                     #
################################################################################
# Environment and executable settings
env = os.environ.copy()
command = 'g16'

# System and calculation settings
link0 = '%Mem = 2GB\n%NProcShared = 8' # "Link 0" commands
charge = 0       # total charge
multiplicity = 1 # spin multiplicity (1 = singlet, 2 = doublet...)
method = 'wB97XD/cc-pVDZ TD(NStates=3,Root=1)' # method, basis set and additional options
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

# Write GAUSSIAN input

keyline = ['#']
keyline.append(method)
keyline.append('NoSymm')
if der > 0:
  keyline.append('Force')

tmpfile = sys.argv[2]+'.tmp'
with open(tmpfile, 'w') as f:
  if link0:
    f.write(f'{link0}\n')
  f.write(' '.join(keyline)+'\n\n')
  f.write('Title\n\n')
  f.write(f'{charge} {multiplicity}\n')
  for l in geometry:
    f.write(f'{l[0]} {l[1]} {l[2]} {l[3]}\n')
  f.write('\n')

# Run GAUSSIAN

subprocess.run([command, tmpfile], env=env)

# Read GAUSSIAN output

eV = 0.0367493221757   # 1 eV in Eh
debye = 0.393430269519 # 1 D in e*a0

nroot = 1
rlxroot = 1
energy = []
dipole = []
gradient = []
tmpoutfile = os.path.splitext(tmpfile)[0]+'.log'
with open(tmpoutfile, 'r') as f:
  for l in f:
    if l.strip().startswith('SCF Done'):
      energy.append(float(l.split()[4]))
    elif re.match(r'Excited State *\d+:', l.strip()):
      nroot += 1
      energy.append(energy[0]+float(l.split()[4])*eV)
    elif l.strip().startswith('This state for optimization'):
      rlxroot = len(energy)
    elif l.strip().startswith('Dipole moment'):
      for i in range(nroot):
        dipole.append([0.0, 0.0, 0.0])
      d = next(f).split()
      dipole[rlxroot-1] = [float(d[i])*debye for i in [1, 3, 5]]
    elif l.strip().endswith('Forces (Hartrees/Bohr)'):
      next(f)
      next(f)
      for i in range(num):
        g = next(f).split()
        gradient.extend([-float(g[i]) for i in [2, 3, 4]])

# Write FALSE output

with open(sys.argv[2], 'w') as f:
  f.write('[ROOTS]\n')
  f.write(f'{nroot}\n')
  f.write('[RELAX ROOT]\n')
  f.write(f'{rlxroot}\n')
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
    f.write(f'{rlxroot}\n')
    for i in range(num):
      f.write(f'{gradient[3*i]} {gradient[3*i+1]} {gradient[3*i+2]}\n')
