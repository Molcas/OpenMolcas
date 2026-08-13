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
# Sample script to run ORCA (https://www.faccts.de/orca) with OpenMolcas FALSE #
#                MAKE SURE TO TEST AND ADAPT TO YOUR NEEDS                     #
################################################################################
# Environment and executable settings
env = os.environ.copy()
command = 'orca'

# System and calculation settings
charge = 0         # total charge
multiplicity = 1   # spin multiplicity (1 = singlet, 2 = doublet...)
method = 'CAM-B3LYP DEF2-SVP TightSCF defgrid3 PAL8' # method, basis set and additional options
opts = '%tddft NRoots 3 IRoot 1 end' # other specific options
der = 1            # requested energy derivatives (0 = energy, 1 = gradient)
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

# Write ORCA input

keyline = ['!']
keyline.append(method)
if der > 0:
  keyline.append('EnGrad')

tmpfile = sys.argv[2]+'.tmp'
with open(tmpfile, 'w') as f:
  f.write(' '.join(keyline)+'\n')
  if opts:
    f.write(opts+'\n')
  f.write(f'* XYZ {charge} {multiplicity}\n')
  for l in geometry:
    f.write(f'{l[0]} {l[1]} {l[2]} {l[3]}\n')
  f.write('*\n')

# Run ORCA

with open(sys.argv[2]+'.log', 'w') as f:
  result = subprocess.run([command, tmpfile], env=env, stdout=f)

# Read ORCA output

nroot = 1
rlxroot = 1
energy = []
dipole = []
gradient = []
tmpoutfile = os.path.splitext(tmpfile)[0]+'.log'
with open(tmpoutfile, 'r') as f:
  for l in f:
    if l.strip().startswith('TOTAL SCF ENERGY'):
      next(f)
      next(f)
      energy.append(float(next(f).split()[3]))
    elif re.match(r'^STATE *\d+:', l.strip()):
      nroot += 1
      energy.append(energy[0]+float(l.split()[3]))
    elif l.strip().startswith('IROOT'):
      rlxroot = int(l.split()[1])+1
    elif l.strip().startswith('Total Dipole Moment'):
      if not dipole:
        for i in range(nroot):
          dipole.append([0.0, 0.0, 0.0])
      dipole[rlxroot-1] = [float(l.split()[i]) for i in [4, 5, 6]]
    elif l.strip().startswith('The cartesian gradient'):
      for i in range(num):
        g = next(f).split()
        gradient.extend([float(g[i]) for i in [3, 4, 5]])
    elif l.strip().startswith('CARTESIAN GRADIENT'):
      next(f)
      next(f)
      for i in range(num):
        g = next(f).split()
        gradient.extend([float(g[i]) for i in [3, 4, 5]])

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
