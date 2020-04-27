#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
# Copyright (C) 2017, Ignacio Fdez. Galván                             *
#***********************************************************************

# This is a convenience script to generate include files for link_blas.f
# The files are generated from a list of procedures from LAPACK/BLAS

from __future__ import (unicode_literals, division, absolute_import, print_function)
import sys
sys.dont_write_bytecode = True

import re
import os.path

# Read the lists of function from blas.fh, lapack.fh and legacy_mod.f

blas_functions = []
with open('../Include/blas.fh', 'r') as f:
  for line in f:
    match = re.match('#include "(.*)\.[fF]"', line)
    if (match):
      blas_functions.append(match.group(1))
    match = re.match('!.* defines (.*)', line)
    if (match):
      blas_functions.extend(match.group(1).split(', '))
blas_functions.sort()

lapack_functions = []
with open('../Include/lapack.fh', 'r') as f:
  for line in f:
    match = re.match('#include "(.*)\.[fF]"', line)
    if (match):
      lapack_functions.append(match.group(1))
    match = re.match('!.* defines (.*)', line)
    if (match):
      lapack_functions.extend(match.group(1).split(', '))
lapack_functions.sort()

legacy_functions = []
with open('legacy_mod.f', 'r') as f:
  for line in f:
    match = re.match('\s*subroutine (.*)\(', line)
    if (match):
      legacy_functions.append(match.group(1))
    match = re.match('!.* defines (.*)', line)
legacy_functions.sort()

# Generate the include files f1.fh, f2.fh, f3.fh, f4.fh, f5.fh

header = '! File generated with {0}\n\n'.format(os.path.basename(__file__))

with open('f1.fh', 'w') as f:
  f.write(header)
  f.write('  use blas_mod')
  for i in blas_functions:
    f.write(', &\n      int_{0}=>{0}'.format(i))
  f.write('\n')
  f.write('  use lapack_mod')
  for i in lapack_functions:
    f.write(', &\n      int_{0}=>{0}'.format(i))
  f.write('\n')
  f.write('  use legacy_mod')
  for i in legacy_functions:
    f.write(', &\n      int_{0}=>{0}'.format(i))
  f.write('\n')

with open('f2.fh', 'w') as f:
  f.write(header)
  f.write('! BLAS procedures\n')
  for i in blas_functions:
    f.write('  procedure(int_{0}), pointer :: lb_{0}\n'.format(i))
  f.write('! LAPACK procedures\n')
  for i in lapack_functions:
    f.write('  procedure(int_{0}), pointer :: lb_{0}\n'.format(i))
  f.write('! Legacy procedures\n')
  for i in legacy_functions:
    f.write('  procedure(int_{0}), pointer :: lb_{0}\n'.format(i))

with open('f3.fh', 'w') as f:
  f.write(header)
  f.write('''!
!     BLAS procedures\n''')
  for i in blas_functions:
    f.write('''!
      funptr=link_func('{0}')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_{0})
      end if\n'''.format(i))
  f.write('''!
!     LAPACK procedures\n''')
  for i in lapack_functions:
    f.write('''!
      funptr=link_func('{0}')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_{0})
      end if\n'''.format(i))
  f.write('''!
!     Legacy procedures\n''')
  for i in legacy_functions:
    f.write('''!
      funptr=link_func('{0}')
      if (c_associated(funptr)) then
        call c_f_procpointer(funptr, lb_{0})
      end if\n'''.format(i))

with open('f4.fh', 'w') as f:
  f.write(header)
  f.write('''      !
      ! BLAS\n''')
  for i in blas_functions:
    f.write('      lb_{0}=>int_{0}\n'.format(i))
  f.write('''      !
      ! LAPACK\n''')
  for i in lapack_functions:
    f.write('      lb_{0}=>int_{0}\n'.format(i))
  f.write('''      !
      ! Legacy\n''')
  for i in legacy_functions:
    f.write('      lb_{0}=>int_{0}\n'.format(i))

with open('f5.fh', 'w') as f:
  f.write(header)
  f.write('''      ! BLAS
      !\n''')
  for i in blas_functions:
    f.write('''      if (DLAddr(c_funloc(lb_{0}),c_loc(info)) /= 0) then
        write(6,*) '{0} from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no {0} found!'
      end if\n'''.format(i))
  f.write('''
      ! LAPACK
      !\n''')
  for i in lapack_functions:
    f.write('''      if (DLAddr(c_funloc(lb_{0}),c_loc(info)) /= 0) then
        write(6,*) '{0} from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no {0} found!'
      end if\n'''.format(i))
  f.write('''
      ! Legacy
      !\n''')
  for i in legacy_functions:
    f.write('''      if (DLAddr(c_funloc(lb_{0}),c_loc(info)) /= 0) then
        write(6,*) '{0} from: ',c_f_string(info%dli_fname)
      else
        write(6,*) 'no {0} found!'
      end if\n'''.format(i))
