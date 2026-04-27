#!/usr/bin/env python3
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
# Copyright (C) 2017,2026, Ignacio Fdez. Galván                        *
#***********************************************************************

# This is a convenience script to generate include files for link_blas.F90
# The files are generated from a list of procedures from LAPACK/BLAS

from __future__ import (unicode_literals, division, absolute_import, print_function)
import sys
sys.dont_write_bytecode = True

import re
import os.path

# Read the lists of function from blas.fh, lapack.fh and legacy_mod.F90

blas_functions = []
with open('../Include/blas.fh', 'r') as f:
  for line in f:
    match = re.match('#include "(.*)\.[fF](90)?"', line)
    if (match):
      blas_functions.append(match.group(1))
    match = re.match('!.* defines (.*)', line)
    if (match):
      blas_functions.extend(match.group(1).split(', '))
blas_functions.sort()

lapack_functions = []
with open('../Include/lapack.fh', 'r') as f:
  for line in f:
    match = re.match('#include "(.*)\.[fF](90)?"', line)
    if (match):
      lapack_functions.append(match.group(1))
    match = re.match('!.* defines (.*)', line)
    if (match):
      lapack_functions.extend(match.group(1).split(', '))
lapack_functions.sort()

legacy_functions = []
with open('legacy_mod.F90', 'r') as f:
  for line in f:
    match = re.match('\s*subroutine (.*)\(', line)
    if (match):
      legacy_functions.append(match.group(1))
    match = re.match('!.* defines (.*)', line)
legacy_functions.sort()

# Generate the include files f1.fh, f2.fh, f3.fh, f4.fh, f5.fh

header = f'''\
!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
! File generated with {os.path.basename(__file__)}\n
'''

with open('f1.fh', 'w') as f:
  f.write(header)
  f.write('use blas_mod, only:')
  for n,i in enumerate(blas_functions):
    if n > 0:
      f.write(',')
    f.write(f' &\n    int_{i}=>{i}')
  f.write('\n')
  f.write('use lapack_mod, only:')
  for n,i in enumerate(lapack_functions):
    if n > 0:
      f.write(',')
    f.write(f' &\n    int_{i}=>{i}')
  f.write('\n')
  f.write('use legacy_mod, only:')
  for n,i in enumerate(legacy_functions):
    if n > 0:
      f.write(',')
    f.write(f' &\n    int_{i}=>{i}')
  f.write('\n')

with open('f2.fh', 'w') as f:
  f.write(header)
  f.write('! BLAS procedures\n')
  for i in blas_functions:
    f.write(f'procedure(int_{i}), pointer, public :: lb_{i}\n')
  f.write('! LAPACK procedures\n')
  for i in lapack_functions:
    f.write(f'procedure(int_{i}), pointer, public :: lb_{i}\n')
  f.write('! Legacy procedures\n')
  for i in legacy_functions:
    f.write(f'procedure(int_{i}), pointer, public :: lb_{i}\n')

with open('f3.fh', 'w') as f:
  f.write(header)
  f.write('\n    ! BLAS procedures\n')
  for i in blas_functions:
    f.write(f'''
    funptr = link_func('{i}')
    if (c_associated(funptr)) call c_f_procpointer(funptr, lb_{i})\n''')
  f.write('\n    ! LAPACK procedures\n')
  for i in lapack_functions:
    f.write(f'''
    funptr = link_func('{i}')
    if (c_associated(funptr)) call c_f_procpointer(funptr, lb_{i})\n''')
  f.write('\n    ! Legacy procedures\n')
  for i in legacy_functions:
    f.write(f'''
    funptr = link_func('{i}')
    if (c_associated(funptr)) call c_f_procpointer(funptr, lb_{i})\n''')

with open('f4.fh', 'w') as f:
  f.write(header)
  f.write('\n    ! BLAS\n')
  for i in blas_functions:
    f.write(f'    lb_{i} => int_{i}\n')
  f.write('\n    ! LAPACK\n')
  for i in lapack_functions:
    f.write(f'    lb_{i} => int_{i}\n')
  f.write('\n    ! Legacy\n')
  for i in legacy_functions:
    f.write(f'    lb_{i} => int_{i}\n')

with open('f5.fh', 'w') as f:
  f.write(header)
  f.write('\n    ! BLAS\n\n')
  for i in blas_functions:
    f.write(f'''    if (DLAddr(c_funloc(lb_{i}),c_loc(info)) /= 0) then
      write(u6,*) '{i} from: ',c_f_string(info%dli_fname)
    else
      write(u6,*) 'no {i} found!'
    end if\n''')
  f.write('\n    ! LAPACK\n\n')
  for i in lapack_functions:
    f.write(f'''    if (DLAddr(c_funloc(lb_{i}),c_loc(info)) /= 0) then
      write(u6,*) '{i} from: ',c_f_string(info%dli_fname)
    else
      write(u6,*) 'no {i} found!'
    end if\n''')
  f.write('\n    ! Legacy\n\n')
  for i in legacy_functions:
    f.write(f'''    if (DLAddr(c_funloc(lb_{i}),c_loc(info)) /= 0) then
      write(u6,*) '{i} from: ',c_f_string(info%dli_fname)
    else
      write(u6,*) 'no {i} found!'
    end if\n''')
