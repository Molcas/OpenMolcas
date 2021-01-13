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
# Copyright (C) 2014, Steven Vancoillie                                *
#               2020, Ignacio Fdez. Galván                             *
#***********************************************************************
#
# modulenames.py:
#
# extract module names from a test input and generate a proper test name for
# use with CMake (no spaces, special chars)
#
# Steven Vancoillie, Lund, 20th of June 2014
# Ignacio Fdez. Galván, 2020: Translated from Perl to Python

import sys
import os
import re

def usage():
  name = os.path.basename(__file__)
  print('''

 Usage: {} <filename>

        Program to generate a test name by concatenating the
        module names used in a test input given by filename.
'''.format(name))
  sys.exit(1)

try:
  filename = sys.argv[1]
except:
  usage()

module = re.compile('&(\w+)')

module_names = []

try:
  with open(filename, 'r') as f:
    for line in f:
      matches = module.findall(line)
      for name in matches:
        lc_name = name.lower()
        if lc_name != 'end':
          module_names.append(lc_name)
except IOError:
  sys.exit('Error: failed to open file {}'.format(filename))

output = '_'.join(sorted(set(module_names)))

print(output)
