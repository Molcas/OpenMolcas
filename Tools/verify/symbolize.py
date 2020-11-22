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
#***********************************************************************

# script to convert stack trace with pc addresses
# to source code files and line numbers

from __future__ import (unicode_literals, division, absolute_import, print_function)

import sys
import re
import subprocess

print('symbolized stack trace:')

for line in sys.stdin:
  match = re.match('(\s+#\d+)\s+0x[0-9a-f]+\s+\(([^\(\)]*)\)', line)
  if match:
    counter = match.group(1)
    address = match.group(2)
    exe, pc = address.split('+')
    source = subprocess.check_output(['addr2line', '-f', '-p', '-e', exe, pc])
    if re.search('\?\?:\?', source):
      print('{} (?) {} {}'.format(counter, exe, pc))
    else:
      print('{} {}'.format(counter, source))
  else:
    match = re.match('(\s+#\d+)\s+0x[0-9a-f]+\s+(in\s+.*)', line)
    if match:
      counter = match.group(1)
      address = match.group(2)
      print('{} {}'.format(counter, address))
