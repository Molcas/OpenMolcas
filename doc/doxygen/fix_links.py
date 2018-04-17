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

from __future__ import print_function
import fileinput
import re

doxygen_re = re.compile(r"^[*Cc!]>")
link_re = re.compile(r"::(\w*)")

# Replace [::FuncName] with [\ref funcname "FuncName"]
reflink = lambda pat: r'\ref {} "{}"'.format(pat.group(1).lower(), pat.group(1))

for line in fileinput.input():
  if (doxygen_re.search(line)):
    line = re.sub(link_re,reflink,line)
    print(line,end='')
  else:
    print(line,end='')
