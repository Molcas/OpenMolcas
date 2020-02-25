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
# Copyright (C) 2020, Ignacio Fdez. Galván                             *
#***********************************************************************
#
# Parse an output file to get/update the list of used keywords.
# Intended for use with the verify script, to get a list of tested
# keywords, e.g.:
#
# pymolcas verify --postproc "\$MOLCAS_SOURCE/Tools/verify/get_keywords.py \$project.out keylist"

from __future__ import (division, print_function)
import sys, os, re, json

# output file
try:
  outfile = sys.argv[1]
except IndexError:
  sys.exit()

# keywords file
try:
  keyfile = sys.argv[2]
except IndexError:
  sys.exit()

used = {}
try:
  with open(keyfile, 'r') as f:
    used.update(json.load(f))
except IOError:
  pass

val = False
with open(outfile, 'r') as f:
  
  for line in f:
    if line.startswith('-- Input validation'):
      val = True
    elif line.strip() == '--':
      val = False
    if val:
      match = re.match(r'Input for: (\w*)', line)
      if match:
        module = match.group(1)
        if module not in used:
          used[module] = []
      match = re.match(r'\s*> entering group (.*)', line)
      if match:
        group = match.group(1)
        if group not in used[module]:
          used[module].append(group)
      match = re.match(r'\s*keyword (.*)', line)
      if match:
        keyword = match.group(1)
        if keyword not in used[module]:
          used[module].append(keyword)

with open(keyfile, 'w') as f:
  json.dump(used, f, sort_keys=True, indent=2)
