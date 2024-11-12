#! /usr/bin/env python
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
# Copyright (C) 2024, Ignacio Fdez. Galván                             *
#***********************************************************************
#
# Parse a set of verify result files to check that every test (in groups 1 and 2)
# has been successfully run (not skipped) in at least one of the jobs.
# It also checks if any test in groups 3 and 4 could be promoted.

import glob
import re
import sys

results = {}
version = None
for fn in glob.glob('result_*'):
  with open(fn, 'r') as f:
    v = f.readline().strip()
    if version is None:
      version = v
    elif v != version:
      print('Versions do not match!')
      sys.exit(1)
    for l in f:
      m = re.match(r'(\d:[^ :]*:[^ :]*) (.*)$', l)
      if m:
        test = m.groups()[0]
        res = m.groups()[1]
        if test not in results:
          results[test] = res
        elif not results[test].startswith('Failed'):
          if res.startswith('Failed'):
            results[test] = res
          elif res == 'OK':
            results[test] = res

skipped = []
passed = []
for t in sorted(results.keys()):
  print(t)
  if re.match(r'[12]:', t):
    if results[t] != 'OK':
      skipped.append(t)
  if re.match(r'[34]:', t):
    if results[t] == 'OK':
      passed.append(t)
print()

if passed:
  print('The following tests did not fail in any job and passed in some:')
  for t in passed:
    print(t)

if skipped:
  print('The following tests did not pass in any job:')
  for t in skipped:
    print(t)
  sys.exit(1)
