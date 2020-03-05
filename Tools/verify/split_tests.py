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
# Copyright (C) 2019, Ignacio Fdez. Galván                             *
#***********************************************************************
#
# Split a list of tests in sets taking roughly the same running time
#
# ./split_tests.py flatlist total_number selected_number

from __future__ import (division, print_function)
import sys, os, re

# default running time for tests with no estimate
default_time = 60

# this file location
path = os.path.dirname(os.path.abspath(__file__))

# flatlist
try:
  flatlist = sys.argv[1]
except IndexError:
  sys.exit()

# total sets
try:
  numsets = int(sys.argv[2])
  assert numsets > 0
except:
  numsets = 1

# this set
try:
  thisnum = int(sys.argv[3])
  assert 0 < thisnum <= numsets
except:
  thisnum = 1

# read running time estimates
timest = {}
try:
  with open(os.path.join(path, 'timing.data'), 'r') as f:
    while True:
      try:
        line = next(f)
        match = re.match(r'--- (\S*)( .)?', line)
        if match:
          test = match.group(1)
          timest[test] = float(next(f))
          if match.group(2):
            timest[test] = max(timest[test], default_time)
      except StopIteration:
        break
except IOError:
  pass

# read test list
testlist = []
with open(flatlist, 'r') as f:
  for l in f:
    testlist.append(l.strip())

# compute total time
total = 0
for i in testlist:
  total += timest.get(i, default_time)

# split in sets
set_time = total / numsets
sets = []
this_set = []
this_time = 0
for i in testlist:
  this_set.append(i)
  this_time += timest.get(i, default_time)
  if this_time > set_time:
    sets.append(this_set[:])
    this_set = []
    this_time = 0
if (len(sets) < numsets):
  sets.append(this_set[:])

# print selected set
try:
  for i in sets[thisnum-1]:
    print(i)
except IndexError:
  pass
