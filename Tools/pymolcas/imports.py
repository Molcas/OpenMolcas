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
# Copyright (C) 2018, Ignacio Fdez. Galván                             *
#***********************************************************************

# These modules should be available in any python installation
standard_modules = [
  '__future__',
  'argparse',
  'ast',
  'base64',
  'binascii',
  'builtins|future.builtins',
  'codecs',
  'contextlib',
  'datetime',
  'decimal',
  'errno',
  'glob',
  'hashlib',
  'io',
  'json',
  'operator',
  'os',
  'os.path',
  'random',
  're',
  'resource',
  'shlex',
  'shutil',
  'signal',
  'stat',
  'subprocess',
  'sys',
  'tempfile',
  'textwrap',
  'threading',
  'types',
  'zlib',
]

# These modules may have to be installed
modules = [
  'pyparsing',
  #'setuptools', (only used in pack.py)
  'six',
]

fail = []

try:
  import importlib
except ImportError:
  fail.append('importlib')
else:
  for item in standard_modules + modules:
    mods = item.split('|')
    result = []
    for m in mods:
      try:
        importlib.import_module(m)
        result.append(True)
      except ImportError:
        result.append(False)
    if (not any(result)):
      fail.append(item)

import sys

if (fail):
  print(' '.join(fail))

sys.exit(len(fail))
