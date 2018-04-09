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
# Copyright (C) 2015,2018, Ignacio Fdez. Galván                        *
#***********************************************************************

'''
This packs pymolcas into a single-file executable python package that can be
copied around. It still uses python and it's trivial to recover the original
files, comments included. Obfuscation is not the goal here, just getting
something that's easy to run, move and distribute.
'''

from __future__ import (unicode_literals, division, absolute_import, print_function)
from builtins import (bytes)

import glob, os.path

from setuptools import setup

files = [os.path.splitext(x)[0] for x in glob.glob('*.py')]
files.remove('export')
files.remove('pack')

setup(
  script_args=['bdist_wheel', '-d', '.'],
  py_modules=files,
)

# Now try to make the package runnable:
# Add the self-running trick code to the top of the file,
# rename it and make it executable

import sys, os, stat

exe_name = 'pymolcas.whl'
magic_code = '''#!/bin/sh
name=`readlink -f "$0"`
exec {0} -c "import sys, os; sys.path.insert(0, '$name'); from pymolcas import main; sys.exit(main(my_name='$name'))" "$@"
'''.format(sys.executable)

wheel = glob.glob('*.whl')[0]
with open(exe_name, 'wb') as new:
  new.write(bytes(magic_code, 'ascii'))
  with open(wheel, 'rb') as original:
    data = True
    while (data):
      data = original.read(4096)
      new.write(data)
os.remove(wheel)
st = os.stat(exe_name)
os.chmod(exe_name, st.st_mode | stat.S_IEXEC)
