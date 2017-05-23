#!/bin/sh
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

# This compiles pymolcas into an executable file that can be copied around. It
# still needs a python installation to run, but it's not possible to recover the
# original files (as far as I know). Thus, this may be a way to obfuscate the
# code into something that's not so easily modifiable.

nuitka --python-version=3.4 --remove-output \
       --recurse-to=abstract_flow \
       --recurse-to=emil_grammar \
       --recurse-to=emil_parse \
       --recurse-to=molcas_aux \
       --recurse-to=molcas_wrapper \
       --recurse-to=python_parse \
       --recurse-to=tee \
       pymolcas.py

rm -rf __pycache__
