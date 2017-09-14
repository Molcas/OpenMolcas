#!/bin/bash

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

# Preprocess file
cat $1 | \
  # Fix case of automatic links
  ./fix_links.py | \
  # Disable #include directives (using "C" to avoid mis-detection as free format)
  sed -e 's/^#include/Cinclude/' -e 's/^[*cC!]/C/' | \
  # Preprocess with cpp
  #cpp -undef -w -fdirectives-only -dU -P -nostdinc
  # Preprocess with gpp
  gpp -n -U "" "" "(" "," ")" "(" ")" "#" "" -M "\n#\w" "\n" " " " " "\n" "" ""
