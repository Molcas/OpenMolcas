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
#                                                                      *
# Copyright (C) 2020, Ignacio Fdez. Galván                             *
#***********************************************************************

# Script to ...
#

if [ -z "$MOLCAS" ] ; then
  MOLCAS=$PWD
fi

export MOLCAS

DEST=$1
if [ -d "$DEST" ] && [ -w "$DEST" ] ; then
  true
else
  echo "$DEST is not writable directory"
  exit 1
fi

for dir in `$MOLCAS/sbin/verify --group | awk -F "= " '/^Physical/,/^Special/ {print $2}'` ; do
  echo "Copying tests from $dir"
  cp -a "$dir" "$DEST"
done
