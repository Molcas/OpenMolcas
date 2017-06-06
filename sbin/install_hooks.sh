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
#                                                                      *
# Copyright (C) 2017, Ignacio Fdez. Galván                             *
#***********************************************************************

if [ -z "$1" ] ; then
  exit 0
else
  MOLCAS="$1"
fi

if [ ! -d "$MOLCAS" ] ; then
  exit 0
fi

copy_hook() {
  echo "Copying hook \"$1\" into \"$2\""
  cp "$1" "$2"
}

gitdir=$(cd $MOLCAS ; (cd $(git rev-parse --git-dir) ; pwd))
if [ -d "$gitdir" ] ; then
  mkdir -p "$gitdir/hooks"
fi
for hook in "pre-commit" ; do
  src="$MOLCAS/sbin/$hook"
  dest="$gitdir/hooks/$hook"
  if [ -f "$src" ] ; then
    if [ ! -f "$dest" ] ; then
      copy_hook "$src" "$dest"
    else
      custom=$(head "$dest" | grep 'DO NOT REPLACE' | wc -l)
      differ=$(diff "$src" "$dest" | wc -l)
      if [[ "$custom" == "0" && "$differ" != 0 ]]; then
        copy_hook "$src" "$dest"
      fi
    fi
  fi
done
