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

# current pymolcas checksum
if [ -z "$1" ] ; then
  pymolcas_file="Tools/pymolcas/pymolcas"
else
  pymolcas_file="$1"
fi
if [ -f $pymolcas_file ] ; then
  current_checksum=`md5sum $pymolcas_file | awk '{print $1}'`
else
  echo "*** $pymolcas_file does not exist"
  exit 1
fi

PYMOLCAS=""

# check location of pymolcas in your PATH
installed=0

orig_IFS=$IFS
IFS=':'
for x in $PATH ; do
  if [ -f "$x/pymolcas" ] ; then
    l=`md5sum "$x/pymolcas" | awk '{print $1}'`
    if [ "$l" != "$current_checksum" ] ; then
      echo "*** Warning! A different pymolcas version was found at: $x"
    else
      echo "A current pymolcas has been found installed at: $x"
      installed=1
    fi
    PYMOLCAS="$x"
    break
  fi
done

# current version already installed, exit
if [ $installed = 1 ] ; then
  exit 0
fi

# try to find proper location for pymolcas
if [ -z "$PYMOLCAS" ] ; then
  # no existing version found, use default PATH
  if [ -z "$PATH" ] ; then
    E_PATH="$HOME/bin"
  else
    E_PATH="$HOME/bin:$PATH"
  fi
else
  # a previous version was found, create it there
  E_PATH="$PYMOLCAS"
fi

# find first writable directory
dir_found=0
for x in $E_PATH ; do
  if [ "$x" = "." ] ; then continue ; fi
  # detect if directory is writable (-w would not work, since it does not handle mount permissions)
  cp "$pymolcas_file" "$x/this_is_not_pymolcas" > /dev/null 2>&1
  if [ "$?" = 0 ] ; then
    rm -f "$x/this_is_not_pymolcas" > /dev/null 2>&1
    PYMOLCAS="$x"
    dir_found=1
    break
  fi
done
IFS=$orig_IFS

# are we in interactive mode?
if [ "$CMAKE_SESSION" = "OpenMolcas" ] ; then
  INTERACTIVE=0
elif ( ! tty -s ) ; then
  INTERACTIVE=0
else
  INTERACTIVE=1
fi

# function to read with timeout in POSIX sh
read_timeout() {
  old=$(stty -g)
  stty -icanon min 0 time 255
  read $1
  stty $old
}

# create a default molcas driver
if [ $dir_found = 0 ] ; then
  echo "*** Warning! Could not find a proper directory to install pymolcas"
  echo ""
  echo "*** Check that there is a directory in your PATH with write access"
  echo "*** (for example $HOME/bin) and restart the installation"
  echo ""
  echo "*** You have to put pymolcas in any directory in your PATH"
  if [ "$INTERACTIVE" = "0" ] ; then
    exit 1
  fi
else
  echo "pymolcas will be installed in $PYMOLCAS"
  echo "Is this OK? [Y/n] (will assume \"Yes\" in 25 seconds)"
  while true ; do
    if [ "$INTERACTIVE" = "0" ] ; then
      echo "Running in non-interactive mode, assuming \"Yes\""
      answer="Yes"
    else
      read_timeout answer
    fi
    case "${answer}_" in
      [Yy]*|_ )
        cp "$pymolcas_file" "$PYMOLCAS/pymolcas"
        chmod +x "$PYMOLCAS/pymolcas"
        # check again the driver was installed
        l=`md5sum "$x/pymolcas" | awk '{print $1}'`
        if [ "$l" = "$current_checksum" -a -x "$PYMOLCAS/pymolcas" ] ; then
          echo "The installation of pymolcas was successful"
        else
          echo "*** The installation of pymolcas failed!"
          echo ""
          echo "*** You have to put pymolcas in any directory in your PATH"
          if [ "$INTERACTIVE" = "0" ] ; then
            exit 2
          fi
        fi
        break ;;
      [Nn]* )
        echo "*** The installation of pymolcas was canceled!"
        echo ""
        echo "*** You have to put pymolcas in any directory in your PATH"
        if [ "$INTERACTIVE" = "0" ] ; then
          exit 3
        fi
        break ;;
      * ) echo "Please answer yes or no"
    esac
  done
fi

exit 0
