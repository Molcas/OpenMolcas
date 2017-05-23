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

# current driver version
if [ -z "$1" ] ; then
  driver_file="sbin/molcas.driver"
else
  driver_file="$1"
fi
if [ -f $driver_file ] ; then
  current_version=`sed -ne 's/\.// ; s/.*Molcas driver shell script version // p' $driver_file`
else
  echo "*** $driver_file does not exist"
  exit 1
fi

MOLCASDRIVER=""

# check location of driver script in your PATH
installed=0

orig_IFS=$IFS
IFS=':'
for x in $PATH; do
  if [ -f "$x/molcas" ]; then
    l=`sed -ne 's/\.// ; s/.*Molcas driver shell script version // p' "$x/molcas"`
    if [ "$l" -lt "$current_version" ]; then
      echo "*** Warning! An old version of the molcas driver script was found at: $x"
    else
      echo "Molcas driver has been found at: $x"
      installed=1
    fi
    MOLCASDRIVER="$x"
    break
  fi
done

# latest version already installed, exit
if [ $installed = 1 ] ; then
  exit 0
fi

# try to find proper location for molcas driver
if [ -z "$MOLCASDRIVER" ]; then
  # no existing version found, use default PATH
  if [ -z "$PATH" ]; then
    E_PATH="$HOME/bin"
  else
    E_PATH="$HOME/bin:$PATH"
  fi
else
  # a previous version was found, create it there
  E_PATH="$MOLCASDRIVER"
fi

# find first writable directory
dir_found=0
for x in $E_PATH ; do
  if [ "$x" = "." ] ; then continue ; fi
  # detect if directory is writable (-w would not work, since it does not handle mount permissions)
  cp "$driver_file" "$x/this_is_not_molcas" > /dev/null 2>&1
  if [ "$?" = 0 ] ; then
    rm -f "$x/this_is_not_molcas" > /dev/null 2>&1
    MOLCASDRIVER="$x"
    dir_found=1
    break
  fi
done

# are we in interactive mode?
if [ "$CMAKE_SESSION" = "OpenMolcas" ] ; then
  INTERACTIVE=0
elif ( ! tty -s ) ; then
  INTERACTIVE=0
else
  INTERACTIVE=1
fi

# create a default molcas driver
if [ $dir_found = 0 ] ; then
  echo "*** Warning! Could not find a proper directory to install the molcas driver"
  echo ""
  echo "*** Check that there is a directory in your PATH with write access"
  echo "*** (for example $HOME/bin) and restart the installation"
  echo ""
  echo "*** You have to put the molcas driver in any directory in your PATH"
  if [ "$INTERACTIVE" = "0" ] ; then
    exit 1
  fi
else
  echo "molcas driver will be installed in $MOLCASDRIVER"
  echo "Is this OK? [Y/n]"
  while true; do
    if [ "$INTERACTIVE" = "0" ] ; then
      echo "Running in non-interactive mode, assuming \"Yes\""
      answer="Yes"
    else
      read answer
    fi
    case "${answer}_" in
      [Yy]*|_ )
        cp "$driver_file" "$MOLCASDRIVER/molcas"
        chmod +x "$MOLCASDRIVER/molcas"
        # check again the driver was installed
        l=`sed -ne 's/\.// ; s/.*Molcas driver shell script version // p' "$x/molcas"`
        if [ "$l" -eq "$current_version" -a -x "$MOLCASDRIVER/molcas" ] ; then
          echo "Driver installation successful"
        else
          echo "*** Driver installation failed!"
          echo ""
          echo "*** You have to put the molcas driver in any directory in your PATH"
          if [ "$INTERACTIVE" = "0" ] ; then
            exit 2
          fi
        fi
        break ;;
      [Nn]* )
        echo "*** Driver installation canceled!"
        echo ""
        echo "*** You have to put the molcas driver in any directory in your PATH"
        if [ "$INTERACTIVE" = "0" ] ; then
          exit 3
        fi
        break ;;
      * ) echo "Please answer yes or no"
    esac
  done
fi

exit 0
