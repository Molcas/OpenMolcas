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
# Copyright (C) 2020-2022, Bruno Tenorio                               *
#***********************************************************************
# Get the spectrum from the processed Auger_OCA.* files. 
# inputs: OCA outputs in a folder
# output: spectrum.out file

# How to use:

#$ bash oca.spectrum.sh auger_outputs/ &

# where 'auger_outputs' is the directory containing the OCA output files Auger_OCA.*
# alternatively, simply use "--spec" option >> auger_main.py --spec
ORIG_DIR=$PWD
OUTDIR=$1

# Collect BE & Intensities into -> spectrum.out
for f in $(ls $OUTDIR/Auger_OCA.*); do grep "Spectrum: BE(eV) and Intensity" $f >> tmp.spectrum.out; done
cat tmp.spectrum.out | cut -b 1-42 --complement >> spectrum.out
rm tmp.spectrum.out

