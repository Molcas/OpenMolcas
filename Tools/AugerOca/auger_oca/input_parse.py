#!/usr/bin/env python3
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
# Copyright (C) 2020-2022, Bruno Tenorio                               *
#***********************************************************************

import argparse,time,subprocess, re, sys, os
from argparse import RawTextHelpFormatter
# input_parse.py
def parseCL():
    d = 'This program will process the 2p-Dyson density files (r2TM_) and provide: \
\nAuger decay rates via one-center approach (OCA).\
\n[Remember to make sure you have the $Project.rassi.h5 file in your current directory.]\
\n \
\n \
usage example (RAES):\
\n \
\n \
$python3 $Tools/auger_oca/auger_main.py --raes -i r2TM_K2V_002_001\
\n \
\n \
usage example (AES for a triplet final state):\
\n \
\n \
$python3 $Tools/auger_oca/auger_main.py --aes --t -i r2TM_K2V_002_001\
\n \
\n \
(recommended) usage based on a directory with a collection of r2TM_ files.\
\n \
This directory can be the scratch directory for the RASSI calculation.\
\n \
\n \
$python3 $MOLCAS/Tools/auger_oca/auger_main.py --raes --spec -d $WorkDir'
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--input",
                        dest="inp",
                        required=False,
                        type=str,
                        help="path of the density file (r2TM_)")
    parser.add_argument("-d", "--directory",
                        dest="directory",
                        required=False,
                        type=str,
                        help="path of a directory with a number of density files (r2TM_)")
    parser.add_argument("--oca",action='store_true',
                        dest="parse_oca",
                        required=False,
                        default=True,
                        help="Perform Auger OCA (default: True)")
    parser.add_argument("--spec",action='store_true',
                        dest="parse_spec",
                        required=False,
                        default=False,
                        help="save auger.spectrum.out (default: False)")
    parser.add_argument("--raes",action='store_true',
                        dest="parse_raes",
                        required=False,
                        default=True,
                        help="Consider resonant Auger [RAES] (default: True if OCA)")
    parser.add_argument("--aes",action='store_true',
                        dest="parse_aes",
                        required=False,
                        default=False,
                        help="Consider non-resonant Auger [AES] (default: False)")
    parser.add_argument("--t",action='store_true',
                        dest="parse_aes_triplet",
                        required=False,
                        default=False,
                        help="For AES, consider triplet final states (default: False)")
    parser.add_argument("--s",action='store_true',
                        dest="parse_aes_singlet",
                        required=False,
                        default=False,
                        help="For AES, consider singlet final states (default: False; True if AES)")
    args = parser.parse_args()
    if args.parse_aes==True:
        args.parse_raes=False

    if args.parse_aes_triplet==args.parse_aes_singlet and args.parse_aes_triplet==True:
        print('Not possible AES for singlet and triple final ion simultaneously. Exit()')
        sys.exit()
    elif args.parse_aes_triplet==True:
        args.parse_aes=True
        args.parse_raes=False
        args.parse_oca=True
    elif args.parse_aes_singlet==True:
        args.parse_aes=True
        args.parse_raes=False
        args.parse_oca=True

    if args.inp and args.directory:
        print('options -i and -d are mutually exclusive. Exit()')
        sys.exit()

    return args

#args = parseCL()
#print(args)
