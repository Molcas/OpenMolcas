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

import numpy as np
import os

# oca_integrals.py

def oca_integrals(OCA_atom):
    # The file OCA.dat has 21 lines of integrals for each element. 
    # The file OCA.3r.dat has 68 lines of integrals for each element.
    # Available elements (currently) are: C, O, N, Ne, S, CL
    # Reads 5 columns: I,J,L,M,G
    #  I,J: 1   2    3    4    5    6    7    8
    #       2s  2pz  2px  2py  3s   3pz  3px  3py
    oca=list()
    nele_oca=84
    oca_dir = str(os.path.dirname(os.path.abspath(__file__)))+'/OCA.dat'  # give the full real path of OCA.dat
    with open(oca_dir,"r") as oc :
        for curline in oc:
            if curline.startswith("#"):
                pass
            else:
                line3=[float(elem) for elem in curline.split()]
                oca.append(line3)
    temp_oca = [x for x in oca if x != []] # 
    oca=np.array(temp_oca)
    #Carbon set
    oca_C=oca[:nele_oca]
    #Oxygen set
    oca_O=oca[nele_oca:nele_oca*2]
    #N set
    oca_N=oca[nele_oca*2:nele_oca*3]
    #F set
    oca_F=oca[nele_oca*3:nele_oca*4]
    #Ne set
    oca_NE=oca[nele_oca*4:nele_oca*5]
    #S set
    oca_S=oca[nele_oca*5:nele_oca*6]
    #Cl set
    oca_Cl=oca[nele_oca*6:nele_oca*7]
    #Mg set
    oca_Mg=oca[nele_oca*7:nele_oca*8]
    #P set
    oca_P=oca[nele_oca*8:nele_oca*9]
    #Ar set
    oca_Ar=oca[nele_oca*9:nele_oca*10]
    #Al set
    oca_Al=oca[nele_oca*10:nele_oca*11]

    # Define One Center Integral according to the OCA_atom variable
    if OCA_atom=='C':
        OCI=oca_C
    elif OCA_atom=='O':
        OCI=oca_O
    elif OCA_atom=='N':
        OCI=oca_N
    elif OCA_atom=='F':
        OCI=oca_F
    elif OCA_atom=='NE':
        OCI=oca_NE
    elif OCA_atom=='S':
        OCI=oca_S
    elif OCA_atom=='CL':
        OCI=oca_Cl
    elif OCA_atom=='MG':
        OCI=oca_Mg
    elif OCA_atom=='P':
        OCI=oca_P
    elif OCA_atom=='AR':
        OCI=oca_Ar
    elif OCA_atom=='AL':
        OCI=oca_Al

    return OCI

def elmij(OCA_atom,OCA_c,c,i,j,l,m):
    OCI = oca_integrals(OCA_atom)
    ver = 0.0
    third = ['S','P','AR','CL','MG','AL']
    if c ==OCA_c:
        if any(c == y for y in third):
        #if c == 'S' or c == 'CL':
            for ez in OCI:
                if i == '2s':
                    ii=1
                elif i == '2pz':
                    ii=2
                elif i == '2px':
                    ii=3
                elif i == '2py':
                    ii=4
                elif i == '3s':
                    ii=5
                elif i == '3pz':
                    ii=6
                elif i == '3px':
                    ii=7
                elif i == '3py':
                    ii=8
                else:
                    break
                if j == '2s':
                    jj=1
                elif j == '2pz':
                    jj=2
                elif j == '2px':
                    jj=3
                elif j == '2py':
                    jj=4
                elif j == '3s':
                    jj=5
                elif j == '3pz':
                    jj=6
                elif j == '3px':
                    jj=7
                elif j == '3py':
                    jj=8
                else:
                    break
                if all([ii,jj,l,m] == ez[:4]) :
                    ver = ez[4]
        else:
        #
            for ez in OCI:
                if i == '2s':
                    ii=1
                elif i == '2pz':
                    ii=2
                elif i == '2px':
                    ii=3
                elif i == '2py':
                    ii=4
                else:
                    break
                if j == '2s':
                    jj=1
                elif j == '2pz':
                    jj=2
                elif j == '2px':
                    jj=3
                elif j == '2py':
                    jj=4
                else:
                    break

                if all([ii,jj,l,m] == ez[:4]) :
                    ver = ez[4]
    return ver
    # print('Eml',elmij('C 1s', 'C 1s','2py','2px',2,-2)) # this test should give: Eml 0.006919730033
# ------------------------------------------
