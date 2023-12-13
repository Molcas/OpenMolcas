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

## symm.py
# Collection of functions to deal with symmetry point group

# identify the symmetry of the dipole operator
def symop(n,symmetry):
# n=1,2,3->x,y,z
    if symmetry==4:
# (a1,b1,a2,b2)->(1,2,3,4)
# b1->x, b2->y, a1->z
        if n==1:
            return 2
        if n==2:
            return 4
        if n==3:
            return 1
    elif symmetry==8:
# (Ag,B2g,B1g,B3g,Au,B2u,B1u,B3u)->(1,2,3,4,5,6,7,8)
# B3u->x, B2u->y, B1u->z
        if n==1:
            return 8
        if n==2:
            return 6
        if n==3:
            return 7
    elif symmetry==2:
# A'->x,y, A"->z
        if n==1:
            return 1
        if n==2:
            return 1
        if n==3:
            return 2
    elif symmetry==1:
        return 1
    else:
        return 'incorrect symmetry group'

def symoprep(n,symmetry):
# given the symmetry number of the operator
# returns the component as x,y z
    if symmetry==4:
        if n==1:
            return 'z'
        if n==2:
            return 'x'
        if n==3:
            return 'no dipole'
        if n==4:
            return 'y'
    elif symmetry==8:
        if n==6:
            return 'y'
        if n==7:
            return 'z'
        if n==8:
            return 'x'
        else:
            return 'no dipole'
    if symmetry==2:
        if n==1:
            return 'x'
        if n==2:
            return 'z'
    elif symmetry==1:
        return 1
    else:
        return 'incorrect symmetry group'

# multiplication table. for now only C2v, D2h, Cs and c1
def mul(i,j,symmetry):
# (a1,b1,a2,b2)->(1,2,3,4)
    if symmetry==4:
        if i and j <= 4:
            if i ==1:
                return i*j
            if i==2:
                if j==1:
                    return 2
                if j==2:
                    return 1
                if j==3:
                    return 4
                if j==4:
                    return 3
            if i ==3:
                if j==1:
                    return 3
                if j==2:
                    return 4
                if j==3:
                    return 1
                if j==4:
                    return 2
            if i==4:
                if j==1:
                    return 4
                if j==2:
                    return 3
                if j==3:
                    return 2
                if j==4:
                    return 1
# (A',A")->(1,2)
    if symmetry==2:
        if i and j <= 2:
            if i ==1:
                return i*j
            if i==2:
                if j==1:
                    return 2
                if j==2:
                    return 1
    elif symmetry==8:
# (Ag,B2g,B1g,B3g,Au,B2u,B1u,B3u)->(1,2,3,4,5,6,7,8)
        if i and j <= 8:
            if i ==1:
                return i*j
            if i==2:
                if j==1:
                    return 2
                if j==2:
                    return 1
                if j==3:
                    return 4
                if j==4:
                    return 3
                if j==5:
                    return 6
                if j==6:
                    return 5
                if j==7:
                    return 8
                if j==8:
                    return 7
            if i ==3:
                if j==1:
                    return 3
                if j==2:
                    return 4
                if j==3:
                    return 1
                if j==4:
                    return 2
                if j==5:
                    return 7
                if j==6:
                    return 8
                if j==7:
                    return 5
                if j==8:
                    return 6
            if i==4:
                if j==1:
                    return 4
                if j==2:
                    return 3
                if j==3:
                    return 2
                if j==4:
                    return 1
                if j==5:
                    return 8
                if j==6:
                    return 7
                if j==7:
                    return 6
                if j==8:
                    return 5
            if i==5:
                if j==1:
                    return 5
                if j==2:
                    return 6
                if j==3:
                    return 7
                if j==4:
                    return 8
                if j==5:
                    return 1
                if j==6:
                    return 2
                if j==7:
                    return 3
                if j==8:
                    return 4
            if i==6:
                if j==1:
                    return 6
                if j==2:
                    return 5
                if j==3:
                    return 8
                if j==4:
                    return 7
                if j==5:
                    return 2
                if j==6:
                    return 1
                if j==7:
                    return 4
                if j==8:
                    return 3
            if i==7:
                if j==1:
                    return 7
                if j==2:
                    return 8
                if j==3:
                    return 5
                if j==4:
                    return 6
                if j==5:
                    return 3
                if j==6:
                    return 4
                if j==7:
                    return 1
                if j==8:
                    return 2
            if i==8:
                if j==1:
                    return 8
                if j==2:
                    return 7
                if j==3:
                    return 6
                if j==4:
                    return 5
                if j==5:
                    return 4
                if j==6:
                    return 3
                if j==7:
                    return 2
                if j==8:
                    return 1
    elif symmetry==1:
        return 1
    else:
        return 'incorrect symmetry group'

#norb catch the number of orbitals in the i symmetry
def norb(i,nmo):
    #returns the number of orbitals in the i symm
    return int(nmo[i-1])

#nbasis catch the number of basis functions in the i symmetry
def nbasis(i,nbasf):
    #returns the number of basis function in the i symm
    return int(nbasf[i-1])

#catch cmo takes the corresponding symmetry CMO for multiplication
def catchcmo(i,nmo,nbasf):
    #x=Num orb in symm
    x=int(nmo[i-1])
    #y=Num basis in symm
    y=int(nbasf[i-1])
#returns a product #orb*#basis and a tuple (#Orb,#Basis)
    return (x*y,(y,x))

