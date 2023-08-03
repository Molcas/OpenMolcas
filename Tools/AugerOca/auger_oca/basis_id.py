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
import sys, h5py

## orbtransf.py
## Collection of MO to AO and 2molden programs
# basis_func_id,phase_lm,tomolden,mo2ao

# basis_func_id is a function to read a list containing (c,n,l,m)
# with c the atomic center, n the principal quantum number,
# l the angular quantum number, and m the projection of l 
# usefull convention is that l=0,1,-1 is z,x,y respectively.
# d0,d+1,d-1,d+2,d-2 -> z2,xz,yz,xy,x2-y2
#
# Works only with full spherical basis

def basis_func_id(cnlm,element):
    c=cnlm[0]
    n=cnlm[1]
    l=cnlm[2]
    m=cnlm[3]
    list_element=element
#
    if l==0:
        func_id='s'
    elif l==1:
        if m==0:
            func_id='pz'
        if m==1:
            func_id='px'
        if m==-1:
            func_id='py'
    elif l==2:
        if m==0:
            func_id='d0'
        if m==1:
            func_id='d+1'
        if m==-1:
            func_id='d-1'
        if m==2:
            func_id='d+2'
        if m==-2:
            func_id='d-2'
    elif l==3:
        if m==0:
            func_id='f0'
        if m==1:
            func_id='f+1'
        if m==-1:
            func_id='f-1'
        if m==2:
            func_id='f+2'
        if m==-2:
            func_id='f-2'
        if m==3:
            func_id='f+3'
        if m==-3:
            func_id='f-3'
    elif l==4:
        if m==0:
            func_id='g0'
        if m==1:
            func_id='g+1'
        if m==-1:
            func_id='g-1'
        if m==2:
            func_id='g+2'
        if m==-2:
            func_id='g-2'
        if m==3:
            func_id='g+3'
        if m==-3:
            func_id='g-3'
        if m==4:
            func_id='g+4'
        if m==-4:
            func_id='g-4'
    else:
        func_id='UKN'
    return '{}{}{}{}'.format(list_element[c-1]," ",n+l, func_id)
# The following order of D, F en G functions is expected:
# 5D: D 0, D+1, D-1, D+2, D-2
# 7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
# 9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4

#==================================

def get_overlap(hd5_file,basis_id_hd5,nbasft,nbasf,symmetry,element):
    from .symm import nbasis
# Get here overlap matrix. Only once.
    h = h5py.File(hd5_file, 'r') # here I open the 'hd5_file' with h5py
    over_length=0
    for x in range(symmetry):
        over_length=over_length+((nbasis(x+1,nbasf))**2)
    ao_overlap=np.reshape( h['AO_OVERLAP_MATRIX'],(over_length))
    h.close()

    # now I make a list len(nbasft) with all basis ids
    basis_id_list=list()
    for i in range(nbasft):
        basis_id_list.append(basis_func_id(basis_id_hd5[i],element))
        #print(ao.basis_func_id(basis_id_hd5[i]))
    basis_id_list_shell=basis_id_list # I need this to define the MBS
    #print(basis_id_list_shell)
    return basis_id_list,ao_overlap
