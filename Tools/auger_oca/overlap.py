#!/usr/bin/env python3
#***********************************************************************
## This file is part of OpenMolcas.                                     *
##                                                                      *
## OpenMolcas is free software; you can redistribute it and/or modify   *
## it under the terms of the GNU Lesser General Public License, v. 2.1. *
## OpenMolcas is distributed in the hope that it will be useful, but it *
## is provided "as is" and without any express or implied warranties.   *
## For more details see the full text of the license in the file        *
## LICENSE or in <http://www.gnu.org/licenses/>.                        *
##                                                                      *
## Copyright (C) 2020, Bruno Tenorio                                    *
##***********************************************************************
import sys, h5py
import numpy as np
import symm as sym
import orbtransf as ao

# overlap.py

def get_overlap(hd5_file,basis_id_hd5,nbasft,nbasf,symmetry,element):
# Get here overlap matrix. Only once.
    ao_dip=list()
    h = h5py.File(hd5_file, 'r') # here I open the 'hd5_file' with h5py
    ao_dipx=np.reshape( h['AO_MLTPL_X'],(nbasft,nbasft))
    ao_dipy=np.reshape( h['AO_MLTPL_Y'],(nbasft,nbasft))
    ao_dipz=np.reshape( h['AO_MLTPL_Z'],(nbasft,nbasft))
    ao_dip=(ao_dipx,ao_dipy,ao_dipz)
    over_length=0
    for x in range(symmetry):
        over_length=over_length+((sym.nbasis(x+1,nbasf))**2)
    ao_overlap=np.reshape( h['AO_OVERLAP_MATRIX'],(over_length))
    h.close()

    # now I make a list len(nbasft) with all basis ids
    basis_id_list=list()
    for i in range(nbasft):
        basis_id_list.append(ao.basis_func_id(basis_id_hd5[i],element))
        #print(ao.basis_func_id(basis_id_hd5[i]))
    basis_id_list_shell=basis_id_list # I need this to define the MBS
    #print(basis_id_list_shell)
    return basis_id_list,ao_dip,ao_overlap
