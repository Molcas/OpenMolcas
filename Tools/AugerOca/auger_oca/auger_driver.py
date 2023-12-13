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

import sys
import numpy as np
from numpy import linalg as LA

# auger_driver.py
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def driver_auger(input_file,RAES,NAES,NAES_T,NAES_S,\
    OCA_atom,OCA_line,benergy,totalSymmetry,symmetry,nbasf,cmotab,nbasft,comtaboff,comtbasoff,comtbasoff2,\
    comtcmoff,cmob,cmoa,nmo,tdmab,dyson,hd5_file,element,basis_id_hd5):

    from .basis_id import get_overlap
    from .orbtransf_sz import mbs
    from .mbs import basis_list_for_oca,init_sz
    from .rt2mzz import rt2mzz
    
    # Get here overlap matrix. Only once.
    orig_stdout = sys.stdout
    twopi=6.283185307179586 # two times Pi.
    au2ev=27.211396132 # To convert to eV.
    basis_id_list,ao_overlap = get_overlap(hd5_file,basis_id_hd5,nbasft,nbasf,symmetry,element)
    #==================================
    
    atombasis_list,shell_basis_sz,shell_py = basis_list_for_oca(basis_id_hd5,nbasft,element,comtbasoff,nbasf,nmo)

    sz_nbasf,sznbasft,ncsz,nbasz,nszorb,nbaszoff,nszoff,ncszoff,norbsztaboff = init_sz(symmetry,nbasf,shell_py,nmo)

    # mbs() brings the projected MOs onto the MBS
    d_sza,d_szb = mbs(ao_overlap,shell_py,symmetry,cmoa,cmob,nbasf,nmo,comtaboff,cmotab,comtbasoff2,ncszoff,ncsz,sz_nbasf,nbaszoff,nbasz)        
    #
    output_file = 'Auger_OCA.'+input_file+'.out' # also define here the output's name
    #open a output file and print all there.
    g = open(output_file, 'w')
    sys.stdout = g
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%           If OCA             %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OCA_c = OCA_line[0]
    # ELM (L,M) contains 9 elements: { [0,0]; [1,-1]; [1,0]; [1,1]; [2,-2]; [2,-1]; [2,0]; [2,1]; [2,2] }
    ELM=list(0 for ij in range(9))
    print('#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('#%         MBS  projection          %')
    print('#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('#The MOs are being projected into the AO basis')
    [print(basis_id_list[ei]) for ei in shell_basis_sz]
    print('#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('# One Center Approximation')
    print('# ')
    print('# Your job looks like:')
    print('# input file:',input_file,',','RAES:',RAES,',',\
    'AES:',NAES,',','AES triplet:',NAES_T,',','AES singlet:',NAES_S)
    print('#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    # printing decay rate from list ELM
    decay_spectrum=[]
    decay_spectrum.append(benergy) # first element is the binding energy in eV
    for c in OCA_line: # loop over scattering centers
        OCA_c = c
        print(' ')
        print('# The scattering center is :', OCA_c )
        ELM,tdmzz_symm,totfrob = rt2mzz(OCA_c,OCA_atom,cmoa,cmob,d_sza,d_szb,\
        tdmab,totalSymmetry,symmetry,nmo,nbasf,nbasft,comtbasoff,cmotab,comtaboff,basis_id_hd5,element)
    
        ELM2=[elimm**2 for elimm in ELM]
        mcguire=2.0000000 # To compensate the continuum normalization in Phys. Rev. 185, 1
        if RAES==True:
            prefactor=2.000000
        else: # then NAES True
            if NAES_S==True:
                prefactor=1.000000
            elif NAES_T==True:
                prefactor=1.500000
        decay_spectrum.append( mcguire*prefactor*twopi*(sum(ELM2)) )
        print('# Auger Width One center approx','(',OCA_c,')' ,'= ', prefactor*twopi*(sum(ELM2)) )
    print(' ')
    print('# Spectrum: BE(eV) and Intensity from OCA:', *decay_spectrum ) #asterisk before list's name with print statement prints only the elements
    print('# The End ')
    sys.stdout = orig_stdout
    g.close()
    #------------------------------------------------
    return
