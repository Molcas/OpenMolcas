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
from numpy import linalg as LA
#==================================
## orbtransf_sz
# orbital projection into atomic "single-zeta" basis set.

# ----------------------------------------------------
def mbs(ao_overlap,shell_py,symmetry,cmoa,cmob,nbasf,nmo,comtaboff,cmotab,comtbasoff2,ncszoff,ncsz,sz_nbasf,nbaszoff,nbasz):
    # mbs() brings the projected MOs onto the MBS 
    from .symm import norb,nbasis,catchcmo
    sz_ovl=list()
    nbsfoff2=0
    for sim in range(symmetry):
        #print(shell_py[sim],nbsfoff2)
        for i in shell_py[sim]:
            for j in shell_py[sim]:
                if i>0 and j>0:
                    #print(i,j)
                    ij_sz=i + (j-1)*nbasf[sim] + nbsfoff2 -1
                    #print(ij_sz,ao_overlap[ij_sz])
                    sz_ovl.append(ao_overlap[ij_sz])
        nbsfoff2=nbsfoff2+nbasf[sim]**2
    sz_b=list()
    nbsfoff2=0
    for sim in range(symmetry):
        for k in shell_py[sim]:
            for l in range(nbasf[sim]+1):
                if k > 0 and l > 0 : #
                    kl_sz=k + (l-1)*nbasf[sim] + nbsfoff2 -1
                    sz_b.append(ao_overlap[kl_sz])
        nbsfoff2=nbsfoff2+nbasf[sim]**2
    
    #print('Length of B', len(sz_b))
    #   Multiplication of inverse of S_sz with B_sz
    d_sza_sym=list()
    d_szb_sym=list()
    #d_A_sza_sym=list()
    #d_A_szb_sym=list()
    for sim in range(symmetry):
        # get Sz for irrep and invert
        sz_reshape=np.array(sz_ovl[ncszoff[sim]:ncszoff[sim]+ncsz[sim]]).reshape((sz_nbasf[sim],sz_nbasf[sim]))
        la_inv_sz=LA.inv(sz_reshape)
        #print('S')
        #with np.printoptions(precision=3, suppress=True):
        #    print(sz_reshape)
        #print('S.S-1')
        #with np.printoptions(precision=3, suppress=True):
        #    print(np.dot(la_inv_sz,sz_reshape))
        #get B
        b_ovl=np.array(sz_b[nbaszoff[sim]:nbaszoff[sim]+nbasz[sim]]).reshape((sz_nbasf[sim],nbasf[sim]))
        #print('B',np.shape(b_ovl))
        #with np.printoptions(precision=4, suppress=True):
        #    print(b_ovl)       
        tb_sz = np.einsum('ik,kj->ij', la_inv_sz,b_ovl)
        #print('T-1.B', np.shape(tb_sz))
        #with np.printoptions(precision=3, suppress=True):
        #    print(tb_sz)
        #
        # C_MO coefficients
        lsym=sim+1
        nol=norb(lsym,nmo)
        nbl=nbasis(lsym,nbasf)
        # D = (T-1.B).C
        q=catchcmo(lsym,nmo,nbasf)
        # cmoa = cmo1 <N-1|
        cmosza=cmoa[comtaboff[sim]:comtaboff[sim]+cmotab[sim]]
        cmosza=np.reshape(cmosza,q[1],order='F')
        #print('CMO1',np.shape(cmosza))        
        #with np.printoptions(precision=4, suppress=True):
        #    print(cmosza)
        orb_sza = np.einsum('ik,kr->ir', tb_sz, cmosza) # D matrix for CMO1
        d_sza_sym.append(orb_sza.tolist())
        # cmob = cmo2 |N>
        cmoszb=cmob[comtaboff[sim]:comtaboff[sim]+cmotab[sim]]
        cmoszb=np.reshape(cmoszb,q[1],order='F')
        #print('CMO2',np.shape(cmoszb))
        #with np.printoptions(precision=4, suppress=True):
        #    print(cmoszb)
        #print(' ----- ' )
        orb_szb = np.einsum('ik,kr->ir', tb_sz, cmoszb) # D matrix for CMO2
        d_szb_sym.append(orb_szb.tolist())
    
        #print('Orbital SZ1', np.shape(orb_sza))
        #with np.printoptions(precision=4, suppress=True):
        #    print(orb_sza)
    
        #print('Orbital SZ2', np.shape(orb_szb))
        #with np.printoptions(precision=4, suppress=True):
        #    print(orb_szb)
        #print(' ----- ' )
    
        #print('Overlap test original MOs')
        s_int_test=ao_overlap[comtbasoff2[sim]:comtbasoff2[sim]+(nbl**2)]
        s_int_test=np.reshape(s_int_test,(nbl,nbl))
        #print('CMO1*t.S.CMO2',np.shape(cmosza.T),np.shape(s_int_test),np.shape(cmoszb) )
        ffg=np.dot(cmosza.T,s_int_test )
        ffh=np.dot(ffg,cmoszb)
        #with np.printoptions(precision=4, suppress=True):
        #    print(ffh)
    
        #print('Overlap test projected MOs')        
        ffr=np.dot(orb_sza.T,sz_reshape )
        ffe=np.dot(ffr,orb_sza)
        #print(' ----- ' )
        #print('Overlap test: < || Psi_v - Psi^sz_v || > :',np.shape(orb_szb.T), np.shape(b_ovl) , np.shape(cmoszb))
        ffw=np.dot(orb_szb.T, b_ovl )
        ffq=np.dot(ffw,cmoszb)
        minimum0=np.add(ffh,ffe)
        minimum=np.add(minimum0,-2*ffq)
        #with np.printoptions(precision=4, suppress=True):
        #    print(minimum)
        #print(' ----- ' )
    temp2_sza = [x for x in d_sza_sym if x != []] # remove null elements from D for CMO1
    temp_sza = [val for sublist in temp2_sza for val in sublist] #Flatten the list 
    d_sza = [val for sublist in temp_sza for val in sublist] #Flatten the list 
    
    temp2_szb = [x for x in d_szb_sym if x != []] # remove null elements from D for CMO2
    temp_szb = [val for sublist in temp2_szb for val in sublist] #Flatten the list 
    d_szb = [val for sublist in temp_szb for val in sublist] #Flatten the list
    
    return d_sza, d_szb
#==================================
