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

# rt2mzz.py
# --------------------------------------------------------------------------------------

def rt2mzz(OCA_c,OCA_atom,cmoa,cmob,d_sza,d_szb,tdmab,totalSymmetry,symmetry,\
           nmo,nbasf,nbasft,comtbasoff,cmotab,comtaboff,basis_id_hd5,element):
    from .symm import mul,norb,nbasis
    from .mbs import nbasis_sz,catchcmo_sz,basis_list_for_oca,init_sz
    from .oca_integrals import elmij

    countorb=0
    countbasis=0
    countbasis_sz=0
    startorb=0
    totfrob=0.0
    tdmzz_symm=0.0
    # ELM (L,M) contains 9 elements: { [0,0]; [1,-1]; [1,0]; [1,1]; [2,-2]; [2,-1]; [2,0]; [2,1]; [2,2] }
    ELM=list(0 for ij in range(9))    
    Auger_width=0.0

    atombasis_list,shell_basis_sz,shell_py = basis_list_for_oca(basis_id_hd5,nbasft,element,comtbasoff,nbasf,nmo)
    sz_nbasf,sznbasft,ncsz,nbasz,nszorb,nbaszoff,nszoff,ncszoff,norbsztaboff = init_sz(symmetry,nbasf,shell_py,nmo)
    for i in range(symmetry):
        isym=i+1
        noi=norb(isym,nmo)
        nbi=nbasis(isym,nbasf)
        nbi_sz=nbasis_sz(isym,sz_nbasf)
        for j in range(symmetry):
            jsym=j+1
            noj=norb(jsym,nmo)
            nbj=nbasis(jsym,nbasf)
            nbj_sz=nbasis_sz(jsym,sz_nbasf)
            for l in range(symmetry):
                lsym=l+1
                nol=norb(lsym,nmo)
                nbl=nbasis(lsym,nbasf)
                nbl_sz=nbasis_sz(lsym,sz_nbasf)
                if noi !=0 and noj !=0 and nol !=0:
                    if mul((mul(isym,jsym,symmetry)),lsym,symmetry)==totalSymmetry:
                        del(tdmzz_symm)
                        norbtotal=noi*noj*nol
                        nbastotal=nbi*nbj*nbl
                        nbastotal_sz=nbi_sz*nbj_sz*nbl_sz
                        countbasis_sz=countbasis_sz + nbastotal_sz
                        countorb=countorb+norbtotal
                        countbasis=countbasis+nbastotal
                        #take symmetry tdmab and reshape
                        rtdmab=tdmab[startorb:countorb]
                        rtdmab=np.reshape(rtdmab,(noi,noj,nol))
                        #catch symmetry cmoa and cmob and reshape
                        #catchcmo(isym)[0] is the product #Orb*#Basis in that symm
                        #catchcmo(isym)[1] is the tuple(#basis,#orb)

                        x=catchcmo_sz(isym,sz_nbasf,nmo)
                        cmo2i=d_sza[norbsztaboff[i]:norbsztaboff[i]+nszorb[i]]
                        cmo2i=np.reshape(cmo2i,x[1],order='C')
                        #
                        y=catchcmo_sz(jsym,sz_nbasf,nmo)
                        cmo2j=d_szb[norbsztaboff[j]:norbsztaboff[j]+nszorb[j]]
                        cmo2j=np.reshape(cmo2j,y[1],order='C')
                        #with np.printoptions(precision=4, suppress=True):
                        #    print(cmo21)
                        z=catchcmo_sz(lsym,sz_nbasf,nmo)
                        cmo2l=d_szb[norbsztaboff[l]:norbsztaboff[l]+nszorb[l]]
                        cmo2l=np.reshape(cmo2l,z[1],order='C')
                        #Block symmetry product isym,jsym,lsym
                        scr1=np.einsum('ijl,lc->ijc', rtdmab,np.transpose(cmo2l))
                        scr2=np.einsum('ijc,jb->ibc',scr1,np.transpose(cmo2j))
                        tdmzz_symm=np.einsum('ai,ibc->abc',cmo2i,scr2)
                        # ---------    
                        #print tdmzz_symm in 1-D
                        fizz=range(np.shape(tdmzz_symm)[0])
                        fjzz=range(np.shape(tdmzz_symm)[1])
                        flzz=range(np.shape(tdmzz_symm)[2])
                        #print('# shape of tdmzz:',np.shape(tdmzz_symm),', num elements:',nbastotal)
                        for izz in fizz:
                            #print('# Frobenius norm of TDMZZ on orb.',izz+1,':',LA.norm(tdmzz_symm[izz],'fro'))
                            totfrob=totfrob+(LA.norm(tdmzz_symm[izz],'fro'))
                            # first element must be the core
                            if atombasis_list[izz+nszoff[i]] == OCA_c :
                                print("# symmetry:",isym,jsym,lsym)
                                #for jzz in fjzz:
                                    #for lzz in flzz:
                                #Accumulate the integral ELM
                                #print('Elm_ij',tdmzz_symm[izz][jzz][lzz])
                                for ll in range(3): # L=0,1,2 
                                    for mm in range(-ll,ll+1): # M = -L,L
                                        for jzz in fjzz:
                                            for lzz in flzz:
                                                if atombasis_list[jzz+nszoff[j]][:2]==OCA_c[:2] and atombasis_list[lzz+nszoff[l]][:2]==OCA_c[:2]:
                                                    gc=atombasis_list[izz+nszoff[i]]
                                                    gi=atombasis_list[jzz+nszoff[j]][2:].replace(" ","")
                                                    gj=atombasis_list[lzz+nszoff[l]][2:].replace(" ","")
                                                    if ll==0 and mm==0:
                                                        ELM[0] = ELM[0] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==1 and mm==-1:
                                                        ELM[1] = ELM[1] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==1 and mm==0:
                                                        ELM[2] = ELM[2] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==1 and mm==1:
                                                        ELM[3] = ELM[3] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==2 and mm==-2:
                                                        ELM[4] = ELM[4] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==2 and mm==-1:
                                                        ELM[5] = ELM[5] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==2 and mm==0:
                                                        ELM[6] = ELM[6] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==2 and mm==1:
                                                        ELM[7] = ELM[7] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    if ll==2 and mm==2:
                                                        ELM[8] = ELM[8] + tdmzz_symm[izz][jzz][lzz]*elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                    #ELM = ELM + tdmzz_symm[izz][jzz][lzz]*elmij(gc,gi,gj,ll,mm)
                                                    if abs(elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm) )>10e-9:
                                                        print(gc,',',atombasis_list[jzz+nszoff[j]],',',atombasis_list[lzz+nszoff[l]],',',\
                                                            "%.12E" % tdmzz_symm[izz][jzz][lzz],';','elm(',ll,mm,')',':',elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm) )
                        startorb=countorb

    return ELM,tdmzz_symm,totfrob

