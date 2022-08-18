#!/usr/bin/env python
import numpy as np
from numpy import linalg as LA
import symm as sym
import mbs as sz
import oca_integrals as ocai

# rt2mzz.py
# --------------------------------------------------------------------------------------

def rt2mzz(fink_projection,conjdys,OCA_c,OCA_atom,cmoa,cmob,d_sza,d_szb,ao_dip,tdmab,totalSymmetry,symmetry,\
           nmo,nbasf,nbasft,comtbasoff,cmotab,comtaboff,basis_id_hd5,element):
#
    countorb=0
    countbasis=0
    countbasis_sz=0
    startorb=0
    totfrob=0.0
    tdmzz_symm=0.0
    norm_cdys=0.0
    norm_cdys_tot=0.0
    # ELM (L,M) contains 9 elements: { [0,0]; [1,-1]; [1,0]; [1,1]; [2,-2]; [2,-1]; [2,0]; [2,1]; [2,2] }
    ELM=list(0 for ij in range(9))    
    Auger_width=0.0

    if fink_projection:
        atombasis_list,shell_basis_sz,shell_py = sz.basis_list_for_oca(basis_id_hd5,nbasft,element,comtbasoff,nbasf,nmo)
        sz_nbasf,sznbasft,ncsz,nbasz,nszorb,nbaszoff,nszoff,ncszoff,norbsztaboff = sz.init_sz(symmetry,nbasf,shell_py,nmo)
    # cdys_r stores the coefficients for the conj. Dys orb in MO basis at each dip. op.
    cdys_x=list(0 for ij in range(symmetry))
    cdys_y=list(0 for ij in range(symmetry))
    cdys_z=list(0 for ij in range(symmetry))
    for ib in range(symmetry):
        no=sym.norb(ib+1,nmo)
        cdys_x[ib]=list(0 for ib in range(no))
        cdys_y[ib]=list(0 for ib in range(no))
        cdys_z[ib]=list(0 for ib in range(no))
    for i in range(symmetry):
        isym=i+1
        noi=sym.norb(isym,nmo)
        nbi=sym.nbasis(isym,nbasf)
        if fink_projection:
            nbi_sz=sz.nbasis_sz(isym,sz_nbasf)
        for j in range(symmetry):
            jsym=j+1
            noj=sym.norb(jsym,nmo)
            nbj=sym.nbasis(jsym,nbasf)
            if fink_projection:
                nbj_sz=sz.nbasis_sz(jsym,sz_nbasf)
            for l in range(symmetry):
                lsym=l+1
                nol=sym.norb(lsym,nmo)
                nbl=sym.nbasis(lsym,nbasf)
                if fink_projection:
                    nbl_sz=sz.nbasis_sz(lsym,sz_nbasf)
                if noi !=0 and noj !=0 and nol !=0:
                    if sym.mul((sym.mul(isym,jsym,symmetry)),lsym,symmetry)==totalSymmetry:
                        del(tdmzz_symm)
                # Compute conjugated Dyson.
                # Start by getting the Dipole matrix in MO basis
                        rmo_logic=False
                        if conjdys:
                            rmo=list()
                            for op in range(3):
                                op=op+1
                                syop=sym.symop(op,symmetry)
                                if sym.mul(isym,jsym,symmetry)==syop:
                                    rmo_logic=True
                                    dip=ao_dip[op-1]
                                    dipao=dip[comtbasoff[i]:comtbasoff[i]+nbasf[i],comtbasoff[j]:comtbasoff[j]+nbasf[j]]
                                    w=sym.catchcmo(isym,nmo,nbasf)
                                    cmo1=cmoa[comtaboff[i]:comtaboff[i]+cmotab[i]]
                                    cmo1=np.reshape(cmo1,w[1],order='F')
                                    q=sym.catchcmo(jsym,nmo,nbasf)
                                    cmo2=cmob[comtaboff[j]:comtaboff[j]+cmotab[j]]
                                    cmo2=np.reshape(cmo2,q[1],order='F')
                                    print("# 1-e Conj. Dys. Orb. symmetry (isym,jsym),(syop):",'(',isym,',',jsym,')','(',syop,')')
                                    #Block symmetry product isym,jsym,lsym
                                    scr11 = np.einsum('ia,ij->aj', cmo1,dipao)
                                    # rmo is the dipole matrix in mo basis
                                    rmo=np.einsum('aj,jb->ab', scr11 ,cmo2)
                                    #print('# shape of rmo:',np.shape(rmo))
                                    opersym=syop
                        norbtotal=noi*noj*nol
                        nbastotal=nbi*nbj*nbl
                        if fink_projection:
                            nbastotal_sz=nbi_sz*nbj_sz*nbl_sz
                            countbasis_sz=countbasis_sz + nbastotal_sz
                        countorb=countorb+norbtotal
                        countbasis=countbasis+nbastotal
                        #take symmetry tdmab and reshape
                        rtdmab=tdmab[startorb:countorb]
                        rtdmab=np.reshape(rtdmab,(noi,noj,nol))
                        # compute conjugated Dyson orbital coeff.
                        if rmo_logic:
                            # make in MO basis
                            #print("# symmetry:",isym,jsym,lsym)
                            #print('# shape of rmo,rtdmab:',np.shape(rmo),np.shape(rtdmab))
                            cdys = np.einsum('ij,ijl->l', rmo, rtdmab)
                            # make in AO basis
                            q=sym.catchcmo(lsym,nmo,nbasf)
                            cmo2=cmob[comtaboff[l]:comtaboff[l]+cmotab[l]]
                            cmo2=np.reshape(cmo2,q[1],order='F')
                            print("# symmetry Conj Dyson Orb:",lsym,', symm Op:',sym.symoprep(opersym,symmetry) )
                            cdys_ao = np.einsum('ij,j->i', cmo2, cdys)
                            if sym.symoprep(opersym,symmetry)=='x':
                                # accumulate the X components in mocdys_x
                                cdys_x[lsym-1]=cdys_x[lsym-1]+cdys
                            elif sym.symoprep(opersym,symmetry)=='y':
                                # accumulate the X components in mocdys_y
                                cdys_y[lsym-1]=cdys_y[lsym-1]+cdys
                            elif sym.symoprep(opersym,symmetry)=='z':
                                # accumulate the X components in mocdys_z
                                cdys_z[lsym-1]=cdys_z[lsym-1]+cdys
                        #catch symmetry cmoa and cmob and reshape
                        #catchcmo(isym)[0] is the product #Orb*#Basis in that symm
                        #catchcmo(isym)[1] is the tuple(#basis,#orb)
                        if not fink_projection:
                            x=sym.catchcmo(isym,nmo,nbasf)
                            cmo1=cmoa[comtaboff[i]:comtaboff[i]+cmotab[i]]
                            cmo1=np.reshape(cmo1,x[1],order='F')
                            #
                            y=sym.catchcmo(jsym,nmo,nbasf)
                            cmo21=cmob[comtaboff[j]:comtaboff[j]+cmotab[j]]
                            cmo21=np.reshape(cmo21,y[1],order='F')
                            #with np.printoptions(precision=4, suppress=True):
                            #    print(cmo21)
                            z=sym.catchcmo(lsym,nmo,nbasf)
                            cmo22=cmob[comtaboff[l]:comtaboff[l]+cmotab[l]]
                            cmo22=np.reshape(cmo22,z[1],order='F')
                            #Block symmetry product isym,jsym,lsym
                            scr1=np.einsum('ai,ijl->ajl', cmo1,rtdmab)
                            scr2=np.einsum('ajl,jb->abl',scr1,np.transpose(cmo21))
                            tdmzz_symm=np.einsum('abl,lc->abc',scr2,np.transpose(cmo22))
                        else:
                            x=sz.catchcmo_sz(isym,sz_nbasf,nmo)
                            cmo1=d_sza[norbsztaboff[i]:norbsztaboff[i]+nszorb[i]]
                            cmo1=np.reshape(cmo1,x[1],order='C')
                            #
                            y=sz.catchcmo_sz(jsym,sz_nbasf,nmo)
                            cmo21=d_szb[norbsztaboff[j]:norbsztaboff[j]+nszorb[j]]
                            cmo21=np.reshape(cmo21,y[1],order='C')
                            #with np.printoptions(precision=4, suppress=True):
                            #    print(cmo21)
                            z=sz.catchcmo_sz(lsym,sz_nbasf,nmo)
                            cmo22=d_szb[norbsztaboff[l]:norbsztaboff[l]+nszorb[l]]
                            cmo22=np.reshape(cmo22,z[1],order='C')
                            #Block symmetry product isym,jsym,lsym
                            scr1=np.einsum('ai,ijl->ajl', cmo1,rtdmab)
                            scr2=np.einsum('ajl,jb->abl',scr1,np.transpose(cmo21))
                            tdmzz_symm=np.einsum('abl,lc->abc',scr2,np.transpose(cmo22))
                        # ---------    
                        #print tdmzz_symm in 1-D
                        if not fink_projection:
                            print("# symmetry:",isym,jsym,lsym)
                            print('# shape of tdmzz:',np.shape(tdmzz_symm),', num elements:',nbastotal)
                            for izz in range(np.shape(tdmzz_symm)[0]):
                                totfrob=totfrob+(LA.norm(tdmzz_symm[izz],'fro')**2)
                                print('# Frobenius norm of TDMZZ on orb.',izz+1,':',LA.norm(tdmzz_symm[izz],'fro'))
                                for jzz in range(np.shape(tdmzz_symm)[1]):
                                    for lzz in range(np.shape(tdmzz_symm)[2]):
                                        print(izz+1,jzz+1,lzz+1, "%.12E" % tdmzz_symm[izz][jzz][lzz])
                        #
                        if fink_projection: # I am doing OCA projection
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
                                                            ELM[0] = ELM[0] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==1 and mm==-1:
                                                            ELM[1] = ELM[1] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==1 and mm==0:
                                                            ELM[2] = ELM[2] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==1 and mm==1:
                                                            ELM[3] = ELM[3] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==2 and mm==-2:
                                                            ELM[4] = ELM[4] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==2 and mm==-1:
                                                            ELM[5] = ELM[5] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==2 and mm==0:
                                                            ELM[6] = ELM[6] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==2 and mm==1:
                                                            ELM[7] = ELM[7] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        if ll==2 and mm==2:
                                                            ELM[8] = ELM[8] + tdmzz_symm[izz][jzz][lzz]*ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm)
                                                        #ELM = ELM + tdmzz_symm[izz][jzz][lzz]*elmij(gc,gi,gj,ll,mm)
                                                        if abs(ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm) )>10e-9:
                                                            print(gc,',',atombasis_list[jzz+nszoff[j]],',',atombasis_list[lzz+nszoff[l]],',',\
                                                                 "%.12E" % tdmzz_symm[izz][jzz][lzz],';','elm(',ll,mm,')',':',ocai.elmij(OCA_atom,OCA_c,gc,gi,gj,ll,mm) )
                        startorb=countorb
    if not fink_projection:
        print('#total number of #orb and #basis:',countorb,countbasis)

    return ELM,tdmzz_symm,totfrob,cdys_x,cdys_y,cdys_z

