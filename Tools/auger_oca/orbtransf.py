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
import symm as sym
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
    if l==1:
        if m==0:
            func_id='pz'
        if m==1:
            func_id='px'
        if m==-1:
            func_id='py'
    if l==2:
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
    if l==3:
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
    return '{}{}{}{}'.format(list_element[c-1]," ",n+l, func_id)
# The following order of D, F en G functions is expected:
# 5D: D 0, D+1, D-1, D+2, D-2
# 7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
# 9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4


#==================================

def phase_lm(cnlm,symmetry):
    phase_coeff=1.0
    # phase_lm reads a list of integers (c,n,l,m) so to multiply
    # func.coeff by -1 if the centers basis change phase with the irrep. transformation.
    # See the "Symmetry adapted Basis Functions" on your molcas.log file.
    if symmetry==4: #C2v
        # py
        if (cnlm[2] ==1 and cnlm[3] ==-1):
            phase_coeff=-1.0
        # d functions: d-1, d-2
        if (cnlm[2] ==2 and cnlm[3] ==-1):
            phase_coeff=-1.0
        if (cnlm[2] ==2 and cnlm[3] ==-2):
            phase_coeff=-1.0
        # f functions: f-1, f-2, f-3
        if (cnlm[2] ==3 and cnlm[3] ==-1):
            phase_coeff=-1.0
        if (cnlm[2] ==3 and cnlm[3] ==-2):
            phase_coeff=-1.0
        if (cnlm[2] ==3 and cnlm[3] ==-3):
            phase_coeff=-1.0
    return phase_coeff

#==================================
" tomolden(orbital) is a program to receive an orbital as symmetrized basis and print as Molden format"
# the orbital input should be a flat list.
# for now I am accepting up to 12 different symmetric equivalalent elements

#==================================

def tomolden(orb,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft):
# search for equal list elements and equalize
        search=basis_id_hd5[:,:4].astype(int)
        orbital=np.array(orb)
        orbital2=np.array(orb)
        #loop over search
        for s in range(len(search)):
                for h in range(len(search)-s-1):
                        if all(search[s]==search[h+s+1]):
                                #this will be true if I have symmetricaly equivalent atoms
                                # use function phase_lm to decide whether orb* + or - 1.
                                phase_coeff=phase_lm(search[s],symmetry)
                                orbital[s]=(orbital2[s]+orbital2[h+s+1])/np.sqrt(2)
                                orbital[h+s+1]=(orbital2[s]-orbital2[h+s+1])*phase_coeff/np.sqrt(2)
                                #here i am saying there is another atom of the same element
                                search[h+s+1][0]=int(search[h+s+1][0]+400)
        #
        lista = np.ndarray((nbasft,5),dtype = object)
        lista[:,:4]=search[:,:4].astype(int)
        lista[:,4]=orbital[:].astype(float)
        #
        tempfile1=[]
        tempfile2=[]
        tempfile3=[]
        tempfile4=[]
        tempfile5=[]
        tempfile6=[]
        tempfile7=[]
        tempfile8=[]
        tempfile9=[]
        tempfile10=[]
        tempfile11=[]
        tempfile12=[]
        tempfile13=[]
        tempfile14=[]
        tempfile15=[]
        tempfile16=[]
        tempfile17=[]
        tempfile18=[]
        tempfile19=[]
        tempfile20=[]
        tempfile21=[]
        tempfile22=[]
        tempfile23=[]
        tempfile24=[]
        tempfile99=[]
        molden=[]
        molden1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        max_ele=n_elements
        # separate the list by different elements
        for j in range(len(lista)):
                if lista[j][0] == 1:
                        tempfile1.append(lista[j].tolist())
                if lista[j][0] == 401:
                        tempfile2.append(lista[j].tolist())
                if lista[j][0] == 2:
                        tempfile3.append(lista[j].tolist())
                if lista[j][0] == 402:
                        tempfile4.append(lista[j].tolist())
                if lista[j][0] == 3:
                        tempfile5.append(lista[j].tolist())
                if lista[j][0] == 403:
                        tempfile6.append(lista[j].tolist())
                if lista[j][0] == 4:
                        tempfile7.append(lista[j].tolist())
                if lista[j][0] == 404:
                        tempfile8.append(lista[j].tolist())
                if lista[j][0] == 5:
                        tempfile9.append(lista[j].tolist())
                if lista[j][0] == 405:
                        tempfile10.append(lista[j].tolist())
                if lista[j][0] == 6:
                        tempfile11.append(lista[j].tolist())
                if lista[j][0] == 406:
                        tempfile12.append(lista[j].tolist())
                if lista[j][0] == 7:
                        tempfile13.append(lista[j].tolist())
                if lista[j][0] == 407:
                        tempfile14.append(lista[j].tolist())
                if lista[j][0] == 8:
                        tempfile15.append(lista[j].tolist())
                if lista[j][0] == 408:
                        tempfile16.append(lista[j].tolist())
                if lista[j][0] == 9:
                        tempfile17.append(lista[j].tolist())
                if lista[j][0] == 409:
                        tempfile18.append(lista[j].tolist())
                if lista[j][0] == 10:
                        tempfile19.append(lista[j].tolist())
                if lista[j][0] == 410:
                        tempfile20.append(lista[j].tolist())
                if lista[j][0] == 11:
                        tempfile21.append(lista[j].tolist())
                if lista[j][0] == 411:
                        tempfile22.append(lista[j].tolist())
                if lista[j][0] == 12:
                        tempfile23.append(lista[j].tolist())
                if lista[j][0] == 412:
                        tempfile24.append(lista[j].tolist())
        tempfile99.append(tempfile1)
        tempfile99.append(tempfile2)
        tempfile99.append(tempfile3)
        tempfile99.append(tempfile4)
        tempfile99.append(tempfile5)
        tempfile99.append(tempfile6)
        tempfile99.append(tempfile7)
        tempfile99.append(tempfile8)
        tempfile99.append(tempfile9)
        tempfile99.append(tempfile10)
        tempfile99.append(tempfile11)
        tempfile99.append(tempfile12)
        tempfile99.append(tempfile13)
        tempfile99.append(tempfile14)
        tempfile99.append(tempfile15)
        tempfile99.append(tempfile16)
        tempfile99.append(tempfile17)
        tempfile99.append(tempfile18)
        tempfile99.append(tempfile19)
        tempfile99.append(tempfile20)
        tempfile99.append(tempfile21)
        tempfile99.append(tempfile22)
        tempfile99.append(tempfile23)
        tempfile99.append(tempfile24)
        list99=tempfile99
        # the next line is a single line command to remove empty lists of a list.
        tempfile99 = [x for x in list99 if x != []]
        for e in range(n_elements_desym):
                xt=tempfile99[e]
                xt=np.array(xt)
                s=[]
                px=[]
                py=[]
                pz=[]
                d0=[]
                dp1=[]
                dm1=[]
                dp2=[]
                dm2=[]
                f0=[]
                fp1=[]
                fm1=[]
                fp2=[]
                fm2=[]
                fp3=[]
                fm3=[]
                p=[]
                d=[]
                f=[]
                xm=list()
                max_ang=max(xt[:,2])
                for j in range(len(xt)):
                        if xt[j][2] == 0:
                                s.append(xt[j])
                        if xt[j][2] == 1:
                                if xt[j][3] == 0:
                                        pz.append(xt[j])
                                if xt[j][3] == 1:
                                        px.append(xt[j])
                                if xt[j][3] == -1:
                                        py.append(xt[j])
                        if xt[j][2] == 2:
                                if xt[j][3] == 0:
                                        d0.append(xt[j])
                                if xt[j][3] == 1:
                                        dp1.append(xt[j])
                                if xt[j][3] == -1:
                                        dm1.append(xt[j])
                                if xt[j][3] == 2:
                                        dp2.append(xt[j])
                                if xt[j][3] == -2:
                                        dm2.append(xt[j])
                        if xt[j][2] == 3:
                                if xt[j][3] == 0:
                                        f0.append(xt[j])
                                if xt[j][3] == 1:
                                        fp1.append(xt[j])
                                if xt[j][3] == -1:
                                        fm1.append(xt[j])
                                if xt[j][3] == 2:
                                        fp2.append(xt[j])
                                if xt[j][3] == -2:
                                        fm2.append(xt[j])
                                if xt[j][3] == 3:
                                        fp3.append(xt[j])
                                if xt[j][3] == -3:
                                        fm3.append(xt[j])
                if max_ang >= 0:
                        s=np.array(s)
                        s=s.tolist()
                        xm.append(s)
                if max_ang >= 1:
                        #p=px+py+pz
                        #py=np.array(py)
                        for m in range(len(py)):
                                p.append(px[m])
                                p.append(py[m])
                                p.append(pz[m])
                        #p=p[p[:, 1].argsort()]
                        #p=p.tolist()
                        xm.append(p)
                if max_ang >= 2:
                        #d=d0+dp1+dm1+dp2+dm2
                        #d=np.array(d)
                        #d=d[d[:, 1].argsort()]
                        #d=d.tolist()
                        for m in range(len(d0)):
                                d.append(d0[m])
                                d.append(dp1[m])
                                d.append(dm1[m])
                                d.append(dp2[m])
                                d.append(dm2[m])
                        xm.append(d)
                if max_ang >= 3:
                        #f=f0+fp1+fm1+fp2+fm2+fp3+fm3
                        #f=np.array(f)
                        #f=f[f[:, 1].argsort()]
                        #f=f.tolist()
                        for m in range(len(f0)):
                                f.append(f0[m])
                                f.append(fp1[m])
                                f.append(fm1[m])
                                f.append(fp2[m])
                                f.append(fm2[m])
                                f.append(fp3[m])
                                f.append(fm3[m])
                        xm.append(f)
                xm=np.concatenate(xm)
                molden1[e]=xm
        for w in range(n_elements_desym):
                molden1[w]=molden1[w].tolist()
                molden.append(molden1[w])
        molden=np.concatenate(molden)
        if max(molden[:,0])>400:
        # this if will be true if I did change the index for some symmetry equivalent atom
                for at in range(len(molden[:,0])):
        # here I changing the extra element index back as in the original basis_id_hd5
                        if molden[at][0] > 400:
                                molden[at][0]=int(molden[at][0]-400)
        #print('MOLDEN FILE Dyson')
        molden_print=molden[:,:4].astype(int)
        molden_ao=np.array(molden[:,4])

        basis_id_list=list()
        for i in range(nbasft):
                basis_id_list.append(basis_func_id(molden_print[i],element))
        #print('# One particle direct Dyson orbital in MOLDEN format')
        for j in range(nbasft):
                print( '{:^3} {:7} {:19.16f}'.format(j+1,basis_id_list[j],molden_ao[j]) )
        return

#==================================

def mo2ao(mo_orb,cmo_ab,symmetry,nmo,nbasf,comtcmoff,comtaboff,cmotab):
    # function mo2ao transforms MO to AO using the transformation matrix cmo_ab.
    # store block symmetry for orbital
    total_orb=list(0 for i in range(symmetry))
    total_orb_ao=list(0 for i in range(symmetry))
    for s in range(symmetry):
        isym=s+1
        total_orb[s]=mo_orb[comtcmoff[s]:comtcmoff[s]+sym.norb(isym,nmo)]
    # Computing orbital AO 
    for l in range(symmetry):
        lsym=l+1
        nol=sym.norb(lsym,nmo)
        nbl=sym.nbasis(lsym,nbasf)
        # direct Dyson AO transformation
        q=sym.catchcmo(lsym,nmo,nbasf)
        cmon=cmo_ab[comtaboff[l]:comtaboff[l]+cmotab[l]]
        cmon=np.reshape(cmon,q[1],order='F')
        orb_mo=total_orb[l]
        orb_ao = np.einsum('ij,j->i', cmon, orb_mo)
#       store symmetry block in direct_dys_ao
        total_orb_ao[l]=orb_ao
    total_orb_ao2=[val for sublist in total_orb_ao for val in sublist] #Flatten the list
    return total_orb_ao2

#==================================

def get_dyson(dyson,symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtcmoff,comtbasoff2,ao_overlap):
    # Computing direct Dyson orbital AO and norm**2
    norm_ddys_tot=0.0
    # store block symmetry for regular (direct) Dyson
    direct_dys=list(0 for i in range(symmetry))
    direct_dys_ao=list(0 for i in range(symmetry))
    for s in range(symmetry):
        isym=s+1
        direct_dys[s]=dyson[comtcmoff[s]:comtcmoff[s]+sym.norb(isym,nmo)]

    for l in range(symmetry):
        lsym=l+1
        nol=sym.norb(lsym,nmo)
        nbl=sym.nbasis(lsym,nbasf)
        # direct Dyson AO transformation
        q=sym.catchcmo(lsym,nmo,nbasf)
        cmo2=cmob[comtaboff[l]:comtaboff[l]+cmotab[l]]
        cmo2=np.reshape(cmo2,q[1],order='F')
        ddys=direct_dys[l]
        ddys_ao = np.einsum('ij,j->i', cmo2, ddys)
    #   store symmetry block in direct_dys_ao
        direct_dys_ao[l]=ddys_ao
        s_int=ao_overlap[comtbasoff2[l]:comtbasoff2[l]+(nbl**2)]
        s_int=np.reshape(s_int,(nbl,nbl))
        norm_ddys=  np.einsum('ij,j->i',s_int,ddys_ao)
        norm_ddys=  np.einsum('j,j->',norm_ddys,ddys_ao)
        # square norm of direct Dyson term
        norm_ddys_tot=norm_ddys_tot+norm_ddys
    # Here I get the direct Dyson Orb. (with symmetrized basis standard)
    # [val for sublist in list_of_lists for val in sublist] is a simple one-line way to flat a list_of_lists
    flatten_dys=[val for sublist in direct_dys_ao for val in sublist]
    direct_dys_ao=flatten_dys    

    return norm_ddys_tot,direct_dys_ao

#==================================

def write_conj_dys(symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtbasoff2,cdys_x,ao_overlap):
    #print('#Molden conj. Dyson orb. component xyz')
    conjugate_dys_ao=list(0 for i in range(symmetry))
    norm_conjdys_tot=0.0
    for l in range(symmetry):
        lsym=l+1
        nol=sym.norb(lsym,nmo)
        nbl=sym.nbasis(lsym,nbasf)
        # conj Dyson AO transformation
        q=sym.catchcmo(lsym,nmo,nbasf)
        cmo2=cmob[comtaboff[l]:comtaboff[l]+cmotab[l]]
        cmo2=np.reshape(cmo2,q[1],order='F')
        conjdys_r=cdys_x[l]
        conjdys_ao = np.einsum('ij,j->i', cmo2, conjdys_r)
    #   store symmetry block in conjugate_dys_ao[l]
        conjugate_dys_ao[l]=conjdys_ao
        s_int=ao_overlap[comtbasoff2[l]:comtbasoff2[l]+(nbl**2)]
        s_int=np.reshape(s_int,(nbl,nbl))
        norm_cdys=  np.einsum('ij,j->i',s_int,conjdys_ao)
        norm_cdys=  np.einsum('j,j->',norm_cdys,conjdys_ao)
        # square norm of direct Dyson term
        norm_conjdys_tot=norm_conjdys_tot+norm_cdys
    flatten_cdys=[val for sublist in conjugate_dys_ao for val in sublist]
    conjugate_dys_ao=flatten_cdys
    return conjugate_dys_ao,norm_conjdys_tot

#==================================
def write_cmo2(symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtbasoff2):
    for x in range(symmetry):
        xsym=x+1
        nox=sym.norb(xsym,nmo)
        if nox !=0:
            for y in range(nox):
                yorb=y+1
                vecmo=list(0 for i in range(symmetry))
                for z in range(symmetry):
                    zsym=z+1
                    vecmo[z]=np.zeros(nmo[z])
                vecmo[x][y]=1.0
                #print('# Symmetry ',xsym,',  Orbital ',yorb)
                cmo2_ao_sym=list(0 for j in range(symmetry))
                for i in range(symmetry):
                    cmo2_ao=list(0 for j in range(symmetry))
                    isym=i+1
                    noi=sym.norb(isym,nmo)
                    nbi=sym.nbasis(isym,nbasf)
                    q=sym.catchcmo(isym,nmo,nbasf)
                    cmo2=cmob[comtaboff[i]:comtaboff[i]+cmotab[i]]
                    cmo2=np.reshape(cmo2,q[1],order='F')
                    scr1 = vecmo[i]
                    #print('cmo2',cmo2)
                    scr2 = np.einsum('ij,j->i', cmo2, scr1)
                    cmo2_ao_sym[i]=scr2
    
                scr3 =[val for sublist in cmo2_ao_sym for val in sublist]
                #cmo2_ao=scr3
    return scr3

#==================================
