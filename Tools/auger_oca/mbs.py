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

## mbs.py
import sys

 ######### definition of MBS via functions shell_basis(lista,mbs_func) and shell_basis(lista)
def shell_basis(lista,func,comtbasoff,nbasf,nbasft,nmo):
# this function returns the MBS from the full GTO based on MBS
    index_e=list()
    MBS_excl=['2s','2px','2py','2pz']
    for e in lista:
        g=e.split()
        if g[1] in func :
            index_e.append(lista.index(e)+1)
            #if g[0][0] == 'H':
            #    if g[1] in MBS_excl:
            #        continue    

    # now check for H 2s 2p functions
    hydro=list() # this collect indexes of 2s 2p of H functions
    for ij in index_e:
        i=ij-1
        ind_dummy=lista[i].split()
        if ind_dummy[0][0] == 'H':
            if ind_dummy[1] in MBS_excl:
                hydro.append(ij)

    index_eh=[x for x in index_e if x not in hydro] # elements of index_e minus hydro
    index_eh= list(dict.fromkeys(index_eh)) # remove repeted basis indexes

    shell= list()
    scr= list()
    prevbasf = comtbasoff[:] + [nbasft]
    for ii in range(len(nbasf)):
        for r in index_eh :
            if prevbasf[ii] < r <= prevbasf[ii+1] :
                scr.append(r)
        shell.append(scr)
        scr=list()
    return shell

def shell_basis_py(lista,nbasf,comtbasoff):
# this function returns ... 
    shell= list()
    scr= list()
    for t in range(len(nbasf)):
        scr = [ (u - comtbasoff[t]) for u in lista[t] ]
        shell.append(scr)
    return shell


def nbasis_sz(i,sz_nbasf):
    #returns the number of basis function in the i symm
    return int(sz_nbasf[i-1])

def catchcmo_sz(i,sz_nbasf,nmo):
    #x=Num orb in symm
    x=int(nmo[i-1])
    #y=Num basis in symm
    y=int(sz_nbasf[i-1])
    #returns a product #orb*#basis and a tuple (#Orb,#Basis)
    return (x*y,(y,x))

#========================================================

# HERE THE MBS IS DEFINED BASED ON basis_id_hd5
def basis_list_for_oca(basis_id_hd5,nbasft,element,comtbasoff,nbasf,nmo):
    import orbtransf as ao

    dunning=True # True -> Automatic selection of MBS. False -> manual. See below

    # now I make a list len(nbasft) with all basis ids
    basis_id_list=list()
    for i in range(nbasft):
        basis_id_list.append(ao.basis_func_id(basis_id_hd5[i],element))
        #print(ao.basis_func_id(basis_id_hd5[i]))
    basis_id_list_shell=basis_id_list # I need this to define the MBS
    #print(basis_id_list_shell)

    if dunning: # define automatic selection of MBS
    # make use of basis_id_list_shell to make MBS
    #print('Automatic basis selection from Dunning basis set')
        MBS=['1s','2s','2px','2py','2pz']
    # shell is the MBS
        shell = shell_basis(basis_id_list_shell,MBS,comtbasoff,nbasf,nbasft,nmo)
    # Else, define the basis manualy below. See example.
    else: # manual selection of MBS with shell.
    # Example: O3 CS cc-pvtz
    # shell=[[1, 2, 5, 8, 21, 22, 25, 28, 41, 42, 45, 48], [61, 71, 81]]
        shell=[]

    shell_py = shell_basis_py(shell,nbasf,comtbasoff)
    shell_basis_sz=[(fg-1) for fg in [val for sublist in shell for val in sublist] ]
    atombasis_list=[basis_id_list[ei] for ei in shell_basis_sz]
    
    return atombasis_list,shell_basis_sz,shell_py

#========================================================

def init_sz(symmetry,nbasf,shell_py,nmo):
    #sz_nbasf = number of basis elements in SZ for each irrep
    sz_nbasf=[len(im) for im in shell_py ]
    #nbasft = total number of basis function
    sznbasft=int(sum(sz_nbasf))

    # ncsz = number of elements in the overlap matrix in SZ basis. 
    # nbasz = number of elements in matrix B_sz
    ncsz= list(sz_nbasf[ind]**2 for ind in range(symmetry))
    nbasz=list(sz_nbasf[ind]*nbasf[ind] for ind in range(symmetry))
    nszorb=list(sz_nbasf[ind]*nmo[ind] for ind in range(symmetry)) # similar to cmotab

    ncszoff=list()  # number of previus elements in S_sz
    off=0
    for ts in ncsz:
        ncszoff.append(off)
        off=off+ts
    nbaszoff=list() # number of previus elements in B_sz
    off=0
    for ts in nbasz:
        nbaszoff.append(off)
        off=off+ts
    nszoff=list() # number of previus basis elements in D_sz
    off=0
    for ts in sz_nbasf:
        nszoff.append(off)
        off=off+ts
    norbsztaboff=list() # number of previus elements in D_sz. Similar to comtaboff
    off=0
    for ts in nszorb:
        norbsztaboff.append(off)
        off=off+ts

    return sz_nbasf,sznbasft,ncsz,nbasz,nszorb,nbaszoff,nszoff,ncszoff,norbsztaboff
#----------------------------------------

