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
import sys
import numpy as np
from numpy import linalg as LA
import symm as sym
import orbtransf as ao
import orbtransf_sz as asz
import overlap as ovlp
import mbs as sz
import rt2mzz as rt2

# auger_driver.py
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def driver_auger(input_file,fink_projection,conjdys,RAES,NAES,NAES_T,NAES_S,print_direct_dys,\
    OCA_atom,OCA_line,benergy,totalSymmetry,symmetry,nbasf,cmotab,nbasft,comtaboff,comtbasoff,comtbasoff2,\
    comtcmoff,cmob,cmoa,nmo,tdmab,dyson,hd5_file,n_elements,n_elements_desym,element,basis_id_hd5):

    # Get here overlap matrix. Only once.
    orig_stdout = sys.stdout
    twopi=6.283185307179586 # two times Pi.
    au2ev=27.211396132 # To convert to eV.
    if True: #if conjdys or print_direct_dys or fink_projection:
        basis_id_list,ao_dip,ao_overlap = ovlp.get_overlap(hd5_file,basis_id_hd5,nbasft,nbasf,symmetry,element)
        print_cmo2 = False # Print CMO2 as Molden file
    #==================================
    
    if fink_projection:
        atombasis_list,shell_basis_sz,shell_py = sz.basis_list_for_oca(basis_id_hd5,nbasft,element,comtbasoff,nbasf,nmo)
        sz_nbasf,sznbasft,ncsz,nbasz,nszorb,nbaszoff,nszoff,ncszoff,norbsztaboff = sz.init_sz(symmetry,nbasf,shell_py,nmo)
        # mbs() brings the projected MOs onto the MBS
        d_sza,d_szb = asz.mbs(ao_overlap,shell_py,symmetry,cmoa,cmob,nbasf,nmo,comtaboff,cmotab,comtbasoff2,ncszoff,ncsz,sz_nbasf,nbaszoff,nbasz)        
        #----------------------------
        output_file = 'Auger_OCA.'+input_file+'.out' # also define here the output's name
    else : 
        output_file = 'r2TM_AO.'+input_file+'.out' # also define here the output's name
    #open a output file and print all there.
    g = open(output_file, 'w')
    sys.stdout = g
    #----------------------------
    # Computing direct Dyson orbital AO and norm**2
    norm_ddys_tot,direct_dys_ao = ao.get_dyson(dyson,symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtcmoff,comtbasoff2,ao_overlap)
    if print_direct_dys:
        print('# One particle direct Dyson orbital with symmetrized basis (non-Molden file)')
        for i in range(nbasft):
            print( '{:^3} {:7} {:19.16f}'.format(i+1,basis_id_list[i],direct_dys_ao[i]) )
        print('# Sqr Norm direct dyson',norm_ddys_tot)
        print(' ')
        print('# DIRECT DYSON ORBITAL AS MOLDEN FILE')
        ao.tomolden(direct_dys_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
        print(' ')
    #------------------------------------------------
    
    if False: # Change to True if you want to see the Dyson.molden in the output
        print('# DIRECT DYSON ORBITAL AS MOLDEN FILE')
    # ao.tomolden() is the function called to print an orbital
        ao.tomolden(direct_dys_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
        print(' ')
    #-------------------
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%           If OCA             %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if fink_projection:
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
        print('# input file:',input_file,',','CDYS:',conjdys,',','OCA:',fink_projection,',','RAES:',RAES,',',\
    'AES:',NAES,',','AES triplet:',NAES_T,',','AES singlet:',NAES_S)
        print('#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        # printing decay rate from list ELM
        decay_spectrum=[]
        decay_spectrum.append(benergy) # first element is the binding energy in eV
        for c in OCA_line: # loop over scattering centers
            OCA_c = c
            print(' ')
            print('# The scattering center is :', OCA_c )
            ELM,tdmzz_symm,totfrob,cdys_x,cdys_y,cdys_z = rt2.rt2mzz(fink_projection,conjdys,OCA_c,OCA_atom,cmoa,cmob,d_sza,d_szb,ao_dip,\
            tdmab,totalSymmetry,symmetry,nmo,nbasf,nbasft,comtbasoff,cmotab,comtaboff,basis_id_hd5,element)
    
            ELM2=[elimm**2 for elimm in ELM]
            if RAES==True:
                prefactor=2.000000
            else: # then NAES True
                if NAES_S==True:
                    prefactor=1.000000
                elif NAES_T==True:
                    prefactor=3.000000
    
            decay_spectrum.append( prefactor*twopi*(sum(ELM2)) )
            print('# Auger Width One center approx','(',OCA_c,')' ,'= ', prefactor*twopi*(sum(ELM2)) )
        print(' ')
        print('# Spectrum: BE(eV) and Intensity from OCA:', *decay_spectrum ) #asterisk before list's name with print statement prints only the elements
    else:
        OCA_c = 'None'
        d_sza=[]
        d_szb=[]
        ELM,tdmzz_symm,totfrob,cdys_x,cdys_y,cdys_z = rt2.rt2mzz(fink_projection,conjdys,OCA_c,OCA_atom,cmoa,cmob,d_sza,d_szb,ao_dip,\
        tdmab,totalSymmetry,symmetry,nmo,nbasf,nbasft,comtbasoff,cmotab,comtaboff,basis_id_hd5,element)
        print('# Total Frobenius norm of TDMZZ',totfrob)
    
    print('# The End ')
    sys.stdout = orig_stdout
    g.close()
    
    #------------------------------------------------
    if print_direct_dys:
        if  norm_ddys_tot > 10**-14:
            ar='direct_dyson.'+input_file+'.out'
            dys_dir_print = open(ar, 'a')
            sys.stdout = dys_dir_print
            print('#1-p DYSON ORBITAL AS MOLDEN FILE:',norm_ddys_tot)
            # ao.tomolden() is the function called to print an orbital
            ao.tomolden(direct_dys_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
            sys.stdout = orig_stdout
            dys_dir_print.close()
    #------------------------------------------------
    if print_cmo2:
        if  True:
            ar='CMO2.'+input_file+'.out'
            cmo2_dir_print = open(ar, 'a')
            sys.stdout = cmo2_dir_print
            print('#CMO2 AS MOLDEN FILE:')
            # ao.tomolden() is the function called to print an orbital
            cmo2_ao = ao.write_cmo2(symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtbasoff2)
            ao.tomolden(cmo2_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
            sys.stdout = orig_stdout
            cmo2_dir_print.close()

    #------------------------------------------------
    if conjdys==True and norm_ddys_tot > 10**-14:
        if sym.norb(sym.mul(totalSymmetry,sym.symop(1,symmetry),symmetry),nmo)!=0:
            ar='conj.x.'+input_file+'.out'
            cx = open(ar, 'a')
            sys.stdout = cx
            #print('#Molden conj. Dyson orb. component x')
            conjugate_dys_ao,norm_conjdys_tot = ao.write_conj_dys(symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtbasoff2,cdys_x,ao_overlap)
            print('# Conj. Dyson symm:',sym.mul(totalSymmetry,sym.symop(1,symmetry),symmetry),'=',totalSymmetry,'.X component; norm:',norm_conjdys_tot)
            ao.tomolden(conjugate_dys_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
            sys.stdout = orig_stdout
            cx.close()
    #-----------------------------------------------
        if sym.norb(sym.mul(totalSymmetry,sym.symop(2,symmetry),symmetry),nmo)!=0:
            ar='conj.y.'+input_file+'.out'
            cy = open(ar, 'a')
            sys.stdout = cy
            #print('#Molden conj. Dyson orb. component y')
            conjugate_dys_ao,norm_conjdys_tot = ao.write_conj_dys(symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtbasoff2,cdys_y,ao_overlap)
            print('# Conj. Dyson symm:',sym.mul(totalSymmetry,sym.symop(2,symmetry),symmetry),'=',totalSymmetry,'.Y component; norm:',norm_conjdys_tot)
            ao.tomolden(conjugate_dys_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
            sys.stdout = orig_stdout
            cy.close()
    #------------------------------------------------
        if sym.norb(sym.mul(totalSymmetry,sym.symop(3,symmetry),symmetry),nmo)!=0:
            ar='conj.z.'+input_file+'.out'
            cz = open(ar, 'a')
            sys.stdout = cz
            #print('#Molden conj. Dyson orb. component z')
            conjugate_dys_ao,norm_conjdys_tot = ao.write_conj_dys(symmetry,cmob,nmo,nbasf,comtaboff,cmotab,comtbasoff2,cdys_z,ao_overlap)
            print('# Conj. Dyson symm:',sym.mul(totalSymmetry,sym.symop(3,symmetry),symmetry),'=',totalSymmetry,'.Z component; norm:',norm_conjdys_tot)
            ao.tomolden(conjugate_dys_ao,basis_id_hd5,symmetry,element,n_elements,n_elements_desym,nbasft)
            sys.stdout = orig_stdout
            cz.close()
    return

