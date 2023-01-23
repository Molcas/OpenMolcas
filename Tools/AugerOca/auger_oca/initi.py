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

import subprocess, logging, re, sys, os, time, h5py
from threading import Thread
import numpy as np

# initi.py
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# READ some information from Density files provided by OpenMolcas:

def call_subprocess(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if err:
        print(err)

def init1():
    from .input_parse import parseCL

    args = parseCL()
 
    # remove the coments ('#') from the 'input' file r2TM.dat and dump into auger_dyson.dat
    #if args.directory:

    RAES=args.parse_raes           # default True if OCA. Controled by input --raes
                                   # RAES defines the prefactor of 2(true) or 1(false).
    # When RAES false, remember to apply the prefactor of 3 on triplets when plotting. The singlet's prefactor is 1.
    NAES=args.parse_aes             # default False. Controled by input --aes. Not to be used with --raes!
    NAES_T = args.parse_aes_triplet
    NAES_S = args.parse_aes_singlet
    
    file_based = args.inp 
    folder_based = args.directory
    return file_based,folder_based,RAES,NAES,NAES_T,NAES_S 

    #==================================

def init2(input_file):

    f=open(input_file,"r")
    lines=list()
    tmp_files=f.readlines()
    for i in tmp_files:
        if not (i.startswith(' #') or i.startswith('#')):
            lines.append(i)

    # 1st line tells me how much scattering centers from OCA I have
    nsc=int(lines[0]) # number of scattering centers
    # first nsc lines, I have the scattering centers
    OCA_line=list()
    for i in range(nsc):
        OCA_line.append(str(lines[i+1]).strip())
    OCA_line=[c.split()[0].upper()+' '+c.split()[1] for c in OCA_line] # make elements upper case letter
    OCA_center= OCA_line[0] #str(lines[1]).strip()

    # The line immediately after the last scattering center gives the binding energy
    benergy=float(lines[nsc+1])
    # nsc+2 line, Total symmetry of WF product <N-1,N>:
    totalSymmetry=int(lines[nsc+2])
    # nsc+3 line, Number of irreps (symmetry number):
    symmetry=int(lines[nsc+3])
    # nsc+4 line, Number of basis functions (per irrep):
    nbasf= [int(num) for num in lines[nsc+4].split()]
    # nsc+5 line, Number of active MOs (per irrep):
    nash= [int(num) for num in lines[nsc+5].split()]
    # nsc+6 line, Number of total MOs (per irrep). NActive + NInactive:
    nmo= [int(num) for num in lines[nsc+6].split()]
    
    npl=int(int(7)+nsc) #npl=number of previus lines. see 'Collecting data from the file'
    #end of reading    

    #define number of orbitals, nmo, and number basis func per symmetry, nbasis, from input.
    cmotab=list(nmo[ind]*nbasf[ind] for ind in range(symmetry))
    tdmtab=list(nmo[ind]**3 for ind in range(symmetry))
    
    #ncmo=total number of MO components
    ncmo=int(sum(cmotab))
    #nbasft = total number of basis function
    nbasft=int(sum(nbasf))
    #nasht=number of active orbitals
    nasht=int(sum(nash))
    #nosht=total number of orbitals
    nosht=int(sum(nmo))
    #comtaboff: number of orb in previus symm
    comtaboff=list()
    comtbasoff=list()
    comtbasoff2=list()
    comtcmoff=list()
    comtnashoff=list()
    off=0
    for x in cmotab:
        comtaboff.append(off)
        off=off+x
    off=0
    for y in nbasf:
        comtbasoff.append(off)
        off=off+y
    off=0
    for v in nbasf:
        comtbasoff2.append(off)
        off=off+v**2
    off=0
    for z in nmo:  # number of previus occupied orbitals
        comtcmoff.append(off)
        off=off+z
    off=0
    for z in nash:  # number of previus active orbitals
        comtnashoff.append(off)
        off=off+z
    
    #==================================
    
    # get the MOs (cmoa,cmob) and 2p-Dyson densities (tdmab) in MO basis.
    cmob=list()
    cmoa=list()
    tdmab=list()
    dyson=list()
    
    #lines=f.readlines()
    # 2st line, Total symmetry of WF product <N-1,N>:
    # CMO1 comes from PSI1 that is the cation.
    # CMO2 comes from PSI2 that is tne neutral
    #totalSymmetry=int(lines[0])
    
    # number previus_lines=npl
    for z in range(ncmo):
        cmoa.append(float(lines[z+npl]))  #cmoa = cmo1
        cmob.append(float(lines[z+ncmo+npl])) #cmob = cmo2
    
    for x in range(nosht):
        val=float(lines[2*ncmo+x+npl].split()[1])
        dyson.append(val)
    
    ntdab=int(lines[2*ncmo+nosht+npl])
    for y in range(ntdab):
        val=float(lines[2*ncmo+nosht+y+npl+1].split()[3])
        tdmab.append(val)
    
    f.close()
    #==================================

    # for equivalent atoms, like in the N2 molecule, where I have N1 and N2: OCA_center='N1 1s' to project on the N1 atom.
    # OCA_center='O 1s'
    #
    OCA_atom=OCA_center[:2] # takes the element from OCA_center
    n1=OCA_atom.replace(" ","") # these lines eliminate spaces and digits of equivalent atoms (if present)
    n1=n1.replace("1","")
    n1=n1.replace("2","")
    n1=n1.replace("3","")
    n1=n1.replace("4","")
    n1=n1.replace("5","")
    n1=n1.replace("6","")
    n1=n1.replace("7","")
    n1=n1.replace("8","")
    n1=n1.replace("9","")
    OCA_atom=n1.upper()

    return OCA_atom,OCA_center,OCA_line,benergy,totalSymmetry,symmetry,nbasf,nash,nmo,cmotab,tdmtab,ncmo,nbasft,nasht,nosht,comtaboff,comtbasoff,comtbasoff2,comtcmoff,comtnashoff,cmob,cmoa,tdmab,dyson

    #==================================

def init3(nbasft,symmetry):
    # The HDF5 file is saved in the project folder as $Project.rassi.h5
    # -----------
    # Here below is a way to get the $Project.rassi.h5 file and pass it to the code.
    # But before check if the file exist.
    logger = logging.getLogger('ftpuploader')
    try:
        pc1 = subprocess.Popen('find . -type f -name *.rassi.h5', stdout=subprocess.PIPE, shell=True)
        file_name=pc1.stdout.readlines()[0].decode(encoding='UTF-8',errors='strict')
        hd5_file=re.split('./|\n',file_name)[1]
        # the file $Project.rassi.h5 was just passed to the string 'hd5_file'.
    except BaseException as e:
        logger.error('The $Project.rassi.h5 file is missing in your current directory. Exit() ')
        sys.exit()
#
    h = h5py.File(hd5_file, 'r') # here I open the 'hd5_file' with h5py
    basis_id_hd5=np.reshape( h['BASIS_FUNCTION_IDS'],(nbasft,4))
    element=h['CENTER_LABELS']
    n_elements=len(element)
    element=np.reshape(element,(n_elements))
    element=np.array(element, dtype =np.str)
    element=np.char.strip(element)
    h.close()
    #print(basis_id_hd5)
    
    return hd5_file,basis_id_hd5,element

#==================================

