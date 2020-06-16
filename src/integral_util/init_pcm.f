************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992,2000, Roland Lindh                                *
************************************************************************
      SubRoutine Init_PCM(NonEq,iCharg)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              Allok2                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*                                                                      *
*             Modified for Langevin polarizabilities, Marsk 2000 (RL)  *
************************************************************************
      use PCM_arrays
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "angstr.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "unixinfo.fh"
      Character*2 Elements(MxAtom*8)
      Logical NonEq
      Integer  iix(2)
      Real*8   rix(2)
#include "periodic_table.fh"
      Real*8, Allocatable:: Coor(:,:), LcCoor(:,:)
      Integer, Allocatable:: ANr(:), LcANr(:)
*
      If (.Not.PCM) Return
*
      iRout=1
      iPrint=nPrint(iRout)
      Call qEnter('Init_PCM')
*
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*                                                                      *
************************************************************************
*                                                                      *
*---- Reinitiate always for gradient calculations
*
      DoDeriv=.False.
cpcm_solvent
c added mckinley for pcm in second derivatives
      If (      ProgName.eq.'alaska'
     &     .or. ProgName.eq.'mckinley'
     &     .or. ProgName.eq.'mclr'    )
     &    DoDeriv=.True.
cpcm_solvent end
      If (DoDeriv) Then
         LcNAtm = ISlPar(42)
         nDeg=3*LcNAtm
         Call mma_allocate(dTes,nTs,lcNAtm,3,Label='dTes')
         Call mma_allocate(dPnt,nTs,lcNAtm,3,3,Label='dPnt')
         Call mma_allocate(dRad,nS ,lcNAtm,3,Label='dRad')
         Call mma_allocate(dCntr,nS ,lcNAtm,3,3,Label='dCntr')
         Call mma_allocate(PCM_SQ,2,nTs,Label='PCM_SQ')
         Call Get_dArray('PCM Charges',PCM_SQ,2*nTs)
         Go To 888
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Check if we have retrievable PCM data
*
      Call Get_iScalar('PCM info length',nPCM_info)
      If (nPCM_info.ne.0) Then
*
*------- If charge or NonEq/Eq status in retrieved data not the same as
*        for request redo the PCM initiation.
*
         If (iCharge_ref.ne.  iCharg) Go To 888
         If (NonEq_ref  .neqv.NonEq ) Go To 888
*
         Call Allocate_Work(ip_Sph,nPCM_info)
         Call Get_dArray('PCM Info',Work(ip_Sph),nPCM_info)
*
*------- Update the pointers
*
         mChunk=0
         Call Init_a_Chunk(ip_Sph,mChunk)
         Call Get_a_Chunk('PCMSph','Real',ip_Sph,4*NS)
         Call Get_a_Chunk('NOrd','Inte',ip_N ,NS)
         Call IZero(iWork(ip_N),NS)
         Call Get_a_Chunk('PCMTess','Real',ip_Tess,4*nTs)
         Call Get_a_Chunk('Vert','Real',ip_Vert,3*MxVert*nTs)
         Call Get_a_Chunk('Centr','Real',ip_Centr,3*MxVert*nTs)
         Call Get_a_Chunk('SSph','Real',ip_SSph,NS)
         Call Get_a_Chunk('ISph','Inte',ip_ISph,nTs)
         Call Get_a_Chunk('NVert','Inte',ip_NVert,nTs)
         Call Get_a_Chunk('IntSph','Inte',ip_IntS,MxVert*nTs)
         Call Get_a_Chunk('NewSph','Inte',ip_NewS,2*NS)
         Call Get_a_Chunk('DM','Real',ip_DM,nTs**2)
         Call nChunk(mChunk)
         If (mChunk.ne.nPCM_info) Then
            Call WarningMessage(2,'Init_PCM: mChunk.ne.nPCM_Info!')
            Write (6,*) 'mChunk=',mChunk
            Write (6,*) 'nPCM_Info=',nPCM_Info
            Call Abend()
         End If
         Go To 999
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Initial processing for PCM
*
 888  Call Get_nAtoms_All(nAtoms)
      Call mma_allocate(Coor,3,nAtoms,Label='Coor')
      Call Get_Coord_All(Coor,nAtoms)
      Call Get_Name_All(Elements)
      Call mma_Allocate(ANr,nAtoms,Label='ANr')
      Do i = 1, nAtoms
         Do j = 0, Num_Elem
            If (PTab(j).eq.Elements(i)) ANr(i)=j
         End Do
      End Do
      Call mma_allocate(LcCoor,3,nAtoms,Label='LcCoor')
      Call mma_allocate(LcANr,nAtoms,Label='LcANr')
*
*---- Initialize PCM model
*
*     iPrint: Print level
*     ICharg: Molecular charge
*     nAtoms: total number of atoms
*     angstr: conversion factor from bohr to Angstrom
*     Coor: Coordinates of atoms
*     ANr: atomic numbers
*     LcCoor: local array for atomic coordinates
*     LcANr: local array for atomic numbers
*     Solvent: string with explicit solvent name
*     Conductor: logical flag to activate conductor approximation
*     aArea: average area of a tessera
*     r_min_sphere: minimum radius of smoothing sphere
*     ip_Ts: pointer to tesserae
*     nTs  : number of tesserae
*
      Call PCM_Init(iPrint,ICharg,nAtoms,angstr,Coor,ANr,LcCoor,
     &              LcANr,nIrrep,NonEq)
      If (iPrint.gt.5) Then
         Write (6,*)
         Write (6,*)
      End If
*
      Call mma_deallocate(LcANr)
      Call mma_deallocate(LcCoor)
      Call mma_deallocate(ANr)
      Call mma_deallocate(Coor)
*                                                                      *
************************************************************************
*                                                                      *
*---- Put the dynamic arrays on COMFILE
*
      Call Save_PCM_Info(cRFStrt,iRFStrt,lRFStrt,rRFStrt)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
 999    Continue
       Call qExit('Init_PCM')
      Return
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine Save_PCM_Info(cRFStrt,iRFStrt,lRFStrt,rRFStrt)
      Use Iso_C_Binding
      Integer, Target :: cRFStrt,iRFStrt,lRFStrt
      Real*8, Target :: rRFStrt
      Integer, Pointer :: p_cRF(:),p_iRF(:),p_lRF(:)
      Real*8, Pointer :: p_rRF(:)
*
      Call Put_iScalar('PCM info length',nPCM_info)
      Call Put_dArray('PCM Info',Work(ip_Sph),nPCM_Info)
      iCharge_ref=iCharg
      NonEq_ref=NonEq
*
*---- Put the reaction field common blocks on disk again!
*
      Len = ilLoc(lRFEnd)-ilLoc(lRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(lRFStrt),p_lRF,[Len])
      Call Put_iArray('RFlInfo',p_lRF,Len)
*
      Len = idLoc(rRFEnd)-idLoc(rRFStrt)
      Len = (Len+nByte_r)/nByte_r
      Call C_F_Pointer(C_Loc(rRFStrt),p_rRF,[Len])
      Call Put_dArray('RFrInfo',p_rRF,Len)
*
      Len = iiLoc(iRFEnd)-iiLoc(iRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(iRFStrt),p_iRF,[Len])
      Call Put_iArray('RFiInfo',p_iRF,Len)
*
      Len = iiLoc(cRFEnd)-iiLoc(cRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(cRFStrt),p_cRF,[Len])
      Call Put_iArray('RFcInfo',p_cRF,Len)
*
      Nullify(p_lRF,p_rRF,p_iRF,p_cRF)
*
      End SubRoutine Save_PCM_Info
*
      End
