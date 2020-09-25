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
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*                                                                      *
*             Modified for Langevin polarizabilities, March 2000 (RL)  *
************************************************************************
      use PCM_arrays
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "angstr.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
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
*        Evolving the new code
*
         Call mma_allocate(PCMSph,4,NS,Label='PCMSph')
         Call mma_allocate(PCMTess,4,nTs,Label='PCMTess')
         Call mma_allocate(Vert,3,MxVert,nTs,Label='Vert')
         Call mma_allocate(Centr,3,MxVert,nTs,Label='Centr')
         Call mma_allocate(SSph,NS,Label='SSph')
         Call mma_allocate(PCMDM,nTs,nTs,Label='PCMDM')
         Call mma_allocate(PCM_N,NS,Label='PCM_N')
         Call mma_allocate(PCMiSph,nTs,Label='PCMiSph')
         Call mma_allocate(NVert,nTs,Label='NVert')
         Call mma_allocate(IntSph,MxVert,nTs,Label='IntSph')
         Call mma_allocate(NewSph,2,NS,Label='NewSph')
*
         Call Get_dArray('PCMSph',PCMSph,4*NS)
         Call Get_dArray('PCMTess',PCMTess,4*nTs)
         Call Get_dArray('Vert',Vert,3*MxVert*nTs)
         Call Get_dArray('Centr',Centr,3*MxVert*nTs)
         Call Get_dArray('SSph',SSph,NS)
         Call Get_dArray('PCMDM',PCMDM,nTs**2)
         Call Get_iArray('PCM_N',PCM_N,NS)
         Call Get_iArray('PCMiSph',PCMiSph,nTs)
         Call Get_iArray('NVert',NVert,nTs)
         Call Get_iArray('IntSph',IntSph,MxVert*nTs)
         Call Get_iArray('NewSph',NewSph,2*NS)

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
*    MxVert*nTs ANr: atomic numbers
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

      Call Put_dArray('PCMSph',PCMSph,4*NS)
      Call Put_dArray('PCMTess',PCMTess,4*nTs)
      Call Put_dArray('Vert',Vert,3*MxVert*nTs)
      Call Put_dArray('Centr',Centr,3*MxVert*nTs)
      Call Put_dArray('SSph',SSph,NS)
      Call Put_dArray('PCMDM',PCMDM,nTs**2)
      Call Put_iArray('PCM_N',PCM_N,NS)
      Call Put_iArray('PCMiSph',PCMiSph,nTs)
      Call Put_iArray('NVert',NVert,nTs)
      Call Put_iArray('IntSph',IntSph,MxVert*nTs)
      Call Put_iArray('NewSph',NewSph,2*NS)
*
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
