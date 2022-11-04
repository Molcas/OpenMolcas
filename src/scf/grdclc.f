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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2016,2017, Roland Lindh                                *
************************************************************************
      SubRoutine GrdClc(Do_All)
************************************************************************
*                                                                      *
*     purpose: Compute gradients and write on disk.                    *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       Do_All  : variable telling what gradients compute: .true. -    *
*                 all gradients, .False. - last gradient               *
*                                                                      *
*     called from: Wfctl_scf                                           *
*                                                                      *
*     calls to:         EGrad                                          *
*               uses SubRoutines and Functions from Module lnklst.f    *
*               -linked list implementation to store series of vectors *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*#define _DEBUGPRINT_
      Use Interfaces_SCF, Only: vOO2OV
      Use InfSCF, only: Iter, Iter_Start, kOV, mOV, nBO,
     &                  nBT, nOO, iUHF
      use LnkLst, only: LLGrad
      use SCF_Arrays, Only: OneHam, CMO_Ref, Ovrlp, FockMO
      use Constants, only: Zero
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
#include "file.fh"
*
      Logical Do_All

! Local variables
      Real*8, Allocatable:: GrdOV(:)
      Integer nD, iOpt, LpStrt
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      nD = iUHF + 1
*
*--- Allocate memory for gradients
      Call mma_allocate(GrdOV,mOV,Label='GrdOV')

*--- Find the beginning of the loop
      If (Do_All) Then
         LpStrt = Iter_Start
         Do_All=.False.
      Else
         LpStrt = Iter
      End If

*--- Compute all gradients / last gradient
*
      Do iOpt = LpStrt, Iter
*
         GrdOV(:)=Zero
*
         Call EGrad(OneHam,Ovrlp,nBT,CMO_Ref,nBO,
     &                 FockMO,nOO,nD,CMO_Ref,iOpt)
*
         Call vOO2OV(FockMO,nOO,GrdOV,mOV,nD,kOV)
*
*------- Write Gradient to linked list
*
         Call PutVec(GrdOV,mOV,iOpt,'OVWR',LLGrad)
*
#ifdef _DEBUGPRINT_
         Write (6,*) 'GrdClc: Put Gradient iteration:',iOpt
         Write (6,*) 'iOpt,mOV=',iOpt,mOV
         Call NrmClc(FockMO,nOO*nD,'GrdClc','FockMO')
         Call NrmClc(GrdOV,mOV,'GrdClc','GrdOV')
#endif
      End Do
*
*     Deallocate memory
*
      Call mma_deallocate(GrdOV)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
