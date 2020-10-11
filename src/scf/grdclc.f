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
      SubRoutine GrdClc(What,QNR)
      use SCF_Arrays
      Implicit Real*8 (a-h,o-z)
      Character What*3
      Logical   QNR
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "WrkSpc.fh"
*
      nD = iUHF + 1
*                                                                      *
************************************************************************
*                                                                      *
*     QNR     : Quasi-Newton in effect
*
      If (QNR) Then
         Call GrdClc_(What,Dens,TwoHam,Vxc,nBT,nDens,nD,OneHam,
     &                CMO   ,nBB,Ovrlp,CMO)
      Else
         Call GrdClc_(What,Dens,TwoHam,Vxc,nBT,nDens,nD,OneHam,
     &                Lowdin,nBB,Ovrlp,CMO)
      End If
*
      Return
      End
      SubRoutine GrdClc_(What,Dens,TwoHam,Vxc,mBT,mDens,nD,OneHam,
     &                   OCMO,mBB,Ovrlp,CMO)
************************************************************************
*                                                                      *
*     purpose: Compute gradients and write on disk.                    *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       What    : variable telling what gradients compute: 'All' -     *
*                 all gradients, 'Lst' - last gradient                 *
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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "stdalloc.fh"
#include "file.fh"
#include "llists.fh"
*
      Real*8 Dens(mBT,nD,mDens), TwoHam(mBT,nD,mDens), CMO(mBB,nD),
     &       OneHam(mBT), OCMO(mBB,nD), Ovrlp(mBT), Vxc(mBT,nD,mDens)
      Real*8, Dimension(:,:), Allocatable:: GrdOO,GrdOV,AuxD,AuxT,AuxV
      Character What*3
#include "interfaces_scf.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*define _DEBUGPRINT_
*
      If (What.ne.'All' .and. What.ne.'Lst') Then
         Write (6,*) 'GrdClc: What.ne."All" .and. What.ne."Lst"'
         Write (6,'(A,A)') 'What=',What
         Call Abend()
      End If

*--- Allocate memory for gradients and gradient contributions
      Call mma_allocate(GrdOO,nOO,nD,Label='GrdOO')
      Call mma_allocate(GrdOV,nOO,nD,Label='GrdOV')

*--- Allocate memory for auxiliary matrices
      Call mma_allocate(AuxD,nBT,nD,Label='AuxD')
      Call mma_allocate(AuxT,nBT,nD,Label='AuxT')
      Call mma_allocate(AuxV,nBT,nD,Label='AuxV')

*--- Find the beginning of the loop
      If (What.eq.'All') Then
         LpStrt = 1
      Else
         LpStrt = kOptim
      End If

*--- Compute all gradients / last gradient
*
      iter_d=iter-iter0
      Do iOpt = LpStrt, kOptim
         iDT = iter_d - kOptim + iOpt
*
         Call dCopy_(nOV*nD,[Zero],0,GrdOV,1)
*
         jDT=MapDns(iDT)
         If (jDT.lt.0) Then
            Call RWDTG(-jDT,AuxD,nBT*nD,'R','DENS  ',iDisk,MxDDsk)
            Call RWDTG(-jDT,AuxT,nBT*nD,'R','TWOHAM',iDisk,MxDDsk)
            Call RWDTG(-jDT,AuxV,nBT*nD,'R','dVxcdR',iDisk,MxDDsk)
*
            Call EGrad(OneHam,AuxT,AuxV,Ovrlp,AuxD,nBT,OCMO,nBO,
     &                 GrdOO,nOO,nD,CMO)
*
         Else
*
            Call EGrad(OneHam,TwoHam(1,1,jDT),Vxc(1,1,jDT),Ovrlp,
     &                 Dens(1,1,jDT),nBT,OCMO,nBO,GrdOO,nOO,nD,CMO)
*
         End If
*
         Call vOO2OV(GrdOO,nOO,GrdOV,nOV,nD)
*
*------- Write Gradient to linked list
*
         Call PutVec(GrdOV,nD*nOV,LuGrd,iDT+iter0,MemRsv,'OVWR',LLGrad)
*
#ifdef _DEBUGPRINT_
         Write (6,*) 'GrdClc: Put Gradient iteration:',iDT+iter0
         Write (6,*) 'iOpt=',iOpt
         Call NrmClc(GrdOO,nOO*nD,'GrdClc','GrdOO')
         Call NrmClc(GrdOV,nOV*nD,'GrdClc','GrdOV')
         Call RecPrt('GrdClc: g(i)',' ',GrdOV,1,nOV*nD)
#endif
      End Do
*
*     Deallocate memory
*
      Call mma_deallocate(AuxD)
      Call mma_deallocate(AuxT)
      Call mma_deallocate(AuxV)
      Call mma_deallocate(GrdOV)
      Call mma_deallocate(GrdOO)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
