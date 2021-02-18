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
* Copyright (C) Giovanni Ghigo                                         *
************************************************************************
      Subroutine Freq1()
      use Slapaf_Info, only: qInt
      use Slapaf_parameters, only: Delta, iter

      nInter=SIZE(qInt,1)
      Call Freq1_Internal(iter,nInter,Delta/2.5d0,qInt)

      Contains
      Subroutine Freq1_Internal(nIter,nInter,Delta,qInt)
************************************************************************
*                                                                      *
* Object: Displacements for Numerical estimation of single rows and    *
*         columns of Hessian                                           *
*                                                                      *
* Called from: RlxCtl when Aloocated(mRowH)=.True.                     *
*                                                                      *
* Author: Giovanni Ghigo, University of Torino, Italy                  *
************************************************************************
      use Slapaf_Info, only: Shift, mRowH
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 qInt(nInter,nIter+1)
*
      iRout = 183
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Write(6,*)' [Freq1] nInter, nIter, Delta =',nInter,nIter,Delta
         Call RecPrt('Initial Shift:','(10F9.6)',Shift,nInter,nIter)
         Call RecPrt('Initial qInt:','(10F9.6)', qInt,nInter,nIter+1)
      EndIf
*
*-----Compute the new shift
*
      call dcopy_(nInter,[Zero],0,Shift(1,nIter),1)
      nRowH=0
      If (Allocated(mRowH)) nRowH=SIZE(mRowH)
      If (nIter.le.nRowH) then
         kInter = mRowH(nIter)
         Shift(kInter,nIter) = Delta
      EndIf
      If (nIter.gt.1) then
         jInter = mRowH(nIter-1)
         Shift(jInter,nIter) = -Delta ! Undo previous displacement
      End If
*
*---- Compute the new parameter set.
*
      call dcopy_(nInter,qInt(1,nIter),1,qInt(1,nIter+1),1)
      Call DaXpY_(nInter,One,Shift(1,nIter),1,qInt(1,nIter+1),1)
*
      If (iPrint.gt.5) Then
         Write (6,*)
     &           ' Accumulate the gradient for yet one parameter set'
         Write (6,*)
      End If
*
      If (iPrint.ge.98) Then
         Write(6,*)' [Freq1] nInter, nIter, Delta =',nInter,nIter,Delta
         Call RecPrt('Final Shift:','(10F9.6)',Shift,nInter,nIter)
         Call RecPrt('Final  q:','(10F9.6)', qInt,nInter,nIter+1)
      EndIf
      End Subroutine Freq1_Internal


      End Subroutine Freq1
