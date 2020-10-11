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
      Subroutine Freq1(nIter,nInter,nRowH,mRowH,Delta,dq,q)
************************************************************************
*                                                                      *
* Object: Displacements for Numerical estimation of single rows and    *
*         columns of Hessian                                           *
*                                                                      *
* Called from: RlxCtl when lRowH=.True.                                *
*                                                                      *
* Author: Giovanni Ghigo, University of Torino, Italy                  *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 q(nInter,nIter+1), dq(nInter,nIter)
      Integer mRowH(10)
*
      iRout = 183
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Write(6,*)' [Freq1] nInter, nIter, Delta =',nInter,nIter,Delta
         Call RecPrt('Initial dq:','(10F9.6)',dq,nInter,nIter)
         Call RecPrt('Initial  q:','(10F9.6)', q,nInter,nIter+1)
      EndIf
*
*-----Compute the new shift
*
      call dcopy_(nInter,[Zero],0,dq(1,nIter),1)
      If (nIter.le.nRowH) then
         kInter = mRowH(nIter)
         dq(kInter,nIter) = Delta
      EndIf
      If (nIter.gt.1) then
         jInter = mRowH(nIter-1)
         dq(jInter,nIter) = -Delta ! Undo previous displacement
      End If
*
*---- Compute the new parameter set.
*
      call dcopy_(nInter,q(1,nIter),1,q(1,nIter+1),1)
      Call DaXpY_(nInter,One,dq(1,nIter),1,q(1,nIter+1),1)
*
      If (iPrint.gt.5) Then
         Write (6,*)
     &           ' Accumulate the gradient for yet one parameter set'
         Write (6,*)
      End If
*
      If (iPrint.ge.98) Then
         Write(6,*)' [Freq1] nInter, nIter, Delta =',nInter,nIter,Delta
         Call RecPrt('Final dq:','(10F9.6)',dq,nInter,nIter)
         Call RecPrt('Final  q:','(10F9.6)', q,nInter,nIter+1)
      EndIf
      Return
      End
