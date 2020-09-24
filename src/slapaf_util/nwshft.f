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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine NwShft(dq,nInter,g,nIter,Delta,q)
************************************************************************
*                                                                      *
* Object: to numerically evaluate the molecular Hessian.               *
*                                                                      *
* Called from: RlxCtl                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             May '92                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 dq(nInter,nIter), g(nInter,nIter),q(nInter,nIter+1)
#include "print.fh"
#include "real.fh"
*
      iRout = 182
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt('NwShft  g',' ', g,nInter,nIter)
         Call RecPrt('NwShft  q',' ', q,nInter,nIter)
         Call RecPrt('NwShft dq',' ',dq,nInter,nIter-1)
      End If
*
*
*-----Compute the new shift
*
*     Write (*,*) ' nIter=',nIter
      If (nIter.lt.2*nInter+1) Then
*
*------- Shifts for the numerical Hessian
*
         jInter = (nIter+1)/2
         call dcopy_(nInter,[Zero],0,dq(1,nIter),1)
         If (Mod(nIter,2).eq.0) Then
            dq(jInter,nIter) = -Two*Delta
         Else
*---------- Undo previous displacement
            If (jInter.gt.1) dq(jInter-1,nIter) = Delta
            dq(jInter,nIter) = Delta
         End If
*
      Else
*
*------- Shifts for the numerical cubic force constants
*
         iCount=(nIter-2*nInter+3)/4
*        Write (*,*) ' iCount=',iCount
         jCount=0
         lInter=0
         Do kInter = 1, nInter
            Do lInter = 1, kInter-1
               jCount = jCount + 1
               If (jCount.eq.iCount) Go To 777
            End Do
         End Do
 777     Continue
         If (lInter.eq.0) Then
            Call WarningMessage(2,'lInter.eq.0')
            Call Abend()
         End If
*        Write (*,*) 'kInter, lInter=',kInter,lInter
         kCount=nIter-2*nInter
         call dcopy_(nInter,[Zero],0,dq(1,nIter),1)
*------- Undo last change for numerical Hessian
         If (iCount.eq.1) dq(nInter,nIter)=Delta
         If (Mod(kCount,4).eq.1) Then
*---------- Undo change due to previous pair
            If (lInter.ne.1) Then
               dq(kInter,nIter)=Delta
               dq(lInter-1,nIter)=Delta
            Else If (lInter.eq.1.and.kInter.ne.2) Then
               dq(kInter-1,nIter)=Delta
               dq(kInter-2,nIter)=Delta
            End If
*---------- +d,+d
*           Write (*,*) ' +d,+d'
            dq(kInter,nIter) = dq(kInter,nIter)+Delta
            dq(lInter,nIter) = dq(lInter,nIter)+Delta
         Else If (Mod(kCount,4).eq.2) Then
*---------- -d,+d
*           Write (*,*) ' -d,+d'
            dq(kInter,nIter) = -Two*Delta
            dq(lInter,nIter) = Zero
         Else If (Mod(kCount,4).eq.3) Then
*---------- +d,-d
*           Write (*,*) ' +d,-d'
            dq(kInter,nIter) = Two*Delta
            dq(lInter,nIter) = -Two*Delta
         Else If (Mod(kCount,4).eq.0) Then
*---------- -d,-d
*           Write (*,*) ' -d,-d'
            dq(kInter,nIter) = -Two*Delta
            dq(lInter,nIter) = Zero
         End If
      End If
*
*---- Compute the new parameter set.
      call dcopy_(nInter,q(1,nIter),1,q(1,nIter+1),1)
      Call DaXpY_(nInter,One,dq(1,nIter),1,q(1,nIter+1),1)
*
      If (iPrint.ge.99) Then
         Call RecPrt('  q',' ', q,nInter,nIter+1)
         Call RecPrt(' dq',' ',dq,nInter,nIter)
      End If
      Return
      End
