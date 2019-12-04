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
* Copyright (C) 1995, Roland Lindh                                     *
************************************************************************
      SubRoutine NewH(nInter,nIter,dq_orig,g,H,iOptH,HUpMet,mIter)
************************************************************************
*                                                                      *
* Object: Driver for inverse Hessian update.                           *
*                                                                      *
* Called from: RlxCtl                                                  *
*                                                                      *
* Calling:     QEnter                                                  *
*              DCopy   (ESSL)                                          *
*              UpdHin                                                  *
*              View                                                    *
*              Minv                                                    *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '95                                              *
************************************************************************
      Use NewH_mod
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer nInter, nIter, mIter, iOptH, i, iPrint, iRout, iSing
      Real*8 dq_orig(nInter,nIter), g(nInter,mIter+1), H(nInter,nInter)
      Character*6 HUpMet
      Logical Test, DoMask
      Real*8, Dimension(:), Allocatable :: dg, gi
      Real*8, Dimension(:,:), Allocatable :: dq
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      Test(i)=iAnd(iOptH,2**(i-1)).eq.2**(i-1)
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('NewH')
      iRout = 112
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Write (6,*)
         Write (6,*) ' NewH: lIter=',nIter
         Call RecPrt(' NewH: dq_orig',' ',dq_orig,nInter,nIter)
         Call RecPrt(' NewH: g',' ',g,nInter,nIter)
         Call RecPrt(' NewH: H(Old)',' ',H,nInter,nInter)
         Write(6,*)' NewH: Test(i)==',(Test(i),i=1,8)
      End If
*
*     Branch out if the first iteration
*
      If (nIter.le.1) Go To 999
*                                                                      *
************************************************************************
*                                                                      *
*-----Get the MM mask
      DoMask=.False.
      If (Allocated(UpdMask)) Then
         If (Size(UpdMask).eq.nInter) DoMask=.True.
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(dg,nInter,label="dg")
      Call mma_allocate(gi,nInter,label="gi")
      Call mma_allocate(dq,nInter,nIter,label="dq")
      call dcopy_(nInter*nIter,dq_orig,1,dq,1)
*                                                                      *
************************************************************************
*                                                                      *
*-----Form the difference between the gradients. Observe the order!
*     This since we are storing the forces rather than the gradients,
*     big mistake!!!
*
      Do i = 1, nInter
         dg(i) = g(i,nIter-1) - g(i,nIter)
         If (DoMask) Then
            If (UpdMask(i).ne.0) Then
               dg(i) = Zero
               dq(i,nIter-1) = Zero
            End If
         End If
      End Do
      If (iPrint.ge.99) Call RecPrt(' NewH: dg',' ',dg,nInter,1)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the update
*
      iSing=0
      If (Test(4)) Then
*
*------- No update
*
         HUpMet=' None '
*
      Else If (Test(1)) Then
*
*------ Fletcher (or Meyer) update
*
         HUpMet='  F   '
         Write (6,*) 'Deleted option in NewH'
         Call Abend()
*
      Else If (Test(2)) Then
*
*------- Broyden-Powel Symmetric Rank-2 update
*
         HUpMet='  BP  '
         Write (6,*) 'Deleted option in NewH'
         Call Abend()
*
      Else If (Test(3)) Then
*
*------- Broyden-Fletcher-Goldfarb-Shanno update
*
         HUpMet=' BFGS '
         Call  DFP(H,nInter,gi,dq(1,nIter-1),dg)
*
      Else If (Test(5)) Then
*
*------- Murtagh-Sargent-Powell update
*
         HUpMet=' MSP  '
         Call dGeMV_('N',nInter,nInter,
     &              -One,H,nInter,
     &              dq(1,nIter-1),1,
     &              One,dg,1)
         If (iPrint.ge.99) Call RecPrt(' NewH: gamma',' ',dg,nInter,1)
         Call MSP(H,gi,dg,dq(1,nIter-1),nInter)
*
      Else If (Test(7)) Then
*
*------- EU update
*
         HUpMet='  EU  '
*
*        Some precalculations:
         Do i=1,nInter
            gi(i) = -g(i,nIter-1)
            If (DoMask) Then
               If (UpdMask(i).ne.0) gi(i) = Zero
            End If
         End Do
*
         Call EU(dq(1,nIter-1),dg,gi,H,nInter)
*
      Else If (Test(8)) Then
*
*------- TS_BFGS update
*
         HUpMet='TSBFGS'
*
*        Some precalculations:
         Do i=1,nInter
            gi(i) = -g(i,nIter-1)
            If (DoMask) Then
               If (UpdMask(i).ne.0) gi(i) = Zero
            End If
         End Do
*
         Call TS_BFGS(dq(1,nIter-1),dg,gi,H,nInter)
*
      Else
*
         Call WarningMessage(2,'Error in NewH')
         Write (6,*) ' Improper value of iOptH:',iOptH
         Call Abend()
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.99) Then
         Call RecPrt(' NewH:  H(New)',' ',H,nInter,nInter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(dg)
      Call mma_deallocate(gi)
      Call mma_deallocate(dq)
*                                                                      *
************************************************************************
*                                                                      *
 999  Continue
      Call QExit('NewH')
      Return
      End Subroutine NewH
