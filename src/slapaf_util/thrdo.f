************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine ThrdO(nInter,g,A,e,Fail)
************************************************************************
*                                                                      *
* Object: to find the error vector given the gradient, and the Hessian *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 g(nInter), A(nInter,nInter), e(nInter,2)
      Logical Fail
*
      i0=1
      i1=2
      iter=0
      iStep=0
      Fail=.True.
*
*---- Compile the error vector, starting in position 1
*     Observe that g is the force and NOT the gradient
*
      call dcopy_(nInter,g,1,e(1,i0),1)
      Call DPOTRS('U',nInter,1,A,nInter,e(1,i0),nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'ThrdO(DPOTRS): iRC=',iRC
         Call Abend()
      End If

      Call RecPrt(' ThrdO: e(0)',' ',e(1,i0),nInter,1)
*
      Thrd=1.D-6
      iterMx=40
*
800   Continue
      iStep=iStep+1
*
900   Continue
*
*
*---- Newton-Raphson scheme
*
      call dcopy_(nInter,g,1,e(1,i1),1)
      Call DPOTRS('U',nInter,1,A,nInter,e(1,i1),nInter,iRC)
      If (iRC.ne.0) Then
         Write (6,*) 'ThrdO(DPOTRS): iRC=',iRC
         Call Abend()
      End If
*
*
*     Call RecPrt(' ThrdO: e',' ',e(1,i1),nInter,1)
      iter=iter+1
*
*---- Check if the error vectors are self consistent.
*
      Test=Zero
      Do i = 1, nInter
         diff = Abs(e(i,i0)-e(i,i1))
         If (diff.gt.Test) Test=diff
      End Do
*     Write (*,*) iter,diff
*
      If (iter.gt.iterMx) Then
         Call WarningMessage(1,'ThrdO: Exceeded max iterations')
         Return
      End If
*
      If (Test.lt.Thrd) Then
*------- Copy converged vectors to slot 1
         If (i1.ne.1) call dcopy_(nInter,e(1,i1),1,e(1,1),1)
         If (iStep.eq.10) Then
            Call RecPrt(' ThrdO: e(Final)',' ',e,nInter,1)
            Fail=.False.
            Return
         Else
            iter=0
            Go To 800
         End If
      Else
*------- If not change vector slot
         itmp=i0
         i0=i1
         i1=itmp
         Go To 900
      End If
      Call WarningMessage(2,'Error in ThrdO')
      Call Abend()
      End
