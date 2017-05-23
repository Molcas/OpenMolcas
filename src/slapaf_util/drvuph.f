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
      Subroutine DrvUpH(nWndw,nIter,H,nInter,
     &                  dq,g,iOptH,HUpMet,nRowH,
     &                  jPrint,IterHess)
      Use NewH_mod
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 H(nInter,nInter), dq(nInter,nIter), g(nInter,nIter+1)
      Logical Found, DoMask
      Character*6 HUpMet
*                                                                      *
************************************************************************
*                                                                      *
      Lu=6
*                                                                      *
************************************************************************
*                                                                      *
      iSt=Max(2,nIter-(nWndw-1))
      Call Qpg_iScalar('HessIter',Found)
      If (Found) Then
         Call Get_iScalar('HessIter',IterHess)
         iSt=Max(iSt,IterHess+1)
      Else
         IterHess=0
      End If
      If (nRowH.GT.0) iSt=Max(iSt,nRowH+2)
      If (jPrint.ge.99) Then
         Write(Lu,*) 'DrvUpH: iSt,kIter=',iSt,nIter
         Call RecPrt('DrvUpH: Initial Hessian',' ',H,nInter,nInter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (jPrint.ge.6) Then
         Write (Lu,*)
         If (nIter.lt.iSt) Then
            Write (Lu,*) 'No update of Hessian on the first iteration'
         Else
            Write (Lu,'(A,30I3)') 'Hessian update from points:',
     &            (lIter,lIter=iSt-1,nIter)
         End If
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      DoMask=.False.
      If (Allocated(UpdMask)) Then
         If (Size(UpdMask).eq.nInter) DoMask=.True.
      End If
      If (DoMask) Then
         Do i=1,nInter
            If (UpdMask(i).ne.0) Then
               Do j=1,nInter
                  H(i,j)=Zero
                  H(j,i)=Zero
               End Do
               H(i,i)=DiagMM
            End If
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Update the Hessian over the window
*
      If (jPrint.ge.99)
     &   Call RecPrt('DrvUpH: Initial Hessian',' ',H,nInter,nInter)
      Do lIter=iSt,nIter
         If (jPrint.ge.99) Write(Lu,*)'DrvUpH: Call NewH, lIter=',lIter
         Call NewH(nInter,lIter,dq,g,H,iOptH,HUpMet,nIter)
      End Do
      If (jPrint.ge.99)
     &   Call RecPrt('DrvUpH: Updated Hessian',' ',H,nInter,nInter)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
