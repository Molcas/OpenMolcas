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
      Subroutine RowHessian(nIter,nInter,nRowH,mRowH,Delta,dq,q,g)
      Implicit Real*8 (A-H,O-Z)
      Real*8 dq(nInter,nIter), g(nInter,nIter),q(nInter,nIter+1)
      Integer mRowH(10)
#include "WrkSpc.fh"
*
      Call Allocate_Work(ipH,nInter**2)
      Call Get_dArray('Hss_Q',Work(ipH),nInter**2)
      Call Put_dArray('Hss_upd',Work(ip_Dummy),0)
*
      Call RowHessian_(nIter,nInter,nRowH,mRowH,Delta,Work(ipH),dq,q,g)
      Write (6,*)
      Write (6,*) ' Partial numerical differentiation is finished!'
*
      Call Put_dArray('Hss_Q',Work(ipH),nInter**2)
      Call Free_Work(ipH)
*
      Return
      End
      Subroutine RowHessian_(nIter,nInter,nRowH,mRowH,Delta,H,dq,q,g)
************************************************************************
*                                                                      *
* Object: Numerical estimation of single rows and columns of Hessian   *
* Called from: RlxCtl when lRowH=.True. & iter.EQ.NmIter               *
* Author: Giovanni Ghigo, University of Torino, Italy                  *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 dq(nInter,nIter), g(nInter,nIter),H(nInter,nInter),
     &       q(nInter,nIter+1)
      Integer mRowH(10)
#include "print.fh"
#include "real.fh"
*
      iRout = 184
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) then
         Write(6,*)
         Write(6,*) 'RowHessian:'
         Call RecPrt('Initial Hessian',' ',H,nInter,nInter)
         Call RecPrt('Displacement dq','(10F9.6)',dq,nInter,nIter)
         Call RecPrt('Coordinates   q:','(10F9.6)', q,nInter,nIter)
         Call RecPrt('Gradient      g:','(10F9.6)', g,nInter,nIter)
         call XFlush (6)
      EndIf
*
* --- Evaluate the Hessian
*
      Do iRowH = 1, nRowH
         iInter = mRowH(iRowH)
         Do jInter = 1, nInter
            H(iInter,jInter) = (g(jInter,1) - g(jInter,iRowH+1)) / Delta
            H(jInter,iInter) = H(iInter,jInter)
         End Do
      EndDo
      If (iPrint.ge.98) then
         Call RecPrt('Final Hessian',' ',H,nInter,nInter)
         call XFlush (6)
      EndIf
*
* --- Symmetrize
*
      Do iInter = 1, nInter
         Do jInter = 1, nInter
            dElement = (H(iInter,jInter)+H(jInter,iInter))/2.0d0
            H(iInter,jInter) = dElement
            H(jInter,iInter) = dElement
         EndDo
      EndDo

      Return
      End
