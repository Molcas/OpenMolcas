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
      Subroutine Hss_q()
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "info_slapaf.fh"
*
      If (lOld) Return
*
      ipdqInt_ref = (iRef-1)*nQQ + ipdqInt
      Call Hss_q_(Degen,nsAtom,nQQ,Smmtrc,Analytic_Hessian,
     &            Work(ipdqInt_ref),nDimBC,Curvilinear)
*
      Return
      End
      Subroutine Hss_q_(Degen,nAtom,nQQ,Smmtrc,Analytic_Hessian,Grad,
     &                 nDim,Curvilinear)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Degen(3*nAtom),Grad(nQQ)
      Logical Smmtrc(3*nAtom), Analytic_Hessian, Curvilinear
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*---- Back-transform from cartesian to internals
*
*     dQ/dx d^2E/dQ^2 dQ/dx + d^2Q/dx^2 dE/dQ = d^2E/dx^2
*
*     Pickup d^2E/dx^2
*
      Call Allocate_Work(ip_Hss_x,nDim**2)
      Call Get_dArray('Hss_X',Work(ip_Hss_x),nDim**2)
      Call Allocate_Work(ip_KtB,nDim*nQQ)
      Call Get_dArray('KtB',Work(ip_KtB),nDim*nQQ)
#ifdef _DEBUGPRINT_
      Call RecPrt('Hss_x',' ',Work(ip_Hss_X),nDim,nDim)
#endif
*
      Call Allocate_Work(ipDegen,nDim)
      i=0
      Do ix = 1, 3*nAtom
         If (Smmtrc(ix)) Then
            Work(ipDegen+i) = Degen(ix)
            i = i + 1
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Work(ipDegen)',' ',Work(ipDegen),nDim,1)
#endif
*
      If (Analytic_Hessian.and.Curvilinear) Then
*
*        Form u^(1/2) (Sum(i) d^2Q_i/dx^2 * dE/dQ_i) u^(1/2)
*
*        and form d^2E/dx^2 - d^2Q/dx^2 dE/dQ
*
         Call dBuu(Work(ipDegen),nQQ,nDim,Grad,Work(ip_Hss_X),.False.)
#ifdef _DEBUGPRINT_
         Call RecPrt('H(X)-BtgQ',' ',Work(ip_Hss_X),nDim,nDim)
#endif
      End If
*
      Call Allocate_Work(ip_Hss_Q,nQQ**2)
      Call Hess_Tra(Work(ip_Hss_X),nDim,Work(ipDegen),
     &              Work(ip_KtB),nQQ,Work(ip_Hss_Q))
*
      Call Put_dArray('Hss_Q',Work(ip_Hss_Q),nQQ**2)
      Call Put_dArray('Hss_upd',Work(ip_Dummy),0)
#ifdef _DEBUGPRINT_
      Call RecPrt('Hss_Q: Hessian',' ',Work(ip_Hss_Q),nQQ,nQQ)
#endif
      Call Free_Work(ip_Hss_Q)
      Call Free_Work(ipDegen)
      Call Free_Work(ip_KtB)
      Call Free_Work(ip_Hss_X)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
