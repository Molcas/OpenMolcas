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
      use Slapaf_Info, only: dqInt
      use Slapaf_Parameters, only: lOld
      Implicit Real*8 (a-h,o-z)
#include "info_slapaf.fh"
*
      If (lOld) Return
*
      Call Hss_q_(nsAtom,nQQ,dqInt(:,iRef))
*
      Return
      End
      Subroutine Hss_q_(nAtom,nQQ,Grad)
      use Slapaf_Info, only: Degen, Smmtrc
      use Slapaf_Parameters, only: Curvilinear, nDimBC, Analytic_Hessian
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Grad(nQQ)
      Real*8 rDum(1)
      Real*8, Allocatable:: Hss_X(:), Degen2(:), Hss_Q(:), KtB(:)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*---- Back-transform from cartesian to internals
*
*     dQ/dx d^2E/dQ^2 dQ/dx + d^2Q/dx^2 dE/dQ = d^2E/dx^2
*
*     Pickup d^2E/dx^2
*
      Call mma_allocate(Hss_x,nDimBC**2,Label='Hss_X')
      Call Get_dArray('Hss_X',Hss_x,nDimBC**2)
      Call mma_allocate(KtB,nDimBC*nQQ,Label='KtB')
      Call Get_dArray('KtB',KtB,nDimBC*nQQ)
#ifdef _DEBUGPRINT_
      Call RecPrt('Hss_x',' ',Hss_X,nDimBC,nDimBC)
#endif
*
      Call mma_allocate(Degen2,nDimBC,Label='Degen2')
      i=0
      Do ix = 1, 3*nAtom
         iAtom = (ix+2)/3
         ixyz = ix - (iAtom-1)*3
         If (Smmtrc(ixyz,iAtom)) Then
            i = i + 1
            Degen2(i) = Degen(ixyz,iAtom)
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Degen2',' ',Degen2,nDimBC,1)
#endif
*
      If (Analytic_Hessian.and.Curvilinear) Then
*
*        Form u^(1/2) (Sum(i) d^2Q_i/dx^2 * dE/dQ_i) u^(1/2)
*
*        and form d^2E/dx^2 - d^2Q/dx^2 dE/dQ
*
         Call dBuu(Degen2,nQQ,nDimBC,Grad,Hss_X,.False.)
#ifdef _DEBUGPRINT_
         Call RecPrt('H(X)-BtgQ',' ',Hss_X,nDimBC,nDimBC)
#endif
      End If
*
      Call mma_allocate(Hss_Q,nQQ**2,Label='Hss_Q')
      Call Hess_Tra(Hss_X,nDimBC,Degen2,KtB,nQQ,Hss_Q)
*
      Call Put_dArray('Hss_Q',Hss_Q,nQQ**2)
      Call Put_dArray('Hss_upd',rDum,0)
#ifdef _DEBUGPRINT_
      Call RecPrt('Hss_Q: Hessian',' ',Hss_Q,nQQ,nQQ)
#endif
      Call mma_deallocate(Hss_Q)
      Call mma_deallocate(KtB)
      Call mma_deallocate(Degen2)
      Call mma_deallocate(Hss_X)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
