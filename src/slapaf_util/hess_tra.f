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
      Subroutine Hess_Tra(Hss_X,nDim,Degen,BMx,nInter,Hss_Q)
      Implicit Real*8 (a-h,o-z)
      Real*8 Hss_X(nDim*nDim), Degen(nDim), BMx(nDim,nInter),
     &       Hss_Q(nInter*nInter)
#include "real.fh"
#include "WrkSpc.fh"
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the Hessian matrix in Internal Cartesian Coordinates
*     by solving:
*
*     d^2E/dx^2 = dQ/dx d^2E/dQ^2 dQ/dx
*
      Do j = 1, nDim
         Do i = 1, nDim
            ij = (j-1)*nDim + i
            Hss_X(ij) = Hss_X(ij) / Sqrt(Degen(i)*Degen(j))
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
      Call RecPrt('BMx',' ',BMx,nDim,nInter)
      Call RecPrt('Hss_X',' ',Hss_X,nDim,nDim)
#endif
      Call Allocate_Work(ipX,nDim*nInter)
      M=nDim
      N=nInter
      NRHS=nDim
      Call Eq_Solver('N',M,N,NRHS,BMx,.FALSE.,Degen,Hss_X,Work(ipX))
*
      Call Allocate_Work(ipXT,nDim*nInter)
      Call TRNSPS(nInter,nDim,Work(ipX),Work(ipXT))
#ifdef _DEBUGPRINT_
      Call RecPrt('Work(ipX)',' ',Work(ipX),nInter,nDim)
      Call RecPrt('Work(ipXT)',' ',Work(ipXT),nDim,nInter)
#endif
*
      M=nDim
      N=nInter
      NRHS=nInter
      Call Eq_Solver('N',M,N,NRHS,BMx,.False.,Degen,Work(ipXT),Hss_Q)
*
      Call Free_Work(ipXT)
      Call Free_Work(ipX)
#ifdef _DEBUGPRINT_
      Call RecPrt('Hss_Q',' ',Hss_Q,nInter,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
