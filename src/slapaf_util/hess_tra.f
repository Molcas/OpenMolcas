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
#include "stdalloc.fh"
      Real*8, Allocatable:: X(:), XT(:)
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
      Call mma_allocate(X,nDim*nInter,Label='X')
      M=nDim
      N=nInter
      NRHS=nDim
      Call Eq_Solver('N',M,N,NRHS,BMx,.FALSE.,Degen,Hss_X,X)
*
      Call mma_allocate(XT,nDim*nInter,Label='XT')
      Call TRNSPS(nInter,nDim,X,XT)
#ifdef _DEBUGPRINT_
      Call RecPrt('X',' ',X,nInter,nDim)
      Call RecPrt('XT',' ',XT,nDim,nInter)
#endif
*
      M=nDim
      N=nInter
      NRHS=nInter
      Call Eq_Solver('N',M,N,NRHS,BMx,.False.,Degen,XT,Hss_Q)
*
      Call mma_deallocate(XT)
      Call mma_deallocate(X)
#ifdef _DEBUGPRINT_
      Call RecPrt('Hss_Q',' ',Hss_Q,nInter,nInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
