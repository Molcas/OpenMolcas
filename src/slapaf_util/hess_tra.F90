!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Hess_Tra(Hss_X,nDim,Degen,BMx,nInter,Hss_Q)

implicit real*8(a-h,o-z)
real*8 Hss_X(nDim*nDim), Degen(nDim), BMx(nDim,nInter), Hss_Q(nInter*nInter)
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: X(:), XT(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the Hessian matrix in Internal Cartesian Coordinates by solving:
!
! d^2E/dx^2 = dQ/dx d^2E/dQ^2 dQ/dx

do j=1,nDim
  do i=1,nDim
    ij = (j-1)*nDim+i
    Hss_X(ij) = Hss_X(ij)/sqrt(Degen(i)*Degen(j))
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('BMx',' ',BMx,nDim,nInter)
call RecPrt('Hss_X',' ',Hss_X,nDim,nDim)
#endif
call mma_allocate(X,nDim*nInter,Label='X')
M = nDim
N = nInter
NRHS = nDim
call Eq_Solver('N',M,N,NRHS,BMx,.false.,Degen,Hss_X,X)

call mma_allocate(XT,nDim*nInter,Label='XT')
call TRNSPS(nInter,nDim,X,XT)
#ifdef _DEBUGPRINT_
call RecPrt('X',' ',X,nInter,nDim)
call RecPrt('XT',' ',XT,nDim,nInter)
#endif

M = nDim
N = nInter
NRHS = nInter
call Eq_Solver('N',M,N,NRHS,BMx,.false.,Degen,XT,Hss_Q)

call mma_deallocate(XT)
call mma_deallocate(X)
#ifdef _DEBUGPRINT_
call RecPrt('Hss_Q',' ',Hss_Q,nInter,nInter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Hess_Tra
