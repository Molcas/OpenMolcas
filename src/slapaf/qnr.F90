!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

subroutine QNR(nInter,nIter,dq,H,g)
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
real*8 dq(nInter,nIter), H(nInter,nInter), g(nInter,nIter+1)
integer, allocatable :: Tmp(:)

! Compute a new independent geometry by relaxation of
! the gradient vector.

call mma_allocate(Tmp,nInter*nInter,Label='Tmp')

Tmp(:) = 0
dq(:,nIter) = Zero

call CG_Solver(nInter,nInter**2,H,Tmp,g(1,nIter),dq(1,nIter),Info)

call mma_deallocate(Tmp)

if (Info < 0) then
  call SysAbendMsg('QNR','Conjugate gradients not converged',' ')
end if

#ifdef _DEBUGPRINT_
write(6,*) 'CG converged in ',Info,' iterations'
call RecPrt(' H ',' ',H,nInter,nInter)
call RecPrt(' dq',' ',dq(1,nIter),1,nInter)
#endif

return

end subroutine QNR
