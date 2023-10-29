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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nInter, nIter
real(kind=wp), intent(inout) :: dq(nInter,nIter)
real(kind=wp), intent(in) :: H(nInter,nInter), g(nInter,nIter+1)
integer(kind=iwp) :: Info
integer(kind=iwp), allocatable :: Tmp(:)

! Compute a new independent geometry by relaxation of
! the gradient vector.

call mma_allocate(Tmp,nInter*nInter,Label='Tmp')

Tmp(:) = 0
dq(:,nIter) = Zero

call CG_Solver(nInter,nInter**2,H,Tmp,g(1,nIter),dq(1,nIter),Info)

call mma_deallocate(Tmp)

if (Info < 0) call SysAbendMsg('QNR','Conjugate gradients not converged',' ')

#ifdef _DEBUGPRINT_
write(u6,*) 'CG converged in ',Info,' iterations'
call RecPrt(' H ',' ',H,nInter,nInter)
call RecPrt(' dq',' ',dq(1,nIter),1,nInter)
#endif

return

end subroutine QNR
