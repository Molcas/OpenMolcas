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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetU_ER(U,R,n)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: compute U = R*[R^T*R]^(-1/2).
!
! (used by ER orbital localisation - hence the _ER)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(out) :: U(n,n)
real(kind=wp), intent(in) :: R(n,n)
integer(kind=iwp) :: iTask, lScr, n2, nn
real(kind=wp), allocatable :: ISqrt(:,:), RTR(:,:), Scr(:), Sqroot(:,:)

if (n < 1) return

! Allocations.
! ------------

nn = n*(n+1)/2
n2 = n**2
lScr = 2*n2+nn
call mma_allocate(RTR,n,n,label='RTR')
call mma_allocate(Sqroot,n,n,label='Sqrt')
call mma_allocate(ISqrt,n,n,label='ISqrt')
call mma_allocate(Scr,lScr,label='Scr')

! Compute R^T*R.
! --------------

call DGEMM_('T','N',n,n,n,One,R,n,R,n,Zero,RTR,n)

! Compute inverse square root of R^T*R.
! -------------------------------------

iTask = 2 ! compute sqrt as well as inverse sqrt
call SqrtMt(RTR,n,iTask,Sqroot,ISqrt,Scr)

! Compute U.
! ----------

call DGEMM_('N','N',n,n,n,One,R,n,ISqrt,n,Zero,U,n)

! De-allocations.
! ---------------

call mma_deallocate(RTR)
call mma_deallocate(Sqroot)
call mma_deallocate(ISqrt)
call mma_deallocate(Scr)

end subroutine GetU_ER
