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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine ShiftHess(Hess,shift,nDim,nDim2)
!  Purpose:
!    Shifts Hessian to make it positive definite.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nDim2
real(kind=wp), intent(inout) :: Hess(nDim,nDim2)
logical(kind=iwp), intent(out) :: shift
integer(kind=iwp) :: i, j, k
real(kind=wp) :: eigen_min, eps
real(kind=wp), allocatable :: Hess_lowT(:), U(:,:)

call mma_allocate(U,nDim,nDim2,label='U')
call mma_allocate(Hess_lowT,nDim*(nDim+1)/2,label='Hess_lowT')

! Initialize.

k = 0
do i=1,nDim
  do j=1,i
    k = k+1
    Hess_lowT(k) = Hess(i,j)
  end do
end do
U(:,:) = Zero
do i=1,nDim
  U(i,i) = One
end do
call Jacob(Hess_lowT,U,nDim,nDim)
call Jacord(Hess_lowT,U,nDim,nDim)
eigen_min = Hess_lowT(1)
shift = eigen_min < Zero
if (shift) then
  eps = 2*eigen_min
  do i=1,nDim
    Hess(i,i) = Hess(i,i)-eps
  end do
end if
call mma_deallocate(U)
call mma_deallocate(Hess_lowT)

end subroutine ShiftHess
