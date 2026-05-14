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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine nize_cvb(c,nnrm,s,n,metr,ierr)
! Normalizes NNRM vectors in C.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nnrm, n, metr
real(kind=wp), intent(inout) :: c(n,nnrm)
real(kind=wp), intent(in) :: s(*)
integer(kind=iwp), intent(inout) :: ierr
integer(kind=iwp) :: i
real(kind=wp) :: cnrm
logical(kind=iwp) :: safe
real(kind=wp), allocatable :: c2(:)
real(kind=wp), parameter :: thresh = 1.0e-8_wp
real(kind=wp), external :: ddot_, dnrm2_

if (metr /= 0) call mma_allocate(c2,n,label='c2')
safe = ierr /= 0
do i=1,nnrm
  if (metr == 0) then
    cnrm = dnrm2_(n,c(:,i),1)
  else
    call saoon_cvb(c(:,i),c2,1,s,n,metr)
    cnrm = sqrt(ddot_(n,c(:,i),1,c2,1))
  end if
  if (safe .and. (cnrm < thresh)) then
    ierr = ierr+1
  else
    c(:,i) = c(:,i)/cnrm
  end if
end do
if (metr /= 0) call mma_deallocate(c2)

return

end subroutine nize_cvb
