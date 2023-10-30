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

subroutine schmidt2_cvb(c,sxc,nvec,sao,n,metr)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec, n, metr
real(kind=wp), intent(inout) :: c(n,nvec), sxc(n,nvec)
real(kind=wp), intent(in) :: sao(*)
integer(kind=iwp) :: i, j
real(kind=wp), allocatable :: cnrm(:)
real(kind=wp), parameter :: thresh = 1.0e-20_wp
real(kind=wp), external :: ddot_

call mma_allocate(cnrm,nvec,label='cnrm')
do i=1,nvec
  do j=1,i-1
    if (cnrm(j) > thresh) c(:,i) = c(:,i)-ddot_(n,c(:,i),1,sxc(:,j),1)/cnrm(j)*c(:,j)
  end do
  if (metr /= 0) call saoon_cvb(c(:,i),sxc(:,i),1,sao,n,metr)
  cnrm(i) = ddot_(n,c(:,i),1,sxc(:,i),1)
end do
call mma_deallocate(cnrm)

return

end subroutine schmidt2_cvb
