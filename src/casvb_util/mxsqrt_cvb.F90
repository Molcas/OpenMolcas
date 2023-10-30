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

subroutine mxsqrt_cvb(a,n,ipow)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, ipow
real(kind=wp), intent(inout) :: a(n,n)
integer(kind=iwp) :: i, ifail, j, k
real(kind=wp), allocatable :: c(:,:), w(:), z(:,:)

call mma_allocate(w,n,label='w')
call mma_allocate(z,n,n,label='z')
call mma_allocate(c,n,n,label='c')
ifail = 0
call rs(n,n,a,w,1,z,ifail)
if (ifail /= 0) then
  write(u6,*) ' Fatal error in diagonalization (MXSQRT) :',ifail
  call abend_cvb()
end if
a(:,:) = Zero
do i=1,n
  a(i,i) = sqrt(w(i))**ipow
end do
call mxatb_cvb(z,a,n,n,n,c)
a(:,:) = Zero
do k=1,n
  do j=1,n
    a(:,j) = a(:,j)+c(:,k)*z(j,k)
  end do
end do
call mma_deallocate(w)
call mma_deallocate(z)
call mma_deallocate(c)

return

end subroutine mxsqrt_cvb
