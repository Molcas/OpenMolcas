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

subroutine dkh_cofu_spec(n,a,m,c)
! Special combination coefficients

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, m
real(kind=wp), intent(in) :: a(n)
real(kind=wp), intent(out) :: c(m)
integer(kind=iwp) :: i
real(kind=wp), allocatable :: b(:)

c(1) = a(m-1)
do i=2,m-1
  c(i) = a(i-1)*a(m-i)*(-1)**(i-1)
end do
c(m) = a(m-1)*(-1)**(m-1)

call mma_allocate(b,m+1,label='b')

b(1) = a(m)
do i=2,m
  b(i) = a(i-1)*a(m-i+1)*(-1)**(i-1)
end do
b(m+1) = a(m)*(-1)**m

do i=1,m
  c(i) = c(i)-b(i)
  b(i+1) = b(i+1)+b(i)
end do
if (abs(b(m+1)) > 1.0e-12_wp) then
  write(u6,*) 'Error in dkh_dkcof_sp ',b(m+1)
  call Abend()
end if

call mma_deallocate(b)

end subroutine dkh_cofu_spec
