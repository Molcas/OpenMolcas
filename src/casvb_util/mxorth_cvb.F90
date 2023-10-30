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

function mxorth_cvb(a,n)
! Returns .TRUE. if A is orthogonal.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: mxorth_cvb
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: a(n,n)
integer(kind=iwp) :: i, j
real(kind=wp) :: tst
real(kind=wp), allocatable :: at(:,:), c(:,:)
real(kind=wp), parameter :: thresh = 1.0e-8_wp

call mma_allocate(at,n,n,label='at')
call mma_allocate(c,n,n,label='c')
! AT <= A transpose
do i=1,n
  at(i,:) = a(:,i)
end do
call mxatb_cvb(at,a,n,n,n,c)
call mma_deallocate(at)
! C identity ??
mxorth_cvb = .true.
outer: do j=1,n
  do i=1,n
    if (i /= j) then
      tst = abs(c(i,j))
    else
      tst = abs(c(i,j)-One)
    end if
    if (tst > thresh) then
      mxorth_cvb = .false.
      exit outer
    end if
  end do
end do outer
call mma_deallocate(c)

return

end function mxorth_cvb
