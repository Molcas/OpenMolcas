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

subroutine svd_cvb(ainp,val,vec,vmat,n1,n2)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n1, n2
real(kind=wp), intent(in) :: ainp(n1,n2)
real(kind=wp), intent(out) :: val(n2), vec(n1,n2), vmat(n2,n2)
integer(kind=iwp) :: i, ierr, n12
integer(kind=iwp), allocatable :: indx(:)
real(kind=wp), allocatable :: a(:,:), rv1(:), u(:,:), v(:,:), w(:)
real(kind=wp), external :: dnrm2_

n12 = max(n1,n2)
call mma_allocate(a,n12,n2,label='a')
call mma_allocate(w,n2,label='w')
call mma_allocate(u,n12,n2,label='u')
call mma_allocate(v,n12,n2,label='v')
call mma_allocate(rv1,n2,label='rv1')
call mma_allocate(indx,n2,label='indx')

if (n12 == n1) then
  a(:,:) = ainp(:,:)
else
  a(1:n1,:) = ainp(:,:)
  a(n1+1:,:) = Zero
end if
ierr = 0
call svd(n12,n1,n2,a,w,.true.,u,.true.,v,ierr,rv1)

call mma_deallocate(rv1)

if (ierr /= 0) then
  write(u6,*) ' Fatal error in SVD_CVB!',ierr
  call abend_cvb()
end if

! Eispack code is broken, in the following u is generated from v:

! First recreate a:
if (n12 == n1) then
  a(:,:) = ainp(:,:)
else
  a(1:n1,:) = ainp(:,:)
  a(n1+1:,:) = Zero
end if

do i=1,n2
  call mxatb_cvb(a,v(:,i),n12,n2,1,u(:,i))
  u(:,i) = u(:,i)/dnrm2_(n12,u(:,i),1)
end do

call mma_deallocate(a)
call mma_allocate(indx,n2,label='indx')

! Sort singular values in ascending order:
call sortindxr_cvb(n2,w,indx)
do i=1,n2
  val(i) = w(indx(i))
  vmat(:,i) = v(1:n2,indx(i))
  vec(:,i) = u(1:n1,indx(i))
end do

call mma_deallocate(a)
call mma_deallocate(w)
call mma_deallocate(u)
call mma_deallocate(v)
call mma_deallocate(indx)

return

end subroutine svd_cvb
