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

subroutine mxinv_cvb(a,n)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n,n)
integer(kind=iwp) :: i, ierr
real(kind=wp) :: rms
integer(kind=iwp), allocatable :: itmp(:)
real(kind=wp), allocatable :: tmp1(:,:), tmp2(:,:)
real(kind=wp), parameter :: thresh = 1.0e-10_wp
real(kind=wp), external :: ddot_

call mma_allocate(tmp1,n,n,label='tmp1')
call mma_allocate(tmp2,n,n,label='tmp2')
call mma_allocate(itmp,n,label='itmp')
ierr = 0
tmp1(:,:) = a(:,:)
call dgetrf_(n,n,tmp1,n,itmp,ierr)
if (ierr /= 0) then
  write(u6,*) ' Error in LU decomposition for inversion:',ierr
  call mxprint_cvb(a,n,n,0)
  call abend_cvb()
end if
call dgetri_(n,tmp1,n,itmp,tmp2,n*n,ierr)
! Check solution
call mxatb_cvb(a,tmp1,n,n,n,tmp2)
do i=1,n
  tmp2(i,i) = tmp2(i,i)-One
end do
rms = sqrt(ddot_(n*n,tmp2,1,tmp2,1)/real(n*n,kind=wp))
if (rms > thresh) then
  write(u6,*) ' Fatal error in matrix inversion - error:',rms
  write(u6,*) ' Singular or near-singular matrix.'
  write(u6,*) ' Matrix :'
  call mxprint_cvb(a,n,n,0)
  write(u6,*) ' Inverted matrix :'
  call mxprint_cvb(tmp1,n,n,0)
  write(u6,*) ' Check :'
  call mxprint_cvb(tmp2,n,n,0)
  call mxdiag_cvb(a,tmp2,n)
  write(u6,*) ' Eigenvalues :'
  call mxprint_cvb(tmp2,1,n,0)
  write(u6,*) ' Eigenvectors :'
  call mxprint_cvb(a,1,n,0)
  call abend_cvb()
end if
a(:,:) = tmp1(:,:)

call mma_deallocate(tmp1)
call mma_deallocate(tmp2)
call mma_deallocate(itmp)

return

end subroutine mxinv_cvb
