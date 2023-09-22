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

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: a(n,n)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1, i2, i3, ierr
real(kind=wp) :: rms
real(kind=wp), parameter :: thresh = 1.0e-10_wp
integer(kind=iwp), external :: mstacki_cvb, mstackr_cvb
real(kind=wp), external :: ddot_

i1 = mstackr_cvb(n*n)
i2 = mstackr_cvb(n*n)
i3 = mstacki_cvb(n)
ierr = 0
call fmove_cvb(a,work(i1),n*n)
call dgetrf_(n,n,work(i1),n,iwork(i3),ierr)
if (ierr /= 0) then
  write(u6,*) ' Error in LU decomposition for inversion:',ierr
  call mxprint_cvb(a,n,n,0)
  call abend_cvb()
end if
call dgetri_(n,work(i1),n,iwork(i3),work(i2),n*n,ierr)
! Check solution
call mxatb_cvb(a,work(i1),n,n,n,work(i2))
do i=1,n
  work(i+(i-1)*n+i2-1) = work(i+(i-1)*n+i2-1)-One
end do
rms = sqrt(ddot_(n*n,work(i2),1,work(i2),1)/real(n*n,kind=wp))
if (rms > thresh) then
  write(u6,*) ' Fatal error in matrix inversion - error:',rms
  write(u6,*) ' Singular or near-singular matrix.'
  write(u6,*) ' Matrix :'
  call mxprint_cvb(a,n,n,0)
  write(u6,*) ' Inverted matrix :'
  call mxprint_cvb(work(i1),n,n,0)
  write(u6,*) ' Check :'
  call mxprint_cvb(work(i2),n,n,0)
  call mxdiag_cvb(a,work(i2),n)
  write(u6,*) ' Eigenvalues :'
  call mxprint_cvb(work(i2),1,n,0)
  write(u6,*) ' Eigenvectors :'
  call mxprint_cvb(a,1,n,0)
  call abend_cvb()
end if
call fmove_cvb(work(i1),a,n*n)
call mfreer_cvb(i1)

return

end subroutine mxinv_cvb
