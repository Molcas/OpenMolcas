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

subroutine svd2_cvb(ainp,val,vec,vmat,n1,n2,n12,a,w,u,v,rv1,indx)

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: n1, n2, n12, indx(n2)
real(kind=wp) :: ainp(n1,n2), val(n2), vec(n1,n2), vmat(n2,n2), a(n12,n2), w(n2), u(n12,n2), v(n12,n2), rv1(n2)
integer(kind=iwp) :: i, ierr
real(kind=wp), external :: dnrm2_

if (n12 == n1) then
  call fmove_cvb(ainp,a,n1*n2)
else
  call fzero(a,n12*n2)
  do i=1,n2
    call fmove_cvb(ainp(1,i),a(1,i),n1)
  end do
end if
ierr = 0
call svd(n12,n1,n2,a,w,.true.,u,.true.,v,ierr,rv1)

if (ierr /= 0) then
  write(u6,*) ' Fatal error in SVD_CVB!',ierr
  call abend_cvb()
end if

! Eispack code is broken, in the following u is generated from v:

! First recreate a:
if (n12 == n1) then
  call fmove_cvb(ainp,a,n1*n2)
else
  call fzero(a,n12*n2)
  do i=1,n2
    call fmove_cvb(ainp(1,i),a(1,i),n1)
  end do
end if

do i=1,n2
  call mxatb_cvb(a,v(1,i),n12,n2,1,u(1,i))
  call dscal_(n12,One/dnrm2_(n12,u(1,i),1),u(1,i),1)
end do

! Sort singular values in ascending order:
call sortindxr_cvb(n2,w,indx)
do i=1,n2
  val(i) = w(indx(i))
  call fmove_cvb(v(1,indx(i)),vmat(1,i),n2)
  call fmove_cvb(u(1,indx(i)),vec(1,i),n1)
end do

return

end subroutine svd2_cvb
