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

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
dimension a(n,n)

i1 = mstackr_cvb(n)
i2 = mstackr_cvb(n*n)
i3 = mstackr_cvb(n)
i4 = mstackr_cvb(n)
i5 = mstackr_cvb(n*n)
ifail = 0
call rs(n,n,a,work(i1),1,work(i2),work(i3),work(i4),ifail)
if (ifail /= 0) then
  write(6,*) ' Fatal error in diagonalization (MXSQRT) :',ifail
  call abend_cvb()
end if
call fzero(a,n*n)
do i=1,n
  a(i,i) = sqrt(work(i+i1-1))**ipow
end do
call mxatb_cvb(work(i2),a,n,n,n,work(i5))
call fzero(a,n*n)
do k=1,n
  do j=1,n
    do i=1,n
      a(i,j) = a(i,j)+work(i+(k-1)*n+i5-1)*work(j+(k-1)*n+i2-1)
    end do
  end do
end do
call mfreer_cvb(i1)

return

end subroutine mxsqrt_cvb
