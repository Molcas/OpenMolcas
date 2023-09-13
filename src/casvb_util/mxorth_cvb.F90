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

logical function mxorth_cvb(a,n)
! Returns .TRUE. if A is orthogonal.

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
dimension a(n,n)
save thresh, one
data thresh/1d-8/,one/1d0/

i1 = mstackr_cvb(n*n)
i2 = mstackr_cvb(n*n)
! Work(I1) <= A transpose
do i=1,n
  do j=1,n
    work(i+(j-1)*n+i1-1) = a(j,i)
  end do
end do
call mxatb_cvb(work(i1),a,n,n,n,work(i2))
! Work(I2) identity ??
mxorth_cvb = .true.
do j=1,n
  do i=1,n
    if (i /= j) then
      tst = abs(work(i+(j-1)*n+i2-1))
      if (tst > thresh) mxorth_cvb = .false.
    else
      tst = abs(work(i+(j-1)*n+i2-1)-one)
      if (tst > thresh) mxorth_cvb = .false.
    end if
  end do
end do

return

end function mxorth_cvb
