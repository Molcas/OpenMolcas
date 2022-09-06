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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

integer function irAmax(n,a,inc)
! Finds the index of element having max. absolute value.
!
! Author:  F. Aquilante

implicit real*8(a-h,o-z)
real*8 a(*)
integer n, inc

iramax = 0
if ((n < 1) .or. (inc <= 0)) return
iramax = 1
if (n == 1) return
if (inc == 1) go to 20

! code for increment not equal to 1

ix = 1
smax = abs(a(1))
ix = ix+inc
do i=2,n
  if (abs(a(ix)) <= smax) go to 5
  iramax = i
  smax = abs(a(ix))
5 ix = ix+inc
end do
return

! code for increment equal to 1

20 smax = abs(a(1))
do i=2,n
  if (abs(a(i)) <= smax) go to 30
  iramax = i
  smax = abs(a(i))
30 continue
end do
return

end function irAmax
