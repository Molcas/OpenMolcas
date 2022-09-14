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

function irAmax(n,a,inc)
! Finds the index of element having max. absolute value.
!
! Author:  F. Aquilante

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: irAmax
integer(kind=iwp), intent(in) :: n, inc
real(kind=wp), intent(in) :: a(*)
integer(kind=iwp) :: i, ix
real(kind=wp) :: smax

iramax = 0
if ((n < 1) .or. (inc <= 0)) return
iramax = 1
if (n == 1) return

if (inc /= 1) then

  ! code for increment not equal to 1

  ix = 1
  smax = abs(a(1))
  ix = ix+inc
  do i=2,n
    if (abs(a(ix)) > smax) then
      iramax = i
      smax = abs(a(ix))
    end if
    ix = ix+inc
  end do

else

  ! code for increment equal to 1

  smax = abs(a(1))
  do i=2,n
    if (abs(a(i)) > smax) then
      iramax = i
      smax = abs(a(i))
    end if
  end do

end if

return

end function irAmax
