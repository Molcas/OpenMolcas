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

subroutine RNFill(length,A,c)
! fill an array with random numbers in interval (-c,c)
! (not true)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: length
real(kind=wp), intent(out) :: A(length)
real(kind=wp), intent(in) :: c
integer(kind=iwp) :: i

#include "macros.fh"
unused_var(c)

do i=1,length
  !A(i) = c*(srand()-Half)
  A(i) = 1.0e-7_wp*i
end do

return

end subroutine RNFill
