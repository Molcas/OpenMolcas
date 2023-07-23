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

subroutine CHO_INVPCK(IJ,I,J,LOW)
!
! Purpose: invert triangular packing index ij to
!          rectangular indices i and j.
!          Flag LOW specifies packing convention:
!          LOW = T: i>=j
!          LOW = F: i<=j

use Constants, only: One, Three, Eight, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IJ, I, J
logical(kind=iwp) :: LOW
integer(kind=iwp) :: ITMP
real(kind=wp) :: XI, XX

if (IJ > 0) then

  XX = Eight*real(IJ,kind=wp)-Three
  XI = (One+sqrt(XX))*Half
  I = int(XI)
  J = IJ-I*(I-1)/2

  if (.not. LOW) then
    ITMP = I
    I = J
    J = ITMP
  end if

else

  I = -1
  J = -2

end if

end subroutine CHO_INVPCK
