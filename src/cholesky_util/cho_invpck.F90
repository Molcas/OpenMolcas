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

use Index_Functions, only: iTri_Rev
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IJ
integer(kind=iwp), intent(out) :: I, J
logical(kind=iwp), intent(in) :: LOW
integer(kind=iwp) :: RES(2)

if (IJ > 0) then

  RES = iTri_Rev(IJ)
  if (LOW) then
    I = RES(1)
    J = RES(2)
  else
    I = RES(2)
    J = RES(1)
  end if

else

  I = -1
  J = -2

end if

end subroutine CHO_INVPCK
