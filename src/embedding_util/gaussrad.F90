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
! Copyright (C) Thomas Dresselhaus                                     *
!***********************************************************************

!***********************************************************************
! Returns the radial part of the value of a GTO with given exponent,   *
! centered at the origin.                                              *
!***********************************************************************
function gaussRad(alpha,r)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: gaussRad
real(kind=wp), intent(in) :: alpha, r(3)
integer(kind=iwp) :: i
real(kind=wp) :: rSquare

rSquare = Zero
do i=1,3
  rSquare = rSquare+r(i)**2
end do

! Radial part
gaussRad = exp(-alpha*rSquare)

return

end function gaussRad
