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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************

function threej(XJ1,XJ2,XJ3,XM1,XM2,XM3)
!  ThreeJ: real Wigner 3-j coefficients. From a modification
!          of Racah's formula for Clebsch-Gordan coeffs.

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: threej
real(kind=wp), intent(in) :: XJ1,XJ2,XJ3,XM1,XM2,XM3
real(kind=wp), external :: DCLEBS
integer(kind=iwp) :: i

threej = DCLEBS(XJ1,XJ2,XJ3,XM1,XM2,-XM3)
if (threej == Zero) return
i = nint(XJ1-XJ2-XM3)
if (i /= (i/2)*2) threej = -threej
threej = threej/sqrt(Two*XJ3+One)

end function threej
