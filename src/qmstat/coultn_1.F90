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

! s-s interaction, normal case.
function CoulTN_1(RA,RB,C,dSepInv,ExpA,ExpB)

use Constants, only: One, Two, Quart
use Definitions, only: wp

implicit none
real(kind=wp) :: CoulTN_1
real(kind=wp), intent(in) :: RA, RB, C, dSepInv, ExpA, ExpB
real(kind=wp) :: T1, T2, TA, TB

T1 = Quart*(Two+C)
T2 = Quart*RA
TA = (One-C)**2*(T1+T2)*ExpA
T1 = Quart*(Two-C)
T2 = Quart*RB
TB = (One+C)**2*(T1+T2)*ExpB
CoulTN_1 = dSepInv*(One-TA-TB)

return

end function CoulTN_1
