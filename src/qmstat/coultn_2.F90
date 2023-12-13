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

! s-p interaction, normal case.
function CoulTN_2(RA,RB,C,dSepInv,ExpA,ExpB)

use Constants, only: One, Two, Three, Five, Ten, Eleven, Half, Quart
use Definitions, only: wp

implicit none
real(kind=wp) :: CoulTN_2
real(kind=wp), intent(in) :: RA, RB, C, dSepInv, ExpA, ExpB
real(kind=wp) :: T1, T2, T3, TA, TB

T1 = (One/16.0_wp)*(Five+Three*C)*(One+Two*RA)
T2 = Quart*RA**2
TA = (One-C)**3*(T1+T2)*ExpA
T1 = (One/16.0_wp)*(Eleven-Ten*C+Three*C**2)*(One+Two*RB)
T2 = Half*(Two-C)*RB**2
T3 = Quart*RB**3
TB = (One+C)**2*(T1+T2+T3)*ExpB
CoulTN_2 = dSepInv**2*(One-TA-TB)

return

end function CoulTN_2
