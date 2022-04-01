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

! p-p (pi), normal case.
function CoulTN_5(RA,RB,C,dSepInv,ExpA,ExpB)

use Constants, only: One, Two, Three, Five, Eight, Nine
use Definitions, only: wp

implicit none
real(kind=wp) :: CoulTN_5
real(kind=wp), intent(in) :: RA, RB, C, dSepInv, ExpA, ExpB
real(kind=wp) :: T1, T2, T3, TA, TB

T1 = (One/16.0_wp)*(Eight+Nine*C+Three*C**2)*(One+Two*RA)
T2 = (One/Eight)*(Five+Three*C)*RA**2
T3 = (One/Eight)*RA**3
TA = (One-C)**3*(T1+T2+T3)*ExpA
T1 = (One/16.0_wp)*(Eight-Nine*C+Three*C**2)*(One+Two*RB)
T2 = (One/Eight)*(Five-Three*C)*RB**2
T3 = (One/Eight)*RB**3
TB = (One+C)**3*(T1+T2+T3)*ExpB
CoulTN_5 = dSepInv**3*(One-TA-TB)

return

end function CoulTN_5
