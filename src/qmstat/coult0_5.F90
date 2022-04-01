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

! p-p (pi), with too small exponent difference.
function CoulT0_5(Rho,dSepInv,Expo)

use Constants, only: One, Two
use Definitions, only: wp

implicit none
real(kind=wp) :: CoulT0_5
real(kind=wp), intent(in) :: Rho, dSepInv, Expo
real(kind=wp) :: T1, T2, T3, T4, T5, T6, T7

T1 = One
T2 = Two*Rho
T3 = Two*Rho**2
T4 = (121.0_wp/96.0_wp)*Rho**3
T5 = (25.0_wp/48.0_wp)*Rho**4
T6 = (Two/15.0_wp)*Rho**5
T7 = (One/60.0_wp)*Rho**6
CoulT0_5 = dSepInv**3*(One-(T1+T2+T3+T4+T5+T6+T7)*Expo)

return

end function CoulT0_5
