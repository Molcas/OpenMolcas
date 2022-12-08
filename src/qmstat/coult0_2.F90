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

! s-p interaction, with too small exponent difference.
function CoulT0_2(Rho,dSepInv,Expo)

use Constants, only: One, Two, Eleven, Twelve
use Definitions, only: wp

implicit none
real(kind=wp) :: CoulT0_2
real(kind=wp), intent(in) :: Rho, dSepInv, Expo
real(kind=wp) :: T1, T2, T3, T4, T5, T6

T1 = One
T2 = Two*Rho
T3 = Two*Rho**2
T4 = (59.0_wp/48.0_wp)*Rho**3
T5 = (Eleven/24.0_wp)*Rho**4
T6 = (One/Twelve)*Rho**5
CoulT0_2 = dSepInv**2*(One-(T1+T2+T3+T4+T5+T6)*Expo)

return

end function CoulT0_2
