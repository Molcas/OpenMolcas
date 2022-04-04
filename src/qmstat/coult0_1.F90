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

! s-s interaction, with too small exponent difference.
function CoulT0_1(Rho,dSepInv,Expo)

use Constants, only: One, Three, Four, Six, Eight, Eleven
use Definitions, only: wp

implicit none
real(kind=wp) :: CoulT0_1
real(kind=wp), intent(in) :: Rho, dSepInv, Expo
real(kind=wp) :: T1, T2, T3, T4

T1 = One
T2 = (Eleven/Eight)*Rho
T3 = (Three/Four)*Rho**2
T4 = (One/Six)*Rho**3
CoulT0_1 = dSepInv*(One-(T1+T2+T3+T4)*Expo)

return

end function CoulT0_1
