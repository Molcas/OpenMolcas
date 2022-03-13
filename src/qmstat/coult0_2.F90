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
real*8 function CoulT0_2(Rho,dSepInv,Expo)

implicit real*8(a-h,o-z)

T1 = 1.0d0
T2 = 2.0d0*Rho
T3 = 2.0d0*Rho**2
T4 = (59.0d0/48.0d0)*Rho**3
T5 = (11.0d0/24.0d0)*Rho**4
T6 = (1.0d0/12.0d0)*Rho**5
CoulT0_2 = dSepInv**2*(1.0d0-(T1+T2+T3+T4+T5+T6)*Expo)

return

end function CoulT0_2
