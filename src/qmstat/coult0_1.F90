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
real*8 function CoulT0_1(Rho,dSepInv,Expo)

implicit real*8(a-h,o-z)

T1 = 1.0d0
T2 = (11.0d0/8.0d0)*Rho
T3 = (3.0d0/4.0d0)*Rho**2
T4 = (1.0d0/6.0d0)*Rho**3
CoulT0_1 = dSepInv*(1.0d0-(T1+T2+T3+T4)*Expo)

return

end function CoulT0_1
