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

! p-p (sigma), with too small exponent difference.
real*8 function CoulT0_4(Rho,dSepInv,Expo)

implicit real*8(a-h,o-z)

T1 = 1.0d0
T2 = 2.0d0*Rho
T3 = 2.0d0*Rho**2
T4 = (263.0d0/192.0d0)*Rho**3
T5 = (71.0d0/96.0d0)*Rho**4
T6 = (77.0d0/240.0d0)*Rho**5
T7 = (1.0d0/10.0d0)*Rho**6
T8 = (1.0d0/60.0d0)*Rho**7
CoulT0_4 = 2.0d0*dSepInv**3*(1.0d0-(T1+T2+T3+T4+T5+T6+T7+T8)*Expo)

return

end function CoulT0_4
