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
!
!-- p-p (pi), with too small exponent difference.
!
      Real*8 Function CoulT0_5(Rho,dSepInv,Expo)
      Implicit Real*8 (a-h,o-z)

      T1=1.0d0
      T2=2.0d0*Rho
      T3=2.0d0*Rho**2
      T4=(121.0d0/96.0d0)*Rho**3
      T5=(25.0d0/48.0d0)*Rho**4
      T6=(2.0d0/15.0d0)*Rho**5
      T7=(1.0d0/60.0d0)*Rho**6
      CoulT0_5=dSepInv**3*(1.0d0-(T1+T2+T3+T4+T5+T6+T7)*Expo)

      Return
      End
