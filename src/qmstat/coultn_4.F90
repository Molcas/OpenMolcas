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
!-- p-p (sigma), normal case.
!
      Real*8 Function CoulTN_4(R,T,RA,RB,C,dSepInv,ExpA,ExpB)
      Implicit Real*8 (a-h,o-z)

      T1=(1.0d0/16.0d0)*(8.0d0+9.0d0*C+3.0d0*C**2)                      &
     &                 *(1.0d0+2.0d0*RA+2.0d0*RA**2)
      T2=(3.0d0/16.0d0)*(3.0d0+2.0d0*C)*RA**3
      T3=(1.0d0/8.0d0)*RA**4
      TA=(1.0d0-C)**3*(T1+T2+T3)*ExpA
      T1=(1.0d0/16.0d0)*(8.0d0-9.0d0*C+3.0d0*C**2)                      &
     &                 *(1.0d0+2.0d0*RB+2.0d0*RB**2)
      T2=(3.0d0/16.0d0)*(3.0d0-2.0d0*C)*RB**3
      T3=(1.0d0/8.0d0)*RB**4
      TB=(1.0d0+C)**3*(T1+T2+T3)*ExpB
      CoulTN_4=2.0d0*dSepInv**3*(1.0d0-TA-TB)

      Return
! Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real(R)
        Call Unused_real(T)
      End If
      End
