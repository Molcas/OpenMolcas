!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
      Real*8 Function Gamma2(m,T)
!***********************************************************************
!                                                                      *
! Object: to compute the auxiliary function in the high argument       *
!         approximation.                                               *
!                                                                      *
!***********************************************************************
      use Constants, only: Zero, One, Two
      Implicit None
      Integer m
      Real*8 T

      Integer i
!
      Gamma2= Sqrt(Two*ACos(Zero)/T)/Two
      Do i = 1, m
         Gamma2 = ((Two*DBLE(i)-One)/(Two*T))*Gamma2
      End Do
      Return
      End Function Gamma2
