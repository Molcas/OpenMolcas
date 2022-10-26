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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
!
! R is an operator
! A is a Cartesian Coordinate in 3D
! RA is the Cartesian Coordinate in 3D after the operator R has operated on A.
!
Subroutine OA(R,A,RA)
use Phase_Info
Implicit None
Real*8, Intent(In):: A(3)
Real*8, Intent(Out):: RA(3)
Integer, Intent(In):: R
RA(1:3) = DBLE(iPhase(1:3,R))*A(1:3)
End Subroutine OA
