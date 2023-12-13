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

! R is an operator
! A is a Cartesian Coordinate in 3D
! RA is the Cartesian Coordinate in 3D after the operator R has operated on A.
subroutine OA(R,A,RA)

use Phase_Info, only: iPhase
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: R
real(kind=wp), intent(in) :: A(3)
real(kind=wp), intent(out) :: RA(3)

RA(:) = real(iPhase(:,R),kind=wp)*A(:)

end subroutine OA
