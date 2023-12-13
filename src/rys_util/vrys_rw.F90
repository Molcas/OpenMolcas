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

module vRys_RW

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MaxRys = 9
integer(kind=iwp) :: iCffR(0:6,MaxRys), iCffW(0:6,MaxRys), iMap(MaxRys), ix0(MaxRys), nMap(MaxRys), nMxRys, nx0(MaxRys)
real(kind=wp), allocatable :: Cff(:), ddx(:), HerR2(:), HerW2(:), TMax(:), x0(:)
integer(kind=iwp), allocatable :: iHerR2(:), iHerW2(:), Map(:)

public :: Cff, ddx, HerR2, HerW2, iCffR, iCffW, iHerR2, iHerW2, iMap, ix0, Map, nMap, nMxRys, nx0, TMax, x0

end module vRys_RW
