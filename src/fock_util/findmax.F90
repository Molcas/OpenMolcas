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

function FindMax(LaJ,numA)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: FindMax
integer(kind=iwp), intent(in) :: numA
real(kind=wp), intent(in) :: LaJ(NumA)
integer(kind=iwp) :: ia
real(kind=wp) :: XMax

XMax = abs(LaJ(1))

do ia=2,numA

  XMax = max(XMax,abs(LaJ(ia)))

end do

FindMax = XMax

return

end function FindMax
