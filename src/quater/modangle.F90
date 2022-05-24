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

function modangle(angle,ref)

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: modangle
real(kind=wp), intent(in) :: angle, ref
integer(kind=iwp) :: n

n = int(Two*angle/ref)

modangle = angle - n*ref

end function modangle
