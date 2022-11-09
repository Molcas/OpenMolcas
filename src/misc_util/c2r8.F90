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

function C2R8(CBuf)

use Definitions, only: wp, RtoB

implicit none
real(kind=wp) :: C2R8
character, intent(in) :: CBuf(RtoB)

C2R8 = transfer(CBuf,C2R8)

end function C2R8
