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

subroutine Store_QLast(Q)

use definitions, only: wp
use TList_Mod, only: QLast

implicit none
real(kind=wp), intent(in) :: Q(2)

QLast(1) = Q(1)
QLast(2) = Q(2)

end subroutine Store_QLast
