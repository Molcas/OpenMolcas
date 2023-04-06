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
! Copyright (C) 2006, Pavel Neogrady                                   *
!***********************************************************************

subroutine extstackhlp1(a,b,dimij,dimb,bb)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimij, dimb, bb
real(kind=wp), intent(out) :: a(dimij)
real(kind=wp), intent(in) :: b(dimij,dimb)

a(:) = b(:,bb)

return

end subroutine extstackhlp1
