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

subroutine cct3_add10(a,b,dimp,fact)
! this routine does:
! B(p) <-- fact * A(p)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp
real(kind=wp), intent(in) :: a(dimp), fact
real(kind=wp), intent(inout) :: b(dimp)

b(:) = b+fact*a

return

end subroutine cct3_add10
