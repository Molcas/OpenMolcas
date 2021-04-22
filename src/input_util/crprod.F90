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

subroutine crprod(a,b,c)

! calculates cross product:   a x b = c

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: a(3), b(3)
real(kind=wp), intent(out) :: c(3)

c(1) = a(2)*b(3)-a(3)*b(2)
c(2) = a(3)*b(1)-a(1)*b(3)
c(3) = a(1)*b(2)-a(2)*b(1)

return

end subroutine crprod
