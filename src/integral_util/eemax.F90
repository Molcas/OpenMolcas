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

subroutine EEMax(n,A,B)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: A(n)
real(kind=wp), intent(in) :: B(n)
integer(kind=iwp) :: i

do i=1,n
  if (B(i) > A(i)) A(i) = B(i)
end do

end subroutine EEMax
