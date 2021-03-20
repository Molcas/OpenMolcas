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

subroutine power(c,n,ipower)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, ipower
real(kind=wp), intent(inout) :: c(n)
integer(kind=iwp) :: i

if (ipower == 1) then
  return
else if (ipower == 2) then
  do i=1,n
    c(i) = c(i)*c(i)
  end do
else if (ipower == -2) then
  do i=1,n
    c(i) = sign(c(i)*c(i),c(i))
  end do
end if

return

end subroutine power
