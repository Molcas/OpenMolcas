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

subroutine UpdateMostNegative(n,X,Val)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: X(n)
real(kind=wp), intent(in) :: Val
integer(kind=iwp) :: i, j

if (Val >= X(n)) return
i = 0
do while (i < n)
  i = i+1
  if (Val < X(i)) then
    do j=n,i+1,-1
      X(j) = X(j-1)
    end do
    X(i) = Val
    i = n+1 ! break while loop
  end if
end do

end subroutine UpdateMostNegative
