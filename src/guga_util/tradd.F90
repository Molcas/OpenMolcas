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

subroutine TRADD(A,B,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(N,N)
real(kind=wp), intent(inout) :: B(*)
integer(kind=iwp) :: I, IIN, J

IIN = 0
do I=1,N
  do J=1,I
    IIN = IIN+1
    B(IIN) = B(IIN)+A(I,J)-A(J,I)
  end do
end do

return

end subroutine TRADD
