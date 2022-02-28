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

subroutine MTRANS(A,B,N,M)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, M
real(kind=wp), intent(in) :: A(M,N)
real(kind=wp), intent(out) :: B(N,M)
integer(kind=iwp) :: I, J

do I=1,N
  do J=1,M
    B(I,J) = A(J,I)
  end do
end do

return

end subroutine MTRANS
