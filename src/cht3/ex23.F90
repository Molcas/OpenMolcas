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

subroutine EX23(Y,X,I1,I2,J1,J2)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: I1, I2, J1, J2
real(kind=wp), intent(in) :: Y(I1,I2,J1,J2)
real(kind=wp), intent(out) :: X(I1,J1,I2,J2)
integer(kind=iwp) :: I, J, L

do J=1,J2
  do L=1,I2
    do I=1,J1
      X(:,I,L,J) = Y(:,L,I,J)
    end do
  end do
end do

return

end subroutine EX23
