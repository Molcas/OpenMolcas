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

implicit none
integer I1, I2, J1, J2, I, J, IJ, L
real*8 Y, X, ONE
parameter(ONE=1.d0)
dimension Y(I1,I2,J1,J2), X(*)

IJ = 1
do J=1,J2
  do L=1,I2
    do I=1,J1
      call DCOPY_(I1,Y(1,L,I,J),1,X(IJ),1)
      IJ = IJ+I1
    end do
  end do
end do

return

end subroutine EX23
