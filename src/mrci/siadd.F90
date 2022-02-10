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

subroutine SIADD(A,B,N)

implicit real*8(A-H,O-Z)
dimension A(N,N), B(*)

IN = 0
do I=1,N
  do J=1,I
    IN = IN+1
    B(IN) = B(IN)+A(I,J)+A(J,I)
  end do
  B(IN) = B(IN)-A(I,I)
end do

return

end subroutine SIADD
