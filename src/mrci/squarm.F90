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

subroutine SQUARM(A,B,N)

implicit real*8(A-H,O-Z)
dimension A(*), B(N,N)

IN = 2
do I=2,N
  call VNEG(A(IN),1,B(1,I),1,I-1)
  IN = IN+I
end do
do I=1,N-1
  call VNEG(B(I,I+1),N,B(I+1,I),1,N-I)
end do
call DCOPY_(N,[0.0d00],0,B,N+1)

return

end subroutine SQUARM
