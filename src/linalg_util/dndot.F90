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

subroutine DNDOT(N,M,S,INCS,ISW,X,INCXI,INCXO,Y,INCYI,INCYO)

! COMPUTE DOT PRODUCT N TIMES

implicit real*8(A-H,O-Z)
real*8 S(1+(N-1)*INCS),X((1+(M-1)*INCXI)*(1+(N-1)*INCXO)),Y((1+(M-1)*INCYI)*(1+(N-1)*INCYO))

if (ISW == 1) then
  do I=1,N
    S(1+(I-1)*INCS) = DDOT_(M,X(1+(I-1)*INCXO),INCXI,Y(1+(I-1)*INCYO),INCYI)
  end do
else if (ISW == 2) then
  do I=1,N
    S(1+(I-1)*INCS) = -DDOT_(M,X(1+(I-1)*INCXO),INCXI,Y(1+(I-1)*INCYO),INCYI)
  end do
else if (ISW == 3) then
  do I=1,N
    S(1+(I-1)*INCS) = S(1+(I-1)*INCS)+DDOT_(M,X(1+(I-1)*INCXO),INCXI,Y(1+(I-1)*INCYO),INCYI)
  end do
else if (ISW == 4) then
  do I=1,N
    S(1+(I-1)*INCS) = S(1+(I-1)*INCS)-DDOT_(M,X(1+(I-1)*INCXO),INCXI,Y(1+(I-1)*INCYO),INCYI)
  end do
else
  call SysAbendMsg('dndot','ISW IS OUT OF RANGE IN DNDOT',' ')
end if

return

end subroutine DNDOT
