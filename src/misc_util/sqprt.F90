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

subroutine SQPRT(A,N)

implicit real*8(A-H,O-Z)
dimension A(N,N)
character*60 FMT

BIG = 0.0d0
do I=1,N
  do J=1,N
    BIG = max(BIG,abs(A(I,J)))
  end do
end do
if ((0.1d0 < BIG) .and. (BIG < 10000.0d0)) then
  FMT = '(8(1X,F12.6))'
else
  FMT = '(8(1X,E12.6))'
end if
do I=1,N
  write(6,FMT) (A(I,J),J=1,N)
end do

return

end subroutine SQPRT
