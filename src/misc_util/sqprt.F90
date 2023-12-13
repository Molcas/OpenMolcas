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

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(N,N)
character(len=60) :: FRMT
integer(kind=iwp) :: I, J
real(kind=wp) :: BIG

BIG = Zero
do I=1,N
  do J=1,N
    BIG = max(BIG,abs(A(I,J)))
  end do
end do
if ((0.1_wp < BIG) .and. (BIG < 1.0e4_wp)) then
  FRMT = '(8(1X,F13.6))'
else
  FRMT = '(8(1X,ES13.6))'
end if
do I=1,N
  write(u6,FRMT) A(I,:)
end do

return

end subroutine SQPRT
