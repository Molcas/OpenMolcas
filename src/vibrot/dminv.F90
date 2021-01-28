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

! This subroutine inverts the matrix a of dimension N,N
Subroutine DMINV(N,NMAX,A)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, NMAX
real(kind=wp), intent(inout) :: A(NMAX,NMAX)
integer(kind=iwp) :: I, J, K
real(kind=wp) :: BIGA, HOLD

do  K=1,N
  BIGA = A(K,K)
  do I=1,N
    if (I /= K) then
      A(I,K) = -A(I,K)/BIGA
    end if
  end do
  do I=1,N
    HOLD = A(I,K)
    do J=1,N
      if ((I /= K).and.(J /= K)) then
        A(I,J) = HOLD*A(K,J)+A(I,J)
      end if
    end do
  end do
  do J=1,N
    if(J /= K) then
      A(K,J) = A(K,J)/BIGA
    end if
  end do
  A(K,K) = One/BIGA
end do

return

end subroutine dminv
