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

subroutine CPLAB(A,B,L,M,N,IA,IB,C,IC,IER)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, M, N, IA, IB, IC
real(kind=wp), intent(in) :: A(IA,M), B(IB,N)
real(kind=wp), intent(inout) :: C(IC,N)
integer(kind=iwp), intent(out) :: IER
integer(kind=iwp) :: I, J, K
real(kind=wp) :: TEMP

if ((IA >= L) .and. (IB >= M) .and. (IC >= L)) then
  IER = 0
  do I=1,L
    do J=1,N
      TEMP = Zero
      do K=1,M
        TEMP = A(I,K)*B(K,J)+TEMP
      end do
      C(I,J) = C(I,J)+TEMP
    end do
  end do
else
  IER = 129
end if

return

end subroutine CPLAB
