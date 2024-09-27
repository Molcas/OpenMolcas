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

subroutine MYDGEMM(DoIt,M,N,K,A,LDA,B,LDB,C,LDC)
!  Purpose
!  =======
!
!  DGEMM for a special case.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, DoIt(N), M, K, LDA, LDB, LDC
real(kind=wp), intent(in) :: A(LDA,K), B(LDB,N)
real(kind=wp), intent(inout) :: C(LDC,N)
integer(kind=iwp) :: J, L

! Form  C := A*B + C.

do J=1,N
  if (DoIt(J) /= 1) cycle
  do L=1,K
    if (B(L,J) /= Zero) C(1:M,J) = C(1:M,J)+B(L,J)*A(1:M,L)
  end do
end do

end subroutine MYDGEMM
