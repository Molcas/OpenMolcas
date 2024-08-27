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

implicit none
integer M, N, K, LDA, LDB, LDC
integer DoIt(*)
real*8 A(LDA,*), B(LDB,*), C(LDC,*)
integer J, L

! Form  C := A*B + C.

do J=1,N
  if (DoIt(J) /= 1) cycle
  do L=1,K
    if (B(L,J) == ZERO) cycle
    call DAxPy_(M,B(L,J),A(:,L),1,C(:,J),1)
  end do
end do

end subroutine MYDGEMM
