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

subroutine TRSMT(A,B,C,N,H,W)
! B*A*BT

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(N*(N+1)/2), B(N,N)
real(kind=wp), intent(out) :: C(N*(N+1)/2), H(N,N), W(N,N)
integer(kind=iwp) :: I, IJ, J, K, L

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    C(IJ) = Zero
    W(I,J) = A(IJ)
    W(J,I) = A(IJ)
    H(I,J) = Zero
    H(J,I) = Zero
  end do
end do
do I=1,N
  do L=1,N
    do K=1,N
      H(I,L) = B(I,K)*W(K,L)+H(I,L)
    end do
  end do
end do
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    do L=1,N
      C(IJ) = H(I,L)*B(J,L)+C(IJ)
    end do
  end do
end do

return

end subroutine TRSMT
