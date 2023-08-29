!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine SlvEqs(N,A,X,B,Error)
!-----------------------------------------------------------------------
! Function : Solve equations (Loesen)
!            Solve AX = B
!-----------------------------------------------------------------------

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: A(40,40), X(40), B(40)
logical(kind=iwp), intent(out) :: Error
integer(kind=iwp) :: I, II, J, K
real(kind=wp) :: P, PInv, S, SInv
real(kind=wp), parameter :: G_Eps = 1.0e-19_wp

Error = .false.
do I=1,N
  do K=I,N
    S = sum(A(K,I:N)**2)
    if (S == Zero) return
    S = sqrt(S)
    SInv = One/S
    B(K) = B(K)*SInv
    A(K,1:N) = A(K,1:N)*SInv
  end do

  P = A(I,I)
  K = I
  do J=I,N
    S = A(J,I)
    if (abs(S) > abs(P)) then
      P = S
      K = J
    end if
  end do
  if (K > I) then
    do J=I,N
      S = A(I,J)
      A(I,J) = A(K,J)
      A(K,J) = S
    end do
    S = B(I)
    B(I) = B(K)
    B(K) = S
  end if
  if ((abs(P) < G_Eps) .and. (P == Zero)) return

  PInv = One/P
  B(I) = B(I)*PInv
  A(I,I+1:N) = A(I,I+1:N)*PInv
  do K=I+1,N
    S = A(K,I)
    if (S /= Zero) then
      B(K) = B(K)-S*B(I)
      A(K,J) = A(K,J)-S*sum(A(I,I+1:N))
    end if
  end do
end do

II = N
do I=1,N
  X(II) = B(II)-sum(A(II,II+1:N)*X(II+1:N))
  II = II-1
end do

Error = .true.

return

end subroutine SlvEqs
