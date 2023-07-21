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

implicit real*8(A-H,O-Z)
parameter(ZERO=0.0D+00,ONE=1.0D+00,G_Eps=1.0D-19)
logical Error
real*8 A(40,40), B(40), X(40)

Error = .false.
do I=1,N
  do K=I,N
    S = ZERO
    do J=I,N
      Temp = A(K,J)*A(K,J)
      S = S+Temp
    end do
    if (S == ZERO) goto 999
    S = sqrt(S)
    SInv = ONE/S
    B(K) = B(K)*SInv
    do J=1,N
      A(K,J) = A(K,J)*SInv
    end do
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
  if ((abs(P) < G_Eps) .and. (P == ZERO)) goto 999

  PInv = ONE/P
  B(I) = B(I)*PInv
  do J=I+1,N
    A(I,J) = A(I,J)*PInv
  end do
  do K=I+1,N
    S = A(K,I)
    if (S /= ZERO) then
      Temp = S*B(I)
      B(K) = B(K)-Temp
      do J=I+1,N
        Temp = S*A(I,J)
        A(K,J) = A(K,J)-Temp
      end do
    end if
  end do
end do

II = N
do I=1,N
  S = B(II)
  do J=II+1,N
    Temp = A(II,J)*X(J)
    S = S-Temp
  end do
  X(II) = S
  II = II-1
end do

Error = .true.

999 continue
return

end subroutine SlvEqs
