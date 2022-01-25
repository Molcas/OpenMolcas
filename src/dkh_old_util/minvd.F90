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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

!#NUMPAC#MINVD               REVISED ON 1984-11-30
subroutine MINVD(A,KA,N,EPS,ILL)

implicit real*8(A-H,O-Z)
dimension A(KA,N)
integer MX(1000)

if ((N < 1) .or. (N > 1000) .or. (N > KA) .or. (EPS <= 0.0)) GO TO 250
! LU DECOMPOSITION
JM1 = 0 ! dummy initialize
M = 0 ! dummy initialize
NM1 = N-1
do J=1,N
  if (J == 1) GO TO 30
  JM1 = J-1
  do I=1,JM1
    M = MX(I)
    S = A(M,J)
    A(M,J) = A(I,J)
    if (I == 1) GO TO 21
    IM1 = I-1
    do K=1,IM1
      S = A(I,K)*A(K,J)+S
    end do
21  A(I,J) = S
  end do
30 AM = 0.
  do I=J,N
    S = A(I,J)
    if (J == 1) GO TO 50
    do K=1,JM1
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S
50  AA = abs(S)
    if (AA <= AM) GO TO 60
    AM = AA
    M = I
60  continue
  end do
  if (AM < EPS) GO TO 240
  MX(J) = M
  if (M == J) GO TO 80
  do K=1,J
    W = A(M,K)
    A(M,K) = A(J,K)
    A(J,K) = W
  end do
80 if (J == N) GO TO 100
  JP1 = J+1
  W = -A(J,J)
  do I=JP1,N
    A(I,J) = A(I,J)/W
  end do
end do
100 if (N <= 2) GO TO 130
! INPLACE INVERSION OF L-COMPONENT
do I=3,N
  IM1 = I-1
  IM2 = I-2
  do J=1,IM2
    S = A(I,J)
    JP1 = J+1
    do K=JP1,IM1
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S
  end do
end do
! INPLACE INVERSION OF U-COMPONENT
130 A(1,1) = 1./A(1,1)
if (N == 1) GO TO 230
do J=2,N
  A(J,J) = 1./A(J,J)
  P = -A(J,J)
  JM1 = J-1
  do I=1,JM1
    S = 0.
    do K=I,JM1
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S*P
  end do
end do
! INPLACE MULTIPLICATION OF L AND U COMPONENT
do J=1,NM1
  JP1 = J+1
  do I=1,J
    S = A(I,J)
    do K=JP1,N
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S
  end do
  do I=JP1,N
    S = 0.
    do K=I,N
      S = A(I,K)*A(K,J)+S
    end do
    A(I,J) = S
  end do
end do
! INTERCHANGE OF COLUMNS
J = NM1
200 M = MX(J)
if (M == J) GO TO 220
do I=1,N
  W = A(I,M)
  A(I,M) = A(I,J)
  A(I,J) = W
end do
220 J = J-1
if (J >= 1) GO TO 200
230 ILL = 0
return
240 ILL = J
return
250 ILL = 30000

return

end subroutine MINVD
