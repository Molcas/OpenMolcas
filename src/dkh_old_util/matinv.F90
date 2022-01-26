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

subroutine MATINV(A,B,N,L,I_DIM)
! IF L=0 RETURNS INVERSE OF A IN A, IF L=1 SOLUTION OF AX=B IN B

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: N, L, I_DIM
real(kind=wp) :: A(I_DIM,I_DIM), B(I_DIM)
integer(kind=iwp) :: I, IC, I_N(I_DIM,2), IP(I_DIM), IR, J, K !IFG
real(kind=wp) :: AMAX, D

IR = 0
IC = 0
D = One
do I=1,N
  IP(I) = 0
end do
do I=1,N
  AMAX = Zero
  do J=1,N
    if (IP(J) > 0) GO TO 3
    if (IP(J) < 0) GO TO 4
    do K=1,N
      if (IP(K) == 1) GO TO 2
      if (IP(K) > 1) GO TO 4
      if (abs(A(J,K)) <= AMAX) GO TO 2
      IR = J
      IC = K
      AMAX = abs(A(J,K))
2     continue
    end do
3   continue
  end do
  IP(IC) = IP(IC)+1
  if (AMAX > 1.0e-30_wp) GO TO 6
4 write(u6,105)
105 format(' * '/' * ',16H SINGULAR MATRIX)
  call Abend()
6 if (IR == IC) GO TO 8
  D = -D
  do K=1,N
    AMAX = A(IR,K)
    A(IR,K) = A(IC,K)
    A(IC,K) = AMAX
  end do
  if (L == 0) GO TO 8
  AMAX = B(IR)
  B(IR) = B(IC)
  B(IC) = AMAX
8 I_N(I,1) = IR
  I_N(I,2) = IC
  AMAX = A(IC,IC)
  D = D*AMAX
  A(IC,IC) = One
  do K=1,N
    A(IC,K) = A(IC,K)/AMAX
  end do
  if (L == 0) GO TO 10
  B(IC) = B(IC)/AMAX
10 do J=1,N
    if (J == IC) GO TO 120
    AMAX = A(J,IC)
    A(J,IC) = Zero
    do K=1,N
      A(J,K) = A(J,K)-A(IC,K)*AMAX
    end do
    if (L == 0) GO TO 120
    B(J) = B(J)-B(IC)*AMAX
120 continue
  end do
end do
if (L == 1) GO TO 15
do I=1,N
  J = N+1-I
  if (I_N(J,1) == I_N(J,2)) GO TO 14
  IR = I_N(J,1)
  IC = I_N(J,2)
  do K=1,N
    AMAX = A(K,IR)
    A(K,IR) = A(K,IC)
    A(K,IC) = AMAX
  end do
14 continue
end do
15 continue

return

end subroutine MATINV
