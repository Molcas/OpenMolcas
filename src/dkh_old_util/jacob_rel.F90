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

subroutine JACOB_REL(A,B,EIG,N,RNG,IC)
! IC=1 EIGENVALUES AND VECTORS REARRANGED
!   =0 LEFT AS THEY ARE

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, IC
real(kind=wp), intent(inout) :: A(N,N)
real(kind=wp), intent(out) :: B(N,N), EIG(N)
real(kind=wp), intent(in) :: RNG
integer(kind=iwp) :: I, IL, IM, IND, J, K, L, LL, M, MM
real(kind=wp) :: COST, COST2, ENUF, ENUI, SINCS, SINT, SINT2, THR, U1, X, XY, Y

ENUI = Zero
U1 = real(N,kind=wp)
do I=1,N
  B(I,I) = One
  EIG(I) = A(I,I)
  do J=1,I-1
    B(I,J) = Zero
    B(J,I) = Zero
    ENUI = ENUI+A(I,J)*A(I,J)
  end do
end do
if (ENUI <= 0) then
  GO TO 200
else
  GO TO 10
end if
10 ENUI = sqrt(Two*ENUI)
ENUF = ENUI*RNG/U1
IND = 0
THR = ENUI
15 THR = THR/U1
20 L = 1
25 M = L+1
30 if (abs(A(M,L))-THR < 0) then
  GO TO 90
else
  GO TO 35
end if
35 IND = 1
X = Half*(EIG(L)-EIG(M))
Y = -A(M,L)/sqrt(A(M,L)*A(M,L)+X*X)
if (X < 0) then
  GO TO 40
else
  GO TO 45
end if
40 Y = -Y
45 if (Y > One) Y = One
if (Y < -One) Y = -One
XY = One-Y*Y
SINT = Y/sqrt(Two*(One+sqrt(XY)))
SINT2 = SINT*SINT
COST2 = One-SINT2
COST = sqrt(COST2)
SINCS = SINT*COST
do I=1,N
  if (I-M < 0) then
    GO TO 50
  else if (I-M == 0) then
    GO TO 80
  else
    GO TO 55
  end if
50 IM = M
  MM = I
  GO TO 60
55 IM = I
  MM = M
60 if (I-L < 0) then
    GO TO 65
  else if (I-L == 0) then
    GO TO 80
  else
    GO TO 70
  end if
65 IL = L
  LL = I
  GO TO 75
70 IL = I
  LL = L
75 X = A(IL,LL)*COST-A(IM,MM)*SINT
  A(IM,MM) = A(IL,LL)*SINT+A(IM,MM)*COST
  A(IL,LL) = X
80 X = B(I,L)*COST-B(I,M)*SINT
  B(I,M) = B(I,L)*SINT+B(I,M)*COST
  B(I,L) = X
end do
X = Two*A(M,L)*SINCS
Y = EIG(L)*COST2+EIG(M)*SINT2-X
X = EIG(L)*SINT2+EIG(M)*COST2+X
A(M,L) = (EIG(L)-EIG(M))*SINCS+A(M,L)*(COST2-SINT2)
EIG(L) = Y
EIG(M) = X
90 if (M-N == 0) then
  GO TO 100
else
  GO TO 95
end if
95 M = M+1
GO TO 30
100 if (L-M+1 == 0) then
  GO TO 110
else
  GO TO 105
end if
105 L = L+1
GO TO 25
110 if (IND-1 == 0) then
  GO TO 115
else
  GO TO 120
end if
115 IND = 0
GO TO 20
120 if (THR-ENUF <= 0) then
  GO TO 200
else
  GO TO 15
end if
200 if (IC == 0) GO TO 230
do I=1,N
  do J=I,N
    if (EIG(I)-EIG(J) <= 0) then
      GO TO 226
    else
      GO TO 210
    end if
210 X = EIG(I)
    EIG(I) = EIG(J)
    EIG(J) = X
    do K=1,N
      Y = B(K,I)
      B(K,I) = B(K,J)
      B(K,J) = Y
    end do
226 continue
  end do
end do
230 continue

return

end subroutine JACOB_REL
