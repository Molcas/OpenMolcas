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
if (ENUI > Zero) then
  ENUI = sqrt(Two*ENUI)
  ENUF = ENUI*RNG/U1
  IND = 0
  THR = ENUI
  outer: do
    THR = THR/U1
    L = 1
    M = L+1
    do
      if (abs(A(M,L))-THR >= Zero) then
        IND = 1
        X = Half*(EIG(L)-EIG(M))
        Y = -A(M,L)/sqrt(A(M,L)*A(M,L)+X*X)
        if (X < Zero) Y = -Y
        if (Y > One) Y = One
        if (Y < -One) Y = -One
        XY = One-Y*Y
        SINT = Y/sqrt(Two*(One+sqrt(XY)))
        SINT2 = SINT*SINT
        COST2 = One-SINT2
        COST = sqrt(COST2)
        SINCS = SINT*COST
        do I=1,N
          if ((I /= M) .and. (I /= L)) then
            IM = max(M,I)
            MM = min(M,I)
            IL = max(L,I)
            LL = min(L,I)
            X = A(IL,LL)*COST-A(IM,MM)*SINT
            A(IM,MM) = A(IL,LL)*SINT+A(IM,MM)*COST
            A(IL,LL) = X
          end if
          X = B(I,L)*COST-B(I,M)*SINT
          B(I,M) = B(I,L)*SINT+B(I,M)*COST
          B(I,L) = X
        end do
        X = Two*A(M,L)*SINCS
        Y = EIG(L)*COST2+EIG(M)*SINT2-X
        X = EIG(L)*SINT2+EIG(M)*COST2+X
        A(M,L) = (EIG(L)-EIG(M))*SINCS+A(M,L)*(COST2-SINT2)
        EIG(L) = Y
        EIG(M) = X
      end if
      if (M /= N) then
        M = M+1
      else if (M /= L+1) then
        L = L+1
        M = L+1
      else if (IND == 1) then
        IND = 0
        L = 1
        M = L+1
      else if (THR-ENUF <= Zero) then
        exit outer
      else
        exit
      end if
    end do
  end do outer
end if
if (IC /= 0) then
  do I=1,N
    do J=I,N
      if (EIG(I)-EIG(J) <= Zero) cycle
      X = EIG(I)
      EIG(I) = EIG(J)
      EIG(J) = X
      do K=1,N
        Y = B(K,I)
        B(K,I) = B(K,J)
        B(K,J) = Y
      end do
    end do
  end do
end if

return

end subroutine JACOB_REL
