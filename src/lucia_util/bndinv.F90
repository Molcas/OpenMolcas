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

subroutine BNDINV(A,EL,N,DETERM,EPSIL,ITEST,NSIZE)
! MATRIX INVERSION SUBROUTINE
! FROM "DLYTAP".

implicit real*8(A-H,O-Z)
dimension A(NSIZE,*), EL(NSIZE,*)

INDSNL = 0
if (N < 2) GO TO 140
ISL2 = 0
!K000FX = 2
if (ISL2 == 0) INDSNL = 2
if (ISL2 == 1) INDSNL = 1
!call SLITET(2,INDSNL)
!call OVERFL(K000FX)
!call DVCHK(K000FX)

! SET EL = IDENTITY MATRIX
do I=1,N
  do J=1,N
    EL(I,J) = 0.0d0
  end do
  EL(I,I) = 1.0d0
end do

! TRIANGULARIZE A, FORM EL

N1 = N-1
M = 2
do J=1,N1
  do I=M,N
    if (A(I,J) == 0.0d0) GO TO 45
    D = sqrt(A(J,J)*A(J,J)+A(I,J)*A(I,J))
    C = A(J,J)/D
    S = A(I,J)/D
    do K=J,N
      D = C*A(J,K)+S*A(I,K)
      A(I,K) = C*A(I,K)-S*A(J,K)
      A(J,K) = D
    end do
    do K=1,N
      D = C*EL(J,K)+S*EL(I,K)
      EL(I,K) = C*EL(I,K)-S*EL(J,K)
      EL(J,K) = D
    end do
45  continue
  end do
  M = M+1
end do
!call OVERFL(K000FX)
!goto (140,51),K000FX

! CALCULATE THE DETERMINANT
DETERP = A(1,1)
do I=2,N
  DETERP = DETERP*A(I,I)
end do
DETERM = DETERP
!call OVERFL(K000FX)
!goto (140,520,520),K000FX

! IS MATRIX SINGULAR
!520 continue
F = A(1,1)
E = A(1,1)
do I=2,N
  if (abs(F) < abs(A(I,I))) F = A(I,I)
  if (abs(E) > abs(A(I,I))) E = A(I,I)
end do
EPSILP = EPSIL
if (EPSILP <= 0.0d0) EPSILP = 1.0D-8
RAT = E/F
if (abs(RAT) < EPSILP) GO TO 130

! INVERT TRIANGULAR MATRIX
J = N
do J1=1,N
  A(J,J) = 1.0d0/A(J,J)
  I = J-1
  do I1=2,J
    D = 0.0d0
    do K=I+1,J
      D = D+A(I,K)*A(K,J)
    end do
    A(I,J) = -D/A(I,I)
    I = I-1
  end do
  J = J-1
end do
!call OVERFL(K000FX)
!goto (140,103,103),K000FX

!103 call DVCHK(K000FX)
!goto (140,105),K000FX

! PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
M = 1
do I=1,N
  do J=1,N
    D = 0.0d0
    do K=M,N
      D = D+A(I,K)*EL(K,J)
    end do
    EL(I,J) = D
  end do
  M = M+1
end do
!call OVERFL(K000FX)
!goto (140,123,123),K000FX

! RECOPY EL TO A
do I=1,N
  do J=1,N
    A(I,J) = EL(I,J)
  end do
end do
ITEST = 0
!126 if (INDSNL == 1) call SLITE(2)
126 if (INDSNL == 1) ISL2 = 1

return

130 ITEST = 1
goto 126
140 ITEST = -1
goto 126

end subroutine BNDINV
