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

subroutine BNDINV(A,EL,N,DETERM,EPSIL,ITEST)
! MATRIX INVERSION SUBROUTINE
! FROM "DLYTAP".

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: A(N,N)
real(kind=wp), intent(out) :: EL(N,N), DETERM
real(kind=wp), intent(in) :: EPSIL
integer(kind=iwp), intent(out) :: ITEST
integer(kind=iwp) :: I, I1, J, J1, K, M, N1
real(kind=wp) :: C, D, E, EPSILP, F, RAT, S

!INDSNL = 0
if (N < 2) then
  ITEST = -1
  return
end if
!ISL2 = 0
!K000FX = 2
!if (ISL2 == 0) INDSNL = 2
!if (ISL2 == 1) INDSNL = 1
!call SLITET(2,INDSNL)
!call OVERFL(K000FX)
!call DVCHK(K000FX)

! SET EL = IDENTITY MATRIX
call unitmat(EL,N)

! TRIANGULARIZE A, FORM EL

N1 = N-1
M = 2
do J=1,N1
  do I=M,N
    if (A(I,J) == Zero) cycle
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
  end do
  M = M+1
end do
!call OVERFL(K000FX)
!if (K000FX == 0) then
!  ITEST = -1
!  goto 126
!end if

! CALCULATE THE DETERMINANT
DETERM = A(1,1)
do I=2,N
  DETERM = DETERM*A(I,I)
end do
!call OVERFL(K000FX)
!if (K000FX == 0) then
!  ITEST = -1
!  goto 126
!end if

! IS MATRIX SINGULAR
F = A(1,1)
E = A(1,1)
do I=2,N
  if (abs(F) < abs(A(I,I))) F = A(I,I)
  if (abs(E) > abs(A(I,I))) E = A(I,I)
end do
EPSILP = EPSIL
if (EPSILP <= Zero) EPSILP = 1.0e-8_wp
RAT = E/F
if (abs(RAT) < EPSILP) then
  ITEST = 1
  return
end if

! INVERT TRIANGULAR MATRIX
J = N
do J1=1,N
  A(J,J) = One/A(J,J)
  I = J-1
  do I1=2,J
    D = sum(A(I,I+1:J)*A(I+1:J,J))
    A(I,J) = -D/A(I,I)
    I = I-1
  end do
  J = J-1
end do
!call OVERFL(K000FX)
!if (K000FX == 0) then
!  ITEST = -1
!  goto 126
!end if

!call DVCHK(K000FX)
!if (K000FX == 0) then
!  ITEST = -1
!  goto 126
!end if

! PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
M = 1
do I=1,N
  do J=1,N
    D = sum(A(I,M:N)*EL(M:N,J))
    EL(I,J) = D
  end do
  M = M+1
end do
!call OVERFL(K000FX)
!if (K000FX == 0) then
!  ITEST = -1
!  goto 126
!end if

! RECOPY EL TO A
A(:,:) = EL(:,:)
ITEST = 0
!126 if (INDSNL == 1) call SLITE(2)
!126 if (INDSNL == 1) ISL2 = 1

end subroutine BNDINV
