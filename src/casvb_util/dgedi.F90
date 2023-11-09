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
!********************************************************
!** Public-domain library routines used by casvb only. **
!********************************************************
!**********************
!** LINPACK ROUTINES **
!**********************

subroutine DGEDI(A,LDA,N,IPVT,DET,W,JOB)
! DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
! USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
!
! ON ENTRY
!
!    A       real*8(LDA, N)
!            THE OUTPUT FROM DGECO OR DGEFA.
!
!    LDA     INTEGER
!            THE LEADING DIMENSION OF THE ARRAY  A .
!
!    N       INTEGER
!            THE ORDER OF THE MATRIX  A .
!
!    IPVT    INTEGER(N)
!            THE PIVOT VECTOR FROM DGECO OR DGEFA.
!
!    W       real*8(N)
!            WORK VECTOR.  CONTENTS DESTROYED.
!
!    JOB     INTEGER
!            = 11   BOTH DETERMINANT AND INVERSE.
!            = 01   INVERSE ONLY.
!            = 10   DETERMINANT ONLY.
!
! ON RETURN
!
!    A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!            OTHERWISE UNCHANGED.
!
!    DET     real*8(2)
!            DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!            OTHERWISE NOT REFERENCED.
!            DETERMINANT = DET(1) * 10.0**DET(2)
!            WITH  1.0 <= abs(DET(1)) < 10.0
!            OR  DET(1) == 0.0 .
!
! ERROR CONDITION
!
!    A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!    A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!    IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!    AND IF DGECO HAS SET RCOND > 0.0 OR DGEFA HAS SET
!    INFO == 0 .
!
! LINPACK. THIS VERSION DATED 08/14/78 .
! CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
! SUBROUTINES AND FUNCTIONS
!
! BLAS DSWAP
! FORTRAN ABS,MOD
!
! Updated to Fortran 90+ (Sep. 2023)

use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LDA, N, IPVT(N), JOB
real(kind=wp), intent(inout) :: A(LDA,N)
real(kind=wp), intent(out) :: DET(2), W(N)
integer(kind=iwp) :: I, J, K, KB, KP1, L, NM1
real(kind=wp) :: T

! COMPUTE DETERMINANT

if (JOB/10 /= 0) then
  DET(1) = One
  DET(2) = Zero
  do I=1,N
    if (IPVT(I) /= I) DET(1) = -DET(1)
    DET(1) = A(I,I)*DET(1)
    ! ...EXIT
    if (DET(1) == Zero) exit
    do
      if (abs(DET(1)) >= One) exit
      DET(1) = Ten*DET(1)
      DET(2) = DET(2)-One
    end do
    do
      if (abs(DET(1)) < Ten) exit
      DET(1) = DET(1)/Ten
      DET(2) = DET(2)+One
    end do
  end do
end if

! COMPUTE INVERSE(U)

if (mod(JOB,10) /= 0) then
  do K=1,N
    A(K,K) = One/A(K,K)
    T = -A(K,K)
    A(1:K-1,K) = T*A(1:K-1,K)
    KP1 = K+1
    if (N < KP1) cycle
    do J=KP1,N
      T = A(K,J)
      A(K,J) = Zero
      A(1:K,J) = A(1:K,J)+T*A(1:K,K)
    end do
  end do

  ! FORM INVERSE(U)*INVERSE(L)

  NM1 = N-1
  if (NM1 >= 1) then
    do KB=1,NM1
      K = N-KB
      KP1 = K+1
      W(KP1:N) = A(KP1:N,K)
      A(KP1:N,K) = Zero
      do J=KP1,N
        T = W(J)
        A(1:N,K) = A(1:N,K)+T*A(1:N,J)
      end do
      L = IPVT(K)
      if (L /= K) call DSWAP_(N,A(1:N,K),1,A(1:N,L),1)
    end do
  end if
end if

return

end subroutine DGEDI
