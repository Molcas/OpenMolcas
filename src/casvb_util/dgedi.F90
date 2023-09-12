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
!  ********************************************************
!  ** Public-domain library routines used by casvb only. **
!  ********************************************************
!  **********************
!  ** LINPACK ROUTINES **
!  **********************
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,W,JOB)
      INTEGER LDA,N,IPVT(N),JOB
      REAL*8 A(LDA,N),DET(2),W(*)
!
!     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
!
!     ON ENTRY
!
!        A       real*8(LDA, N)
!                THE OUTPUT FROM DGECO OR DGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.
!
!        W       real*8(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     real*8(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. abs(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
!        INFO .EQ. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,DSWAP
!     FORTRAN ABS,MOD
!
!     INTERNAL VARIABLES
!
      REAL*8 T
      REAL*8 TEN
      INTEGER I,J,K,KB,KP1,L,NM1
!
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
!        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (abs(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (abs(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(U)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL_(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY_(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM INVERSE(U)*INVERSE(L)
!
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               W(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = W(J)
               CALL DAXPY_(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP_(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
