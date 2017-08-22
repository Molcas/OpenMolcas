************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE DGEICD(A,LDA,N,IOPT,RCOND,DET,AUX,NAUX)
C
C     COMPUTE THE INVERSE, THE RECIPROCAL CONDITION NUMBER AND THE
C     DETERMINANT OF MATRIX A
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
      DIMENSION A(LDA,N),AUX(NAUX),DET(2)
      DIMENSION IPVT(16),Z(16),TEMP2(2)
C
      IF ( N.LT.0 .OR. N.GT.LDA ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') 'ORDER OF MATRIX A IS OUT OF BOUNDS'
         WRITE(6,*)
         Call Abend()
      END IF
      IF ( IOPT.LT.0 .OR. IOPT.GT.3 ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') '    OPTION KEY IS OUT OF BOUNDS'
         WRITE(6,*)
         Call Abend()
      END IF
      IF ( N.GE.16 .AND. NAUX.LT.2*N ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') '      WORK AREA IS TO SMALL'
         WRITE(6,*)
         Call Abend()
      END IF
C
      IF ( N.LT.16 ) THEN
         CALL DGECO(A,LDA,N,IPVT,TEMP1,Z)
      ELSE
C        CALL DGECO(A,LDA,N,AUX(1),TEMP1,AUX(N+1))
*
*------- This trick to have the correct data type for the 4th argument.
*
         ipAux=ip_of_iWork(Aux(1))
         CALL DGECO(A,LDA,N,iWork(ipAux),TEMP1,AUX(N+1))
      END IF
      IF ( (1.0D0+TEMP1).EQ.1.0D0 ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') '      THIS A SINGULAR MATRIX'
         WRITE(6,*)
         Call Abend()
      END IF
      IF ( IOPT.EQ.1 .OR. IOPT.EQ.3 ) RCOND=TEMP1
      JOB=11
      IF ( N.LT.16 ) THEN
         CALL DGEDI(A,LDA,N,IPVT,TEMP2,Z,JOB)
      ELSE
C        CALL DGEDI(A,LDA,N,AUX(1),TEMP2,AUX(N+1),JOB)
*
*------- This trick to have the correct data type for the 4th argument.
*
         ipAux=ip_of_iWork(Aux(1))
         CALL DGEDI(A,LDA,N,iWork(ipAux),TEMP2,AUX(N+1),JOB)
      END IF
      IF ( IOPT.EQ.2 .OR. IOPT.EQ.3 ) THEN
         DET(1)=TEMP2(1)
         DET(2)=TEMP2(2)
      END IF
C
      RETURN
      END
************************************************************************
* This file from LINPACK:                                              *
*   http://www.netlib.org/linpack/                                     *
*                                                                      *
* To the best of our knowledge, the routines in LINPACK are public     *
* domain or freely distributable and modifiable.                       *
************************************************************************

      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(N)
      REAL*8 A(LDA,N),Z(N)
      REAL*8 RCOND
C
C     DGECO FACTORS A real*8 MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       real*8(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   real*8
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       real*8(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN ABS,MAX,SIGN
C
C     INTERNAL VARIABLES
C
      REAL*8 DDOT_,EK,T,WK,WKM
      REAL*8 ANORM,S,DASUM_,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = MAX(ANORM,DASUM_(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = sign(EK,-Z(K))
         IF (abs(EK-Z(K)) .LE. abs(A(K,K))) GO TO 30
            S = abs(A(K,K))/abs(EK-Z(K))
            CALL DSCAL_(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = abs(WK)
         SM = abs(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + abs(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + abs(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM_(N,Z,1)
      CALL DSCAL_(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT_(N-K,A(K+1,K),1,Z(K+1),1)
         IF (abs(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/abs(Z(K))
            CALL DSCAL_(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM_(N,Z,1)
      CALL DSCAL_(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY_(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (abs(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/abs(Z(K))
            CALL DSCAL_(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM_(N,Z,1)
      CALL DSCAL_(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (abs(Z(K)) .LE. abs(A(K,K))) GO TO 150
            S = abs(A(K,K))/abs(Z(K))
            CALL DSCAL_(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY_(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM_(N,Z,1)
      CALL DSCAL_(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
************************************************************************
* This file from LINPACK:                                              *
*   http://www.netlib.org/linpack/                                     *
*                                                                      *
* To the best of our knowledge, the routines in LINPACK are public     *
* domain or freely distributable and modifiable.                       *
************************************************************************

      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(N),JOB
      REAL*8 A(LDA,N),DET(2),WORK(*)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       real*8(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    real*8(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     real*8(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. abs(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN ABS,MOD
C
C     INTERNAL VARIABLES
C
      REAL*8 T
      REAL*8 TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
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
C
C     COMPUTE INVERSE(U)
C
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
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY_(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP_(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
