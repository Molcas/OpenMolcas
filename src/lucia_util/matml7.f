************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2003, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE MATML7(       C,       A,       B,   NCROW,   NCCOL,
     &                     NAROW,   NACOL,   NBROW,   NBCOL, FACTORC,
     &                  FACTORAB, ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(T) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(T) FOR ITRNSP = 2
C
C     C =  FACTORC*C + FACTORAB* A(T) * B(T) FOR ITRNSP = 3
*
* Warning ITRNSP = 3 should only be used for small matrices,
* as this path involves notunit strides in inner loops.
* As Lasse points out, it is better to calculate C(T) = BA
* and then transpose C
*
C... JEPPE OLSEN,
*
* ITRNSP = 3 added, march 2003
C
*. Notice : If the summation index has dimension zero nothing
*           is performed
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)

      IF(NAROW.eq.0.or.NACOL.eq.0.or.NBROW.eq.0.or.
     &   NBCOL.eq.0.or.NCROW.eq.0.or.NCCOL .EQ. 0) THEN
        IZERO = 1
      ELSE
        IZERO = 0
      END IF
*
      IF(IZERO.EQ.1.AND.NCROW*NCCOL.NE.0) THEN
        IF(FACTORC.NE.0.0D0) THEN
         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
        ELSE
         ZERO = 0.0D0
         CALL SETVEC(C,ZERO,NCROW*NCCOL)
        END IF
      END IF
*
      IF (IZERO.EQ.0 ) THEN
*. DGEMM from CONVEX/ESSL  lib
        LDA = MAX(1,NAROW)
        LDB = MAX(1,NBROW)
*
        LDC = MAX(1,NCROW)
        IF(ITRNSP.EQ.0) THEN
        CALL DGEMM_(      'N',      'N',    NAROW,    NBCOL,    NACOL,
     &               FACTORAB,        A,      LDA,        B,      LDB,
     &                FACTORC,        C,      LDC)
        ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM_(      'T',      'N',    NACOL,    NBCOL,    NAROW,
     &               FACTORAB,        A,      LDA,        B,      LDB,
     &                FACTORC,        C,      LDC)
        ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM_(      'N',      'T',    NAROW,    NBROW,    NACOL,
     &               FACTORAB,        A,      LDA,        B,      LDB,
     &                FACTORC,        C,      LDC)
        END IF

c here if iZERO == 1
      ELSE
* Use Jeppes version ( it should be working )
        IF( ITRNSP .EQ. 0 ) THEN
* ======
* C=A*B
* ======
*
C         CALL SCALVE(C,FACTORC,NCROW*NCCOL)
          DO J =1, NCCOL
*. Initialize with FACTORC*C(I,J) + FACTORAB*A(I,1)*B(1,J)
            IF(NBROW.GE.1) THEN
              B1J = FACTORAB*B(1,J)
              DO I = 1, NCROW
                C(I,J) = FACTORC*C(I,J) + B1J*A(I,1)
              END DO
            END IF
*. and the major part
            DO K =2, NBROW
              BKJ = FACTORAB*B(K,J)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BKJ*A(I,K)
              END DO
            END DO
          END DO
*
        END IF
        IF ( ITRNSP .EQ. 1 ) THEN
*
* =========
* C=A(T)*B
* =========
*
          DO J = 1, NCCOL
            DO I = 1, NCROW
              T = 0.0D0
              DO K = 1, NBROW
                T = T  + A(K,I)*B(K,J)
              END DO
              C(I,J) = FACTORC*C(I,J) + FACTORAB*T
            END DO
          END DO
        END IF
C
        IF ( ITRNSP .EQ. 2 ) THEN
C ===========
C. C = A*B(T)
C ===========
          DO J = 1,NCCOL
*. Initialization
            IF(NBCOL.GE.1) THEN
              BJ1 = FACTORAB*B(J,1)
              DO I = 1, NCROW
                C(I,J) = FACTORC*C(I,J) + BJ1*A(I,1)
              END DO
            END IF
*. And the rest
            DO K = 2,NBCOL
              BJK = FACTORAB*B(J,K)
              DO I = 1, NCROW
                C(I,J) = C(I,J) + BJK*A(I,K)
              END DO
            END DO
          END DO
        END IF
      END IF
*
c end of iZero ==1

      IF(ITRNSP.EQ.3) THEN
C ================
C. C = A(T)*B(T)
C ================
* C(I,J) = FACTORC*C(I,J) + FACTORAB*sum(K) A(K,I)*B(J,K)
        CALL SCALVE(C,FACTORC,NCROW*NCCOL)
        DO I = 1, NCROW
          DO K = 1, NAROW
            AKI = FACTORAB*A(K,I)
            DO J = 1,NBROW
              C(I,J) = C(I,J) + AKI*B(J,K)
            END DO
          END DO
        END DO
      END IF
C
c      IF ( NTEST .NE. 0 ) THEN
c      WRITE(6,*)
c      WRITE(6,*) ' C MATRIX FROM MATML7 '
c      WRITE(6,*)
c      CALL WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
c      END IF
C
C
      RETURN
      END
