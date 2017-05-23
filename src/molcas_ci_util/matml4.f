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
      SUBROUTINE MATML4(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     *                  NBROW,NBCOL,ITRNSP )
C
C     MULTIPLY A AND B TO GIVE C
C
C     C = A * B             FOR ITRNSP = 0
C     C = A(TRANSPOSED) * B FOR ITRNSP = 1
C     C = A * B(TRANSPOSED) FOR ITRNSP = 2
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION A(NAROW,NACOL)
      DIMENSION B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
C
      IZERO = 0
      IF ( (NAROW*NACOL*NBROW*NBCOL*NCROW*NCCOL).EQ.0 ) IZERO = 1
C
      IF ( ITRNSP.EQ.0 ) THEN
        IF ( IZERO.EQ.1 ) THEN
          CALL DCOPY_(NCROW*NCCOL,0.0D0,0,C,1)
          DO J = 1,NCCOL
            DO K = 1,NBROW
              BKJ = B(K,J)
              CALL DAXPY_(NCROW,BKJ,A(1,K),1,C(1,J),1)
            END DO
          END DO
        ELSE
          CALL DGEMM_('N','N',
     &                NCROW,NCCOL,NACOL,
     &                1.0d0,A,NAROW,
     &                B,NBROW,
     &                0.0d0,C,NCROW)
        END IF
      END  IF
C
      IF ( ITRNSP.EQ.1 ) THEN
        IF ( IZERO.EQ.1 ) THEN
          DO J = 1, NCCOL
            DO I = 1, NCROW
              C(I,J) = DDOT_(NBROW,A(1,I),1,B(1,J),1)
            END DO
          END DO
        ELSE
          CALL DGEMM_('T','N',
     &                NCROW,NCCOL,NAROW,
     &                1.0d0,A,NAROW,
     &                B,NBROW,
     &                0.0d0,C,NCROW)
        END IF
      END IF
C
      IF ( ITRNSP.EQ.2 ) THEN
        IF ( IZERO.EQ.1 ) THEN
          CALL DCOPY_(NCROW*NCCOL,0.0D0,0,C,1)
          DO J = 1,NCCOL
            DO K = 1,NBCOL
              BJK = B(J,K)
              CALL DAXPY_(NCROW,BJK,A(1,K),1,C(1,J),1)
            END DO
          END DO
        ELSE
          CALL DGEMM_('N','T',
     &                NCROW,NCCOL,NACOL,
     &                1.0d0,A,NAROW,
     &                B,NBROW,
     &                0.0d0,C,NCROW)
        END IF
      END IF
C
      RETURN
      END
