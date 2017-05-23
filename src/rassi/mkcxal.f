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
C NOTE: THE FOLLOWING ROUTINE MAY NOT BE VECTORIZED.
C THERE IS A COMPILER BUG ON FORTRAN VERSION 2.2.0 (JUNE 1987).
C WRONG RESULTS PRODUCED EVEN ON VECTORIZED LEVEL 1.
      SUBROUTINE MKCXAL(NDIMEN,TRAL,CXAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TRAL(NDIMEN,NDIMEN),CXAL(NDIMEN,NDIMEN)
      DO I=1,NDIMEN
        DO J=I,NDIMEN
          CXAL(I,J)=0.0D00
        END DO
        CXAL(I,I)=1.0D00
      END DO
      DO K=1,NDIMEN
        DO I=1,K-1
          SUMMA=0.0D00
          DO J=1,K-1
            SUMMA=SUMMA+CXAL(I,J)*TRAL(J,K)
          END DO
          CXAL(I,K)=-(SUMMA/TRAL(K,K))
        END DO
        DO I=K,NDIMEN
          SUMMA=TRAL(I,K)
          IF (I.EQ.K) SUMMA=-1.0D00
          DO J=1,K-1
            SUMMA=SUMMA+CXAL(I,J)*TRAL(J,K)
          END DO
          CXAL(I,K)=-(SUMMA/TRAL(K,K))
        END DO
        CONTINUE
      END DO
      RETURN
      END
