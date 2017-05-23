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

      SUBROUTINE DGMM2 (AOUT,AIN,DIAG,IWAY,NRDIM,NCDIM)
C
C PRODUCT OF DIAGONAL MATRIX AND MATRIX :
C
C     IWAY = 1 : AOUT(I,J) = DIAG(I)*AIN(I,J)
C     IWAY = 2 : AOUT(I,J) = DIAG(J)*AIN(I,J)
C
      IMPLICIT REAL*8          (A-H,O-Z)
      DIMENSION AIN(NRDIM,NCDIM),DIAG(*)
      DIMENSION AOUT(NRDIM,NCDIM)
C
      IF ( IWAY .EQ. 1 ) THEN
         DO 100 J = 1, NCDIM
           CALL VVTOV(AIN(1,J),DIAG(1),AOUT(1,J),NRDIM)
  100    CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        DO 200 J = 1, NCDIM
          FACTOR = DIAG(J)
          CALL VECSUM(AOUT(1,J),AOUT(1,J),AIN(1,J),0.0D0,
     &                FACTOR,NRDIM)
  200   CONTINUE
      END IF
C
      NTEST = 00
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AIN DIAG AOUT  FROM DGMTMT '
        CALL WRTMAT(AIN ,NRDIM,NCDIM,NRDIM,NCDIM)
        IF(IWAY.EQ.1) THEN
        CALL WRTMAT(DIAG,1   ,NRDIM,1,NRDIM)
        ELSE
        CALL WRTMAT(DIAG,1   ,NCDIM,1,NCDIM)
        END IF
        CALL WRTMAT(AOUT,NRDIM,NCDIM,NRDIM,NCDIM)
      END IF
C
      RETURN
      END
