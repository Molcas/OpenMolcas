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
      SUBROUTINE MKDAW(NVERT,IDOWN,IDAW,IPRINT)
C
C     PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
C
C
C     BEGIN TO CONSTRUCT DOWN CHAIN TABLE
C
      DO 10 IC=0,3
10      IDAW(NVERT,IC)=0
      IDAW(NVERT,4)=1
      DO 30 IV=NVERT-1,1,-1
        ISUM=0
        DO 20 IC=0,3
          IDAW(IV,IC)=0
          IDWN=IDOWN(IV,IC)
          IF(IDWN.EQ.0) GOTO 20
          IDAW(IV,IC)=ISUM
          ISUM=ISUM+IDAW(IDWN,4)
20      CONTINUE
        IDAW(IV,4)=ISUM
30    CONTINUE
C
      IF(IPRINT.GT.5) THEN
        WRITE(6,*)
        WRITE(6,*)' DIRECT ARC WEIGHTS:'
        DO 50 IV=1,NVERT
          WRITE(6,'(1X,I4,5X,5(1X,I6))') IV,(IDAW(IV,IC),IC=0,4)
50      CONTINUE
        WRITE(6,*)
      ENDIF
C
C
C     EXIT
C
      RETURN
      END
