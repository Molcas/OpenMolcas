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
      SUBROUTINE MKDAW(IDOWN,IDAW,IPRINT)
C     PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "gugx.fh"
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='MKDAW   ')
C
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
C
C     BEGIN TO CONSTRUCT DOWN CHAIN TABLE
C
      DO IC=0,3
       IDAW(NVERT,IC)=0
      END DO
      IDAW(NVERT,4)=1
      DO IV=NVERT-1,1,-1
        ISUM=0
        DO IC=0,3
          IDAW(IV,IC)=0
          IDWN=IDOWN(IV,IC)
          IF(IDWN.EQ.0) GOTO 20
          IDAW(IV,IC)=ISUM
          ISUM=ISUM+IDAW(IDWN,4)
  20      CONTINUE
        END DO
        IDAW(IV,4)=ISUM
      END DO
C
      IF(IPRINT.GT.5) THEN
        Write(LF,*)
        Write(LF,*)' DIRECT ARC WEIGHTS:'
        DO IV=1,NVERT
          Write(LF,'(1X,I4,5X,5(1X,I6))') IV,(IDAW(IV,IC),IC=0,4)
        END DO
        Write(LF,*)
      ENDIF
      RETURN
      END
