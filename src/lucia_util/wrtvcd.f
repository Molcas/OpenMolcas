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
      SUBROUTINE WRTVCD(SEGMNT,LU,IREW,LBLK)
C
C PRINT VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
      DIMENSION SEGMNT(*),IDUMMY(1)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          IDISK(LU)=0
        ELSE
          IDISK(LU)=0
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          CALL IDAFILE(LU,2,IDUMMY,1,IDISK(LU))
          LBL=IDUMMY(1)
        ELSE
          CALL IDAFILE(LU,2,IDUMMY,1,IDISK(LU))
          LBL=IDUMMY(1)
          CALL IDAFILE(LU,2,IDUMMY,1,IDISK(LU))
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(  SEGMNT,    LBL ,    KBLK,      LU,  IMZERO,
     &                  IAMPACK)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMAT(SEGMNT,1,LBL,1,LBL)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
