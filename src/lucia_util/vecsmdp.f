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
      SUBROUTINE VECSMDP(    VEC1,    VEC2,    FAC1,    FAC2,     LU1,
     &                        LU2,     LU3,    IREW,    LBLK)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
*
* Packed version, May 1996
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        IDISK(LU1)=0
        IDISK(LU2)=0
        IDISK(LU3)=0
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          CALL IDAFILE(LU1,2,NBL1,1,IDISK(LU1))
          CALL IDAFILE(LU2,2,NBL2,1,IDISK(LU2))
          CALL IDAFILE(LU3,1,NBL1,1,IDISK(LU3))
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IDAFILE(LU1,2,NBL1,1,IDISK(LU1))
          CALL IDAFILE(LU1,2,IDUMMY,1,IDISK(LU1))
          CALL IDAFILE(LU2,2,NBL2,1,IDISK(LU2))
          CALL IDAFILE(LU2,2,IDUMMY,1,IDISK(LU2))
          CALL IDAFILE(LU3,1,NBL1,1,IDISK(LU3))
          CALL IDAFILE(LU3,1,-1,1,IDISK(LU3))
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &  NBL1,NBL2
*       STOP ' INCOMPATIBLE BLOCKSIZES IN VECSMF '
        CALL SYSABENDMSG('lucia_util/vecsmf','Different block sizes',
     &                   ' ')
      END IF
C
      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        NO_ZEROING = 1
        CALL FRMDSC2(     VEC1,     NBL1,     KBLK,      LU1,  IMZERO1,
     &                 IAMPACK,NO_ZEROING)
        CALL FRMDSC2(     VEC2,     NBL1,     KBLK,      LU2,  IMZERO2,
     &                 IAMPACK,NO_ZEROING)
        IF( NBL1 .GT. 0 ) THEN
          IF(IMZERO1.EQ.1.AND.IMZERO2.EQ.1) THEN
*. Simple zero record
            CALL ZERORC(NBL1,LU3,IAMPACK)
          ELSE
*. Nonvanishing record
            ZERO = 0.0D0
            IF(IMZERO1.EQ.1) THEN
              CALL VECSUM(    VEC1,    VEC1,    VEC2,    ZERO,    FAC2,
     &                        NBL1)
            ELSE IF(IMZERO2.EQ.1) THEN
              CALL VECSUM(    VEC1,    VEC1,    VEC2,    FAC1,    ZERO,
     &                        NBL1)
            ELSE
              CALL VECSUM(    VEC1,    VEC1,    VEC2,    FAC1,    FAC2,
     &                        NBL1)
            END IF
            CALL TODSCP(VEC1,NBL1,KBLK,LU3)
          END IF
        ELSE IF (NBL1.EQ.0) THEN
          CALL TODSCP(VEC1,NBL1,KBLK,LU3)
        END IF
      END IF
C
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
