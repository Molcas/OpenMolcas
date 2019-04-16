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
      SUBROUTINE COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILE
*
* Structure of output file is inherited by output file,
* if input file is packed, so is output file
*
*
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
      DIMENSION SEGMNT(*),LBL(1),IDUMMY(1)
C
      IF( IREW .NE. 0 ) THEN
        IDISK(LUIN)=0
        IDISK(LUOUT)=0
      END IF

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL(1) = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          CALL IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
          CALL IDAFILE(LUOUT,1,LBL,1,IDISK(LUOUT))
C?        write(6,*) ' COPVCD LBL : ', LBL(1)
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IDAFILE(LUIN,2,LBL,1,IDISK(LUIN))
          CALL IDAFILE(LUIN,2,IDUMMY,1,IDISK(LUIN))
          CALL IDAFILE(LUOUT,1,LBL,1,IDISK(LUOUT))
          CALL IDAFILE(LUOUT,1,[-1],1,IDISK(LUOUT))
        END IF
        IF( LBL(1) .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL(1)
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL(1),KBLK
          NO_ZEROING = 1
          CALL FRMDSC2(  SEGMNT,  LBL(1),    KBLK,    LUIN,  IMZERO,
     &                  IAMPACK,NO_ZEROING)
          IF(IAMPACK.NE.0) THEN
C?          WRITE(6,*) ' COPVCD, IAMPACK,FILE = ', IAMPACK,LUIN
          END IF
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC (SEGMNT,LBL(1),KBLK,LUOUT)
            ELSE
              CALL TODSCP(SEGMNT,LBL(1),KBLK,LUOUT)
            END IF
          ELSE
            CALL ZERORC(LBL(1),LUOUT,IAMPACK)
          END IF
        END IF
      IF( LBL(1) .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
