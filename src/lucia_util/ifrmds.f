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
      SUBROUTINE IFRMDS(IARRAY,NDIM,MBLOCK,IFILE)
C
C     TRANSFER INTEGER ARRAY FROM DISC FILE IFILE
C
C NBLOCK .LT. 0 INDICATES USE OF FASTIO
C
C If nblock .eq. 0 NBLOCK = NDIM
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
      DIMENSION IARRAY(*)
C
      NBLOCK = MBLOCK

C       DO NOT USE FASTIO
        IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
          IF(IREST.GT.NBLOCK) THEN
            CALL IDAFILE(IFILE,2,IARRAY(IBASE+1),NBLOCK,IDISK(IFILE))
            IBASE=IBASE+NBLOCK
            IREST=IREST-NBLOCK
          ELSE
            CALL IDAFILE(IFILE,2,IARRAY(IBASE+1),IREST,IDISK(IFILE))
            IREST=0
          END IF
          CALL IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
        IF( IREST .GT. 0 ) GOTO 100
cvv      ELSE
C       USE FAST IO
cvv        CALL SQFILE(IFILE,2,IARRAY,NDIM)
cvv      END IF
      RETURN
      END
