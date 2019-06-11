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
      SUBROUTINE SKPRCD2(NDIM,MBLOCK,IFILE)
C
C     Skip record in file IFILE
C
*. Version allowing zero and packed blocks
*
* Dos not work with FASTIO - I expect
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
*
      DIMENSION ISCR(2)
c      PARAMETER(LPBLK=50000)

C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
        IF(IMZERO.EQ.1) THEN
          GOTO 1001
        END IF
      END IF
*
      IF(I_AM_PACKED.EQ.1) THEN
*. Loop over packed records of dimension LPBLK
*. The next LPBLK elements
  999   CONTINUE
*. Read next batch
          CALL IDAFILE(IFILE,2,ISCR,1,IDISK(IFILE))
          LBATCH=ISCR(1)
          IF(LBATCH.GT.0) THEN
            CALL IDAFILE(IFILE,0,[0],LBATCH,IDISK(IFILE))
            CALL DDAFILE(IFILE,0,[0.0d0],LBATCH,IDISK(IFILE))
          END IF
          CALL IDAFILE(IFILE,2,ISCR,1,IDISK(IFILE))
          ISTOP=ISCR(1)
        IF(ISTOP.EQ.0) GOTO 999
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
cvv        IF(.true.) THEN
          NBLOCK = MBLOCK
          IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
          IREST=NDIM
          IBASE=0
  100     CONTINUE
          IF(IREST.GT.NBLOCK) THEN
           CALL DDAFILE(IFILE,0,[0.0d0],NBLOCK,IDISK(IFILE))
           IBASE=IBASE+NBLOCK
           IREST=IREST-NBLOCK
          ELSE
           CALL DDAFILE(IFILE,0,[0.0d0],IREST,IDISK(IFILE))
           IREST=0
          END IF
          CALL IDAFILE(IFILE,0,[0],1,IDISK(IFILE))
        IF( IREST .GT. 0 ) GOTO 100
cvv        END IF
C
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
