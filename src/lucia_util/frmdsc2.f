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
      SUBROUTINE FRMDSC2(   ARRAY,    NDIM,  MBLOCK,   IFILE,  IMZERO,
     &                   I_AM_PACKED,NO_ZEROING)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
*. Version allowing zero and packed blocks
*
* If NO_ZEROING = 1, the elements of zero blocks
*    are not set to zero, the routine just returns with
*    IMZERO = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
      DIMENSION ARRAY(*)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)

      IMZERO = 0
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
         MMBLOCK = MBLOCK
C       IF(MMBLOCK.GE.2) MMBLOCK = 2
C       CALL IFRMDS(ISCR,2,MMBLOCK,IFILE)
         CALL IFRMDS(ISCR,2,2,IFILE)
         IMZERO=ISCR(1)
         I_AM_PACKED=ISCR(2)
         IF(IMZERO.EQ.1) THEN
            IF(NO_ZEROING.EQ.0) THEN
               ZERO = 0.0D0
               CALL SETVEC(ARRAY,ZERO,NDIM)
            END IF
            GOTO 1001
         END IF
      END IF
*
      IF(I_AM_PACKED.EQ.1) THEN
         ZERO = 0.0D0
         CALL SETVEC(ARRAY,ZERO,NDIM)
*. Loop over packed records of dimension LPBLK
         NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements
         LBATCH=-2**30
 999     CONTINUE
         NBATCH = NBATCH + 1
         IF(NBATCH.NE.1) THEN
            LBATCHP = LBATCH
         END IF
*. Read next batch
         CALL IDAFILE(IFILE,2,LBATCH,1,IDISK(IFILE))
         IF(LBATCH.GT.0) THEN
           CALL IDAFILE(IFILE,2,IPAK,LBATCH,IDISK(IFILE))
           CALL DDAFILE(IFILE,2,XPAK,LBATCH,IDISK(IFILE))
         END IF
         CALL IDAFILE(IFILE,2,ISTOP,1,IDISK(IFILE))
         DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
               WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT
               WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
               WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE
               IF(NBATCH.EQ.1) THEN
                  WRITE(6,*) ' NBATCH = 1 '
               ELSE
                  WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
               END IF
               WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
*              STOP ' problem in FRMDSC '
               CALL SYSABENDMSG('lucia_util/frmdsc','Internal error',
     &                          ' ')
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
         END DO
         IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
            NBLOCK = MBLOCK
            IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
            IREST=NDIM
            IBASE=0
 100        CONTINUE
            IF(IREST.GT.NBLOCK) THEN
               CALL DDAFILE(IFILE,2,ARRAY(IBASE+1),NBLOCK,IDISK(IFILE))
               IBASE=IBASE+NBLOCK
               IREST=IREST-NBLOCK
            ELSE
               CALL DDAFILE(IFILE,2,ARRAY(IBASE+1),IREST,IDISK(IFILE))
               IREST=0
            END IF
            CALL IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
            IF( IREST .GT. 0 ) GOTO 100
C
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
