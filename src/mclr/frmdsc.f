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
      SUBROUTINE FRMDSC_MCLR(ARRAY,NDIM,MBLOCK,IFILE,IMZERO)
C
C     TRANSFER ARRAY FROM DISC FILE IFILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY(*)
      DIMENSION IDUM(1)
C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(IDUM,1,MBLOCK,IFILE)
        IMZERO=IDUM(1)
        IF(IMZERO.EQ.1) THEN
          ZERO = 0.0D0
*         CALL SETVEC(ARRAY,ZERO,NDIM)
          call dcopy_(NDIM,[ZERO],0,ARRAY,1)
          GOTO 1001
        END IF
      END IF
*
      ICRAY = 1

      IF( MBLOCK .GE. 0 .OR.ICRAY.EQ.1) THEN
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      IREST=NDIM
      IBASE=0
  100 CONTINUE
       IF(IREST.GT.NBLOCK) THEN
        READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
        IBASE=IBASE+NBLOCK
        IREST=IREST-NBLOCK
       ELSE
        READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
        IREST=0
       END IF
      IF( IREST .GT. 0 ) GOTO 100
      END IF
C
      IF( MBLOCK.LT.0.AND.NDIM.GT.0.AND.ICRAY.EQ.0 ) THEN
*      CALL SQFILE(IFILE,2,ARRAY,2*NDIM)
       Call SysHalt('frmdsc')
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
