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
      SUBROUTINE TODSC_MCLR(A,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      INTEGER START,STOP
*
C?    write(6,*) ' entering TODSC '
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = ddot_(nDim,A,1,A,1)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
C?      WRITE(6,*) ' I am going to call ITODS'
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.1) MMBLOCK = 1
        CALL ITODS(IMZERO,1,MMBLOCK,IFIL)
C?      WRITE(6,*) ' back from ITODS '
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      ICRAY = 1
      IF( MBLOCK .GE.0 .OR.ICRAY .EQ. 1 ) THEN
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (A(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
      END IF
C
      IF( ICRAY.EQ.0.AND.MBLOCK.LT.0.AND.NDIM.GT.0) THEN
*      CALL SQFILE(IFIL,1,A,2*NDIM)
       Call SysHalt('todsc')
      END IF
*
 1001 CONTINUE
C
C?    write(6,*) ' leaving TODSC '
      RETURN
      END
