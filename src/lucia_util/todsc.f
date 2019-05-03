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
      SUBROUTINE TODSC(A,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY REAL*8  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8 (A-H,O-Z)
#include "io_util.fh"
      DIMENSION A(NDIM)
      INTEGER START,ISTOP
      REAL*8 INPROD
      INTEGER ISCR(2)
*
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. No packing
        ISCR(2) = 0
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
C
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      ISTOP=0
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
       START=ISTOP+1
       ISTOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       CALL DDAFILE(IFIL,1,A(START),ISTOP-START+1,IDISK(IFIL))
       CALL IDAFILE(IFIL,1,[NLABEL],1,IDISK(IFIL))
      IF(NBACK.NE.0) GOTO 100
*
 1001 CONTINUE
C
C?    write(6,*) ' leaving TODSC '
      RETURN
      END
