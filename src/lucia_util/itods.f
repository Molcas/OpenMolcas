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
      SUBROUTINE ITODS(IA,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY INTEGER IA(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
      DIMENSION IA(*)
      INTEGER START,ISTOP
*
      NBLOCK = MBLOCK
C
      IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
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
       CALL IDAFILE(IFIL,1,IA(START),ISTOP-START+1,IDISK(IFIL))
       CALL IDAFILE(IFIL,1,NLABEL,1,IDISK(IFIL))
      IF(NBACK.NE.0) GOTO 100
C
      RETURN
      END
