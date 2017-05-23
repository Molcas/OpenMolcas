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
      SUBROUTINE TODSCN(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSC(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,0)
        END IF
      END DO
*
      RETURN
      END

*
