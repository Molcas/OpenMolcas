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
      SUBROUTINE FRMDSCN(VEC,NREC,LBLK,LU)
*
* Read  VEC as multiple record file, NREC records read
*
      IMPLICIT REAL*8(A-H,O-Z)
*. OUtput
      DIMENSION VEC(*)
*
      IOFF = 1
      DO IREC = 1, NREC
        CALL IFRMDS(LREC,1,LBLK,LU)
        CALL FRMDSC(VEC(IOFF),     LREC,     LBLK,       LU,   IMZERO,
     &                IAMPACK)
        IOFF = IOFF + LREC
      END DO
*
      RETURN
      END
C !!! End trace !!!
