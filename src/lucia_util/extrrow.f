************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE EXTRROW(INMAT,IROW,NROW,NCOL,IOUTVEC)
*
* Extract row IROW from integer matrix INMAT
*
* Jeppe Olsen, Winter 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION INMAT(NROW,NCOL)
      DIMENSION IOUTVEC(NCOL)
*
      DO ICOL = 1, NCOL
        IOUTVEC(ICOL) = INMAT(IROW,ICOL)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output vector from EXTRROW '
        WRITE(6,*) ' Extracted ROW ', IROW
        CALL IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
      END IF
*
      RETURN
      END
C !!! Start trace !!!
