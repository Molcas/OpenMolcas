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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE SMOST_MCLR(NSMST,NSMCI,MXPCSM,ISMOST)
*
* ISMOST(ISYM,ITOTSM) : Symmetry of an internal state if ITOTSM
*                       if symmetry of 1 string is ISYM, the
*                       symmetry of the other string is
*                       ISMOST(ISYM,ITOTSM)
*
* Jeppe Olsen , Spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISMOST(MXPCSM,MXPCSM)
*
      DO 1000 ITOTSM = 1, NSMCI
        DO 900 ISTSM  = 1, NSMST
          JSTSM = 1 + IEOR(ISTSM-1,ITOTSM-1)
          ISMOST(ISTSM,ITOTSM) = JSTSM
900     CONTINUE
1000  CONTINUE
*
      RETURN
      END
