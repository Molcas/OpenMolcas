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
      SUBROUTINE SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
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
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM(2,1,ISTSM,JSTSM,ITOTSM)
        ISMOST(ISTSM,ITOTSM) = JSTSM
  900  CONTINUE
 1000 CONTINUE
*
      NTEST = 0
      IF( NTEST.NE. 0 ) THEN
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' Info from SMOST '
        WRITE(6,*) ' ==============='
        DO 1010 ITOTSM = 1, NSMCI
          WRITE(6,*) ' ISMOST array for ITOTSM = ', ITOTSM
          CALL IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
 1010   CONTINUE
      END IF
*
      RETURN
      END
