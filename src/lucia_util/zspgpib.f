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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ZSPGPIB(NSTSO,ISTSO,NSPGP,NSMST)
*
* Offset for supergroups of strings with given sym.
*. Each supergroup is given start address 1
*
* Jeppe Olsen, Still summer of 95
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER NSTSO(NSMST,NSPGP)
*. Output
      INTEGER ISTSO(NSMST,NSPGP)
*
      DO ISPGP = 1, NSPGP
        ISTSO(1,ISPGP) = 1
        DO ISMST = 2, NSMST
          ISTSO(ISMST,ISPGP) = ISTSO(ISMST-1,ISPGP) + NSTSO(ISMST,ISPGP)
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ZSPGPIB '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        CALL IWRTMA(ISTSO,NSMST,NSPGP,NSMST,NSPGP)
      END IF
*
      RETURN
      END
