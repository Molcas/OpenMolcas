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
* Copyright (C) 1997, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE NEXT_SYM_DISTR(  NGAS, MINVAL, MAXVAL,   ISYM,ISYM_TOT,
     &                           IFIRST,  NONEW)
*
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION MINVAL(NGAS),MAXVAL(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*
*. Symmetry of first NGAS -1 spaces
*
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 1, NGAS-1
          ISYM(IGAS) = MINVAL(IGAS)
        END DO
        NONEW = 0
      END IF
 1001 CONTINUE
      IF(IFIRST.EQ.0) CALL NXTNUM3(ISYM,NGAS-1,MINVAL,MAXVAL,NONEW)
      IFIRST = 0
*
*. Symmetry of last space
*
      IF(NONEW.EQ.0) THEN
C       JSYM = 1
C       DO IGAS = 1, NGAS-1
C         CALL SYMCOM(3,0,JSYM,ISYM(IGAS),KSYM)
C         JSYM = KSYM
C       END DO
        JSYM = ISYMSTR(ISYM,NGAS-1)
        CALL SYMCOM(2,0,JSYM,ISYM(NGAS),ISYM_TOT)
*
        IF(MINVAL(NGAS).GT.ISYM(NGAS).OR.
     &     MAXVAL(NGAS).LT.ISYM(NGAS)    )GOTO 1001
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
*
