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
* Copyright (C) 1987, Jeppe Olsen                                      *
*               1989, Markus P. Fuelscher                              *
************************************************************************
      SUBROUTINE DETSTR2(IDET,IASTR,IBSTR,NEL,NAEL,NBEL,
     *                  ISIGN,IWORK,IPRINT)
C
C     AUTHOR:        J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
C     MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
C                    M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
C
C     PURPOSE:
C
C    A DETERMINANT,IDET,IS GIVEN AS A SET OF OCCUPIED SPIN ORBITALS,
C    POSITIVE NUMBER INDICATES ALPHA ORBITAL AND NEGATIVE NUMBER
C    INDICATES BETA ORBITAL .
C    FIND CORRESPONDING ALPHA STRING AND BETA STRING
C
C     SUBROUTINE CALLS:
C
C     IWRTMA,ORDSTR,ICOPY
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDET(NEL)
      DIMENSION IASTR(NAEL),IBSTR(NBEL)
      DIMENSION IWORK(*)
C required length of IWORK : NEL
C
C FIRST REORDER SPIN ORBITALS IN ASCENDING SEQUENCE
C THIS WILL AUTOMATICALLY SPLIT ALPHA AND BETASTRING
C
      CALL ORDSTR(IDET,IWORK,NEL,ISIGN,IPRINT)
C
C ALPHA STRING IS LAST NAEL ORBITALS
      CALL ICOPY(NAEL,IWORK(NBEL+1),1,IASTR,1)
C
C BETA  STRING MUST BE COMPLETELY TURNED AROUND
      DO 10 IBEL = 1, NBEL
        IBSTR(IBEL) = -IWORK(NBEL+1-IBEL)
10    CONTINUE
C SIGN CHANGE FOR SWITCH OF BETA ORBITALS
      IEXPO = (NBEL*NBEL+NBEL)/2
      ISIGN = ISIGN * (-1) ** IEXPO
C
      IF( IPRINT.GT.30 ) THEN
        Write(6,*) ' INPUT DETERMINANT '
        CALL IWRTMA(IDET,1,NEL,1,NEL)
        Write(6,*) ' CORRESPONDING ALPHA STRING '
        CALL IWRTMA(IASTR,1,NAEL,1,NAEL)
        Write(6,*) ' CORRESPONDING BETA STRING '
        CALL IWRTMA(IBSTR,1,NBEL,1,NBEL)
        Write(6,*) ' ISIGN FOR SWITCH ', ISIGN
      END IF
C
      RETURN
      END
