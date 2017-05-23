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
      SUBROUTINE SHTOOB(  NSHPIR,  NIRREP,  MXPOBS,   NSMOB,  NOSPIR,
     &                    IOSPIR,   NOBPS,     NOB)
*
* Number of shells per irrep => Number of orbitals per symmetry
*
* =====
* Input
* =====
*
*  NSHPIR : Number of shells per irrep
*  NIRREP : Number of irreps
*  MXPOBS : Largest allowed number of orbitals symmetries
*  NSMOB  : Number of orbital symmetries
*  NOSPIR : Number of orbital symmetries per irrep
*  IOSPIR : Orbital symmetries per irrep
*
* ======
* Output
* ======
*  NOBPS  : Number of orbitals per symmetry
*  NOB    : Number of orbitals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NSHPIR(*),NOSPIR(*),IOSPIR(MXPOBS,*)
*. Output
      DIMENSION NOBPS(*)
      CALL ISETVC(NOBPS,0,NSMOB)
      NOB = 0
      DO 100 IRREP = 1, NIRREP
        DO 90 ISM = 1, NOSPIR(IRREP)
          IISM = IOSPIR(ISM,IRREP)
          NOBPS(IISM) = NOBPS(IISM) + NSHPIR(IRREP)
          NOB = NOB + NSHPIR(IRREP)
   90   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' SHTOOB Speaking '
         WRITE(6,*) ' =============== '
         WRITE(6,*) ' Number of orbitals obtained ', NOB
         WRITE(6,*) ' Number of orbitals per symmetry '
         CALL IWRTMA(NOBPS,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
