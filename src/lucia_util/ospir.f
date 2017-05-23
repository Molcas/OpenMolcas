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
      SUBROUTINE OSPIR(  NOSPIR,  IOSPIR,  PNTGRP,  NIRREP,  MXPIRR,
     &                   MXPOBS,   IPRNT)
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   PNTGRP  : type of pointgroup
*         = 1 => D2h or a subgroup of D2H
*         = 2 => C inf v
*         = 3 => D inf h
*         = 4 => O 3
*   NIRREP : Number of irreducible representations per point group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPOBS : Largest allowed number of orbital symmetries
*
* ======
* Output
* ======
*
*   NOSPIR : Number of orbital symmetries per irrep
*   IOSPIR : Orbital symmetries corresponding to a given irrep
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER PNTGRP
*. Output
      DIMENSION NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR)
*
      IF(PNTGRP.EQ.1) THEN
*=====
*.D2h
*=====
        NSMOB = 0
        DO 10 IRREP = 1, 8
          NOSPIR(IRREP) = 1
          IOSPIR(1,IRREP) = IRREP
   10   CONTINUE
      ELSE
        WRITE(6,*) ' Sorry  PNTGRP out of range , PNTGRP = ', PNTGRP
        WRITE(6,*) ' OSPIR fatally wounded '
*        STOP 5
        CALL SYSABENDMSG('lucia_util/ospir','Internal error',' ')
      END IF
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' OSPIR speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Number of orbitals per irrep '
        CALL IWRTMA(NOSPIR,1,NIRREP,1,NIRREP)
        WRITE(6,*) ' Orbital symmetries per irrep '
        DO 100 IRREP = 1, NIRREP
          CALL IWRTMA(IOSPIR(1,IRREP),1,NOSPIR(IRREP),1,NOSPIR(IRREP))
  100   CONTINUE
      END IF
*
      RETURN
      END
