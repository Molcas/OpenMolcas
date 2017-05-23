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
* Copyright (C) 1989, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GETCNF(KCNF,KTYP,K,ICONF,IREFSM,NEL,NTEST)
*
* Obtain configuration number K .
* Occupation in KCNF
* Type in KTYP
*
* Jeppe Olsen , summer of 89
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "detdim.fh"
#include "spinfo_mclr.fh"
*. General input
      DIMENSION KCNF(*)
*. Output
      INTEGER ICONF(*)
*
*
      ICNFB1 = 1
      ICNFB2 = 1
      KTYP   = 0
      DO 100 JTYP = 1, NTYP
        JOP = JTYP - 1 + MINOP
        JCL = (NEL-JOP)/2
        JOCC = JOP + JCL
*
        NJCNF = NCNATS(JTYP,IREFSM)
        IF(K.GE.ICNFB1.AND.K.LE.ICNFB1+NJCNF-1) THEN
          KREL = K - ICNFB1 + 1
          KADD = (KREL-1)*JOCC
          KTYP = JTYP
          CALL ICOPVE(ICONF(ICNFB2+KADD),KCNF,JOCC)
        END IF
*
        ICNFB1 = ICNFB1 + NJCNF
        ICNFB2 = ICNFB2 + NJCNF*JOCC
  100 CONTINUE
*
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NTEST)
      END
