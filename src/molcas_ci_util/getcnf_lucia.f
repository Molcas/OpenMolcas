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
* Copyright (C) 2002, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GETCNF_LUCIA(KCNF,KTYP,K,ICONF,IREFSM,NEL)
*
* Obtain configuration number K .
* Occupation in KCNF in form of old RASSCF ( doubly occ orbs first)
* Type in KTYP
*
* Adapted for LUCIA Jeppe Olsen, summer of 02
*
      IMPLICIT REAL*8 (A-H,O-Z)

#include "spinfo.fh"
#include "ciinfo.fh"
*
      DIMENSION KCNF(*),ICONF(*)
*. Configuration list is assumed to be in the form used
*. in LUCIA, i.e. doubly occupied orbitals are flagged by
*. a minus
*
      ICNFB1 = 1
      ICNFB2 = 1
      KTYP   = 0
      DO 100 JTYP = 1, NTYP
        JOP = JTYP - 1 + MINOP
        JCL = (NEL-JOP)/2
        JOCC = JOP + JCL
*
        NJCNF = NCNFTP(JTYP,IREFSM)
        IF(K.GE.ICNFB1.AND.K.LE.ICNFB1+NJCNF-1) THEN
          KREL = K - ICNFB1 + 1
          KADD = (KREL-1)*JOCC
          KTYP = JTYP
*. Outdated ....
C         CALL ICOPY(JOCC,ICONF(ICNFB2+KADD),1,KCNF,1)
*
*. Obtain configuration in standard RASSCF form
          IIBOP = 1
          IIBCL = 1
          DO KOCC = 1, JOCC
            KORB = ICONF(ICNFB2+KADD-1+KOCC)
            IF(KORB.LT.0) THEN
*. Doubly occupied orbitals
              KCNF(IIBCL) = ABS(KORB)
              IIBCL = IIBCL + 1
            ELSE
*. Singly occupied orbital
              KCNF(JCL+IIBOP) = KORB
              IIBOP = IIBOP + 1
            END IF
          END DO
        END IF
*
        ICNFB1 = ICNFB1 + NJCNF
        ICNFB2 = ICNFB2 + NJCNF*JOCC
100   CONTINUE
*
      NTEST = 0
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Output from GETCNF '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Input configuration number : ', K
        WRITE(6,*) ' Corresponding type : ', KTYP
        WRITE(6,*) ' Occupation : '
        NOCC = (NEL+KTYP-1+MINOP)/2
        CALL IWRTMA(KCNF,1,NOCC,1,NOCC)
      END IF
*
      RETURN
      END
