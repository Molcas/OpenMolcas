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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE SDCI_CPF(H,iH,LIC0)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
c      DIMENSION H(LIC0), iH(RtoI*LIC0)
      DIMENSION H(*), iH(*)
      CALL QENTER('SDCI')
      CALL SDCI_CPF_INTERNAL(H)
      CALL QEXIT('SDCI')
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE SDCI_CPF_INTERNAL(H)
      USE ISO_C_BINDING
      REAL*8, TARGET :: H(*)
      INTEGER, POINTER :: iH1(:),iH2(:),iH3(:)
      LIC=LIC0
      IPRINT=5
      IDENS=0
C INPUT, SORTING AND DIAGONAL ELEMENTS
      CALL READIN_CPF(H,iH)
      CALL DIAGCT_CPF(H)
      ITER=1
      IF(IREST.EQ.1)ITER=ITER+1
      ITPUL=1
      IF(IREST.EQ.0)CALL START_CPF(H(LW(26)),JSC(4),IREF0)
      IF(IREST.EQ.1)CALL RESTART_CPFMCPF(H(LW(26)),JSC(4))
      IF(ICPF.EQ.0.AND.INCPF.EQ.0.AND.ISDCI.EQ.0) THEN
        CALL C_F_POINTER(C_LOC(H(LW(1))),iH1,[1])
        CALL THETSET(iH1,H(LW(29)),IRC(4))
        NULLIFY(iH1)
      END IF
100   CALL C_F_POINTER(C_LOC(H(LW(1))),iH1,[1])
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL NPSET(iH2,iH3,H(LW(26)),H(LW(30)),H(LW(31)),
     *H(LW(72)),H(LW(27)),H(LW(28)),H(LW(32)),iH1)
      NULLIFY(iH1,iH2,iH3)
      CALL TWOCT(H)
      CALL ONECT(H)
      CALL CPFCTL(H)
      ITER=ITER+1
      ITPUL=ITPUL+1
      IF(ITER.GT.MAXIT.OR.ICONV.EQ.1)GO TO 395
      GO TO 100
C     FIRST ORDER DENSITY MATRIX
395   IDENS=1
*
      CALL DENSCT_CPF(H,LIC0)
*
      IF(NREF.GT.1) THEN
         WRITE(6,*) ' This is a single reference program, but more than'
         WRITE(6,*) ' one reference state has been specified in the'
         WRITE(6,*) ' GUGA program. Change input to GUGA and run again.'
      CALL XFLUSH(6)
         CALL QEXIT('SDCI')
         RETURN
      END IF
      RETURN
      END SUBROUTINE SDCI_CPF_INTERNAL
*
      END
