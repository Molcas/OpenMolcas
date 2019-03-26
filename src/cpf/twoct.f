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
      SUBROUTINE TWOCT(H)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION H(*)
      CALL QENTER('TWOCT')
      CALL TWOCT_INTERNAL(H)
      CALL QEXIT('TWOCT')
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE TWOCT_INTERNAL(H)
      USE ISO_C_BINDING
      REAL*8, TARGET :: H(*)
      INTEGER, POINTER :: iH2(:),iH3(:),iH4(:),iH39(:),iH47(:),iH51(:)
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      IF(ISDCI.EQ.0.AND.ICPF.EQ.0.AND.INCPF.EQ.0)GO TO 30
C CPF, ACPF AND SDCI
      IF(ITER.EQ.1)GO TO 25
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL DIAGC(iH2,H(LW(26)),H(LW(27)))
      NULLIFY(iH2)
      IF(IFIRST.NE.0)GO TO 15
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(51))),iH51,[1])
      CALL ABCI(iH2,iH3,H(LW(26)),H(LW(27)),
     *H(LW(50)),iH51,H(LW(52)),H(LW(53)),H(LW(54)))
      NULLIFY(iH2,iH3,iH51)
15    CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(47))),iH47,[1])
      CALL IJKL(iH2,iH3,H(LW(26)),H(LW(27)),
     *H(LW(46)),H(LW(47)),iH47,H(LW(31)),H(LW(32)))
      NULLIFY(iH2,iH3,iH47)
      IF(IFIRST.NE.0)GO TO 25
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(4))),iH4,[1])
      CALL ABCD(iH2,iH3,iH4,H(LW(26)),H(LW(27)),
     *H(LW(57)),H(LW(58)),H(LW(59)))
      NULLIFY(iH2,iH3,iH4)
25    CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(39))),iH39,[1])
      CALL FAIBJ(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(36)),
     *H(LW(37)),H(LW(38)),H(LW(39)),iH39,H(LW(40)),H(LW(41)),
     *H(LW(42)),H(LW(43)),H(LW(31)),H(LW(32)))
      NULLIFY(iH2,iH3,iH39)
      GO TO 50
C     MCPF
30    IF(ITER.EQ.1)GO TO 45
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL MDIAGC(iH2,H(LW(26)),H(LW(27)),H(LW(28)),H(LW(29)),
     *H(LW(31)),IRC(ILIM))
      NULLIFY(iH2)
      IF(IFIRST.NE.0)GO TO 35
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(51))),iH51,[1])
      CALL MABCI(iH2,iH3,H(LW(26)),H(LW(27)),
     *H(LW(50)),iH51,H(LW(52)),H(LW(53)),H(LW(54)),
     *H(LW(28)),H(LW(29)),H(LW(31)),IRC(ILIM))
      NULLIFY(iH2,iH3,iH51)
35    CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(47))),iH47,[1])
      CALL MIJKL(iH2,iH3,H(LW(26)),H(LW(27)),
     *H(LW(46)),H(LW(47)),iH47,H(LW(28)),H(LW(29)),H(LW(31)),
     *H(LW(32)),IRC(ILIM))
      NULLIFY(iH2,iH3,iH47)
      IF(IFIRST.NE.0)GO TO 45
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(4))),iH4,[1])
      CALL MABCD(iH2,iH3,iH4,H(LW(26)),H(LW(27)),
     *H(LW(57)),H(LW(58)),H(LW(59)),H(LW(28)),H(LW(29)),
     *H(LW(31)),IRC(ILIM))
      NULLIFY(iH2,iH3,iH4)
45    CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(39))),iH39,[1])
      CALL MFAIBJ(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(36)),
     *H(LW(37)),H(LW(38)),H(LW(39)),iH39,H(LW(40)),H(LW(41)),
     *H(LW(42)),H(LW(43)),H(LW(28)),H(LW(29)),H(LW(31)),H(LW(32)),
     *IRC(ILIM))
      NULLIFY(iH2,iH3,iH39)
50      Continue
      RETURN
      END SUBROUTINE TWOCT_INTERNAL
*
      END
