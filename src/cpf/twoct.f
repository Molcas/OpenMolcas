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
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      IF(ISDCI.EQ.0.AND.ICPF.EQ.0.AND.INCPF.EQ.0)GO TO 30
C CPF, ACPF AND SDCI
      IF(ITER.EQ.1)GO TO 25
      CALL DIAGC(H(LW(2)),H(LW(26)),H(LW(27)))
      IF(IFIRST.NE.0)GO TO 15
      CALL ABCI(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),
     *H(LW(50)),H(LW(51)),H(LW(52)),H(LW(53)),H(LW(54)))
15    CALL IJKL(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),
     *H(LW(46)),H(LW(47)),H(LW(47)),H(LW(31)),H(LW(32)))
      IF(IFIRST.NE.0)GO TO 25
      CALL ABCD(H(LW(2)),H(LW(3)),H(LW(4)),H(LW(26)),H(LW(27)),
     *H(LW(57)),H(LW(58)),H(LW(59)))
25    CALL FAIBJ(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),H(LW(36)),
     *H(LW(37)),H(LW(38)),H(LW(39)),H(LW(39)),H(LW(40)),H(LW(41)),
     *H(LW(42)),H(LW(43)),H(LW(31)),H(LW(32)))
      GO TO 50
C     MCPF
30    IF(ITER.EQ.1)GO TO 45
      CALL MDIAGC(H(LW(2)),H(LW(26)),H(LW(27)),H(LW(28)),H(LW(29)),
     *H(LW(31)),IRC(ILIM))
      IF(IFIRST.NE.0)GO TO 35
      CALL MABCI(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),
     *H(LW(50)),H(LW(51)),H(LW(52)),H(LW(53)),H(LW(54)),
     *H(LW(28)),H(LW(29)),H(LW(31)),IRC(ILIM))
35    CALL MIJKL(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),
     *H(LW(46)),H(LW(47)),H(LW(47)),H(LW(28)),H(LW(29)),H(LW(31)),
     *H(LW(32)),IRC(ILIM))
      IF(IFIRST.NE.0)GO TO 45
      CALL MABCD(H(LW(2)),H(LW(3)),H(LW(4)),H(LW(26)),H(LW(27)),
     *H(LW(57)),H(LW(58)),H(LW(59)),H(LW(28)),H(LW(29)),
     *H(LW(31)),IRC(ILIM))
45    CALL MFAIBJ(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),H(LW(36)),
     *H(LW(37)),H(LW(38)),H(LW(39)),H(LW(39)),H(LW(40)),H(LW(41)),
     *H(LW(42)),H(LW(43)),H(LW(28)),H(LW(29)),H(LW(31)),H(LW(32)),
     *IRC(ILIM))
50      Continue
       CALL QEXIT('TWOCT')
      RETURN
      END
