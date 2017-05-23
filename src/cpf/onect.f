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
      SUBROUTINE ONECT(H)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION H(*)
      CALL QENTER('ONECT')
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      IF(ICPF.EQ.0.AND.ISDCI.EQ.0.AND.INCPF.EQ.0)GO TO 15
C     CPF AND SDCI
      IF(IDENS.EQ.1)GO TO 10
C     (AI/JK) INTEGRALS
      CALL AI_CPF(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),H(LW(62)),
     *H(LW(63)),H(LW(63)),H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),
     *H(LW(31)),H(LW(32)),1)
10    CALL FIJ(H(LW(1)),H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),
     *H(LW(62)),H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),
     *H(LW(31)),H(LW(32)))
      GO TO 20
C     MCPF
15    IF(IDENS.EQ.1)GO TO 5
C     (AI/JK) INTEGRALS
      CALL MAI(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),H(LW(62)),
     *H(LW(63)),H(LW(63)),H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),
     *H(LW(28)),H(LW(29)),H(LW(31)),H(LW(32)),IRC(ILIM),1)
5     CALL MFIJ(H(LW(1)),H(LW(2)),H(LW(3)),H(LW(26)),H(LW(27)),
     *H(LW(62)),H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),H(LW(28)),
     *H(LW(29)),H(LW(31)),H(LW(32)),IRC(ILIM))
20      Continue
       CALL QEXIT('ONECT')
      RETURN
      END
