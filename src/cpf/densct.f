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
      SUBROUTINE DENSCT(H,LIC0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(LIC0)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      CALL QENTER('DENSCT')
C
      CALL DENS(H(LW(26)),H(LW(62)),H(LW(1)),A)
C
C     MULTIPLY C BY MP
C
      CALL NPSET(H(LW(2)),H(LW(3)),H(LW(26)),H(LW(30)),H(LW(31)),
     *H(LW(72)),H(LW(27)),H(LW(28)),H(LW(32)),H(LW(34)))
*
      CALL ONECT(H)
      IF(A.GT.D1) THEN
        WRITE(6,*)'DENSCT Error: A>1.0D0 (See code.)'
      END IF
      CALL NATCT(H,LIC0)
*
      CALL QEXIT('DENSCT')
      RETURN
      END
