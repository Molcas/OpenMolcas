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
      SUBROUTINE DIAGCT(H)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION H(*)
      CALL QENTER('DIAGCT')
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      NCONF=JSC(ILIM)
* Initialize sorting buffer, so that automatic detection of
* uninitialized variables does not give false alarms.
      CALL DCOPY_(LW(21)-LW(20),0.0D0,0,H(LW(20)),1)
* Now H(LW(20)) and up are filled with zeroes.
* Similar before SORTB and SORT_CPF.
      CALL SORTA(H(LW(20)),H(LW(20)),H(LW(21)),H(LW(22)),H(LW(10)),
     &   H(LW(4)),H(LW(25)),H(LW(25)),H(LW(23)),H(LW(24)),NINTGR)
      IF(IFIRST.EQ.0) THEN
        CALL DCOPY_(LW(18)-LW(17),0.0D0,0,H(LW(17)),1)
        CALL SORTB(H(LW(17)),H(LW(17)),H(LW(18)),
     &   H(LW(19)),H(LW(10)),H(LW(94)),H(LW(95)),H(LW(4)),H(LW(96)))
      END IF
      CALL DCOPY_(LW(12)-LW(11),0.0D0,0,H(LW(11)),1)
      CALL SORT_CPF(H(LW(11)),H(LW(11)),H(LW(12)),H(LW(13)),H(LW(14)),
     &   H(LW(15)),H(LW(16)),H(LW(10)))
      CALL DIAG_CPF(H(LW(1)),H(LW(2)),H(LW(11)),H(LW(14)),H(LW(15)),
     &   H(LW(16)))
      CALL QEXIT('DIAGCT')
      RETURN
      END
