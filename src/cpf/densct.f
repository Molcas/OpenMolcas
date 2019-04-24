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
      SUBROUTINE DENSCT_CPF(H,LIC0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(LIC0)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      CALL QENTER('DENSCT_CPF')
      CALL DENSCT_INTERNAL(H)
      CALL QEXIT('DENSCT_CPF')
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE DENSCT_INTERNAL(H)
      USE ISO_C_BINDING
      REAL*8, TARGET :: H(*)
      INTEGER, POINTER :: iH1(:),iH2(:),iH3(:),iH34(:)
C
      CALL C_F_POINTER(C_LOC(H(LW(1))),iH1,[1])
      CALL DENS_CPF(H(LW(26)),H(LW(62)),iH1,A)
      NULLIFY(iH1)
C
C     MULTIPLY C BY MP
C
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(34))),iH34,[1])
      CALL NPSET(iH2,iH3,H(LW(26)),H(LW(30)),H(LW(31)),
     *H(LW(72)),H(LW(27)),H(LW(28)),H(LW(32)),iH34)
      NULLIFY(iH2,iH3,iH34)
*
      CALL ONECT(H)
      IF(A.GT.D1) THEN
        WRITE(6,*)'DENSCT_CPF Error: A>1.0D0 (See code.)'
      END IF
      CALL NATCT(H,LIC0)
*
      RETURN
      END SUBROUTINE DENSCT_INTERNAL
      END
