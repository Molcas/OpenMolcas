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
      SUBROUTINE DIAGCT_CPF(H)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION H(*)
      CALL DIAGCT_INTERNAL(H)
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE DIAGCT_INTERNAL(H)
      USE ISO_C_BINDING
      REAL*8, TARGET :: H(*)
      INTEGER, POINTER :: iH1(:),iH2(:),iH4(:),iH11(:),iH12(:),iH13(:),
     &iH17(:),iH18(:),iH19(:),iH20(:),iH21(:),iH22(:),iH25(:)
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
      NCONF=JSC(ILIM)
* Initialize sorting buffer, so that automatic detection of
* uninitialized variables does not give false alarms.
      CALL DCOPY_(LW(21)-LW(20),[0.0D0],0,H(LW(20)),1)
* Now H(LW(20)) and up are filled with zeroes.
* Similar before SORTB_CPF and SORT_CPF.
      CALL C_F_POINTER(C_LOC(H(LW(4))),iH4,[1])
      CALL C_F_POINTER(C_LOC(H(LW(20))),iH20,[1])
      CALL C_F_POINTER(C_LOC(H(LW(21))),iH21,[1])
      CALL C_F_POINTER(C_LOC(H(LW(22))),iH22,[1])
      CALL C_F_POINTER(C_LOC(H(LW(25))),iH25,[1])
      CALL SORTA_CPF(H(LW(20)),iH20,iH21,iH22,H(LW(10)),
     &   iH4,H(LW(25)),iH25,H(LW(23)),H(LW(24)),NINTGR)
      NULLIFY(iH4,iH20,iH21,iH22,iH25)
      IF(IFIRST.EQ.0) THEN
        CALL DCOPY_(LW(18)-LW(17),[0.0D0],0,H(LW(17)),1)
        CALL C_F_POINTER(C_LOC(H(LW(4))),iH4,[1])
        CALL C_F_POINTER(C_LOC(H(LW(17))),iH17,[1])
        CALL C_F_POINTER(C_LOC(H(LW(18))),iH18,[1])
        CALL C_F_POINTER(C_LOC(H(LW(19))),iH19,[1])
        CALL SORTB_CPF(H(LW(17)),iH17,iH18,
     &   iH19,H(LW(10)),H(LW(94)),H(LW(95)),iH4,H(LW(96)))
        NULLIFY(iH4,iH17,iH18,iH19)
      END IF
      CALL DCOPY_(LW(12)-LW(11),[0.0D0],0,H(LW(11)),1)
      CALL C_F_POINTER(C_LOC(H(LW(11))),iH11,[1])
      CALL C_F_POINTER(C_LOC(H(LW(12))),iH12,[1])
      CALL C_F_POINTER(C_LOC(H(LW(13))),iH13,[1])
      CALL SORT_CPF(H(LW(11)),iH11,iH12,iH13,H(LW(14)),
     &   H(LW(15)),H(LW(16)),H(LW(10)))
      NULLIFY(iH11,iH12,iH13)
      CALL C_F_POINTER(C_LOC(H(LW(1))),iH1,[1])
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL DIAG_CPF(iH1,iH2,H(LW(11)),H(LW(14)),H(LW(15)),
     &   H(LW(16)))
      NULLIFY(iH1,iH2)
      RETURN
      END SUBROUTINE DIAGCT_INTERNAL
*
      END
