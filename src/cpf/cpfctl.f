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
      SUBROUTINE CPFCTL(H)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(*)

#include "SysDef.fh"

#include "cpfmcpf.fh"
C
      CALL QENTER('CPFCTL')
      CALL CPFCTL_INTERNAL(H)
      CALL QEXIT('CPFCTL')
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE CPFCTL_INTERNAL(H)
      USE ISO_C_BINDING
      REAL*8, TARGET :: H(*)
      INTEGER, POINTER :: iH1(:),iH2(:),iH3(:)
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL EPSBIS(iH2,iH3,H(LW(26)),H(LW(28)),H(LW(75)))
      NULLIFY(iH2,iH3)
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL EPSPRIM(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(32)))
      NULLIFY(iH2,iH3)
      IP=IRC(4)
      CALL VECSUM_CPFMCPF(H(LW(32)),ENER,IP)
      ETOTT=ENER+POTNUC
      DELE=ETOTT-ETOT
      ETOT=ETOTT
      IF(ITER.EQ.1) THEN
        WRITE(6,'(1X,A)')' ITER      TOTAL ENERGY'
     &       //'          CORR ENERGY           DECREASE'
        CALL XFLUSH(6)
      END IF
      WRITE(6,'(1X,I3,3(5X,F16.8))')ITER,ETOT,ENER,DELE
      CALL XFLUSH(6)
      IF(ABS(DELE).LT.ETHRE.AND.ITPUL.NE.1)ICONV=1
      IF(ICONV.EQ.0.AND.ITER.NE.MAXIT)GO TO 20
C If more iterations should be done, goto 20.

      IF(ICONV.EQ.1)WRITE(6,37)
37    FORMAT(/,5X,'CALCULATION CONVERGED')
      IF(ICONV.EQ.0)WRITE(6,38)
38    FORMAT(/,5X,'CALCULATION NOT COMPLETELY CONVERGED')
      IF(ISDCI.EQ.1) WRITE(6,30)ETOT
30    FORMAT(/,5X,'FINAL CI ENERGY',6X,F17.8)
      IF(ICPF.EQ.1)WRITE(6,35)ETOT
35    FORMAT(/,5X,'FINAL CPF ENERGY',5X,F17.8)
      IF(INCPF.EQ.1)WRITE(6,39)ETOT
39    FORMAT(/,5X,'FINAL ACPF ENERGY',4X,F17.8)
      IF(ISDCI.EQ.0.AND.ICPF.EQ.0.AND.INCPF.EQ.0)WRITE(6,36)ETOT
36    FORMAT(/,5X,'FINAL MCPF ENERGY',5X,F17.8)
      WRITE(6,31)ENER,POTNUC
      CALL XFLUSH(6)
31    FORMAT(5X,'FINAL CORRELATION ENERGY',F14.8,'  REFERENCE ENERGY',
     *F17.8)
      If (ISDCI.EQ.1) Call Add_Info('E_SDCI',[ETOT],1,8)
      If (ICPF.EQ.1)  Call Add_Info('E_CPF',[ETOT],1,8)
      If (INCPF.EQ.1) Call Add_Info('E_ACPF',[ETOT],1,8)
      If (ISDCI.EQ.0.AND.ICPF.EQ.0.AND.INCPF.EQ.0)
     &                Call Add_Info('E_MCPF',[ETOT],1,8)
      CALL XFLUSH(6)
      IF(ISDCI.EQ.0)GO TO 21
      EENP=H(LW(31)+IRC(4)-1)
      C0=D1/SQRT(EENP)
      DECORR=ENER*(EENP-D1)
      DETOT=ETOT+DECORR
      WRITE(6,32)DETOT
32    FORMAT(5X,'DAVIDSON CORR. ENERGY',F17.8)
      WRITE(6,33)DECORR,C0
33    FORMAT(5X,'DAVIDSON CORRECTION',F19.8,'  C0 = ',F12.6)
      CALL XFLUSH(6)

21    CONTINUE
      IF(IPRINT.GT.5) THEN
        ISTA=LW(31)
        IEND=ISTA+IRC(4)-1
        IF(IPRINT.GT.5)WRITE(6,34)(H(I),I=ISTA,IEND)
34      FORMAT(/,(5X,'ENP',5F10.6))
      END IF

      RETURN

20    CONTINUE
C Here if ICONV.EQ.0 and  ITER.NE.MAXIT (More iterations to do).
      IDIIS=0
      IF(ITPUL.EQ.MAXITP)IDIIS=1
      CALL C_F_POINTER(C_LOC(H(LW(1))),iH1,[1])
      CALL  APPRIM(H(LW(32)),H(LW(75)),H(LW(30)),H(LW(76)),H(LW(31)),
     &             H(LW(79)),H(LW(80)),iH1)
      NULLIFY(iH1)
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL CUPDATE(iH2,iH3,H(LW(26)),H(LW(27)),
     &             H(LW(76)),H(LW(33)),H(LW(80)),H(LW(31)))
      NULLIFY(iH2,iH3)
      ITP=ITPUL+1
      CALL    DIIS_CPF(H(LW(26)),H(LW(27)),H(LW(33)),
     &             MAXIT,H(LW(77)),ITP,H(LW(78)))
      RETURN
      END SUBROUTINE CPFCTL_INTERNAL
*
      END
