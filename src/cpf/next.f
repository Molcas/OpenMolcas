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
      SUBROUTINE NEXT(P,DPS,CN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(*),DPS(*),CN(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
C
      IAD=IADDP(1)
      CALL dDAFILE(Lu_CI,2,P,NCONF,IAD)
      ITM=ITPUL-1
      DO 5 I=1,ITM
        IN=I+1
        CTOT=0.0D00
        DO 6 J=IN,ITPUL
          CTOT=CTOT+CN(J)
6       CONTINUE
        IAD=IADDP(I+1)
        CALL dDAFILE(Lu_CI,2,DPS,NCONF,IAD)
        CALL VSMA(DPS,1,CTOT,P,1,P,1,NCONF)
5     CONTINUE
      IF(IPRINT.GE.15)WRITE(6,19)(P(I),I=1,NCONF)
19    FORMAT(6X,'C(NEXT)',5F10.6)
C
      IADC(ITPUL+2)=IAD
      RETURN
      END
