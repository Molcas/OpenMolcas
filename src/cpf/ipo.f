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
      SUBROUTINE IPO(IPOA,NVIR,MUL,NSYM,KLS,IFT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPOA(*),NVIR(*),MUL(8,8)
      NSUM=0
      DO 10 N=1,NSYM
      IPOA(N)=NSUM
      M=MUL(N,KLS)
      IF(IFT.GE.0)GO TO 20
      NSUM=NSUM+NVIR(N)*NVIR(M)
      GO TO 10
20    IF (N-M.LT.0) THEN
        GO TO 10
      ELSE IF (N-M.EQ.0) THEN
        GO TO 11
      ELSE
        GO TO 12
      END IF
11    NSUM=NSUM+NVIR(N)*(NVIR(N)+1)/2
      GO TO 10
12    NSUM=NSUM+NVIR(N)*NVIR(M)
10    CONTINUE
      IPOA(NSYM+1)=NSUM
      RETURN
      END
