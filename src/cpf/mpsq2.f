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
      SUBROUTINE MPSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC3,LSYM,
     *NVIRT,SQ2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),S(*),W(*),MUL(8,8),INDEX(*),JSY(*),NDIAG(*)
CPAM97      INTEGER UNPACK
CPAM97      EXTERNAL UNPACK
CRL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
CPAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP_CPF(JSY,L)
      DO 10 I=1,INUM
        II1=IRC3+I
        NS1=JSYM(II1)
        NS1L=MUL(NS1,LSYM)
        IF(NS1L.NE.1)GO TO 10
        NA=INDEX(II1)
        DO 20 MA=1,NVIRT
          C(NA+NDIAG(MA))=SQ2*C(NA+NDIAG(MA))
          S(NA+NDIAG(MA))=S(NA+NDIAG(MA))/SQ2
          W(NA+NDIAG(MA))=W(NA+NDIAG(MA))/SQ2
20      CONTINUE
10    CONTINUE
      RETURN
      END
