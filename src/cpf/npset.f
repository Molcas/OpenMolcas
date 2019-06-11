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
      SUBROUTINE NPSET(JSY,INDEX,C,TPQ,ENP,T,S,W,EPP,ICASE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JSY(*),INDEX(*),C(*),TPQ(*),ENP(*),T(*),S(*)
      DIMENSION W(*),EPP(*),ICASE(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
*
      JSYM(L)=JSUNP_CPF(JSY,L)
C
      IF(IDENS.EQ.1)GO TO 65
      IF(ITPUL.NE.1)GO TO 60
      IAD=0
      IADDP(1)=0
      CALL dDAFILE(Lu_CI,1,C,NCONF,IAD)
      IADDP(2)=IAD
C
C     VALENCE
C
60    IQ=IRC(1)
      DO 30 I=1,IQ
        T(I)=C(I)*C(I)
30    CONTINUE
C
C     SINGLES
C
      IQ=IRC(2)-IRC(1)
      IND=IRC(1)
      IN=IRC(1)
      DO 5 I=1,IQ
        IND=IND+1
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NVIR(NSIL)
        IST=INDEX(IN+I)+1
        T(IND)=DDOT_(INUM,C(IST),1,C(IST),1)
5     CONTINUE
C
C     DOUBLES
C
      IQ=IRC(4)-IRC(2)
      IN=IRC(2)
      DO 10 I=1,IQ
        IND=IND+1
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NNS(NSIL)
        IST=INDEX(IN+I)+1
        T(IND)=DDOT_(INUM,C(IST),1,C(IST),1)
10    CONTINUE
      IP=IRC(4)
      DO 15 I=1,IP
        CALL TPQSET(ICASE,TPQ,I)
        ENP(I)=DDOT_(IP,TPQ,1,T,1)
        ENP(I)=ENP(I)+D1
15    CONTINUE
      IP=IRC(4)
      IF(IPRINT.GT.5)WRITE(6,12)(ENP(I),I=1,IP)
12    FORMAT(6X,'ENP ',5F14.8)
C
C     VALENCE
C
65    IQ=IRC(1)
      DO 6 I=1,IQ
        IF(IDENS.EQ.0)EMPI=D1/SQRT(ENP(I))
        IF(IDENS.EQ.1)EMPI=SQRT(ENP(I))
        C(I)=C(I)*EMPI
6     CONTINUE
C
C     SINGLES
C
      IQ=IRC(2)-IRC(1)
      IN=IRC(1)
      DO 16 I=1,IQ
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NVIR(NSIL)
        IST=INDEX(IN+I)+1
        IF(IDENS.EQ.0)EMPI=D1/SQRT(ENP(IN+I))
        IF(IDENS.EQ.1)EMPI=SQRT(ENP(IN+I))
        CALL VSMUL(C(IST),1,EMPI,C(IST),1,INUM)
16    CONTINUE
C
C     DOUBLES
C
      IQ=IRC(4)-IRC(2)
      IN=IRC(2)
      DO 11 I=1,IQ
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NNS(NSIL)
        IST=INDEX(IN+I)+1
        IF(IDENS.EQ.0)EMPI=D1/SQRT(ENP(IN+I))
        IF(IDENS.EQ.1)EMPI=SQRT(ENP(IN+I))
        CALL VSMUL(C(IST),1,EMPI,C(IST),1,INUM)
11    CONTINUE
      IF(IPRINT.GE.15)WRITE(6,13)(C(I),I=1,NCONF)
13    FORMAT(6X,'C(NP)',5F10.6)
      IF(IDENS.EQ.1)RETURN
C
      CALL SETZ(EPP,IRC(4))
      CALL SETZ(S,JSC(4))
      IF(ICPF.NE.1.AND.ISDCI.NE.1.AND.INCPF.NE.1)CALL SETZ(W,JSC(4))
      RETURN
      END
