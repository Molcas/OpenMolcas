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
      SUBROUTINE EPSPRIM(JSY,INDEX,C,S,EPP)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),S(*),EPP(*)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CRL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
CPAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(JSY,L)
C
C     VALENCE
C
      IP=IRC(1)
      DO 15 I=1,IP
      EPP(I)=EPP(I)+C(I)*S(I)
15    CONTINUE
C
C
C     SINGLES
C
      IP=IRC(2)-IRC(1)
      IN=IRC(1)
      DO 5 I=1,IP
      NS1=JSYM(IN+I)
      NSIL=MUL(NS1,LSYM)
      INUM=NVIR(NSIL)
      IST=INDEX(IN+I)+1
CRL   CALL DOTPR(C(IST),1,S(IST),1,T,INUM)
      T=DDOT_(INUM,C(IST),1,S(IST),1)
      EPP(IN+I)=EPP(IN+I)+T
5     CONTINUE
C
C     DOUBLES
C
      IP=IRC(4)-IRC(2)
      IN=IRC(2)
      DO 10 I=1,IP
      NS1=JSYM(IN+I)
      NSIL=MUL(NS1,LSYM)
      INUM=NNS(NSIL)
      IST=INDEX(IN+I)+1
CRL   CALL DOTPR(C(IST),1,S(IST),1,T,INUM)
      T=DDOT_(INUM,C(IST),1,S(IST),1)
      EPP(IN+I)=EPP(IN+I)+T
10    CONTINUE
C
      IP=IRC(4)
      IF(IPRINT.GT.5)WRITE(6,998)(EPP(I),I=1,IP)
998   FORMAT(6X,'EPP ',5F10.6)
      RETURN
      END
