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
      SUBROUTINE EPSBIS(JSY,INDEX,C,W,EPB)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),W(*),EPB(*)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CRL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
CPAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(JSY,L)
C
      CALL SETZ(EPB,IRC(4))
      IF(ICPF.EQ.1.OR.ISDCI.EQ.1.OR.INCPF.EQ.1)RETURN
C
C     VALENCE
C
      IP=IRC(1)
      DO 15 I=1,IP
        EPB(I)=C(I)*W(I)
15    CONTINUE
C
C     SINGLES
C
      IP=IRC(2)-IRC(1)
      IN=IRC(1)
      DO 5 I=1,IP
CFUE    IND=IND+1
        NS1=JSYM(IN+I)
        NSIL=MUL(NS1,LSYM)
        INUM=NVIR(NSIL)
        IST=INDEX(IN+I)+1
CRL     CALL DOTPR(C(IST),1,W(IST),1,EPB(IN+I),INUM)
        EPB(IN+I)=DDOT_(INUM,C(IST),1,W(IST),1)
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
CRL     CALL DOTPR(C(IST),1,W(IST),1,EPB(IN+I),INUM)
        EPB(IN+I)=DDOT_(INUM,C(IST),1,W(IST),1)
10    CONTINUE
C
      IP=IRC(4)
      IF(IPRINT.GT.5)WRITE(6,998)(EPB(I),I=1,IP)
998   FORMAT(6X,'EPB ',5F10.6)
C
      RETURN
      END
