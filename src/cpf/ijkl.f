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
cpgi$g opt=1
      SUBROUTINE IJKL_CPF(JSY,INDEX,C,S,FIJKL,BUFIN,IBUFIN,
     *ENP,EPP)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),S(*),FIJKL(*),BUFIN(*),
     &           IBUFIN(*),ENP(*),EPP(*)
      PARAMETER(IPOW8=2**8,IPOW16=2**16,IPOW24=2**24)
      PARAMETER(IPOW6=2**6,IPOW13=2**13,IPOW19=2**19)
*
      JSYM(L)=JSUNP_CPF(JSY,L)
*
      FINI=0.0D0 ! dummy initialize
      NCONF=JSC(4)
      ICHK=0
      NIJ=IROW(LN+1)
      NIJKL=NIJ*(NIJ+1)/2
      DO 5 I=1,NIJKL
      FIJKL(I)=D0
5     CONTINUE
      KKBUF0=(RTOI*(KBUFF1+2)-2)/(RTOI+1)
      KKBUF1=RTOI*KKBUF0+KKBUF0+1
      KKBUF2=KKBUF1+1
      IADR=LASTAD(1)
201   CALL iDAFILE(Lu_TiABCI,2,IBUFIN,KKBUF2,IADR)
      LENGTH=IBUFIN(KKBUF1)
      IADR=IBUFIN(KKBUF2)
      IF(LENGTH.EQ.0)GO TO 209
      CALL SCATTER(LENGTH,FIJKL,IBUFIN(RTOI*KKBUF0+1),BUFIN)
209   IF(IADR.NE.-1) GO TO 201
      IADD10=IAD10(5)
100   CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0)GO TO 200
      DO 10 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0)GO TO 460
      IF(IND.NE.0)GO TO 22
      ICHK=1
      GO TO 10
460   ICHK=0
      INDI=IND
*      IP=MOD(INDI,IPOW8)
*      JP=MOD(INDI/IPOW8,IPOW8)
*      KP=MOD(INDI/IPOW16,IPOW8)
*      LP=MOD(INDI/IPOW24,IPOW8)
      IP=IBITS(INDI,0,8)
      JP=IBITS(INDI,8,8)
      KP=IBITS(INDI,16,8)
      LP=IBITS(INDI,24,8)
      NIJ=IROW(IP)+JP
      NKL=IROW(KP)+LP
      IND=NIJ*(NIJ-1)/2+NKL
      FINI=FIJKL(IND)
      GO TO 10
22    IF(ABS(FINI).LT.1.d-06)GO TO 10
CPAM97      IVL=IAND(IND,63)
CPAM97      IC2=IAND(ISHFT(IND,-6),8191)
CPAM97      IC1=IAND(ISHFT(IND,-19),8191)
*      IVL=MOD(IND,IPOW6)
*      IC2=MOD(IND/IPOW6,IPOW13)
*      IC1=MOD(IND/IPOW19,IPOW13)
      IVL=IBITS(IND, 0,6)
      IC2=IBITS(IND,6,13)
      IC1=IBITS(IND,19,13)
      COPI=COP(IN)*FINI
      IF(IVL.NE.0)GO TO 13
      IF(IC1.NE.IREF0)GO TO 16
      COPI=COPI/SQRT(ENP(IC2))
      S(IC2)=S(IC2)+COPI
      IF(ITER.EQ.1)GO TO 10
      EPP(IC2)=EPP(IC2)+COPI*C(IC2)
      GO TO 10
16    IF(IC2.NE.IREF0)GO TO 17
      COPI=COPI/SQRT(ENP(IC1))
      S(IC1)=S(IC1)+COPI
      IF(ITER.EQ.1)GO TO 10
      EPP(IC1)=EPP(IC1)+COPI*C(IC1)
      GO TO 10
17    FACS=D1
      S(IC1)=S(IC1)+FACS*COPI*C(IC2)
      S(IC2)=S(IC2)+FACS*COPI*C(IC1)
      GO TO 10
13    INDA=IRC(IVL)+IC1
      INDB=IRC(IVL)+IC2
      FACS=D1
      NA=INDEX(INDA)
      NB=INDEX(INDB)
      NS1=JSYM(INDA)
      NS1L=MUL(NS1,LSYM)
      INUM=NVIR(NS1L)
      IF(IVL.GE.2)INUM=NNS(NS1L)
      CALL DAXPY_(INUM,COPI*FACS,C(NB+1),1,S(NA+1),1)
      CALL DAXPY_(INUM,COPI*FACS,C(NA+1),1,S(NB+1),1)
10    CONTINUE
      GO TO 100
200   RETURN
      END
