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
      SUBROUTINE MABCI(JSY,INDEX,C,S,BMN,IBMN,BIAC,BICA,BUFIN,
     &                 W,THET,ENP,NII)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),S(*),BMN(*),IBMN(*),
     &          BIAC(*),BICA(*),BUFIN(*),W(*),THET(NII,NII),ENP(*)
      PARAMETER (IPOW6=2**6,IPOW13=2**13,IPOW19=2**19)
*
      JSYM(L)=JSUNP(JSY,L)
*
      INUM=IRC(4)-IRC(3)
      CALL MPSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
      ICHK=0
      INSIN=KBUFF1
      IAD15=IADABCI
      IADD10=IAD10(4)
      CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IN=2
      NSAVE=ICOP1(IN)
100   NI=NSAVE
      IOUT=0
110   IN=IN+1
      IF(IN.LE.LEN)GO TO 15
      CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.LE.0)GO TO 5
      IN=1
15    IF(ICHK.NE.0)GO TO 460
      IF(ICOP1(IN).EQ.0)GO TO 10
      IOUT=IOUT+1
      BMN(IOUT)=COP(IN)
      IBMN(IOUT)=ICOP1(IN)
      GO TO 110
10    ICHK=1
      GO TO 110
460   ICHK=0
      NSAVE=ICOP1(IN)
5     CONTINUE
      DO 20 NB=1,NVIRT
        NSIB=MUL(NSM(LN+NB),NSM(NI))
        NSLB=MUL(NSM(LN+NB),LSYM)
        LB=NB-NSYS(NSM(LN+NB))
        INS=NNS(NSIB)
        ILOOP=0
72      CONTINUE
        DO 75 I=1,INS
          IF ( INSIN.GE.KBUFF1 ) THEN
             CALL dDAFILE(Lu_TiABCI,2,BUFIN,KBUFF1,IAD15)
             INSIN=0
          END IF
          INSIN=INSIN+1
          IF(ILOOP.EQ.0)BIAC(I)=BUFIN(INSIN)
          IF(ILOOP.EQ.1)BICA(I)=BUFIN(INSIN)
75      CONTINUE
        ILOOP=ILOOP+1
        IF(ILOOP.EQ.1)GO TO 72
        DO 25 IT=1,IOUT
          IND=IBMN(IT)
CPAM97          ICP1=IAND(ISHFT(IND,-19),8191)
*          ICP1=MOD(IND/IPOW19,IPOW13)
          ICP1=IBITS(IND,19,13)
          INDA=IRC(1)+ICP1
          IF(JSYM(INDA).NE.NSLB)GO TO 25
          MA=INDEX(INDA)+LB
CPAM97          ICP2=IAND(ISHFT(IND,-6),8191)
*          ICP2=MOD(IND/IPOW6,IPOW13)
CPAM97          ITYP=IAND(IND,63)
*          ITYP=MOD(IND,IPOW6)
          ICP2=IBITS(IND,6,13)
          ITYP=IBITS(IND, 0,6)
          IF(INS.EQ.0)GO TO 25
          COPL=BMN(IT)*C(MA)
          INDB=IRC(ITYP)+ICP2
          D1=1.0D0
          D2=2.0D0
          XXX=THET(INDA,INDB)/2.0D0
          ENPQ=(D1-XXX)*(ENP(INDA)+ENP(INDB)-D1)+XXX
          FACS=SQRT(ENP(INDA))*SQRT(ENP(INDB))/ENPQ
          FACW=FACS*(D2-THET(INDA,INDB))/ENPQ
          FACWA=FACW*ENP(INDA)-FACS
          FACWB=FACW*ENP(INDB)-FACS
          ICCB=INDEX(INDB)+1
          IF ( ITYP.EQ.3 ) THEN
            TERM=DDOT_(INS,C(ICCB),1,BIAC,1)
            CALL DAXPY_(INS,COPL*FACS,BIAC,1,S(ICCB),1)
            CALL DAXPY_(INS,COPL*FACWB,BIAC,1,W(ICCB),1)
          ELSE
            TERM=DDOT_(INS,C(ICCB),1,BICA,1)
            CALL DAXPY_(INS,COPL*FACS,BICA,1,S(ICCB),1)
            CALL DAXPY_(INS,COPL*FACWB,BICA,1,W(ICCB),1)
          END IF
          S(MA)=S(MA)+BMN(IT)*FACS*TERM
          W(MA)=W(MA)+BMN(IT)*FACWA*TERM
25      CONTINUE
20    CONTINUE
      IF(LEN.GE.0)GO TO 100
      CALL MDSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),
     *LSYM,NVIRT,SQ2)
      RETURN
      END
