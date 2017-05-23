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
      SUBROUTINE MAB(ICASE,JSY,INDEX,C,S,FC,A,B,F,W,THET,ENP,NII)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION JSY(*),INDEX(*),C(*),S(*),FC(*),A(*),B(*),
     &          F(*),W(*),THET(NII,NII),ENP(*)
      DIMENSION ICASE(*)
      DIMENSION IPOA(9),IPOF(9)
      DIMENSION IOC(55)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CRL   JO(L)=IAND(ISHFT(QOCC((L+29)/30),-2*((L+29)/30*30-L)),3)
CPAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICASE,L)
CRL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
CPAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
      JSYM(L)=JSUNP(JSY,L)
      NAB=0 ! dummy initialize
      NOB2=IROW(NORBT+1)
      IF(IPRINT.GE.15) THEN
        WRITE(6,'(A,/,(10F12.6))')' S,AB',(S(I),I=1,JSC(4))
      CALL XFLUSH(6)
        WRITE(6,'(A,/,(10F12.6))')' W,AB',(W(I),I=1,JSC(4))
      CALL XFLUSH(6)
        IF(IDENS.EQ.1)WRITE(6,'(A,/,(10F12.6))')
     &                      ' FC,AB',(FC(I),I=1,NOB2)
      END IF
      INUM=IRC(4)-IRC(3)
      CALL MPSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),
     *LSYM,NVIRT,SQ2)
      NCLIM=4
      IF(IFIRST.NE.0)NCLIM=2
C     MOVE FOCK (DENSITY) MATRIX TO F IN SYMMETRY BLOCKS
      CALL IPO(IPOF,NVIR,MUL,NSYM,1,-1)
      ITURN=0
90    DO 10 IASYM=1,NSYM
      IAB=IPOF(IASYM)
      NA1=NSYS(IASYM)+1
      NA2=NSYS(IASYM+1)
      IF(NA2.LT.NA1)GO TO 10
      DO 15 NA=NA1,NA2
      DO 20 NB=NA1,NA2
      IAB=IAB+1
      IF(NA.GE.NB)NAB=IROW(LN+NA)+LN+NB
      IF(NB.GT.NA)NAB=IROW(LN+NB)+LN+NA
      IF(ITURN.EQ.1)GO TO 320
      IF(IDENS.EQ.0)F(IAB)=D0
      IF(IDENS.EQ.1)F(IAB)=FC(NAB)
      IF(NA.NE.NB)F(IAB)=FC(NAB)
      GO TO 20
320   IF(NA.LT.NB)FC(NAB)=F(IAB)
20    CONTINUE
15    CONTINUE
10    CONTINUE
      IF(ITURN.EQ.0)GO TO 11
      TR=D0
      IJ=0
      DO 510 I=1,NORBT
      IJ=IJ+I
      TR=TR+FC(IJ)
510   CONTINUE
      If (iPrint.ge.15) WRITE(6,310)TR
310   FORMAT(/,6X,'TRACE OF DENSITY MATRIX',F16.8)
      GO TO 300
11    II1=0
      ITAIL=IRC(NCLIM)
      DO 40 INDA=1,ITAIL
      IF(IDENS.EQ.0)GO TO 111
      DO 110 I=1,LN
      II1=II1+1
      JOJ=JO(II1)
      IF(JOJ.GT.1)JOJ=JOJ-1
      IOC(I)=JOJ
110   CONTINUE
111   IF(INDA.GT.IRC(1))GO TO 120
      IF(IDENS.EQ.0.OR.INDA.EQ.IREF0)GO TO 40
      ENPQ=(D1-THET(INDA,INDA)/D2)*(ENP(INDA)+ENP(INDA)-D1)+
     *THET(INDA,INDA)/D2
      TSUM=C(INDA)*C(INDA)/ENPQ
      GO TO 106
120   MYSYM=JSYM(INDA)
      MYL=MUL(MYSYM,LSYM)
      INMY=INDEX(INDA)+1
      ENPQ=(D1-THET(INDA,INDA)/D2)*(ENP(INDA)+ENP(INDA)-D1)+
     *THET(INDA,INDA)/D2
      FACS=SQRT(ENP(INDA))*SQRT(ENP(INDA))/ENPQ
      FACW=(FACS*(D2-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
      IF(INDA.GT.IRC(2))GO TO 25
C     DOUBLET-DOUBLET INTERACTIONS
      IF(NVIR(MYL).EQ.0)GO TO 40
      IF(IDENS.EQ.1)GO TO 65
      CALL SETZ(A,NVIR(MYL))
      CALL FMMM(F(IPOF(MYL)+1),C(INMY),A,NVIR(MYL),1,NVIR(MYL))
      CALL DAXPY_(NVIR(MYL),FACS,A,1,S(INMY),1)
      CALL DAXPY_(NVIR(MYL),FACW,A,1,W(INMY),1)
      GO TO 40
65    CALL FMUL2(C(INMY),C(INMY),A,NVIR(MYL),NVIR(MYL),1)
      IPF=IPOF(MYL)+1
      IN=IPOF(MYL+1)-IPOF(MYL)
      ENPQ=(D1-THET(INDA,INDA)/D2)*(ENP(INDA)+ENP(INDA)-D1)+
     *THET(INDA,INDA)/D2
      COPI=D1/ENPQ
      CALL VSMA(A,1,COPI,F(IPF),1,F(IPF),1,IN)
      NVIRA=NVIR(MYL)
      LNA=LN+NSYS(MYL)
      IIA=IROW(LNA+1)
      TSUM=D0
      DO 130 I=1,NVIRA
      SUM=COPI*C(INMY)*C(INMY)
      INMY=INMY+1
      TSUM=TSUM+SUM
      IIA=IIA+LNA+I
      FC(IIA)=FC(IIA)+SUM
130   CONTINUE
      GO TO 106
C     TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
25    IFT=1
      IF(INDA.GT.IRC(3))IFT=0
      CALL IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
      IN=0
      TSUM=D0
      DO 70 IASYM=1,NSYM
      IAB=IPOF(IASYM+1)-IPOF(IASYM)
      IF(IAB.EQ.0)GO TO 70
      ICSYM=MUL(MYL,IASYM)
      IF(NVIR(ICSYM).EQ.0)GO TO 70
      IF(IDENS.EQ.1)GO TO 75
      IF(MYL.NE.1)GO TO 30
      IF(IFT.EQ.0)CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
C      IF(IFT.EQ.1)CALL SQUARN(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      IF(IFT.EQ.1)CALL SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      NAA=NVIR(IASYM)*NVIR(IASYM)
      CALL SETZ(B,NAA)
      CALL FMMM(F(IPOF(IASYM)+1),A,B,NVIR(IASYM),NVIR(IASYM),
     *NVIR(IASYM))
      CALL SETZ(A,NAA)
      CALL DAXPY_(NAA,FACS,B,1,A,1)
      IF(IFT.EQ.1)GO TO 230
      CALL SIADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
      CALL SETZ(A,NAA)
      CALL DAXPY_(NAA,FACW,B,1,A,1)
      CALL SIADD(A,W(INMY+IPOA(IASYM)),NVIR(IASYM))
      GO TO 70
230   CALL TRADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
      CALL SETZ(A,NAA)
      CALL DAXPY_(NAA,FACW,B,1,A,1)
      CALL TRADD(A,W(INMY+IPOA(IASYM)),NVIR(IASYM))
      GO TO 70
30    NAC=NVIR(IASYM)*NVIR(ICSYM)
      CALL SETZ(A,NAC)
      IF(IASYM.GT.ICSYM)GO TO 31
      CALL FMMM(F(IPOF(IASYM)+1),C(INMY+IPOA(ICSYM)),A,
     *NVIR(IASYM),NVIR(ICSYM),NVIR(IASYM))
      CALL DAXPY_(NAC,FACS,A,1,S(INMY+IPOA(ICSYM)),1)
      CALL DAXPY_(NAC,FACW,A,1,W(INMY+IPOA(ICSYM)),1)
      GO TO 70
31    CALL FMMM(C(INMY+IPOA(IASYM)),F(IPOF(IASYM)+1),A,
     *NVIR(ICSYM),NVIR(IASYM),NVIR(IASYM))
      CALL DAXPY_(NAC,FACS,A,1,S(INMY+IPOA(IASYM)),1)
      CALL DAXPY_(NAC,FACW,A,1,W(INMY+IPOA(IASYM)),1)
      GO TO 70
75    IF(MYL.NE.1)GO TO 330
      IF(IFT.EQ.0)CALL SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      IF(IFT.EQ.1)CALL SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
      GO TO 255
330   IF(IASYM.GT.ICSYM)GO TO 231
      NAC=NVIR(IASYM)*NVIR(ICSYM)
      IF(IFT.EQ.0)CALL DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
      IF(IFT.EQ.1)CALL VNEG(C(INMY+IPOA(ICSYM)),1,A,1,NAC)
      GO TO 255
231   CALL MTRANS(C(INMY+IPOA(IASYM)),A,NVIR(IASYM),NVIR(ICSYM))
255   CALL FMUL2(A,A,B,NVIR(IASYM),NVIR(IASYM),NVIR(ICSYM))
      IPF=IPOF(IASYM)+1
      ENPQ=(D1-THET(INDA,INDA)/D2)*(ENP(INDA)+ENP(INDA)-D1)+
     *THET(INDA,INDA)/D2
      COPI=D1/ENPQ
      CALL VSMA(B,1,COPI,F(IPF),1,F(IPF),1,IAB)
      NVIRA=NVIR(IASYM)
      NVIRC=NVIR(ICSYM)
      INN=1
      LNC=LN+NSYS(ICSYM)
      IIC=IROW(LNC+1)
      DO 105 I=1,NVIRC
      SUM=DDOT_(NVIRA,A(INN),1,A(INN),1)
      SUM=COPI*SUM
      TSUM=TSUM+SUM
      IIC=IIC+LNC+I
      FC(IIC)=FC(IIC)+SUM
      INN=INN+NVIRA
105   CONTINUE
70    CONTINUE
      IF(IDENS.EQ.0)GO TO 40
      TSUM=TSUM/D2
106   IJ=0
      DO 107 I=1,LN
      IJ=IJ+I
      FC(IJ)=FC(IJ)+IOC(I)*TSUM
107   CONTINUE
40    CONTINUE
      ITURN=1
      IF(IDENS.EQ.1)GO TO 90
300   CALL MDSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),
     *LSYM,NVIRT,SQ2)
      IF(IPRINT.GE.15) THEN
        WRITE(6,'(A,/,(10F12.6))')' S,AB',(S(I),I=1,JSC(4))
      CALL XFLUSH(6)
        WRITE(6,'(A,/,(10F12.6))')' W,AB',(W(I),I=1,JSC(4))
      CALL XFLUSH(6)
        IF(IDENS.EQ.1)WRITE(6,'(A,/,(10F12.6))')
     &                      ' FC,AB',(FC(I),I=1,NOB2)
      END IF
      RETURN
      END
