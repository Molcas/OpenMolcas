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
************************************************************************
      SUBROUTINE AIBJ(L0,L1,L2,L3,ITAI)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
      DIMENSION L0(*),L1(*),L2(*),L3(*),ITAI(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
      DIMENSION NUMM(7)
*
      JO(L)=ICUNP(ICASE,L)
*
      IC1=0    ! dummy initialize
      IC2=0    ! dummy initialize
      COPLA=D0 ! dummy initialize
      DO 5 I=1,7
      NUMM(I)=0
5     CONTINUE
      IOUT=0
      NMAT=0
      DO 10 NI=1,LN
      DO 20 NJ=1,NI
      I=ICH(NI)
      J=ICH(NJ)
      IF(I.GT.J)GO TO 19
      I=ICH(NJ)
      J=ICH(NI)
19    LTYP=0
      IOUT=IOUT+1
      ICOP1(IOUT)=0
      IF(IOUT.LT.NBUF)GO TO 460
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
460   IOUT=IOUT+1
CPAM96      ICOP1(IOUT)=IOR(I,ISHFT(J,10))
*      ICOP1(IOUT)=I+2**10*J
      ICOP1(IOUT)=IOR(I,ISHFT(J,10))
      IF(IOUT.LT.NBUF)GO TO 11
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
11    IJS=IJ(I+1)+1
      IJM=IJ(I)
C     FIRST ORDER INTERACTION
C     TRIPLET-VALENCE INTERACTIONS
      JTURN=1
      ITURN=0
      ITT1=2
      ITT2=0
150   IT1=ITT1*MXVERT
culf      IT2=ITT2*300
      IT2=ITT2*MXVERT
      II=0
      IID=0
      IF(ITT2.EQ.0)GO TO 149
      II=IRC(ITT2)
      IID=JRC(ITT2)
149   JJ=IRC(ITT1)
      JJD=JRC(ITT1)
      ITYP=JTURN
      IF(ITURN.NE.0)ITYP=JTURN-2
      DO 30 IJJ=IJS,IJM
      ITAIL=IX(IT2+IJJ)
      IF(IT1.NE.IT2)CALL TAIL(I,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(I)=1
32    KM=I
      J2(KM+1)=IJJ
      J1(KM+1)=IJJ
      IF(I.EQ.J)GO TO 51
      CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 30
41    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.J)GO TO 51
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 41
      KM=KM+1
      IF(KM.EQ.I)GO TO 32
      GO TO 42
51    IWAY(J)=1
52    KM=J
      JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      IABIJ=0
      IF(I.EQ.J)GO TO 12
      IF(ITURN.EQ.0)CALL LOOP10(KM,ISTOP,IT1,IT2)
      IFAI=1
      IF(ITURN.EQ.1)CALL LOOP13(KM,ISTOP,IFAI,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 14
      IF(J1(KM).NE.J2(KM))GO TO 13
      IABIJ=1
      IC1=ICOUP(KM)
      IC2=ICOUP1(KM)
      KM1=KM+1
      IDIF=IA(J1(KM1))-IA(J2(KM1))
      IF(IWAY(J).EQ.2)COPLA=COUP(J+1)
      IF(IDIF.EQ.0.AND.IWAY(J).EQ.5)COPLA=COUP(J+1)*BS4(IB(J2(J+1))+1)
      IF(IDIF.EQ.1.AND.IWAY(J).EQ.4)COPLA=COUP(J+1)*BS3(IB(J2(J+1))+1)
      GO TO 53
12    IF(ITURN.EQ.0)CALL LOOP7(KM,ISTOP,IT1,IT2)
      IFAI=1
      IF(ITURN.EQ.1)CALL LOOP26(KM,ISTOP,IFAI,IT1,IT2)
      LTYP=1
      IF(IWAY(I).EQ.5)LTYP=0
13    IF(ISTOP.EQ.0)GO TO 53
14    IF(I.EQ.J)GO TO 30
      KM=KM+1
      IF(KM.EQ.I)GO TO 32
      GO TO 42
53    KM=KM-1
      IF(KM.EQ.0)GO TO 61
      IWAY(KM)=1
62    JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      IF(ITURN.EQ.0)CALL LOOP17(KM,ISTOP,IT1,IT2)
      IFAI=0
      IF(J1(KM+1).NE.J2(KM+1))GO TO 165
      IF(I.EQ.J)GO TO 162
      IF(IABIJ.EQ.0)GO TO 165
      IC11=ICOUP(KM+1)-IC1
      IC22=ICOUP1(KM+1)-IC2
      IF(IC11.EQ.IC22)IFAI=1
      GO TO 165
162   IF(ICOUP(KM+1).NE.ICOUP1(KM+1))GO TO 165
      IFAI=1
165   IF(ITURN.EQ.1)CALL LOOP23(KM,ISTOP,IFAI,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 53
      KM=KM+1
      IF(KM.EQ.J)GO TO 52
      GO TO 62
61    COPL=COUP(1)
      IF(JTURN.GE.6.OR.JTURN.EQ.3)COPL=COUP1(1)
      IFAB=0
      IF(JTURN.NE.5)GO TO 63
      IF(I.EQ.J.AND.LTYP.EQ.1)GO TO 72
63    IF(ITT1.NE.ITT2)GO TO 70
      IF(I.NE.J)GO TO 65
      IF(LTYP.EQ.1.AND.ICOUP(1).GT.ICOUP1(1))GO TO 72
65    IF(IABIJ.EQ.0)GO TO 70
      IC11=ICOUP(1)-IC1
      IC22=ICOUP1(1)-IC2
      IF(IC11.EQ.IC22)IFAB=1
70    COPLA0=COPLA
      IF(IFAB.EQ.0)COPLA0=D0
      DO 80 IN=1,ITAIL
      ICP1=ICOUP(1)+IN
      JND1=JNDX(II+ICP1)
      IF(JND1.EQ.0)GO TO 80
      ICP1=JND1-IID
      IF(ITT1.NE.ITT2)GO TO 289
      IN2=IN
      GO TO 288
289   IN2=ITAI(IN)
      IF(IN2.EQ.0)GO TO 80
288   ICP2=ICOUP1(1)+IN2
      JND2=JNDX(JJ+ICP2)
      IF(JND2.EQ.0)GO TO 80
      ICP2=JND2-JJD
      IF(ITT1.NE.ITT2.OR.ICP1.NE.ICP2)GO TO 100
      JJ1=(JJD+ICP1-1)*LN+I
      JOJ=JO(JJ1)
      IF(JOJ.GT.1)JOJ=JOJ-1
      COPLA0=JOJ-2
      IFAB=1
100   IOUT=IOUT+1
      NUMM(JTURN)=NUMM(JTURN)+1
      COP(IOUT)=COPL
CPAM96      IND1=IOR(IFAB,ISHFT(ITURN,1))
CPAM96      IND2=IOR(IND1,ISHFT(ITYP,2))
CPAM96      IND3=IOR(IND2,ISHFT(ICP1,5))
CPAM96      ICOP1(IOUT)=IOR(IND3,ISHFT(ICP2,18))
*      IND1=IFAB+2**1*ITURN
*      IND2=IND1+2**2*ITYP
*      IND3=IND2+2**5*ICP1
*      ICOP1(IOUT)=IND3+2**18*ICP2
      IND1=IOR(IFAB,ISHFT(ITURN,1))
      IND2=IOR(IND1,ISHFT(ITYP,2))
      IND3=IOR(IND2,ISHFT(ICP1,5))
      ICOP1(IOUT)=IOR(IND3,ISHFT(ICP2,18))
      IF(IOUT.LT.NBUF)GO TO 71
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
71    IF(IFAB.EQ.0)GO TO 80
      IOUT=IOUT+1
      NUMM(JTURN)=NUMM(JTURN)+1
      COP(IOUT)=COPLA0
      ICOP1(IOUT)=1
      IF(IOUT.LT.NBUF)GO TO 80
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
80    CONTINUE
72    KM=KM+1
      IF(KM.EQ.J)GO TO 52
      GO TO 62
30    CONTINUE
      GO TO (151,152,153,154,155,156,20),JTURN
C     SINGLET-VALENCE INTERACTIONS
151   JTURN=2
      ITT1=3
      GO TO 150
C     TRIPLET-TRIPLET INTERACTIONS
152   JTURN=3
      ITURN=1
      ITT1=2
      ITT2=2
      GO TO 150
C     SINGLET-SINGLET INTERACTIONS
153   JTURN=4
      ITT1=3
      ITT2=3
      GO TO 150
C     TRIPLET-SINGLET INTERACTIONS
154   JTURN=5
      ITT1=2
      ITT2=3
      GO TO 150
C     SINGLET-TRIPLET INTERACTIONS
155   JTURN=6
      ITT1=3
      ITT2=2
      GO TO 150
C     DOUBLET-DOUBLET INTERACTIONS
156   JTURN=7
      ITT1=1
      ITT2=1
      GO TO 150
20    CONTINUE
10    CONTINUE
      ICOP1(nCOP+1)=IOUT
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+IOUT
      ICOP1(nCOP+1)=-1
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      WRITE(IW,600)NMAT
600   FORMAT(/,6X,'COEFFICIENTS FOR AIBJ',I9)
      WRITE(IW,610)(NUMM(I),I=1,7)
610   FORMAT(6X,'DIFFERENT TYPES',7I9)
      RETURN
      END
