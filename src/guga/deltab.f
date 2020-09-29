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
      SUBROUTINE DELTAB(NREF,IOCR,L0,L1,L2,L3,INTNUM,LV,IFCORE,
     *ICOR,NONE,JONE,K00,K11,K22,K33,L00,L11,L22,L33)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IOCR(*),L0(*),L1(*),L2(*),L3(*),ICOR(*),JONE(*),
     *K00(*),K11(*),K22(*),K33(*),L00(*),L11(*),L22(*),L33(*)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      DIMENSION IOC(55),ISP(55)
      DIMENSION K0M(MXVERT),K1M(MXVERT),K2M(MXVERT),K3M(MXVERT)
      DIMENSION L0M(MXVERT),L1M(MXVERT),L2M(MXVERT),L3M(MXVERT)
      DO 2 I=1,4*MXVERT
      K0(I)=0
      K1(I)=0
      K2(I)=0
      K3(I)=0
      L0(I)=0
      L1(I)=0
      L2(I)=0
      L3(I)=0
2     CONTINUE
      IBS=0
      IEL=2
      LSYM=1
      LNS=NIORB+LV+1
      IRR=0
      DO 5 I=LNS,LN
      IRR=IRR+1
      IF(IOCR(IRR).EQ.1)LSYM=MUL(LSYM,NSM(I))
5     CONTINUE
      DO 10 IIJ=1,ILIM
      ISTA=(IIJ-1)*MXVERT
      DO 3 I=1,IV0
      K0M(I)=0
      K1M(I)=0
      K2M(I)=0
      K3M(I)=0
      L0M(I)=0
      L1M(I)=0
      L2M(I)=0
      L3M(I)=0
3     CONTINUE
      IJJ=IV0+1-IIJ
      KM=1
      J2(KM)=IJJ
11    KM=KM+1
      IWAY(KM)=0
12    KM1=KM-1
      IF(L00(J2(KM1)).EQ.0.OR.IWAY(KM).GE.1)GO TO 14
      J2(KM)=L00(J2(KM1))
      IWAY(KM)=1
      IOC(KM1)=0
      ISP(KM1)=0
      GO TO 20
14    IF(L11(J2(KM1)).EQ.0.OR.IWAY(KM).GE.2)GO TO 15
      J2(KM)=L11(J2(KM1))
      IWAY(KM)=2
      IOC(KM1)=1
      ISP(KM1)=1
      GO TO 20
15    IF(L22(J2(KM1)).EQ.0.OR.IWAY(KM).GE.3)GO TO 16
      J2(KM)=L22(J2(KM1))
      IWAY(KM)=3
      IOC(KM1)=1
      ISP(KM1)=2
      GO TO 20
16    IF(L33(J2(KM1)).EQ.0.OR.IWAY(KM).GE.4)GO TO 17
      J2(KM)=L33(J2(KM1))
      IWAY(KM)=4
      IOC(KM1)=2
      ISP(KM1)=3
      GO TO 20
17    KM=KM-1
      IF(KM.EQ.1)GO TO 210
      GO TO 12
20    IF(KM1.EQ.NIORB+LV)IBS=IB(J2(KM))
      IF(KM.NE.LN+1)GO TO 11
      NSJ=1
      INHOLE=0
      DO 110 I=1,LN
      IF(IOC(I).EQ.1)NSJ=MUL(NSJ,NSM(I))
      IF(I.LE.NIORB+LV.AND.I.GT.LV)INHOLE=INHOLE+2-IOC(I)
110   CONTINUE
C     STRIKE OUT INTERNAL CONFIGURATIONS
      IPART=0
      IF(IIJ.GT.1)IPART=IPART+1
      IF(IIJ.GT.2)IPART=IPART+1
      JJ1=0
      DO 111 IREF=1,NREF
      JHOLE=0
      JPART=IPART
      DO 112 I=1,LN
      IF(I.GT.LV)GO TO 250
      IDIF=IOC(I)
      GO TO 251
250   IF(I.GT.NIORB+LV)GO TO 252
      IDIF=IOC(I)-2
      GO TO 251
252   JJ1=JJ1+1
      IF(IOC(I).EQ.IOCR(JJ1))GO TO 112
      IDIF=IOC(I)-IOCR(JJ1)
251   IF(IDIF.GT.0)GO TO 114
      JHOLE=JHOLE-IDIF
      GO TO 112
114   JPART=JPART+IDIF
112   CONTINUE
      If (JPART.NE.JHOLE) Then
         Write (6,*) 'DeltaB: JPART.NE.JHOLE'
         Write (6,*) 'JPART,JHOLE=',JPART,JHOLE
         Write (6,*) 'iREF=',iREF
         Call QTrace
         Call Abend
      End If
      IF(JPART.LE.IEL)GO TO 113
111   CONTINUE
      GO TO 12
113   IF(IPART.EQ.0.AND.NSJ.NE.LSYM)GO TO 12
      IF(IPART.NE.2.OR.INTNUM.EQ.0)GO TO 115
C     INTERACTING SPACE
      IF(INHOLE.EQ.2.AND.IBS.NE.0)GO TO 12
C     NO CORE-CORE CORRELATION
115   IF(IFCORE.EQ.0)GO TO 116
      NCORR=0
      DO 117 I=1,LN
      IF(ICOR(I).EQ.0)GO TO 117
      NCORR=NCORR+2-IOC(I)
117   CONTINUE
      IF(NCORR.GT.1)GO TO 12
C     SINGLY OCCUPIED ORBITALS
116   IF(NONE.EQ.0)GO TO 118
      DO 119 I=1,NONE
      IF(IOC(JONE(I)).NE.1)GO TO 12
119   CONTINUE
118   DO 130 K=1,LN
      IF (ISP(K)-1.LT.0) THEN
         GO TO 131
      ELSE IF (ISP(K)-1.EQ.0) THEN
         GO TO 132
      ELSE
         GO TO 133
      END IF
131   K0M(J2(K+1))=K00(J2(K+1))
      L0M(J2(K))=L00(J2(K))
      GO TO 130
132   K1M(J2(K+1))=K11(J2(K+1))
      L1M(J2(K))=L11(J2(K))
      GO TO 130
133   IF(ISP(K).EQ.3)GO TO 134
      K2M(J2(K+1))=K22(J2(K+1))
      L2M(J2(K))=L22(J2(K))
      GO TO 130
134   K3M(J2(K+1))=K33(J2(K+1))
      L3M(J2(K))=L33(J2(K))
130   CONTINUE
      GO TO 12
210   DO 1 I=1,IV0
      K0(ISTA+I)=K0M(I)
      K1(ISTA+I)=K1M(I)
      K2(ISTA+I)=K2M(I)
      K3(ISTA+I)=K3M(I)
      L0(ISTA+I)=L0M(I)
      L1(ISTA+I)=L1M(I)
      L2(ISTA+I)=L2M(I)
      L3(ISTA+I)=L3M(I)
1     CONTINUE
10    CONTINUE
*
      Return
      End
