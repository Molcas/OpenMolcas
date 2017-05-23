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
      SUBROUTINE LOOP21(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('LOOP21')
      ISTOP=0
      KM1=KM+1
      IDIF=IA(J1(KM1))-IA(J2(KM1))
      IF(IDIF.LT.-1.OR.IDIF.GT.1)GO TO 55
      IF (IDIF.LT.0) THEN
        GO TO 51
      ELSE IF (IDIF.EQ.0) THEN
        GO TO 52
      ELSE
        GO TO 53
      END IF
51    IWAYKM=IWAY(KM)
      GO TO (39,41,42,43,44,55),IWAYKM
39    IWAY(KM)=2
C     (J+R,Q+O)
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 41
      IF(K1F(JM(KM1)).EQ.0)GO TO 141
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM)=-BL1(IB(J2(KM1))+3)*COUP(KM1)/(IB(J2(KM1))+2)
      JM1(KM)=K1F(JM(KM1))
      IF(K2F(JM(KM1)).EQ.0)GO TO 40
      GO TO 241
141   IF(K2F(JM(KM1)).EQ.0)GO TO 41
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
241   COUP(KM)=BL2(IB(J2(KM1))+1)*COUP(KM1)/(IB(J2(KM1))+2)
      JM(KM)=K2F(JM(KM1))
      GO TO 40
C     I+M
41    IWAY(KM)=3
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 42
      IF(K0F(JM(KM1)).EQ.0)GO TO 42
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
      COUP(KM)=COUP(KM1)
      JM(KM)=K0F(JM(KM1))
      GO TO 40
C     J+N
42    IWAY(KM)=4
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 43
      IF(K1F(JM(KM1)).EQ.0)GO TO 43
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM)=-BL1(IB(J2(KM1))+3)*COUP(KM1)
      JM(KM)=K1F(JM(KM1))
      GO TO 40
C     K+O
43    IWAY(KM)=5
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 44
      IF(K2F(JM(KM1)).EQ.0)GO TO 44
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM)=-BL2(IB(J2(KM1))+1)*COUP(KM1)
      JM(KM)=K2F(JM(KM1))
      GO TO 40
C     L+P
44    IWAY(KM)=6
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM(KM1)).EQ.0)GO TO 55
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=COUP(KM1)
      JM(KM)=K3F(JM(KM1))
      GO TO 40
52    IWAYKM=IWAY(KM)
      GO TO (59,61,62,63,64,65,55),IWAYKM
59    IWAY(KM)=2
C     (J+J,Q+Q,N+N)
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 61
      WMP=D0
      WPP=D0
      IF(K1F(JM1(KM1)).EQ.0)GO TO 161
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM)=COUP1(KM1)*BL1(IB(J2(KM1))+1)**2
      JM1(KM)=K1F(JM1(KM1))
      IF(K2F(JM1(KM1)).NE.0)GO TO 162
      IF(K1F(JM(KM1)).NE.0)GO TO 163
      GO TO 40
161   IF(K2F(JM1(KM1)).EQ.0)GO TO 164
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
162   WMP=D1/(IB(J2(KM1))**2)
      JM(KM)=K2F(JM1(KM1))
      IF(K1F(JM(KM1)).EQ.0)GO TO 165
      GO TO 163
164   IF(K1F(JM(KM1)).EQ.0)GO TO 61
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
163   WPP=D1
      JM(KM)=K1F(JM(KM1))
165   COUP(KM)=WMP*COUP1(KM1)+WPP*COUP(KM1)
      GO TO 40
C     (O+O,K+K,R+R)
61    IWAY(KM)=3
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 62
      WMM=D0
      WPM=D0
      IF(K2F(JM(KM1)).EQ.0)GO TO 261
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM)=COUP(KM1)*BL2(IB(J2(KM1))+1)**2
      JM(KM)=K2F(JM(KM1))
      IF(K2F(JM1(KM1)).NE.0)GO TO 262
      IF(K1F(JM(KM1)).NE.0)GO TO 263
      GO TO 40
261   IF(K2F(JM1(KM1)).EQ.0)GO TO 264
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
262   WMM=D1
      JM1(KM)=K2F(JM1(KM1))
      IF(K1F(JM(KM1)).EQ.0)GO TO 265
      GO TO 263
264   IF(K1F(JM(KM1)).EQ.0)GO TO 62
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
263   WPM=D1/((IB(J2(KM1))+2)**2)
      JM1(KM)=K1F(JM(KM1))
265   COUP1(KM)=WMM*COUP1(KM1)+WPM*COUP(KM1)
      GO TO 40
C     (I+I,M+M)
62    IWAY(KM)=4
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 63
      IF(K0F(JM1(KM1)).EQ.0)GO TO 361
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
      COUP1(KM)=COUP1(KM1)
      JM1(KM)=K0F(JM1(KM1))
      IF(K0F(JM(KM1)).EQ.0)GO TO 40
      GO TO 362
361   IF(K0F(JM(KM1)).EQ.0)GO TO 63
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
362   COUP(KM)=COUP(KM1)
      JM(KM)=K0F(JM(KM1))
      GO TO 40
C     (L+L,P+P)
63    IWAY(KM)=5
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 64
      IF(K3F(JM1(KM1)).EQ.0)GO TO 461
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM)=COUP1(KM1)
      JM1(KM)=K3F(JM1(KM1))
      IF(K3F(JM(KM1)).EQ.0)GO TO 40
      GO TO 462
461   IF(K3F(JM(KM1)).EQ.0)GO TO 64
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
462   COUP(KM)=COUP(KM1)
      JM(KM)=K3F(JM(KM1))
      GO TO 40
C     (K+Q,R+N)
64    IWAY(KM)=6
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 65
      WMP=D0
      WPP=D0
      IF(K2F(JM1(KM1)).EQ.0)GO TO 561
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      WMP=-D1/IB(J2(KM1))
      JM(KM)=K2F(JM1(KM1))
      IF(K1F(JM(KM1)).EQ.0)GO TO 562
      GO TO 563
561   IF(K1F(JM(KM1)).EQ.0)GO TO 65
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
563   WPP=D1/(IB(J2(KM1))+2)
      JM(KM)=K1F(JM(KM1))
562   COUP(KM)=WMP*COUP1(KM1)+WPP*COUP(KM1)
      GO TO 40
C     (Q+K,N+R)
65    IWAY(KM)=7
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 55
      WMM=D0
      WPM=D0
      IF(K2F(JM1(KM1)).EQ.0)GO TO 661
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      WMM=-D1/IB(J2(KM1))
      JM1(KM)=K2F(JM1(KM1))
      IF(K1F(JM(KM1)).EQ.0)GO TO 662
      GO TO 663
661   IF(K1F(JM(KM1)).EQ.0)GO TO 55
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
663   WPM=D1/(IB(J2(KM1))+2)
      JM1(KM)=K1F(JM(KM1))
662   COUP1(KM)=WMM*COUP1(KM1)+WPM*COUP(KM1)
      GO TO 40
53    IWAYKM=IWAY(KM)
      GO TO (69,71,72,73,74,55),IWAYKM
69    IWAY(KM)=2
C     (R+J,O+Q)
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 71
      IF(K1F(JM1(KM1)).EQ.0)GO TO 171
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM)=-BL1(IB(J2(KM1))+1)*COUP1(KM1)/IB(J2(KM1))
      JM1(KM)=K1F(JM1(KM1))
      IF(K2F(JM1(KM1)).EQ.0)GO TO 40
      GO TO 172
171   IF(K2F(JM1(KM1)).EQ.0)GO TO 71
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
172   COUP(KM)=BL2(IB(J2(KM1))-1)*COUP1(KM1)/IB(J2(KM1))
      JM(KM)=K2F(JM1(KM1))
      GO TO 40
C     M+I
71    IWAY(KM)=3
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 72
      IF(K0F(JM1(KM1)).EQ.0)GO TO 72
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
      COUP1(KM)=COUP1(KM1)
      JM1(KM)=K0F(JM1(KM1))
      GO TO 40
C     N+J
72    IWAY(KM)=4
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 73
      IF(K1F(JM1(KM1)).EQ.0)GO TO 73
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM)=-BL1(IB(J2(KM1))+1)*COUP1(KM1)
      JM1(KM)=K1F(JM1(KM1))
      GO TO 40
C     O+K
73    IWAY(KM)=5
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 74
      IF(K2F(JM1(KM1)).EQ.0)GO TO 74
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM)=-BL2(IB(J2(KM1))-1)*COUP1(KM1)
      JM1(KM)=K2F(JM1(KM1))
      GO TO 40
C     P+L
74    IWAY(KM)=6
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM1(KM1)).EQ.0)GO TO 55
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM)=COUP1(KM1)
      JM1(KM)=K3F(JM1(KM1))
      GO TO 40
55    ISTOP=1
40      Continue
       CALL QEXIT('LOOP21')
      RETURN
      END
