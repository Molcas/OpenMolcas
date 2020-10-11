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
      SUBROUTINE LOOP19(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
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
      GO TO (39,41,55),IWAYKM
39    IWAY(KM)=2
C     E+O
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 41
      IF(K2F(JM(KM1)).EQ.0)GO TO 41
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM)=BL2(IB(J2(KM1))+1)*COUP(KM1)
      GO TO 40
C     F+P
41    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM(KM1)).EQ.0)GO TO 55
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=-BS3(IB(J2(KM1))+3)*COUP(KM1)
      GO TO 40
52    IWAYKM=IWAY(KM)
      GO TO (59,61,62,63,55),IWAYKM
59    IWAY(KM)=2
C     (E+K,G+R)
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 61
      WM0=D0
      WP0=D0
      IF(K2F(JM1(KM1)).EQ.0)GO TO 161
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      WM0=-D1
      IF(K1F(JM(KM1)).EQ.0)GO TO 162
      GO TO 163
161   IF(K1F(JM(KM1)).EQ.0)GO TO 61
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
163   WP0=-D1/(IB(J2(KM1))+2)
162   COUP(KM)=WM0*COUP1(KM1)+WP0*COUP(KM1)
      GO TO 40
C     (E+Q,G+N)
61    IWAY(KM)=3
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 62
      WM0=D0
      WP0=D0
      IF(K2F(JM1(KM1)).EQ.0)GO TO 261
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      WM0=D1/IB(J2(KM1))
      IF(K1F(JM(KM1)).EQ.0)GO TO 262
      GO TO 263
261   IF(K1F(JM(KM1)).EQ.0)GO TO 62
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
263   WP0=-D1
262   COUP(KM)=WM0*COUP1(KM1)+WP0*COUP(KM1)
      GO TO 40
C     F+L
62    IWAY(KM)=4
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 63
      IF(K3F(JM1(KM1)).EQ.0)GO TO 63
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=-BS3(IB(J2(KM1))+1)*COUP1(KM1)
      GO TO 40
C     H+P
63    IWAY(KM)=5
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM(KM1)).EQ.0)GO TO 55
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=-BS4(IB(J2(KM1))+1)*COUP(KM1)
      GO TO 40
53    IWAYKM=IWAY(KM)
      GO TO (69,71,55),IWAYKM
69    IWAY(KM)=2
C     G+J
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 71
      IF(K1F(JM1(KM1)).EQ.0)GO TO 71
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM)=BL1(IB(J2(KM1))+1)*COUP1(KM1)
      GO TO 40
C     H+L
71    IWAY(KM)=3
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM1(KM1)).EQ.0)GO TO 55
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=-BS4(IB(J2(KM1))-1)*COUP1(KM1)
      GO TO 40
55    ISTOP=1
40      Continue
      RETURN
      END
