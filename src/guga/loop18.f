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
      SUBROUTINE LOOP18(KM,ISTOP,IT1,IT2)
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
51    IF(IWAY(KM).EQ.2)GO TO 55
      IWAY(KM)=2
C     F+H
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM(KM1)).EQ.0)GO TO 55
      J2(KM)=K2(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM)=BS3(IB(J2(KM1))+3)*BS4(IB(J2(KM1))+1)*COUP(KM1)
      GO TO 40
52    IWAYKM=IWAY(KM)
      GO TO (59,61,62,55),IWAYKM
59    IWAY(KM)=2
C     (E+E,G+G)
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 61
      WM0=D0
      WP0=D0
      IF(K2F(JM1(KM1)).EQ.0)GO TO 161
      J2(KM)=K0(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
      WM0=D1
      IF(K1F(JM(KM1)).EQ.0)GO TO 162
      GO TO 163
161   IF(K1F(JM(KM1)).EQ.0)GO TO 61
      J2(KM)=K0(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
163   WP0=D1
162   COUP(KM)=WM0*COUP1(KM1)+WP0*COUP(KM1)
      GO TO 40
C     F+F
61    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 62
      IF(K3F(JM1(KM1)).EQ.0)GO TO 62
      J2(KM)=K1(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM)=COUP1(KM1)*BS3(IB(J2(KM1))+1)**2
      GO TO 40
C     H+H
62    IWAY(KM)=4
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM(KM1)).EQ.0)GO TO 55
      J2(KM)=K2(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP(KM)=COUP(KM1)*BS4(IB(J2(KM1))+1)**2
      GO TO 40
53    IF(IWAY(KM).EQ.2)GO TO 55
      IWAY(KM)=2
C     H+F
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(JM1(KM1)).EQ.0)GO TO 55
      J2(KM)=K1(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM)=BS4(IB(J2(KM1))-1)*BS3(IB(J2(KM1))+1)*COUP1(KM1)
      GO TO 40
55    ISTOP=1
40      Continue
      RETURN
      END
