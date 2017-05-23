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
      SUBROUTINE LOOP26(KM,ISTOP,IFAI,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('LOOP26')
      ISTOP=0
      KM1=KM+1
      J2F=IPO(J2(KM1))
C     FOUR POSSIBLE STARTS
      IWAYKM=IWAY(KM)
      GO TO (34,35,36,37,55),IWAYKM
C     (CC+,AA+)
34    IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 35
      IF(K1F(J2F).EQ.0)GO TO 135
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=0
      COUP1(KM)=D1
      JM1(KM)=K1F(J2F)
      IF(K2F(J2F).EQ.0)GO TO 40
      GO TO 136
135   IF(K2F(J2F).EQ.0)GO TO 137
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=0
136   COUP(KM)=D1
      JM(KM)=K2F(J2F)
      GO TO 40
137   IF(IFAI.EQ.0)GO TO 35
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=0
      COUP1(KM)=D0
      COUP(KM)=D0
      GO TO 40
C     BB+
35    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 36
      IF(K3F(J2F).EQ.0)GO TO 138
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      JM(KM)=K3F(J2F)
      ICOUP1(KM)=IY(IT1+J1(KM1),1)
      ICOUP(KM)=IY(IT2+J2(KM1),1)
      COUP(KM)=BS1(IB(J2(KM1))+1)**2
      GO TO 40
138   IF(IFAI.EQ.0)GO TO 36
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=IY(IT1+J1(KM1),1)
      ICOUP(KM)=IY(IT2+J2(KM1),1)
      COUP1(KM)=D0
      COUP(KM)=D0
      GO TO 40
C     DD+
36    IWAY(KM)=4
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 37
      IF(K3F(J2F).EQ.0)GO TO 139
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      JM1(KM)=K3F(J2F)
      ICOUP1(KM)=IY(IT1+J1(KM1),2)
      ICOUP(KM)=IY(IT2+J2(KM1),2)
      COUP1(KM)=BS2(IB(J2(KM1))+1)**2
      GO TO 40
139   IF(IFAI.EQ.0)GO TO 37
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=IY(IT1+J1(KM1),2)
      ICOUP(KM)=IY(IT2+J2(KM1),2)
      COUP1(KM)=D0
      COUP(KM)=D0
      GO TO 40
C     BD+
37    IWAY(KM)=5
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(J2F).EQ.0)GO TO 55
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      JM1(KM)=K3F(J2F)
      ICOUP1(KM)=IY(IT1+J1(KM1),1)
      ICOUP(KM)=IY(IT2+J2(KM1),2)
      COUP1(KM)=BS1(IB(J2(KM1))+1)*BS2(IB(J2(KM1))+1)
      GO TO 40
55    ISTOP=1
40      Continue
       CALL QEXIT('LOOP26')
      RETURN
      END
