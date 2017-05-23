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
      SUBROUTINE LOOP6(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('LOOP6')
      ISTOP=0
      KM1=KM+1
      IDIF=IA(J2(KM1))-IA(J1(KM1))
      IF(IDIF.LT.0.OR.IDIF.GT.1)GO TO 55
      IF(IDIF.EQ.0)GO TO 45
      IWAYKM=IWAY(KM)
      GO TO (39,41,42,43,44,55),IWAYKM
C     CASE I-L AND Q
39    IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 41
      COUP(KM)=COUP(KM1)
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
      GO TO 40
41    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 42
      COUP(KM)=BL1(IB(J1(KM1))+1)*COUP(KM1)
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      GO TO 40
42    IWAY(KM)=4
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 43
      COUP(KM)=-COUP(KM1)
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      GO TO 40
43    IWAY(KM)=5
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 44
      COUP(KM)=-COUP(KM1)
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      GO TO 40
44    IWAY(KM)=6
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 55
      COUP(KM)=COUP(KM1)/IB(J1(KM1))
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      GO TO 40
C     CASE M-P AND R
45    IWAYKM=IWAY(KM)
      GO TO (59,61,62,63,64,55),IWAYKM
59    IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 61
      COUP(KM)=COUP(KM1)
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)
      GO TO 40
61    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 62
      COUP(KM)=-COUP(KM1)
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      GO TO 40
62    IWAY(KM)=4
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 63
      COUP(KM)=BL2(IB(J1(KM1))+1)*COUP(KM1)
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      GO TO 40
63    IWAY(KM)=5
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 64
      COUP(KM)=-COUP(KM1)
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      GO TO 40
64    IWAY(KM)=6
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 55
      COUP(KM)=-COUP(KM1)/(IB(J1(KM1))+2)
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      GO TO 40
55    ISTOP=1
40      Continue
       CALL QEXIT('LOOP6')
      RETURN
      END
