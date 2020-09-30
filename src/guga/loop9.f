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
      SUBROUTINE LOOP9(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      ISTOP=0
      KM1=KM+1
      J1F=IPO(J1(KM1))
      IDIF=IA(J1(KM1))-IA(J2(KM1))
      IF(IDIF.LT.0.OR.IDIF.GT.1)GO TO 55
      IF(IDIF.EQ.1)GO TO 41
      IF(IWAY(KM).EQ.2)GO TO 55
      IWAY(KM)=2
C     B+G
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K1F(J1F).EQ.0)GO TO 55
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)
      COUP(KM)=BS1(IB(J2(KM1))+2)*COUP(KM1)
      GO TO 40
C     D+E
41    IF(IWAY(KM).EQ.2)GO TO 55
      IWAY(KM)=2
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K2F(J1F).EQ.0)GO TO 55
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=K0(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)
      COUP(KM)=BS2(IB(J2(KM1)))*COUP(KM1)
      GO TO 40
55    ISTOP=1
40      Continue
      RETURN
      END
