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
      SUBROUTINE LOOP1(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      ISTOP=0
      KM1=KM+1
C     FOUR POSSIBLE STARTS
      IWAYKM=IWAY(KM)
      GO TO (34,35,36,37,30),IWAYKM
34    IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 35
      COUP(KM)=D1
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=IY(IT2+J2(KM1),2)
      GO TO 40
35    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 36
      COUP(KM)=BS1(IB(J2(KM1))+1)
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=IY(IT1+J1(KM1),1)
      ICOUP(KM)=IY(IT2+J2(KM1),3)
      GO TO 40
36    IWAY(KM)=4
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 37
      COUP(KM)=D1
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=IY(IT2+J2(KM1),1)
      GO TO 40
37    IWAY(KM)=5
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 30
      COUP(KM)=BS2(IB(J2(KM1))+1)
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=IY(IT1+J1(KM1),2)
      ICOUP(KM)=IY(IT2+J2(KM1),3)
      GO TO 40
30    ISTOP=1
40      Continue
      RETURN
      END
