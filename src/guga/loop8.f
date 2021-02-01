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
      SUBROUTINE LOOP8(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      ISTOP=0
      KM1=KM+1
      J2F=IPO(J2(KM1))
      IWAYKM=IWAY(KM)
      GO TO (39,41,42,43,55),IWAYKM
39    IWAY(KM)=2
C     (B+B,D+D)
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 41
      IF(K1F(J2F).EQ.0)GO TO 141
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=J1(KM)
      ICOUP1(KM)=IY(IT1+J1(KM1),3)
      ICOUP(KM)=IY(IT2+J2(KM1),3)
      COUP1(KM)=BS1(IB(J2(KM1))+1)**2
      JM1(KM)=K1F(J2F)
      IF(K2F(J2F).EQ.0)GO TO 40
      GO TO 142
141   IF(K2F(J2F).EQ.0)GO TO 41
      J1(KM)=K3(IT1+J1(KM1))
      J2(KM)=J1(KM)
      ICOUP1(KM)=IY(IT1+J1(KM1),3)
      ICOUP(KM)=IY(IT2+J2(KM1),3)
142   COUP(KM)=BS2(IB(J2(KM1))+1)**2
      JM(KM)=K2F(J2F)
      GO TO 40
C     A+A
41    IWAY(KM)=3
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 42
      IF(K0F(J2F).EQ.0)GO TO 42
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=J1(KM)
      ICOUP1(KM)=IY(IT1+J1(KM1),2)
      ICOUP(KM)=IY(IT2+J2(KM1),2)
      COUP1(KM)=D1
      JM1(KM)=K0F(J2F)
      GO TO 40
C     C+C
42    IWAY(KM)=4
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 43
      IF(K0F(J2F).EQ.0)GO TO 43
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=J1(KM)
      ICOUP1(KM)=IY(IT1+J1(KM1),1)
      ICOUP(KM)=IY(IT2+J2(KM1),1)
      COUP(KM)=D1
      JM(KM)=K0F(J2F)
      GO TO 40
C     C+A
43    IWAY(KM)=5
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K0F(J2F).EQ.0)GO TO 55
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=IY(IT1+J1(KM1),1)
      ICOUP(KM)=IY(IT2+J2(KM1),2)
      COUP1(KM)=D1
      JM1(KM)=K0F(J2F)
      GO TO 40
55    ISTOP=1
40      Continue
      RETURN
      END
