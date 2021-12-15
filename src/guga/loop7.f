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
      SUBROUTINE LOOP7(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      ISTOP=0
      KM1=KM+1
      J2F=IPO(J2(KM1))
      IF(IWAY(KM).EQ.2)GO TO 55
      IWAY(KM)=2
C     (CB,AD)
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K1F(J2F).EQ.0)GO TO 141
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=IY(IT2+J2(KM1),3)
      COUP1(KM)=BS1(IB(J2(KM1))+1)
      JM1(KM)=K1F(J2F)
      IF(K2F(J2F).EQ.0)GO TO 40
      GO TO 143
141   IF(K2F(J2F).EQ.0)GO TO 55
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=0
      ICOUP(KM)=IY(IT2+J2(KM1),3)
143   COUP(KM)=BS2(IB(J2(KM1))+1)
      JM(KM)=K2F(J2F)
      GO TO 40
55    ISTOP=1
40      Continue
      RETURN
      END
