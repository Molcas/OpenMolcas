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
      SUBROUTINE PATH(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      ISTOP=0
      KM1=KM+1
      IWAYKM=IWAY(KM)
      GO TO (132,133,134,135,131),IWAYKM
132   IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 133
      J2(KM)=K0(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)
      ICOUP1(KM)=ICOUP1(KM1)
      GO TO 140
133   IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 134
      J2(KM)=K1(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      GO TO 140
134   IWAY(KM)=4
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 135
      J2(KM)=K2(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      GO TO 140
135   IWAY(KM)=5
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 131
      J2(KM)=K3(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      GO TO 140
131   ISTOP=1
140     Continue
      RETURN
      END
