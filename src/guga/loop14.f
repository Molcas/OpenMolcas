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
      SUBROUTINE LOOP14(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('LOOP14')
      ISTOP=0
      KM1=KM+1
      IDIF=IA(J1(KM1))-IA(J2(KM1))
      IF(IDIF.NE.1)GO TO 55
      IF(IWAY(KM).EQ.2)GO TO 55
      IWAY(KM)=2
C     (HE,FG)
      IF(K3(IT1+J1(KM1)).EQ.0.OR.K0(IT2+J2(KM1)).EQ.0)GO TO 55
      WM0=D0
      WP0=D0
      IF(K2F(JM1(KM1)).EQ.0)GO TO 141
      J2(KM)=K0(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)
      WM0=BS4(IB(J2(KM1)))
      IF(K1F(JM(KM1)).EQ.0)GO TO 142
      GO TO 143
141   IF(K1F(JM(KM1)).EQ.0)GO TO 55
      J2(KM)=K0(IT2+J2(KM1))
      J1(KM)=J2(KM)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),3)
      ICOUP(KM)=ICOUP(KM1)
143   WP0=BS3(IB(J2(KM1))+2)
142   COUP(KM)=WM0*COUP1(KM1)+WP0*COUP(KM1)
      GO TO 40
55    ISTOP=1
40      Continue
       CALL QEXIT('LOOP14')
      RETURN
      END
