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
      SUBROUTINE LOOP4(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('LOOP4')
      ISTOP=0
C     STOP THE LOOP
      KM1=KM+1
      IDIF=IA(J2(KM1))-IA(J1(KM1))
      IF(IDIF.LT.0.OR.IDIF.GT.1)GO TO 52
      IF(IDIF.EQ.0)GO TO 60
C     CASE E-F
      IWAYKM=IWAY(KM)
      GO TO (49,51,52),IWAYKM
49    IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 51
      COUP(KM)=COUP(KM1)
      J2(KM)=K0(IT1+J1(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      ICOUP1(KM)=ICOUP1(KM1)
      GO TO 40
51    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 52
      COUP(KM)=COUP(KM1)*BS3(IB(J1(KM1))+1)
      J2(KM)=K1(IT1+J1(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      GO TO 40
C     CASE G-H
60    IWAYKM=IWAY(KM)
      GO TO (64,65,52),IWAYKM
64    IWAY(KM)=2
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 65
      COUP(KM)=COUP(KM1)
      J2(KM)=K0(IT1+J1(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      ICOUP1(KM)=ICOUP1(KM1)
      GO TO 40
65    IWAY(KM)=3
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 52
      COUP(KM)=COUP(KM1)*BS4(IB(J1(KM1))+1)
      J2(KM)=K2(IT1+J1(KM1))
      J1(KM)=J2(KM)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      GO TO 40
52    ISTOP=1
40      Continue
       CALL QEXIT('LOOP4')
      RETURN
      END
