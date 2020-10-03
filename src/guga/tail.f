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
      SUBROUTINE TAIL(LL,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
#include "integ.fh"
       L=LL+1
      IF(ITAIL.EQ.0)THEN
        RETURN
      END IF
      DO 5 I=1,ITAIL
      ITAI(I)=0
5     CONTINUE
      IF(L.EQ.LN+1)ITAI(1)=1
      KM=L
      J2(KM)=IJJ
      ICOUP(L)=1
      ICOUP1(L)=1
11    KM=KM+1
      IWAY(KM)=0
12    KM1=KM-1
      IF(IWAY(KM).GE.1)GO TO 14
      IF(L0(IT1+J2(KM1)).EQ.0.OR.L0(IT2+J2(KM1)).EQ.0)GO TO 14
      J2(KM)=L0(IT2+J2(KM1))
      IWAY(KM)=1
      ICOUP(KM)=ICOUP(KM1)
      ICOUP1(KM)=ICOUP1(KM1)
      GO TO 20
14    IF(IWAY(KM).GE.2)GO TO 15
      IF(L1(IT1+J2(KM1)).EQ.0.OR.L1(IT2+J2(KM1)).EQ.0)GO TO 15
      J2(KM)=L1(IT2+J2(KM1))
      IWAY(KM)=2
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM),1)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J2(KM),1)
      GO TO 20
15    IF(IWAY(KM).GE.3)GO TO 16
      IF(L2(IT1+J2(KM1)).EQ.0.OR.L2(IT2+J2(KM1)).EQ.0)GO TO 16
      J2(KM)=L2(IT2+J2(KM1))
      IWAY(KM)=3
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM),2)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J2(KM),2)
      GO TO 20
16    IF(IWAY(KM).GE.4)GO TO 17
      IF(L3(IT1+J2(KM1)).EQ.0.OR.L3(IT2+J2(KM1)).EQ.0)GO TO 17
      J2(KM)=L3(IT2+J2(KM1))
      IWAY(KM)=4
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM),3)
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J2(KM),3)
      GO TO 20
17    KM=KM-1
      IF(KM.EQ.L)GO TO 10
      GO TO 12
20    IF(KM.NE.LN+1)GO TO 11
      ITAI(ICOUP(LN+1))=ICOUP1(LN+1)
      GO TO 12
10      Continue
      RETURN
      END
