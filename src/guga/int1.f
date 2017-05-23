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
      SUBROUTINE INT1(I,J,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,
     *L0,L1,L2,L3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
C     I.LT.J.LT.K.LT.L
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('INT1')
      LJS=IJ(L+1)+1
      LJM=IJ(L)
      ITYP=0
      DO 10 LJ=LJS,LJM
      ITAIL=IX(IT2+LJ)
      IF(IT1.NE.IT2)CALL TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(L)=1
32    KM=L
      J2(KM+1)=LJ
      J1(KM+1)=LJ
      CALL LOOP1(L,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 10
41    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.K)GO TO 51
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 41
      KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
51    IWAY(K)=1
52    KM=K
      CALL LOOP3(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 30
      KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
30    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.J)GO TO 133
62    CALL PATH(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 30
      KM=KM+1
      IF(KM.EQ.K)GO TO 52
      GO TO 62
133   IWAY(J)=1
132   KM=J
      IF(IT1.GE.IT2)CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(IT2.GT.IT1)CALL LOOP2(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 141
      KM=KM+1
      IF(KM.EQ.K)GO TO 52
      GO TO 62
141   KM=KM-1
      IF(KM.EQ.0)GO TO 20
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 151
142   IF(IT1.GE.IT2)CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(IT2.GT.IT1)CALL LOOP6(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 141
      KM=KM+1
      IF(KM.EQ.J)GO TO 132
      GO TO 142
151   IWAY(I)=1
152   KM=I
      IF(IT1.GE.IT2)CALL LOOP3(KM,ISTOP,IT1,IT2)
      IF(IT2.GT.IT1)CALL LOOP4(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 153
      COUP(I)=COUP(I)*COUP(K)
      IF(ABS(COUP(I)).LT.1.D-06)GO TO 152
      ICPI=ICOUP(I)
      ICP1I=ICOUP1(I)
      ICOUP(I)=ICOUP(J+1)+ICPI
      ICOUP1(I)=ICOUP1(J+1)+ICP1I
      CALL COMP(I,LJ,ITYP,I,IT1,IT2)
      ICOUP(I)=ICOUP(J+1)+ICP1I
      ICOUP1(I)=ICOUP1(J+1)+ICPI
      CALL COMP(I,LJ,ITYP,I,IT1,IT2)
      GO TO 152
153   KM=KM+1
      IF(KM.EQ.J)GO TO 132
      GO TO 142
20    COUP(1)=COUP(1)*COUP(K)
      IF(ABS(COUP(1)).LT.1.D-06)GO TO 154
      ICOUP(1)=ICOUP(J+1)+ICOUP(1)
      ICOUP1(1)=ICOUP1(J+1)+ICOUP1(1)
      CALL COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
154   KM=1
      IF(KM.EQ.J)GO TO 132
      GO TO 142
10    CONTINUE
      CALL QEXIT('INT1')
      RETURN
      END
