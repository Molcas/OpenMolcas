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
      SUBROUTINE INT61(I,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,
     *L0,L1,L2,L3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
C     I.LT.K.LT.L  J.EQ.L
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('INT61')
      ITYP=0
      LJS=IJ(L+1)+1
      LJM=IJ(L)
      DO 10 LJ=LJS,LJM
      ITAIL=IX(IT2+LJ)
      IF(IT1.NE.IT2)CALL TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(L)=1
32    KM=L
      J2(KM+1)=LJ
      J1(KM+1)=LJ
      JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      CALL LOOP7(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 10
53    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.K)GO TO 61
62    JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      CALL LOOP17(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 53
      KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 62
61    IWAY(K)=1
72    KM=K
      CALL LOOP16(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 63
      KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 62
63    KM=KM-1
      IF(KM.EQ.0)GO TO 20
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 74
82    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 63
      KM=KM+1
      IF(KM.EQ.K)GO TO 72
      GO TO 82
74    IWAY(I)=1
71    KM=I
      CALL LOOP3(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 73
      IF(ABS(COUP(I)).LT.1.D-06)GO TO 71
      CALL COMP(I,LJ,ITYP,I,IT1,IT2)
      GO TO 71
73    KM=KM+1
      IF(KM.EQ.K)GO TO 72
      GO TO 82
20    CALL COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
      KM=1
      IF(KM.EQ.K)GO TO 72
      GO TO 82
10    CONTINUE
      CALL QEXIT('INT61')
      RETURN
      END
