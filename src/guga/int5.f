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
      SUBROUTINE INT5(I,J,L)
C     I.LT.J.LT.L   I.EQ.K
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      ITYP=0
      LJS=IJ(L+1)+1
      LJM=IJ(L)
      DO 100 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      DO 10 LJ=LJS,LJM
      IWAY(L)=1
32    KM=L
      J2(KM+1)=LJ
      J1(KM+1)=LJ
      CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 10
41    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.J)GO TO 51
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 41
      KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
51    ITURN=0
55    IWAY(J)=1
52    KM=J
      JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      IF(ITURN.EQ.0)CALL LOOP10(KM,ISTOP,IT1,IT2)
      IFAI=0
      IF(ITURN.EQ.1)CALL LOOP13(KM,ISTOP,IFAI,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 53
      IF(ITURN.EQ.1)GO TO 54
      ITURN=1
      GO TO 55
54    KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
53    KM=KM-1
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 71
62    JM(KM)=IVF0+1
      JM1(KM)=IVF0+1
      IF(ITURN.EQ.0)CALL LOOP17(KM,ISTOP,IT1,IT2)
      IFAI=0
      IF(ITURN.EQ.1)CALL LOOP23(KM,ISTOP,IFAI,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 53
      KM=KM+1
      IF(KM.EQ.J)GO TO 52
      GO TO 62
71    KM=I
      IF(ITURN.EQ.0)CALL LOOP14(KM,ISTOP,IT1,IT2)
      IF(ITURN.EQ.1)CALL LOOP22(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 73
      IF(ABS(COUP(I)).LT.1.D-06)GO TO 71
      CALL COMP(I,LJ,ITYP,I,IT1,IT2)
      GO TO 71
73    KM=KM+1
      IF(KM.EQ.J)GO TO 52
      GO TO 62
10    CONTINUE
100   CONTINUE
      RETURN
      END
