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
      SUBROUTINE INT9(I,J,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,
     *L0,L1,L2,L3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
C     I.LT.L . CASE 1 J=K=L , CASE 2 I=J=K .
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
c      COMMON/ADDR/IADD10,IAD10(9),IADD11,IDUM,COP(600),ICOP1(601)
      ITYP=0
      FAC=D1
      LJS=IJ(L+1)+1
      LJM=IJ(L)
      DO 10 LJ=LJS,LJM
      ITAIL=IX(IT2+LJ)
      IF(IT1.NE.IT2)CALL TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(L)=1
32    KM=L
      J2(KM+1)=LJ
      J1(KM+1)=LJ
      CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 10
      IF(I.EQ.J)GO TO 41
      FAC=D1
      IF(IWAY(KM).EQ.2.OR.IWAY(KM).EQ.4)FAC=D0
41    KM=KM-1
      IF(KM.EQ.0)GO TO 20
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 51
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.0)GO TO 41
      KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
51    IWAY(I)=1
52    KM=I
      CALL LOOP3(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 53
      IF(I.NE.J)GO TO 54
      FAC=D1
      IF(IWAY(KM).EQ.2)FAC=D0
54    IF(FAC.EQ.D0)GO TO 52
      IF(ABS(COUP(I)).LT.1.D-06)GO TO 52
      CALL COMP(I,LJ,ITYP,I,IT1,IT2)
      GO TO 52
53    KM=KM+1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
20    IF(FAC.EQ.D0)GO TO 33
      CALL COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
33    KM=1
      IF(KM.EQ.L)GO TO 32
      GO TO 42
10    CONTINUE
      RETURN
      END
