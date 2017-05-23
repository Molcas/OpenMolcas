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
      SUBROUTINE INT8(I,J,L,IT1,IT2,II,IID,JJT,JJD,JTYP,ITAI,
     *L0,L1,L2,L3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ITAI(*),L0(*),L1(*),L2(*),L3(*)
C     K.EQ.L , I.LT.J
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
c      COMMON/ADDR/IADD10,IAD10(9),IADD11,IDUM,COP(600),ICOP1(601)
      CALL QENTER('INT8')
      ITYP=0
      IF(L.LT.I.OR.L.GT.J)ITYP=1
      IF(I.EQ.0)ITYP=0
      IF(I.EQ.0.AND.L.GT.J)ITYP=1
      FAC=D1
      JJS=IJ(J+1)+1
      JJM=IJ(J)
      DO 10 JJ=JJS,JJM
      ITAIL=IX(IT2+JJ)
      IF(IT1.NE.IT2)CALL TAIL(J,JJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(J)=1
32    KM=J
      J2(KM+1)=JJ
      J1(KM+1)=JJ
      CALL LOOP1(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 10
41    KM=KM-1
      IF(KM.EQ.0)GO TO 20
      IWAY(KM)=1
      IF(KM.EQ.I)GO TO 51
42    CALL LOOP5(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 43
      IF(KM.NE.L)GO TO 41
      IF(IWAY(KM).EQ.2)GO TO 42
      FAC=D1
      IF(IWAY(KM).EQ.5)FAC=D2
      GO TO 41
43    KM=KM+1
      IF(KM.EQ.J)GO TO 32
      GO TO 42
51    IWAY(I)=1
52    KM=I
      CALL LOOP3(KM,ISTOP,IT1,IT2)
      IF(ISTOP.EQ.1)GO TO 53
      IF(ABS(COUP(I)).LT.1.D-06)GO TO 52
      COUP(I)=FAC*COUP(I)
      CALL COMP(I,JJ,ITYP,L,IT1,IT2)
      GO TO 52
53    KM=KM+1
      IF(KM.EQ.J)GO TO 32
      GO TO 42
20    COUP(1)=FAC*COUP(1)
      CALL COMP1(JJ,ITYP,L,IT2,II,IID,JJT,JJD,JTYP,ITAI)
      KM=1
      IF(KM.EQ.J)GO TO 32
      GO TO 42
10    CONTINUE
      CALL QEXIT('INT8')
      RETURN
      END
