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
      SUBROUTINE CIALL(LSYM,NREF,IOCR,nIOCR,L0,L1,L2,L3,LV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IOCR(nIOCR),L0(*),L1(*),L2(*),L3(*)
#include "integ.fh"
      DIMENSION IOC(55)
*
      IIN=0
      NREF=0
      IJJ=IV0
      KM=1
      J2(KM)=IJJ
11    KM=KM+1
      IWAY(KM)=0
12    KM1=KM-1
      IF(L0(J2(KM1)).EQ.0.OR.IWAY(KM).GE.1)GO TO 14
      J2(KM)=L0(J2(KM1))
      IWAY(KM)=1
      IOC(KM1)=0
      GO TO 20
14    IF(L1(J2(KM1)).EQ.0.OR.IWAY(KM).GE.2)GO TO 15
      J2(KM)=L1(J2(KM1))
      IWAY(KM)=2
      IOC(KM1)=1
      GO TO 20
15    IF(L2(J2(KM1)).EQ.0.OR.IWAY(KM).GE.3)GO TO 16
      J2(KM)=L2(J2(KM1))
      IWAY(KM)=3
      IOC(KM1)=1
      GO TO 20
16    IF(L3(J2(KM1)).EQ.0.OR.IWAY(KM).GE.4)GO TO 17
      J2(KM)=L3(J2(KM1))
      IWAY(KM)=4
      IOC(KM1)=2
      GO TO 20
17    KM=KM-1
      IF(KM.EQ.1)GO TO 10
      GO TO 12
20    IF(KM.NE.LN+1)GO TO 11
      NSJ=1
      DO 110 I=1,LN
      IF(I.GT.LV)GO TO 113
      IF(IOC(I).NE.0)GO TO 12
      GO TO 110
113   IF(I.GT.NIORB+LV)GO TO 114
      IF(IOC(I).NE.2)GO TO 12
      GO TO 110
114   IF(IOC(I).EQ.1)NSJ=MUL(NSJ,NSM(I))
110   CONTINUE
      IF(NSJ.NE.LSYM)GO TO 12
      NREF=NREF+1
      DO 111 I=1,LN
         IF(I.LE.NIORB+LV)GO TO 111
         IIN=IIN+1
         IF (IIN.GT.nIOCR) Then
            Write (6,*) 'CIall: IIN.GT.nIOCR'
            Write (6,*) 'IIN=',IIN
            Write (6,*) 'nIOCR=',nIOCR
            Call Abend()
         End If
         IOCR(IIN)=IOC(I)
111   CONTINUE
      GO TO 12
*
10      Continue
      RETURN
      End
