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
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE TPQSET(ICASE,TPQ,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TPQ(*),IOCR(100)
      DIMENSION ICASE(*)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      LOGICAL LWSP
      COMMON /SPIN/ LWSP
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CRL   JO(L)=IAND(ISHFT(QOCC((L+29)/30),-2*((L+29)/30*30-L)),3)
CPAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
      JO(L)=ICUNP(ICASE,L)
C
      D0=0.0D0
      D1=1.0D0
      D2=2.0D0
C
      IOR=0
      II1=(IREF0-1)*LN
      DO 35 I=1,LN
      JOJ=JO(II1+I)
      IOR=IOR+1
      IOCR(IOR)=JOJ
35    CONTINUE
C
      IINT=IRC(4)
      DO 7 IQ=1,IINT
      TPQ(IQ)=D1
      IF(INCPF.EQ.1)TPQ(IQ)=D2/N
      IF(IQ.EQ.IREF0.OR.IP.EQ.IREF0)TPQ(IQ)=D0
7     CONTINUE
      IF(ISDCI.EQ.1.OR.INCPF.EQ.1.OR.IP.EQ.IREF0)RETURN
C
      II=0
      IJ=0
      DO 15 I=1,LN
      JJ=(IP-1)*LN+I
      IF(JO(JJ).EQ.IOCR(I).OR.JO(JJ).EQ.3)GO TO 15
      IF(LWSP.AND.JO(JJ)*IOCR(I).EQ.2) GO TO 15
      IF(II.NE.0)GO TO 16
      II=I
16    IJ=I
15    CONTINUE
      NI=IOCR(II)
      IF(NI.GT.1)NI=NI-1
      NJ=IOCR(IJ)
      IF(NJ.GT.1)NJ=NJ-1
      DO 20 IQ=1,IINT
      IK=0
      IL=0
      DO 25 I=1,LN
      JJ=(IQ-1)*LN+I
      IF(JO(JJ).EQ.IOCR(I).OR.JO(JJ).EQ.3)GO TO 25
      IF(LWSP.AND.JO(JJ)*IOCR(I).EQ.2) GO TO 25
      IF(IK.NE.0)GO TO 26
      IK=I
26    IL=I
25    CONTINUE
      DIK=D0
      DIL=D0
      DJK=D0
      DJL=D0
      IF(II.EQ.IK)DIK=D1
      IF(II.EQ.IL)DIL=D1
      IF(IJ.EQ.IK)DJK=D1
      IF(IJ.EQ.IL)DJL=D1
      TPQ(IQ)=(DIK+DIL)/(D2*NI)+(DJK+DJL)/(D2*NJ)
      IF(IQ.EQ.IREF0)TPQ(IQ)=D0
20    CONTINUE
      IF(IPRINT.LT.15)RETURN
      IF(IPRINT.GT.5)WRITE(6,11)(TPQ(IQ),IQ=1,IINT)
11    FORMAT(5X,'TPQ',10F5.2)
      RETURN
      END
