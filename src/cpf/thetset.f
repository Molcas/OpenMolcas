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
      SUBROUTINE THETSET(ICASE,THE,NII)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION THE(NII,NII),IOCR(100)
      DIMENSION ICASE(*)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      LOGICAL LWSP
      COMMON /SPIN/ LWSP
      JO(L)=ICUNP(ICASE,L)
CPAM97      EXTERNAL UNPACK
CPAM97      INTEGER UNPACK
CPAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
*PAM06 This routine is called from SDCI if this is an MCPF calculation
C
      IOR=0
      II1=(IREF0-1)*LN
      DO 35 I=1,LN
         JOJ=JO(II1+I)
         IOR=IOR+1
         IOCR(IOR)=JOJ
35    CONTINUE
      IF(IPRINT.GT.5)WRITE(6,888)IREF0,(IOCR(I),I=1,LN)
888   FORMAT(5X,'IREF0=',I3/5X,'IOCR=',10I5)
C
      IINT=IRC(4)
      DO 8 IP=1,IINT
         DO 7 IQ=1,IINT
            THE(IQ,IP)=D1
7        CONTINUE
8     CONTINUE
C
      DO 6 IP=1,IINT
         DO 5 IQ=1,IINT
            THE(IQ,IP)=0.0D0
5        CONTINUE
6     CONTINUE
      DO 10 IP=1,IINT
         II=0
         IJ=0
         DO 15 I=1,LN
            JJ=(IP-1)*LN+I
            IF(JO(JJ).EQ.IOCR(I).OR.JO(JJ).EQ.3)GO TO 15
            IF(LWSP.AND.JO(JJ)*IOCR(I).EQ.2) GO TO 15
            IF(II.NE.0)GO TO 16
            II=I
16          IJ=I
15       CONTINUE
*PAM06 BUG: What if we come down here with II.eq.0 still??
* the IOCR will be accessed below first element. Provisional fix:
         IF(II.EQ.0) GOTO 10
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
26             IL=I
25          CONTINUE
            IF(II.EQ.IK.AND.IJ.EQ.IL)THE(IQ,IP)=1.0D0
20       CONTINUE
10    CONTINUE
      RETURN
      END
