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
      SUBROUTINE TAB2F(IVER,LV)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      DIMENSION IORB(MXVERT)
      nijj=0
      IEL=2
      IF(IFIRST.NE.0)IEL=1
      IEL=IEL+1
      IUT=0
      IBF(1)=INT(2*S)
      IAF(1)=INT(N-2*S)/2
      IJF(LN+1)=0
      IJF(LN)=1
      NIJ=1
      IJR=1
      IJS=2
      IJRL=IJR
      IORB(1)=0
      DO 10 II=1,LN
      IIM=LN-II+1-LV
C     S=0
16    NIJ=NIJ+1
c      IF(NIJ.GT.IVER)
      nijj=max(nij,nijj)
      IAF(NIJ)=IAF(IJR)
      IBF(NIJ)=IBF(IJR)
      IORB(NIJ)=IORB(IJR)+2
      IF(IIM.LE.0)GO TO 11
      CALL CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
      IF(ISTOP.EQ.1)NIJ=NIJ-1
11    IF(IBF(IJR).EQ.0)GO TO 12
C     S=1
      NIJ=NIJ+1
c      IF(NIJ.GT.IVER)
      nijj=max(nij,nijj)
      IAF(NIJ)=IAF(IJR)
      IBF(NIJ)=IBF(IJR)-1
      IORB(NIJ)=IORB(IJR)+1
      IF(IIM.LE.0)GO TO 12
      CALL CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
      IF(ISTOP.EQ.1)NIJ=NIJ-1
12    IF(IAF(IJR).EQ.0)GO TO 13
C     S=2
      NIJ=NIJ+1
c      IF(NIJ.GT.IVER)
      nijj=max(nij,nijj)
      IAF(NIJ)=IAF(IJR)-1
      IBF(NIJ)=IBF(IJR)+1
      IORB(NIJ)=IORB(IJR)+1
      IF(IIM.LE.0)GO TO 13
      CALL CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
      IF(ISTOP.EQ.1)NIJ=NIJ-1
13    IF(IAF(IJR).EQ.0)GO TO 14
C     S=3
      NIJ=NIJ+1
c      IF(NIJ.GT.IVER)
      nijj=max(nij,nijj)
      IAF(NIJ)=IAF(IJR)-1
      IBF(NIJ)=IBF(IJR)
      IORB(NIJ)=IORB(IJR)
      IF(IIM.LE.0)GO TO 14
      CALL CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
      IF(ISTOP.EQ.1)NIJ=NIJ-1
14    IF(IJR.EQ.IJRL)GO TO 15
      IJR=IJR+1
      GO TO 16
C     DELETE VERTICES
15    CONTINUE
      NIJ1=NIJ-1
      IN=IJS
      IUT=IJS
      IF(NIJ1.LT.IJS)GO TO 21
      DO 20 IJD=IJS,NIJ1
         JJ1=NIJ-IJD+IJS-1
         J=JJ1+1
         DO 25 K=IJS,JJ1
            IF(IAF(J).NE.IAF(K))GO TO 25
            IF(IBF(J).NE.IBF(K))GO TO 25
            GO TO 26
25       CONTINUE
         GO TO 20
26       IAF(J)=-1
         IBF(J)=-1
20    CONTINUE
C     PACK VERTICES
      IJS1=IJS+1
      DO 30 J=IJS1,NIJ
         IF(IAF(J).NE.-1)GO TO 31
         IF(IBF(J).NE.-1)GO TO 31
         IN=IN+1
         GO TO 30
31       IN=IN+1
         IUT=IUT+1
         IAF(IUT)=IAF(IN)
         IBF(IUT)=IBF(IN)
         IORB(IUT)=IORB(IN)
30    CONTINUE
C     ORDER VERTICES
      IUT1=IUT-1
      IF(IUT1.LT.IJS)GO TO 21
      DO 41 J=IJS,IUT1
         J11=J+1
         DO 42 K=J11,IUT
            IF (IAF(J)-IAF(K).LT.0) THEN
               GO TO 43
            ELSE IF (IAF(J)-IAF(K).EQ.0) THEN
               GO TO 44
            ELSE
               GO TO 42
            END IF
44          IF(IBF(J).GT.IBF(K))GO TO 42
43          IAT=IAF(J)
            IBT=IBF(J)
            IAF(J)=IAF(K)
            IBF(J)=IBF(K)
            IAF(K)=IAT
            IBF(K)=IBT
42       CONTINUE
41    CONTINUE
21    IF(II.NE.LN)IJF(LN-II)=IUT
      IJR=IJS
      IJS=IUT+1
      IJRL=IUT
      NIJ=IUT
10    CONTINUE
      JJ2=0
      DO 40 II=1,LN
      I=LN-II+1
      JJ1=IJF(I+1)+1
      JJ2=IJF(I)
      J3=JJ2+1
      IF(I.NE.1)J4=IJF(I-1)
      IF(I.EQ.1)J4=IUT
C     DETERMINE CASE DOWN
      DO 50 J=JJ1,JJ2
      IA1=IAF(J)
      IB1=IBF(J)
      K0F(J)=0
      K1F(J)=0
      K2F(J)=0
      K3F(J)=0
      DO 60 JJ=J3,J4
      IF(IA1.EQ.IAF(JJ))GO TO 61
      IF((IA1-IAF(JJ)).NE.1)GO TO 60
      IF(IB1.EQ.IBF(JJ))GO TO 62
      IF((IBF(JJ)-IB1).NE.1)GO TO 60
      K2F(J)=JJ
      GO TO 60
62    K3F(J)=JJ
      GO TO 60
61    IF(IB1.EQ.IBF(JJ))GO TO 63
      IF((IB1-IBF(JJ)).NE.1)GO TO 60
      K1F(J)=JJ
      GO TO 60
63    K0F(J)=JJ
60    CONTINUE
50    CONTINUE
40    CONTINUE
      IVF0=IUT
      IVF1=IUT-1
      IVF2=IUT-2
      IVF3=IUT-3
      K0F(IUT)=0
      K1F(IUT)=0
      K2F(IUT)=0
      K3F(IUT)=0
      K0F(IUT+1)=0
      K1F(IUT+1)=0
      K2F(IUT+1)=0
      K3F(IUT+1)=0
      IF(IPRINT.GE.5)WRITE(IW,101)
101   FORMAT(///,6X,'TAB2F',//,8X,'J',7X,'AF',2X,'BF',6X,
     *'K0F',1X,'K1F',1X,'K2F',1X,'K3F',/)
      IF(IPRINT.GE.5)WRITE(IW,100)(J,IAF(J),IBF(J),K0F(J),
     *K1F(J),K2F(J),K3F(J),J=1,IUT)
100   FORMAT(6X,I3,5X,2I4,5X,4I4)
      write(6,*) ' Number of vertices', nijj, iut
      if (nijj.gt.iver) go to 300
      RETURN
300   WRITE(IW,310)IVER
310   FORMAT(/,6X,'NUMBER OF VERTICES EXCEEDS',I7)
      CALL Abend
      END
