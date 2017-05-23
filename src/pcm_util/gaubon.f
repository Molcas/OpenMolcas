************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GauBon(MaxT,XE,YE,ZE,RE,IntSph,
     +           NV,NS,PTS,CCC,PP,AREA,IPRINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      Parameter(MxVert=20)
      Dimension XE(*),YE(*),ZE(*),RE(*),IntSph(MxVert,*)
      DIMENSION P1(3),P2(3),P3(3),U1(3),U2(3)
      DIMENSION PTS(3,MxVert),CCC(3,MxVert),PP(3)
      Save Zero, One
      DATA ZERO/0.00D0/, One/1.0d0/
      PI = 4.0D+00*ATan(1.0D+00)
C
C     Sfrutta il teorema di Gauss-Bonnet per calcolare l'area
C     della tessera con vertici PTS(3,NV). Consideriamo sempre
C     che il lato N della tessera e' quello compreso tra i vertici
C     N e N+1 (oppure NV e 1). In CCC(3,NV) sono le posizioni dei
C     centri degli archi che sottendono i vari lati della tessera.
C     La formula di Gauss-Bonet per le sfere e':
C            Area=R^2[2pi+S(Phi(N)cosT(N))-S(Beta(N))]
C     dove Phi(N) e' la lunghezza d'arco (in radianti) del lato N,
C     T(N) e' l'angolo polare del lato N, Beta(N) l'angolo esterno
C     relativo al vertice N.
C
C     Questa SUBROUTINE viene chiamata anche da VOLCHG (in calcoli con
C     dielettrico non omogeneo) per calcolare l'area delle tessere
C     proiettate sulle relative "sfere del centro di carica": in questo caso
C     i vertici della tessera sono in PTS, i centri degli archi (tutti
C     coincidenti con il centro di carica) sono in CCC e il punto
C     rappresentativo e' in PP. La cosa piu' importante e' che NS < 0.
C
      Small = 1.d-35
      TPI=2*PI
      N0 = 0 ! dummy initialize
      N2 = 0 ! dummy initialize
C
C     Calcola la prima sommatoria
      SUM1 = 0.D0
      DO 100 N = 1, NV
      X1 = PTS(1,N) - CCC(1,N)
      Y1 = PTS(2,N) - CCC(2,N)
      Z1 = PTS(3,N) - CCC(3,N)
      IF(N.LT.NV) THEN
        X2 = PTS(1,N+1) - CCC(1,N)
        Y2 = PTS(2,N+1) - CCC(2,N)
        Z2 = PTS(3,N+1) - CCC(3,N)
      ELSE
        X2 = PTS(1,1) - CCC(1,N)
        Y2 = PTS(2,1) - CCC(2,N)
        Z2 = PTS(3,1) - CCC(3,N)
      ENDIF
      DNORM1 = X1*X1 + Y1*Y1 + Z1*Z1
      DNORM2 = X2*X2 + Y2*Y2 + Z2*Z2
      SCAL = X1*X2 + Y1*Y2 + Z1*Z2
      COSPHIN = SCAL / (sqrt(DNORM1*DNORM2))
      IF(COSPHIN.GT.1.D0) COSPHIN = 1.D0
      PHIN = ACOS(COSPHIN)
C
C     Normalmente NS>0, ma se GAUBON e' stata chiamata da VOLCHG NS<0
C
      IF(NS.GT.0) THEN
C     NSFE1 e' la sfera con cui la sfera NS si interseca (eventualmente)
        NSFE1 = INTSPH(N,MaxT)
        X1 = XE(NSFE1) - XE(NS)
        Y1 = YE(NSFE1) - YE(NS)
        Z1 = ZE(NSFE1) - ZE(NS)
      ELSE
        X1 = 0
        Y1 = 0
        Z1 = 0
      ENDIF
      DNORM1 = sqrt(X1*X1 + Y1*Y1 + Z1*Z1)
      IF(DNORM1.EQ.ZERO) DNORM1 = 1.D0
      IF(NS.GT.0) THEN
        X2 = PTS(1,N) - XE(NS)
        Y2 = PTS(2,N) - YE(NS)
        Z2 = PTS(3,N) - ZE(NS)
      ELSE
        X2 = PTS(1,N) - CCC(1,1)
        Y2 = PTS(2,N) - CCC(2,1)
        Z2 = PTS(3,N) - CCC(3,1)
      ENDIF
      DNORM2 = sqrt(X2*X2 + Y2*Y2 + Z2*Z2)
      COSTN = (X1*X2+Y1*Y2+Z1*Z2)/(DNORM1*DNORM2)
      SUM1 = SUM1 + PHIN * COSTN
 100  CONTINUE
C
C     Calcola la seconda sommatoria: l'angolo esterno Beta(N) e'
C     definito usando i versori (u(N-1),u(N)) tangenti alla sfera
C     nel vertice N lungo le direzioni dei lati N-1 e N:
C                cos( Pi-Beta(N) )=u(N-1)*u(N)
C            u(N-1) = [V(N) x (V(N) x V(N-1))]/NORM
C            u(N) = [V(N) x (V(N) x V(N+1))]/NORM
C     dove V(I) e' il vettore posizione del vertice I RISPETTO AL
C     CENTRO DELL'ARCO CHE SI STA CONSIDERANDO.
C
      SUM2 = 0.D0
C     Loop sui vertici
      DO 200 N = 1, NV
      DO 2000 JJ = 1, 3
      P1(JJ) = 0.D0
      P2(JJ) = 0.D0
      P3(JJ) = 0.D0
 2000 continue
      N1 = N
      IF(N.GT.1) N0 = N - 1
      IF(N.EQ.1) N0 = NV
      IF(N.LT.NV) N2 = N + 1
      IF(N.EQ.NV) N2 = 1
C     Trova i vettori posizione rispetto ai centri corrispondenti
C     e i versori tangenti
C
C     Lato N0-N1:
      DO 2010 JJ = 1, 3
      P1(JJ) = PTS(JJ,N1) - CCC(JJ,N0)
      P2(JJ) = PTS(JJ,N0) - CCC(JJ,N0)
 2010 continue
C
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2020 JJ = 1, 3
      P2(JJ) = P3(JJ)
 2020 continue
      CALL VECP(P1,P2,P3,DNORM3)
      If(DNorm3.lt.Small) DNorm3 = One
      DO 2030 JJ = 1, 3
      U1(JJ) = P3(JJ)/DNORM3
 2030 continue
C
C     Lato N1-N2:
      DO 2040 JJ = 1, 3
      P1(JJ) = PTS(JJ,N1) - CCC(JJ,N1)
      P2(JJ) = PTS(JJ,N2) - CCC(JJ,N1)
 2040 continue
C
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2050 JJ = 1, 3
      P2(JJ) = P3(JJ)
 2050 continue
      CALL VECP(P1,P2,P3,DNORM3)
      If(DNorm3.lt.Small) DNorm3 = One
      DO 2060 JJ = 1, 3
      U2(JJ) = P3(JJ)/DNORM3
 2060 continue
C
      BETAN = ACOS(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
      SUM2 = SUM2 + (PI - BETAN)
 200  CONTINUE
C     Calcola l'area della tessera
      IF(NS.GT.0) THEN
        AREA = RE(NS)*RE(NS)*(TPI + SUM1 - SUM2)
      ELSE
        RR2 = (PTS(1,1)-CCC(1,1))**2+(PTS(2,1)-CCC(2,1))**2+
     *        (PTS(3,1)-CCC(3,1))**2
        AREA = RR2*(TPI + SUM1 - SUM2)
      ENDIF
      IF(NS.GT.0) THEN
C     Trova il punto rappresentativo (come media dei vertici)
      DO 2070 JJ = 1, 3
      PP(JJ) = 0.D0
 2070 continue
      DO 2080 I = 1, NV
      PP(1) = PP(1) + (PTS(1,I)-XE(NS))
      PP(2) = PP(2) + (PTS(2,I)-YE(NS))
      PP(3) = PP(3) + (PTS(3,I)-ZE(NS))
 2080 continue
      DNORM = 0.D0
      DO 2090 JJ = 1, 3
      DNORM = DNORM + PP(JJ)*PP(JJ)
 2090 continue
      PP(1) = XE(NS) + PP(1) * RE(NS) / sqrt(DNORM)
      PP(2) = YE(NS) + PP(2) * RE(NS) / sqrt(DNORM)
      PP(3) = ZE(NS) + PP(3) * RE(NS) / sqrt(DNORM)
      ENDIF
C
C     A causa delle approssimazioni numeriche, l'area di alcune piccole
C     tessere puo' risultare negativa, e viene in questo caso trascurata
      IF(AREA.LT.0.D0) THEN
        AREA = 0.D0
        IF(IPRINT.GE.1) WRITE(6,1000)NS
 1000   FORMAT(/,'ATTENTION: THE SURFACE OF A TESSERA IN SPHERE ',I3,
     *  ' IS NEGLECTED')
      ENDIF
      RETURN
      END
