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
      Subroutine Tessera(IPRINT,MaxT,Nesf,NS,NV,XE,YE,ZE,RE,
     +           IntSph,PTS,CCC,PP,AREA)
      IMPLICIT REAL*8 (A-H,O-Z)
      Parameter(MxVert=20)
      Dimension XE(*),YE(*),ZE(*),RE(*),IntSph(MxVert,*)

      DIMENSION CCC(3,MxVert),P1(3),P2(3),P3(3),P4(3),PP(3),
     &          PTS(3,MxVert)
      DIMENSION PSCR(3,MxVert),CCCP(3,MxVert),IND(MxVert)
      DIMENSION LTYP(MxVert),POINT(3),POINTL(3,MxVert)
      DIMENSION INTSCR(MxVert)
C     Coord. del centro che sottende l`arco tra i vertici
C     n e n+1 (per i primi tre vertici e' sicuramente il centro della
C     sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
C     appartiene alla sfera originaria INTSPH(N,MxTs)=NS)
C
      Small = 1.0d-12
      AREA = 0.D0
      DO 2000 J=1, 3
      CCC(1,J) = XE(NS)
      CCC(2,J) = YE(NS)
      CCC(3,J) = ZE(NS)
 2000 continue
C     INTSPH viene riferito alla tessera MxTs, e in seguito riceve il
C     numero corretto.
      DO 2010 N = 1, 3
        INTSPH(N,Maxt) = NS
 2010 continue
C     Loop sulle altre sfere
      DO 150 NSFE1=1,NESF
      IF(NSFE1.EQ.NS)GOTO 150
C     Memorizza i vertici e i centri che sottendono gli archi
      DO 2020 J =1, NV
      INTSCR(J) = INTSPH(J,MaxT)
      DO 2030 I = 1,3
      PSCR(I,J) = PTS(I,J)
      CCCP(I,J) = CCC(I,J)
 2030 continue
 2020 continue
      ICOP = 0
      DO 2040 J =1, MxVert
      IND(J) = 0
      LTYP(J) = 0
 2040 continue
C     Loop sui vertici della tessera considerata
      DO 100 I=1,NV
        DELR2=(PTS(1,I)-XE(NSFE1))**2+(PTS(2,I)-YE(NSFE1))**2+
     *  (PTS(3,I)-ZE(NSFE1))**2
        DELR=sqrt(DELR2)
        IF(DELR.LT.(RE(NSFE1)-Small))THEN
          IND(I) = 1
          ICOP = ICOP+1
        ENDIF
 100  CONTINUE
C     Se la tessera e' completamente coperta, la trascura
      IF(ICOP.EQ.NV) RETURN
C     Controlla e classifica i lati della tessera: LTYP = 0 (coperto),
C     1 (tagliato con il II vertice coperto), 2 (tagliato con il I
C     vertice coperto), 3 (bitagliato), 4 (libero)
C     Loop sui lati
      DO 2050 L = 1, NV
        IV1 = L
        IV2 = L+1
        IF(L.EQ.NV) IV2 = 1
        IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.1) THEN
        LTYP(L) = 0
        ELSEIF(IND(IV1).EQ.0.AND.IND(IV2).EQ.1) THEN
        LTYP(L) = 1
        ELSEIF(IND(IV1).EQ.1.AND.IND(IV2).EQ.0) THEN
        LTYP(L) = 2
        ELSEIF(IND(IV1).EQ.0.AND.IND(IV2).EQ.0) THEN
        LTYP(L) = 4
        RC2 = (CCC(1,L)-PTS(1,L))**2 + (CCC(2,L)-PTS(2,L))**2 +
     &          (CCC(3,L)-PTS(3,L))**2
        RC = sqrt(RC2)
C     Su ogni lato si definiscono 11 punti equispaziati, che vengono
C     controllati
        TOL = - 1.D-10
        DO 2060 II = 1, 11
          POINT(1) = PTS(1,IV1) + II * (PTS(1,IV2)-PTS(1,IV1)) / 11
          POINT(2) = PTS(2,IV1) + II * (PTS(2,IV2)-PTS(2,IV1)) / 11
          POINT(3) = PTS(3,IV1) + II * (PTS(3,IV2)-PTS(3,IV1)) / 11
        POINT(1) = POINT(1) - CCC(1,L)
        POINT(2) = POINT(2) - CCC(2,L)
        POINT(3) = POINT(3) - CCC(3,L)
        DNORM = sqrt(POINT(1)**2 + POINT(2)**2 + POINT(3)**2)
        POINT(1) = POINT(1) * RC / DNORM + CCC(1,L)
        POINT(2) = POINT(2) * RC / DNORM + CCC(2,L)
        POINT(3) = POINT(3) * RC / DNORM + CCC(3,L)
        DIST = sqrt( (POINT(1)-XE(NSFE1))**2 +
     &    (POINT(2)-YE(NSFE1))**2 + (POINT(3)-ZE(NSFE1))**2 )
          IF((DIST - RE(NSFE1)) .LT. TOL) THEN
c         IF(DIST.LT.RE(NSFE1))THEN
          LTYP(L) = 3
          DO 2070 JJ = 1, 3
            POINTL(JJ,L) = POINT(JJ)
 2070 continue
          GOTO 160
          ENDIF
 2060 continue
        ENDIF
 160    CONTINUE
 2050 continue
C     Se la tessera e' spezzata in due o piu' tronconi, la trascura
      ICUT = 0
      DO 2080 L = 1, NV
      IF(LTYP(L).EQ.1.OR.LTYP(L).EQ.2) ICUT = ICUT + 1
      IF(LTYP(L).EQ.3) ICUT = ICUT + 2
 2080 continue
      ICUT = ICUT / 2
      IF(ICUT.GT.1) RETURN
C     Creazione dei nuovi vertici e lati della tessera
C     Loop sui lati
      N = 1
      DO 300 L = 1, NV
C     Se il lato L e' coperto:
        IF(LTYP(L).EQ.0) GOTO 300
      IV1 = L
      IV2 = L+1
      IF(L.EQ.NV) IV2 = 1
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato L e' tagliato (con il I vertice scoperto):
      IF(LTYP(L).EQ.1) THEN
        DO 2090 JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
 2090 continue
        INTSPH(N,MaxT) = INTSCR(IV1)
        N = N+1
C     Trova l'intersezione tra i due vertici del lato L
C
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
        DO 2100 JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
 2100 continue

        CALL INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,0,
     +             IPRINT)
C     Aggiorna i vertici della tessera e il centro dell'arco
        DO 2110 JJ = 1,3
          PTS(JJ,N) = P4(JJ)
 2110 continue
C
C     Il nuovo arco sara' sotteso tra questo e il prossimo punto
C     di intersezione: il centro che lo sottende
C     sara' il centro del cerchio di intersezione tra la sfera NS
C     e la sfera NSFE1.
C
        DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+
     *        (ZE(NSFE1)-ZE(NS))**2
        CCC(1,N)=XE(NS)+(XE(NSFE1)-XE(NS))*
     *           (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.D0*DE2)
        CCC(2,N)=YE(NS)+(YE(NSFE1)-YE(NS))*
     *           (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.D0*DE2)
        CCC(3,N)=ZE(NS)+(ZE(NSFE1)-ZE(NS))*
     *           (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.D0*DE2)
        INTSPH(N,MaxT) = NSFE1
        N = N+1
      ENDIF
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato L e' tagliato (con il II vertice scoperto):
      IF(LTYP(L).EQ.2) THEN
C     Trova l'intersezione tra i due vertici del lato L
C
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
        DO 2120 JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
 2120 continue

        CALL INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,1,
     +             IPRINT)
C     Aggiorna i vertici della tessera e il centro dell'arco
        DO 2130 JJ = 1,3
          PTS(JJ,N) = P4(JJ)
          CCC(JJ,N) = CCCP(JJ,IV1)
 2130 continue
        INTSPH(N,MaxT) = INTSCR(IV1)
        N = N+1
        ENDIF
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato e' intersecato due volte:
      IF(LTYP(L).EQ.3) THEN
        DO 2140 JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
 2140 continue
        INTSPH(N,MaxT) = INTSCR(IV1)
        N = N+1
C     Trova l'intersezione tra il primo vertice e un punto intermedio
C     coperto
C
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
        DO 2150 JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = POINTL(JJ,L)
          P3(JJ) = CCCP(JJ,IV1)
 2150 continue
        CALL INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,0,
     +            IPRINT)
C     Aggiorna i vertici della tessera e il centro dell'arco
        DO 2160 JJ = 1,3
          PTS(JJ,N) = P4(JJ)
 2160 continue
C
C     Il nuovo arco sara' sotteso tra questo e il prossimo punto
C     di intersezione: il centro che lo sottende
C     sara' il centro del cerchio di intersezione tra la sfera NS
C     e la sfera NSFE1.
C
        DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+
     *        (ZE(NSFE1)-ZE(NS))**2
        CCC(1,N)=XE(NS)+(XE(NSFE1)-XE(NS))*
     *           (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.D0*DE2)
        CCC(2,N)=YE(NS)+(YE(NSFE1)-YE(NS))*
     *           (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.D0*DE2)
        CCC(3,N)=ZE(NS)+(ZE(NSFE1)-ZE(NS))*
     *           (RE(NS)**2-RE(NSFE1)**2+DE2)/(2.D0*DE2)
        INTSPH(N,MaxT) = NSFE1
        N = N+1
C
C     Trova l'intersezione tra un punto intermedio coperto e il
C     secondo vertice
C
C     P1 = coord. del primo vertice
C     P2 = coord. del secondo vertice
C     P3 = coord. del centro dell`arco sotteso
C     P4 = coord. dell'intersezione
C
        DO 2170 JJ = 1, 3
          P1(JJ) = POINTL(JJ,L)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
 2170 continue

        CALL INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,1,
     +           IPRINT)
C     Aggiorna il vertice e il centro dell'arco
        DO 2180 JJ = 1,3
          PTS(JJ,N) = P4(JJ)
          CCC(JJ,N) = CCCP(JJ,IV1)
 2180 continue
      INTSPH(N,MaxT) = INTSCR(IV1)
      N = N + 1
        ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Se il lato e' scoperto:
      IF(LTYP(L).EQ.4) THEN
        DO 2190 JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
 2190 continue
        INTSPH(N,MaxT) = INTSCR(IV1)
        N = N+1
        ENDIF
C     Controlla che il numero di vertici creati non sia eccessivo
      If(N.GT.11) then
        Write(6,'(/,'' TESSERA: too many vertices in a tessera'')')
        Call Abend()
      EndIf
 300  CONTINUE
C
      NV = N - 1
 150  CONTINUE
C
C     Se la tessera non e' stata scartata, a questo punto ne troviamo
C     l'area e il punto rappresentativo
C
      CALL GAUBON(MaxT,XE,YE,ZE,RE,IntSph,NV,NS,
     +            PTS,CCC,PP,AREA,IPRINT)
      RETURN
      END
