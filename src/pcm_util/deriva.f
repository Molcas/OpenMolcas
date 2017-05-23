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
      Subroutine Deriva(IPrint,ToAng,NAt,nTs,NEsf,NEsfP,RSolv,Tessera,
     $                  Vert,Centr,Sphere,ISphe,IntSph,NOrd,NVert,
     $                  NewSph,DerTes,DerPunt,DerRad,DerCentr)
      Implicit Real*8 (a-h,o-z)
      Parameter (MxVert=20)
      Integer AlGe(63), Casca(MxVert)
      Dimension DISTKI(3), DERP(3), DERPT(3)
      Dimension Tessera(4,*),Vert(3,MxVert,*),Centr(3,MxVert,*)
      Dimension Sphere(4,*),ISphe(*)
      Dimension IntSph(MxVert,*),NewSph(2,NEsf),NOrd(*),NVert(*)
      Dimension DerTes(nTs,NAt,3),DerPunt(nTs,NAt,3,3)
      Dimension DerRad(NEsf,NAt,3),DerCentr(NEsf,NAt,3,3)
      Save Zero, One
      DATA ZERO, One/0.0D0,1.0d0/
*
*     Compute the derivatives of area and of representative point position
*     for each tessera
*
C     Le derivate contengono termini dovuti direttamente allo spostamento
C     del centro della sfera NSJ, e termini "mediati" dagli
C     spostamenti del centro e dal cambiamento del raggio delle sfere
C     "aggiunte" (create da PEDRA, oltre a quelle originarie).
C     Chiamiamo ITS la tessera, I la sfera a cui ITS appartiene, J la sfera
C     che si muove. Le derivate che cerchiamo sono dS(ITS)/dx(J) e
C     dQ(ITS)/dx(J) (indicando con x(J) una qualsiasi coordinata di J.
C     Per prima cosa calcoliamo le derivate dell'area e del punto
C     rappresentativo delle tessere DIRETTAMENTE tagliate da J, poi di quelle
C     che appartengono a J e sono tagliate da altre sfere.
C     Poi consideriamo le "sfere aggiunte" che vengono modificate dal
C     movimento di J: per ogni sfera aggiunta si esaminano tutte le
C     tessere e caso per caso si calcolano le derivate.
C
C     In ogni caso le derivate totali dS/dx(J) e dQ/dx(J) si scrivono come una
C     combinazione delle seguenti derivate parziali:
C          dS/dR(I)  ;  dS/dz(I)  ;  dS/dR(K)  ;  dS/dz(K)
C          dQ/dR(I)  ;  dQ/dz(I)  ;  dQ/dR(K)  ;  dQ/dz(K)
C     dove, ancora una volta, z indica una qualunque coordinata e K un'altra
C     sfera (la stessa J o una sfera aggiunta modificata dal movimento
C     di J) che taglia S.
C
C     NESFJ e' la sfera che sta attorno all'atomo NSJ: se NSJ non ha
C     nessuna sfera, NESFJ = 0
*
*     Initiate the derivative arrays
      NDeg = 3 * NAt
      Call FZero(DerTes(1,1,1),nTs*NDeg)
      Call FZero(DerPunt(1,1,1,1),3*nTs*NDeg)
      Call FZero(DerRad(1,1,1),NEsf*NDeg)
      Call FZero(DerCentr(1,1,1,1),3*NEsf*NDeg)
*     The geometric quantities are expected to be in Angstrom
      Call DScal_(4*nTs,ToAng,Tessera(1,1),1)
      Call DScal_(  nTs,ToAng,Tessera(4,1),4)
      Call DScal_(4*NEsf,ToAng,Sphere(1,1),1)
      Call DScal_(3*MxVert*nTs,ToAng,Vert(1,1,1),1)
      Call DScal_(3*MxVert*nTs,ToAng,Centr(1,1,1),1)
*
      Do 1000 NSJ = 1, NAt
      Do 1010 ICoord = 1, 3
*
      NESFJ = 0
      DO 2000 I = 1, NESFP
        IF( NSJ.EQ.NORD(I) ) NESFJ = I
 2000 continue
C
C     Memorizza in DERRAD(NS,NSJ,ICOORD) la derivata del raggio di ogni
C     sfera NS e in DERCENTR(NS,NSJ,ICOORD,3) le derivate delle
C     coordinate del centro di NS rispetto alla coord. ICOORD della
C     sfera NSJ.
C
C     Se NS e' una sfera originaria queste derivate sono 0, tranne
C     DERCENTR(NESFJ,NSJ,ICOORD,ICOORD)=1 :
C
      DO 2010 NSFE = 1, NESFP
        DERRAD(NSFE,NSJ,ICOORD) = ZERO
        DO 2020 JJ = 1, 3
          DERCENTR(NSFE,NSJ,ICOORD,JJ) = ZERO
 2020 continue
 2010 continue
      IF( NESFJ.NE.0 )
     $ DERCENTR(NESFJ,NSJ,ICOORD,ICOORD) = 1.0D0
C
C     Se non c'e' nessuna sfera attorno all'atomo NSJ ...
      IF(NESFJ.EQ.0) Goto 1000
C
C     1) Effetti diretti.
C     Loop sulle tessere.
      DO 100 ITS= 1, NTS
      DERTS = ZERO
      DO 2030 JJ = 1, 3
        DERPT(JJ) = ZERO
 2030 continue
      NV = NVERT(ITS)
      NSI = ISPHE(ITS)
      IF(NSI.EQ.NESFJ) GOTO 200
C
C     Derivate nel caso I non = J
C
C               dS/dx(J), dQ/dx(J)
C     ITS ha un lato su J ?
      ICONT = 0
      DO 2040 N = 1, NV
        IF(INTSPH(N,ITs).EQ.NESFJ) ICONT = ICONT + 1
 2040 continue
      IF(ICONT.GE.1) THEN
      CALL dSd(0,ITS,ICOORD,NESFJ,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
      DERTS = DERTS + DERS
      DO 2050 JJ = 1, 3
        DERPT(JJ) = DERPT(JJ) + DERP(JJ)
 2050 continue
      ENDIF
      GOTO 150
 200  CONTINUE
C     Derivate nel caso I = J
C
C     Loop sulle sfere K che intersecano I
      DO 300 NSK = 1, NESF
      IF(NSK.EQ.NSI) GOTO 300
C     ITS ha un lato su K ?
      ICONT = 0
      DO 2060 N = 1, NV
      IF(INTSPH(N,ITs).EQ.NSK) ICONT = ICONT + 1
 2060 continue
      IF(ICONT.EQ.0) GOTO 300
      CALL dSd(0,ITS,ICOORD,NSK,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
      DERTS = DERTS - DERS
      DO 2070 JJ = 1, 3
        DERPT(JJ) = DERPT(JJ) - DERP(JJ)
 2070 continue
 300  CONTINUE
 150  DERTES(ITS,NSJ,ICOORD) = DERTS
      DO 2080 JJ = 1, 3
      DERPUNT(ITS,NSJ,ICOORD,JJ) = DERPT(JJ)
 2080 continue
 100  CONTINUE
C
C     2) Effetti indiretti.
C     Loop sulle sfere aggiunte
      DO 500 NSA = NESFP+1, NESF
      DO 2090 II = 1, 63
      ALGE(II) = 0
 2090 continue
C     Costruiamo l'"albero genealogico" della sfera NSA
      ALGE(1) = NSA
      ALGE(2) = ABS(NEWSPH(1,NSA))
      ALGE(3) = ABS(NEWSPH(2,NSA))
      LIVEL = 3
      NUMBER = 2
 510  NSUB = 1
      DO 2100 II = LIVEL-NUMBER+1, LIVEL
      IF(ALGE(II).GT.NESFP) THEN
        ALGE(LIVEL+NSUB)   = ABS(NEWSPH(1,ALGE(II)))
        ALGE(LIVEL+NSUB+1) = ABS(NEWSPH(2,ALGE(II)))
        ENDIF
        NSUB = NSUB + 2
 2100 continue
      NUMBER = NUMBER * 2
      LIVEL = LIVEL + NUMBER
      IF(NUMBER.LT.32) GOTO 510
C     Si accerta che nell'ultimo livello ci siano solo sfere originarie
      DO 2110 II = 34, 63
      IF(ALGE(II).GT.NESFP) THEN
        WRITE(6,10) NSA, ALGE(II)
 10       FORMAT(/,'SFERA ',I3,/,
     %    'AL LIVELLO 5 LA SFERA',I3,' E" AGGIUNTA')
      ENDIF
 2110 continue
C     Quando un elemento di ALGE e' = NESFJ, costruisce la corrispondente
C     "cascata" di sfere aggiunte che collega NESFJ a NSA
      DO 600 LIVEL = 2, 6
      MIN = 2**(LIVEL-1)
      MAX = 2**(LIVEL) - 1
      DO 700 II = MIN, MAX
      IF(ALGE(II).NE.NESFJ) GOTO 700
      DO 2120 K = 1, MxVert
      CASCA(K) = 0
 2120 continue
      CASCA(1) = NESFJ
      INDEX = II
      K = 2
      DO 2130 LL = LIVEL, 2, -1
        FACT = (INDEX - 2**(LL-1)) / 2.D0
        INDEX = 2**(LL-2) + iRToInt(FACT)
        CASCA(K) = ALGE(INDEX)
        K = K + 1
 2130 continue
C     Contiamo gli elementi diversi da 0 in CASCA
      ICONT = 0
      DO 2140 K = 1, MxVert
      IF(CASCA(K).NE.0) ICONT = ICONT + 1
 2140 continue
C
C     Costruiamo le derivate composte del raggio e delle coordinate di
C     NSA (ultimo elemento di CASCA)
C     rispetto alla coordinata ICOORD di NESFJ (primo elemento di CASCA)
      NS1 = CASCA(1)
      NS2 = CASCA(2)
      CALL dRdC(NS2,ICOORD,NS1,DR,RSolv,Sphere,NewSph)
      CALL dCdC(1,NS2,ICOORD,NS1,DX,Sphere,NewSph)
      CALL dCdC(2,NS2,ICOORD,NS1,DY,Sphere,NewSph)
      CALL dCdC(3,NS2,ICOORD,NS1,DZ,Sphere,NewSph)
      DO 800 I = 3, ICONT
      DDR = ZERO
      DDX = ZERO
      DDY = ZERO
      DDZ = ZERO
      NS1 = CASCA(I-1)
      NS2 = CASCA(I)
C
C     Derivata del raggio dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
      CALL dRdR(NS2,NS1,DER,RSolv,Sphere,NewSph)
      DDR = DER * DR
      CALL dRdC(NS2,1,NS1,DER,RSolv,Sphere,NewSph)
      DDR = DDR + DER * DX
      CALL dRdC(NS2,2,NS1,DER,RSolv,Sphere,NewSph)
      DDR = DDR + DER * DY
      CALL dRdC(NS2,3,NS1,DER,RSolv,Sphere,NewSph)
      DDR = DDR + DER * DZ
C
C     Derivata della coord. X dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
      CALL dCdR(1,NS2,NS1,DER,Sphere,NewSph)
      DDX = DER * DR
      CALL dCdC(1,NS2,1,NS1,DER,Sphere,NewSph)
      DDX = DDX + DER * DX
      CALL dCdC(1,NS2,2,NS1,DER,Sphere,NewSph)
      DDX = DDX + DER * DY
      CALL dCdC(1,NS2,3,NS1,DER,Sphere,NewSph)
      DDX = DDX + DER * DZ
C
C     Derivata della coord. Y dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
      CALL dCdR(2,NS2,NS1,DER,Sphere,NewSph)
      DDY = DER * DR
      CALL dCdC(2,NS2,1,NS1,DER,Sphere,NewSph)
      DDY = DDY + DER * DX
      CALL dCdC(2,NS2,2,NS1,DER,Sphere,NewSph)
      DDY = DDY + DER * DY
      CALL dCdC(2,NS2,3,NS1,DER,Sphere,NewSph)
      DDY = DDY + DER * DZ
C
C     Derivata della coord. Z dell'elemento I di CASCA rispetto
C     alla coord. ICOORD dell'elemento 1 di CASCA
      CALL dCdR(3,NS2,NS1,DER,Sphere,NewSph)
      DDZ = DER * DR
      CALL dCdC(3,NS2,1,NS1,DER,Sphere,NewSph)
      DDZ = DDZ + DER * DX
      CALL dCdC(3,NS2,2,NS1,DER,Sphere,NewSph)
      DDZ = DDZ + DER * DY
      CALL dCdC(3,NS2,3,NS1,DER,Sphere,NewSph)
      DDZ = DDZ + DER * DZ
C
      DR = DDR
      DX = DDX
      DY = DDY
      DZ = DDZ
 800  CONTINUE
C
C     Se NS e' una sfera aggiunta, memorizza le derivate del raggio
C     e delle coordinate del centro:
      DERRAD(NSA,NSJ,ICOORD) = DR
      DERCENTR(NSA,NSJ,ICOORD,1) = DX
      DERCENTR(NSA,NSJ,ICOORD,2) = DY
      DERCENTR(NSA,NSJ,ICOORD,3) = DZ
C
C     Per ogni sfera aggiunta connessa a NESFJ, fa passare tutte le
C     tessere, calcolando, se e' il caso, le derivate
      DO 900 ITS = 1, NTS
      DERTS = ZERO
      DO 2150 JJ = 1, 3
      DERPT(JJ) = ZERO
 2150 continue
      NV = NVERT(ITS)
      NSI = ISPHE(ITS)
C     Derivate delle tessere appartenenti a NSA
      IF(NSI.EQ.NSA) THEN
C       Derivate relative al raggio di NSI: 1) scaling
        DERTS = 2.D0 * Tessera(4,ITS) * DR / Sphere(4,NSI)
      DERPT(1) = (Tessera(1,ITS) - Sphere(1,NSI)) * DR / Sphere(4,NSI)
      DERPT(2) = (Tessera(2,ITS) - Sphere(2,NSI)) * DR / Sphere(4,NSI)
      DERPT(3) = (Tessera(3,ITS) - Sphere(3,NSI)) * DR / Sphere(4,NSI)
C       Loop sulle altre sfere K che tagliano ITS
        DO 950 NSK = 1, NESF
          IF(NSK.EQ.NSI) GOTO 950
          ICONT = 0
          DO 2160 N = 1, NV
          IF(INTSPH(N,ITs).EQ.NSK) ICONT = ICONT + 1
 2160 continue
          IF(ICONT.EQ.0) GOTO 950
C         Derivate relative al raggio di NSI: 2) raggio delle altre sfere
          CALL dSd(1,ITS,IJunk,NSK,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
          DERTS = DERTS - DERS * Sphere(4,NSK) * DR / Sphere(4,NSI)
        DO 2170 JJ = 1, 3
          DERPT(JJ) = DERPT(JJ) -
     $                DERP(JJ) * Sphere(4,NSK) * DR / Sphere(4,NSI)
 2170 continue
C         Derivate relative al raggio di NSI: 3) coord. delle altre sfere
          DISTKI(1) = Sphere(1,NSK) - Sphere(1,NSI)
          DISTKI(2) = Sphere(2,NSK) - Sphere(2,NSI)
          DISTKI(3) = Sphere(3,NSK) - Sphere(3,NSI)
          FAC = DR / Sphere(4,NSI)
          DO 2180 JJ = 1, 3
            CALL dSd(0,ITS,JJ,NSK,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
                   IF(DISTKI(JJ).NE.ZERO) THEN
            DERTS = DERTS - FAC * DERS * DISTKI(JJ)
            DO 2190 LL = 1, 3
              DERPT(LL) = DERPT(LL) - FAC * DERP(LL) * DISTKI(JJ)
 2190 continue
          ENDIF
 2180 continue
C         Derivate relative alle coordinate di NSI
          CALL dSd(0,ITS,1,NSK,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
          DERTS = DERTS - DX * DERS
        DO 2200 JJ = 1, 3
          DERPT(JJ) = DERPT(JJ) - DX * DERP(JJ)
 2200 continue
          CALL dSd(0,ITS,2,NSK,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
          DERTS = DERTS - DY * DERS
        DO 2210 JJ = 1, 3
          DERPT(JJ) = DERPT(JJ) - DY * DERP(JJ)
 2210 continue
          CALL dSd(0,ITS,3,NSK,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
          DERTS = DERTS - DZ * DERS
        DO 2220 JJ = 1, 3
          DERPT(JJ) = DERPT(JJ) - DZ * DERP(JJ)
 2220 continue
 950    CONTINUE
        DERTES(ITS,NSJ,ICOORD) = DERTES(ITS,NSJ,ICOORD) + DERTS
      DO 2230 JJ = 1, 3
        DERPUNT(ITS,NSJ,ICOORD,JJ) = DERPUNT(ITS,NSJ,ICOORD,JJ) +
     *    DERPT(JJ)
 2230 continue
      ELSE
C       Derivate delle tessere tagliate da NSA
        ICONT = 0
        DO 2240 N = 1, NV
          IF(INTSPH(N,ITs).EQ.NSA) ICONT = ICONT + 1
 2240 continue
        IF(ICONT.EQ.0) GOTO 900
C       Derivate relative al raggio di NSA
      CALL dSd(1,ITS,IJunk,NSA,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
      DERTS = DERS * DR
      DO 2250 JJ = 1, 3
        DERPT(JJ) = DERP(JJ) * DR
 2250 continue
C       Derivate relative alle coordinate di NSA
        CALL dSd(0,ITS,1,NSA,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
        DERTS = DERTS + DX * DERS
      DO 2260 JJ = 1, 3
        DERPT(JJ) = DERPT(JJ) + DX * DERP(JJ)
 2260 continue
        CALL dSd(0,ITS,2,NSA,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
        DERTS = DERTS + DY * DERS
      DO 2270 JJ = 1, 3
        DERPT(JJ) = DERPT(JJ) + DY * DERP(JJ)
 2270 continue
        CALL dSd(0,ITS,3,NSA,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
        DERTS = DERTS + DZ * DERS
      DO 2280 JJ = 1, 3
        DERPT(JJ) = DERPT(JJ) + DZ * DERP(JJ)
 2280 continue
C
        DERTES(ITS,NSJ,ICOORD) = DERTES(ITS,NSJ,ICOORD) + DERTS
      DO 2290 JJ = 1, 3
        DERPUNT(ITS,NSJ,ICOORD,JJ) = DERPUNT(ITS,NSJ,ICOORD,JJ) +
     *    DERPT(JJ)
 2290 continue
      ENDIF
 900  CONTINUE
 700  CONTINUE
 600  CONTINUE
 500  CONTINUE
C
C     Controlla che nessuna derivata sia eccessiva (per imprecisioni numeriche)
      DO 2300 ITS = 1, NTS
      IF( abs(DERTES(ITS,NSJ,ICOORD)) .GT. 10.D0 ) THEN
        IF(IPrint.ge.2) THEN
          WRITE(6,*)'DERIVATA GEOMETRICA ECCESSIVA TRASCURATA'
          WRITE(6,*)'TESSERA, SFERA, COORDINATA',ITS,NESFJ,ICOORD
        ENDIF
        DERTES(ITS,NSJ,ICOORD) = 0.D0
        ENDIF
 2300 continue
*
*     End of loop on atoms and coordinates
*
 1010 Continue
 1000 Continue
*     Put the geometric quantities in Bohr again
      ToBohr = One / ToAng
      Call DScal_(4*nTs,ToBohr,Tessera(1,1),1)
      Call DScal_(  nTs,ToBohr,Tessera(4,1),4)
      Call DScal_(4*NEsf,ToBohr,Sphere(1,1),1)
      Call DScal_(30*nTs,ToBohr,Vert(1,1,1),1)
      Call DScal_(30*nTs,ToBohr,Centr(1,1,1),1)
*
      RETURN
      END
*
      SUBROUTINE DCDC(JJ,NSI,ICOORD,NESFJ,DC,Sphere,NewSph)
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension Sphere(4,*),NewSph(2,*)
      DIMENSION COORDJ(3), COORDK(3)
C
C     Trova la derivata della coordinata JJ del centro della sfera
C     NSI rispetto alla coordinata ICOORD di NSJ, che interseca NSI.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
C     dipende dalle due sfere "precedenti" NESFJ e NSK
C
C     Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(1,NSI).LT.0.OR.NEWSPH(2,NSI).LT.0) GOTO 100
      K = NEWSPH(1,NSI)
      IF(K.EQ.NESFJ) K = NEWSPH(2,NSI)
      COORDJ(1) = Sphere(1,NESFJ)
      COORDJ(2) = Sphere(2,NESFJ)
      COORDJ(3) = Sphere(3,NESFJ)
      COORDK(1) = Sphere(1,K)
      COORDK(2) = Sphere(2,K)
      COORDK(3) = Sphere(3,K)
      D2 = (Sphere(1,NESFJ)-Sphere(1,K))**2
     $   + (Sphere(2,NESFJ)-Sphere(2,K))**2
     $   + (Sphere(3,NESFJ)-Sphere(3,K))**2
      D = sqrt(D2)
      DC = (Sphere(4,NESFJ) - Sphere(4,K))
     $   * (COORDJ(ICOORD) - COORDK(ICOORD))
     $   * (COORDJ(JJ) - COORDK(JJ)) / (2.D0 * D**3)
      IF(JJ.EQ.ICOORD)DC = DC + 0.5D0 - (Sphere(4,NESFJ)-Sphere(4,K))
     $                   / (2.D0*D)
      GOTO 200
C
 100  CONTINUE
      NSK = NEWSPH(1,NSI)
      IF(ABS(NSK).EQ.NESFJ) NSK = NEWSPH(2,NSI)
      COORDJ(1) = Sphere(1,NESFJ)
      COORDJ(2) = Sphere(2,NESFJ)
      COORDJ(3) = Sphere(3,NESFJ)
      COORDK(1) = Sphere(1,ABS(NSK))
      COORDK(2) = Sphere(2,ABS(NSK))
      COORDK(3) = Sphere(3,ABS(NSK))
      D2 = (COORDJ(1)-COORDK(1))**2 + (COORDJ(2)-COORDK(2))**2 +
     &     (COORDJ(3)-COORDK(3))**2
      D = sqrt(D2)
      IF(NSK.GT.0) THEN
      DC = Sphere(4,NESFJ) * (COORDJ(JJ) - COORDK(JJ))
     $   * (COORDJ(ICOORD) - COORDK(ICOORD)) / D**3
        IF(ICOORD.EQ.JJ) DC = DC + 1.D0 - Sphere(4,NESFJ) / D
      ELSE
      DC = - Sphere(4,ABS(NSK)) * (COORDK(JJ) - COORDJ(JJ))
     $     * (COORDK(ICOORD) - COORDJ(ICOORD)) / D**3
        IF(ICOORD.EQ.JJ) DC = DC + Sphere(4,ABS(NSK)) / D
      ENDIF
C
 200  CONTINUE
      RETURN
      END
*
      SUBROUTINE dCdR(JJ,NSI,NESFJ,DC,Sphere,NewSph)
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension Sphere(4,*),NewSph(2,*)
      DIMENSION COORDJ(3), COORDK(3)
      DATA ZERO/0.0D0/
C
C     Trova la derivata della coordinata JJ del centro della sfera
C     NSI rispetto al raggio dellla sfera NSJ.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
C     dipende dalle due sfere "precedenti" NESFJ e NSK
C
C     Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(1,NSI).LT.0.OR.NEWSPH(2,NSI).LT.0) GOTO 100
      NSK = NEWSPH(1,NSI)
      IF(NSK.EQ.NESFJ) NSK = NEWSPH(2,NSI)
      COORDJ(1) = Sphere(1,NESFJ)
      COORDJ(2) = Sphere(2,NESFJ)
      COORDJ(3) = Sphere(3,NESFJ)
      COORDK(1) = Sphere(1,NSK)
      COORDK(2) = Sphere(2,NSK)
      COORDK(3) = Sphere(3,NSK)
      D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2
     $   + (Sphere(2,NESFJ)-Sphere(2,NSK))**2
     $   + (Sphere(3,NESFJ)-Sphere(3,NSK))**2
      D = sqrt(D2)
      DC = - (COORDJ(JJ) - COORDK(JJ)) / (2.D0 * D)
      GOTO 200
C
 100  CONTINUE
      NSK = NEWSPH(1,NSI)
      IF(ABS(NSK).EQ.NESFJ) NSK = NEWSPH(2,NSI)
      DC = ZERO
      IF(NSK.LT.ZERO) GOTO 200
      COORDJ(1) = Sphere(1,NESFJ)
      COORDJ(2) = Sphere(2,NESFJ)
      COORDJ(3) = Sphere(3,NESFJ)
      COORDK(1) = Sphere(1,NSK)
      COORDK(2) = Sphere(2,NSK)
      COORDK(3) = Sphere(3,NSK)
      D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2
     $   + (Sphere(2,NESFJ)-Sphere(2,NSK))**2
     $   + (Sphere(3,NESFJ)-Sphere(3,NSK))**2
      D = sqrt(D2)
      DC = - ( COORDJ(JJ) - COORDK(JJ) ) / D
C
 200  CONTINUE
      RETURN
      END
*
      SUBROUTINE dRdC(NSI,ICOORD,NESFJ,DR,RSolv,Sphere,NewSph)
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension Sphere(4,*),NewSph(2,*)
      DIMENSION COORDJ(3), COORDK(3)
C
C     Trova la derivata del raggio della sfera NSI rispetto alla
C     coordinata ICOORD (1=X, 2=Y, 3=Z) della sfera NSJ, che interseca
C     NSI.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
C     dipende dalle due sfere "precedenti" NESFJ e K
C
C     Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(1,NSI).LT.0.OR.NEWSPH(2,NSI).LT.0) GOTO 100
      K = NEWSPH(1,NSI)
      IF(K.EQ.NESFJ) K = NEWSPH(2,NSI)
      COORDJ(1) = Sphere(1,NESFJ)
      COORDJ(2) = Sphere(2,NESFJ)
      COORDJ(3) = Sphere(3,NESFJ)
      COORDK(1) = Sphere(1,K)
      COORDK(2) = Sphere(2,K)
      COORDK(3) = Sphere(3,K)
      D2 = (Sphere(1,NESFJ)-Sphere(1,K))**2
     $   + (Sphere(2,NESFJ)-Sphere(2,K))**2
     $   + (Sphere(3,NESFJ)-Sphere(3,K))**2
      D = sqrt(D2)
      B = 0.5D0 * (D + Sphere(4,NESFJ) - Sphere(4,K))
      RS = RSOLV
      A = ((Sphere(4,NESFJ)+RS)**2 + D2 - (Sphere(4,K)+RS)**2) / D
      DR = (2.D0*A*B - 2.D0*B*D - A*D)
     $   * (COORDJ(ICOORD)-COORDK(ICOORD))
     $   / (4.D0*D2*(Sphere(4,NSI)+RS))
      GOTO 200
C
 100  CONTINUE
      NSK = NEWSPH(1,NSI)
      IF(ABS(NSK).EQ.NESFJ) NSK = NEWSPH(2,NSI)
      COORDJ(1) = Sphere(1,NESFJ)
      COORDJ(2) = Sphere(2,NESFJ)
      COORDJ(3) = Sphere(3,NESFJ)
      COORDK(1) = Sphere(1,ABS(NSK))
      COORDK(2) = Sphere(2,ABS(NSK))
      COORDK(3) = Sphere(3,ABS(NSK))
      RI = Sphere(4,NSI) + RSOLV
      RJ = Sphere(4,NESFJ) + RSOLV
      RK = Sphere(4,ABS(NSK)) + RSOLV
      DIFF = COORDJ(ICOORD) - COORDK(ICOORD)
      D2 = (COORDJ(1)-COORDK(1))**2 + (COORDJ(2)-COORDK(2))**2 +
     &     (COORDJ(3)-COORDK(3))**2
      D = sqrt(D2)
      FAC = Sphere(4,NESFJ) * ( RJ*RJ - D*D - RK*RK )
      IF(NSK.LT.0) FAC = Sphere(4,ABS(NSK)) * (RK*RK - D*D - RJ*RJ )
      DR = DIFF * FAC / ( 2.D0 * D**3 * RI )
C
 200  CONTINUE
      RETURN
      END
*
      SUBROUTINE dRdR(NSI,NESFJ,DR,RSolv,Sphere,NewSph)
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension Sphere(4,*),NewSph(2,*)
C
C     Trova la derivata del raggio della sfera NSI rispetto al raggio
C     della sfera NSJ.
C
C     La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
C     dipende dalle due sfere "precedenti" NESFJ e NSK
C     Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
C     e la generatrice "principale" corrisponde al label negativo
C     (cfr. JCC 11, 1047 (1990))
C
      IF(NEWSPH(1,NSI).LT.0.OR.NEWSPH(2,NSI).LT.0) GOTO 100
      NSK = NEWSPH(1,NSI)
      IF(NSK.EQ.NESFJ) NSK = NEWSPH(2,NSI)
      RS = RSOLV
      RJ = Sphere(4,NESFJ) + RS
      RK = Sphere(4,NSK) + RS
      RI = Sphere(4,NSI) + RS
      D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2
     $   + (Sphere(2,NESFJ)-Sphere(2,NSK))**2
     $   + (Sphere(3,NESFJ)-Sphere(3,NSK))**2
      D = sqrt(D2)
      DR = (-3.D0*RJ*RJ + RK*RK + 2.D0*RJ*RK + 3.D0*D*RJ - D*RK) /
     %      (4.D0*D*RI)
      GOTO 200
C
 100  CONTINUE
      NSK = NEWSPH(1,NSI)
      IF(ABS(NSK).EQ.NESFJ) NSK = NEWSPH(2,NSI)
C
      IF(NSK.GT.0) THEN
        RS = RSOLV
        RJ = Sphere(4,NESFJ) + RS
        RK = Sphere(4,NSK) + RS
        RI = Sphere(4,NSI) + RS
        D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2
     $     + (Sphere(2,NESFJ)-Sphere(2,NSK))**2
     $     + (Sphere(3,NESFJ)-Sphere(3,NSK))**2
        D = sqrt(D2)
        DR = ( 2.D0*D*RJ + 2.D0*D*Sphere(4,NESFJ)
     $     - 2.D0*RJ*Sphere(4,NESFJ)
     $     + D*D - RJ*RJ - RK*RK ) / (2.D0*D*RI)
      ELSE
        RS = RSOLV
        RJ = Sphere(4,NESFJ) + RS
        RI = Sphere(4,NSI) + RS
        D2 = (Sphere(1,NESFJ)-Sphere(1,ABS(NSK)))**2 +
     %       (Sphere(2,NESFJ)-Sphere(2,ABS(NSK)))**2 +
     %       (Sphere(3,NESFJ)-Sphere(3,ABS(NSK)))**2
        D = sqrt(D2)
        DR = ( Sphere(4,ABS(NSK)) * RJ ) / ( D*RI)
      ENDIF
 200  CONTINUE
      RETURN
      END
*
      SUBROUTINE dSd(IOpt,ITS,ICOORD,NESFJ,DERS,DERP,Tessera,
     $         Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
      IMPLICIT REAL*8 (A-H,O-Z)
      Parameter (MxVert=20)
      Dimension Sphere(4,*),ISphe(*),NewSph(2,*)
      Dimension IntSph(MxVert,*),NVert(*)
      Dimension Tessera(4,*),Vert(3,MxVert,*),Centr(3,MxVert,*)
      DIMENSION DP(MxVert,3), DERP(3), SUMDER(3), SUMVERT(3),  PT(3)
      DATA ZERO/0.0D0/
C
C     Trova le derivate dell'area e del punto rappresentativo della
C     tessera ITS rispetto
C
C       IOpt = 0 : alla coordinata ICOORD(1=X, 2=Y, 3=Z) dell'atomo NSJ.
C                  (ex dsdx)
C       IOpt = 1 : al raggio della sfera NSJ. (ex dsdr)
C
C     NESFJ e' la sfera attorno all'atomo NSJ.
C     ITS appartiene alla sfera NSI (diversa da NESFJ, e intersecata
C     da NESFJ).
C
C
      DERTS = ZERO
      NSI = ISPHE(ITS)
      NV = NVERT(ITS)
      DO 2000 L = 1, NV
      DP(L,1) = ZERO
      DP(L,2) = ZERO
      DP(L,3) = ZERO
 2000 continue
      DO 100 L = 1, NV
      IF(INTSPH(L,ITs).NE.NESFJ) GOTO 100
C     L1 e' il lato di ITS che sta sulla sfera NESFJ
      L1 = L
      L0 = L - 1
      IF(L1.EQ.1) L0 = NV
      L2 = L1 + 1
      IF(L1.EQ.NV) L2 = 1
      L3 = L2 + 1
      IF(L2.EQ.NV) L3 = 1
C     Trova la derivata del vertice L1: DP(L1)
      CALL DVer(IOpt,ICOORD,ITS,L0,L1,L2,DP(L1,1),DP(L1,2),DP(L1,3),
     $          Vert,Centr,nTs,Sphere,IntSph)
C     Trova la derivata del vertice L2: DP(L2)
      CALL DVer(IOpt,ICOORD,ITS,L1,-L2,L3,DP(L2,1),DP(L2,2),DP(L2,3),
     $          Vert,Centr,nTs,Sphere,IntSph)
C
C     1) Derivata dell'area della tessera.
C
C     Calcola i contributi dovuti ai lati L0-L1, L1-L2, L2-L3
      IF(INTSPH(L0,ITs).NE.NSI) THEN
        CALL DerPhi(IOpt,ICOORD,NESFJ,ITS,L0,L1,DP,DA,
     $              Vert,Centr,nTs,Sphere,IntSph,ISphe)
        DERTS = DERTS + DA
      ENDIF
      CALL DerPhi(IOpt,ICOORD,NESFJ,ITS,L1,L2,DP,DA,
     $              Vert,Centr,nTs,Sphere,IntSph,ISphe)
      DERTS = DERTS + DA
      IF(INTSPH(L2,ITs).NE.ISPHE(ITS)) THEN
        CALL DerPhi(IOpt,ICOORD,NESFJ,ITS,L2,L3,DP,DA,
     $              Vert,Centr,nTs,Sphere,IntSph,ISphe)
        DERTS = DERTS + DA
      ENDIF
C
C     Calcola i contributi dovuti ai vertici L1 e L2
      CALL DerBet(IOpt,ICOORD,NESFJ,ITS,L0,L1,L2,DP,DA,
     $            Vert,Centr,nTs,Sphere,IntSph,ISphe)
      DERTS = DERTS - DA
      CALL DerBet(IOpt,ICOORD,NESFJ,ITS,L1,L2,L3,DP,DA,
     $            Vert,Centr,nTs,Sphere,IntSph,ISphe)
      DERTS = DERTS - DA
 100  CONTINUE
      DERS = DERTS
C
C     2) Derivata del punto rappresentativo.
C
C     Trova il punto rappresentativo riferito al centro della sfera
      PT(1) = Tessera(1,ITS) - Sphere(1,NSI)
      PT(2) = Tessera(2,ITS) - Sphere(2,NSI)
      PT(3) = Tessera(3,ITS) - Sphere(3,NSI)
C     Trova il modulo della somma dei vertici della tessera
      DO 2010 JJ = 1, 3
      SUMVERT(JJ) = ZERO
 2010 continue
      DO 2020 L = 1, NV
      SUMVERT(1) = SUMVERT(1) + (VERT(1,L,ITs) - Sphere(1,NSI))
      SUMVERT(2) = SUMVERT(2) + (VERT(2,L,ITs) - Sphere(2,NSI))
      SUMVERT(3) = SUMVERT(3) + (VERT(3,L,ITs) - Sphere(3,NSI))
 2020 continue
      SUM = ZERO
      DO 2030 JJ = 1, 3
      SUM = SUM + SUMVERT(JJ) * SUMVERT(JJ)
 2030 continue
      SUM = sqrt(SUM)
C     Trova la somma delle derivate dei vertici
      DO 2040 JJ = 1, 3
      SUMDER(JJ) = DP(L1,JJ) + DP(L2,JJ)
 2040 continue
C
      PROD = ZERO
      DO 2050 JJ = 1, 3
        PROD = PROD + PT(JJ) * SUMDER(JJ)
 2050 continue
      DO 2060 JJ = 1, 3
      DERP(JJ) = ( SUMDER(JJ) * Sphere(4,NSI) / SUM ) -
     *             ( PT(JJ) * PROD / (Sphere(4,NSI)*SUM) )
 2060 continue
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(NewSph)
      END
*
      SUBROUTINE DVer(IOpt,IC,ITS,L0,L,L2,DX,DY,DZ,
     $                Vert,Centr,nTs,Sphere,IntSph)
      IMPLICIT REAL*8 (A-H,O-Z)
      Parameter (MxVert=20)
      Dimension Sphere(4,*),IntSph(MxVert,*)
      Dimension Vert(3,MxVert,*),Centr(3,MxVert,*)
      DIMENSION P(3),P1(3),P2(3),P3(3)
C
C     Trova la derivata della posizione del vertice L della tessera ITS
C
C     IOpt = 0 : rispetto alla coordinata IC della sfera che, intersecando
C                ISPHE(ITS), forma il lato in esame. (ex derver)
C     IOpt = 1 : rispetto al raggio della sfera aggiunta NSJ che,
C                intersecando la tessera ITS, crea il lato considerato.
C                (ex derver1, IC not referenced)
C
C     Se stiamo considerando il primo vertice del lato : L > 0,
C     il lato che si muove e' L-L2, e il vettore derivata e' diretto
C     lungo la tangente al lato precedente (L0-L).
C     Se stiamo considerando il secondo vertice del lato : L < 0,
C     il lato che si muove e' L0-L, e il vettore derivata e' tangente
C     al lato seguente (L-L2).
C
      ZERO = 0.d0
      IF(L.GT.0) THEN
      L1 = L
        NSJ = INTSPH(L1,ITs)
      ELSE
        L1 = - L
        NSJ = INTSPH(L0,ITs)
      ENDIF
C
C     Il vettore P indica la posizione del vertice rispetto al centro
C     della sfera NSJ
      P(1) = VERT(1,L1,ITs) - Sphere(1,NSJ)
      P(2) = VERT(2,L1,ITs) - Sphere(2,NSJ)
      P(3) = VERT(3,L1,ITs) - Sphere(3,NSJ)
C
C     Trova il versore tangente al lato L0-L1 se stiamo considerando il
C     primo vertice, L2-L1 se consideriamo il secondo
      IF(L.GT.0) THEN
      DO 2000 JJ = 1, 3
        P1(JJ) = VERT(JJ,L1,ITs) - CENTR(JJ,L0,ITs)
        P2(JJ) = VERT(JJ,L0,ITs) - CENTR(JJ,L0,ITs)
 2000 continue
      ELSE
      DO 2010 JJ = 1, 3
        P1(JJ) = VERT(JJ,L1,ITs) - CENTR(JJ,L1,ITs)
        P2(JJ) = VERT(JJ,L2,ITs) - CENTR(JJ,L1,ITs)
 2010 continue
      ENDIF
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2020 JJ = 1, 3
        P2(JJ) = P3(JJ)
 2020 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2030 JJ = 1, 3
        P3(JJ) = P3(JJ) / DNORM3
 2030 continue
C     Trova la derivata
      PROD = P(1)*P3(1) + P(2)*P3(2) + P(3)*P3(3)
      If(IOpt.eq.0) Then
        IF(PROD.EQ.ZERO.AND.P(IC).NE.ZERO) then
          write(6,'("Stop in DVer.")')
          Call Abend()
        EndIf
        IF(PROD.EQ.ZERO) PROD = 1.D0
        FACT = P(IC) / PROD
      Else If(IOpt.eq.1) Then
        IF(PROD.EQ.ZERO) then
          write(6,'("Stop in DVer.")')
          Call Abend()
        EndIf
        FACT = Sphere(4,NSJ) / PROD
      Else
        FACT = 0.0D0
        write(6,'("Illegal IOpt in DVer.")')
        Call Abend()
      EndIf
      DX = FACT * P3(1)
      DY = FACT * P3(2)
      DZ = FACT * P3(3)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nTs)
      END
*
      SUBROUTINE DerPhi(IOpt,IC,NESFJ,ITS,L1,L2,DP,DA,
     $                  Vert,Centr,nTs,Sphere,IntSph,ISphe)
      IMPLICIT REAL*8 (A-H,O-Z)
      Parameter (MxVert=20)
      Dimension Sphere(4,*),ISphe(*),IntSph(MxVert,*)
      Dimension Vert(3,MxVert,*),Centr(3,MxVert,*)
      DIMENSION DP(MxVert,3),V1(3),V2(3),T12(3)
      DIMENSION VEC1(3),VEC2(3),VEC3(3),VEC4(3)
      Save Zero, One, Small1
      Data Zero/0.0d0/, One/1.0d0/, Small1/1.d-6/
C
C     Find the derivative of: Phi(L1) Cos[Theta(L1)]
C
C     IOpt = 0 : refers to the L1 side of tessera ITS, knowing the derivatives
C                of the position of the vertices by which is delimited
C
C     IOpt = 1 : refers to the L1 side of tessera ITS wrt the radius of the
C                NSJ sphere knowing the derivatives of the positions of the
C                vertices by which is delimited
C
C     NS1 is the sphere to which the tessera belongs, NS2 is the sphere that
C     creates the side L1 by intersecting NS1
C
      Small = 1.0d-12
      NS1 = ISPHE(ITS)
      NS2 = INTSPH(L1,ITs)
C
C     Finds the coordinates of the vertices wrt the center of the circle on
C     which is the arc L1 and the radius of the circle
C
      DO 2000 JJ = 1, 3
        V1(JJ) = VERT(JJ,L1,ITs) - CENTR(JJ,L1,ITs)
        V2(JJ) = VERT(JJ,L2,ITs) - CENTR(JJ,L1,ITs)
 2000 continue
      RC2 = V1(1)**2 + V1(2)**2 + V1(3)**2
      PROD = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)
      COSPHI = PROD / RC2
C
C     If the tessera is intersected very near its vertex the numerical
C     approximation could lead to PROD>RC2 and COSPHI>1. Avoid this by:
C
      Delta = 1.d-12
      If(Abs(CosPhi).gt.One) CosPhi = Sign(One-Delta,CosPhi)
      SENPHI = SQRT(One - COSPHI*COSPHI)
C
C     Compute the derivative of Phi(L1)
C
      DO 2010 JJ = 1, 3
      VEC1(JJ) = V1(JJ) - COSPHI*V2(JJ)
      VEC2(JJ) = DP(L2,JJ)
      VEC3(JJ) = V2(JJ) - COSPHI*V1(JJ)
      VEC4(JJ) = DP(L1,JJ)
 2010 continue
C
C     If the side is just that created by the moving sphere some components
C     are corrected
C     IOpt = 0: due to the derivative of the distance between sphere centers
C     IOpt = 1: due to the derivative of the radius of sphere NESFJ
C
      IF(NS2.EQ.NESFJ) THEN
        T12(1) = Sphere(1,NESFJ) - Sphere(1,NS1)
        T12(2) = Sphere(2,NESFJ) - Sphere(2,NS1)
        T12(3) = Sphere(3,NESFJ) - Sphere(3,NS1)
        DIST2 = T12(1)*T12(1) + T12(2)*T12(2) + T12(3)*T12(3)
        If(IOpt.eq.0) Then
          FACT = (Sphere(4,NS1)**2 - Sphere(4,NESFJ)**2 + DIST2)
     $         / (2.d0 * DIST2)
          VEC2(IC) = VEC2(IC) - FACT
          VEC4(IC) = VEC4(IC) - FACT
        Else If(IOpt.eq.1) Then
            DO 2020 JJ = 1, 3
              VEC2(JJ) = VEC2(JJ) + Sphere(4,NESFJ) * T12(JJ) / DIST2
              VEC4(JJ) = VEC4(JJ) + Sphere(4,NESFJ) * T12(JJ) / DIST2
 2020 continue
        Else
          write(6,'("Illegal IOpt in DerPhi.")')
          Call Abend()
        EndIf
      ENDIF
C
      DPHI = 0.D0
      DO 2030 JJ = 1, 3
        DPHI = DPHI - (VEC1(JJ)*VEC2(JJ)+VEC3(JJ)*VEC4(JJ))
 2030 continue
      If(Abs(SenPhi).lt.Small) then
        If(Abs(DPhi).gt.Small1) then
          write(6,'("SenPhi small but not DPhi in DerPhi.")')
          Call Abend()
        EndIf
        DPhi = Zero
      else
        DPHI = DPHI / (RC2*SENPHI)
        endIf
C
C     Compute the cosine of the polar angle
C
      DNORM1 = 0.D0
      DNORM2 = 0.D0
      V1(1) = VERT(1,L1,ITs) - Sphere(1,NS1)
      V1(2) = VERT(2,L1,ITs) - Sphere(2,NS1)
      V1(3) = VERT(3,L1,ITs) - Sphere(3,NS1)
      V2(1) = Sphere(1,NS2) - Sphere(1,NS1)
      V2(2) = Sphere(2,NS2) - Sphere(2,NS1)
      V2(3) = Sphere(3,NS2) - Sphere(3,NS1)
      DO 2040 JJ = 1, 3
        DNORM1 = DNORM1 + V1(JJ)*V1(JJ)
        DNORM2 = DNORM2 + V2(JJ)*V2(JJ)
 2040 continue
      DNORM1 = SQRT(DNORM1)
      DNORM2 = SQRT(DNORM2)
      COSTH=(V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)) / (DNORM1*DNORM2)
C
C     If the side is not that formed by the moving sphere the derivative of
C     the polar angle vanishes
C
      DCOS = 0.D0
      IF(NS2.NE.NESFJ) GOTO 100
C
C     Otherwise:
C
      DO 2050 JJ = 1, 3
      DCOS = DCOS + V2(JJ) * DP(L1,JJ)
 2050 continue
      If(IOpt.eq.0) DCOS = DCOS + V1(IC) -
     $                     Sphere(4,NS1)*COSTH*V2(IC)/DNORM2
      DCOS = DCOS / ( Sphere(4,NS1) * DNORM2 )
 100  CONTINUE
      DA = COSTH * DPHI + ACOS(COSPHI) * DCOS
      DA = Sphere(4,NS1) * Sphere(4,NS1) * DA
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nTs)
      END
*
      SUBROUTINE DerBet(IOpt,IC,NSFE,ITS,L0,L1,L2,DP,DA,
     $                  Vert,Centr,nTs,Sphere,IntSph,ISphe)
      IMPLICIT REAL*8 (A-H,O-Z)
      Parameter (MxVert=20)
      Dimension Sphere(4,*),ISphe(*),IntSph(MxVert,*)
      Dimension Vert(3,MxVert,*),Centr(3,MxVert,*)
      DIMENSION DP(MxVert,3)
      DIMENSION V1(3),V2(3),V3(3),V4(3),DV1(3),DV2(3),DV3(3),DV4(3)
      DIMENSION P1(3),P2(3),P3(3),T0(3),T1(3),T12(3),DT0(3),DT1(3)
C
C     Trova la derivata dell'angolo esterno del vertice L1 della tessera
C     ITS, formato dai lati L0-L1 e L1-L2,
C
C     IOpt = 0 : rispetto alla coordinata IC (ex derbeta)
C     IOpt = 1 : rispetto al raggio (ex derbeta1)
C
C     della sfera NSFE, conoscendo la derivata della posizione del vertice P1.
C
      PI   = 4.0d0*ATan(1.0d0)
      NS1 = ISPHE(ITS)
      NS2 = INTSPH(L0,ITs)
      NS3 = INTSPH(L1,ITs)
C
C     Trova i vettori posizione dei vertici L0, L1 rispetto al centro
C     dell'arco L0, e dei vertici L1, L2 rispetto al centro dell'arco L1
      DO 2000 JJ = 1, 3
        V1(JJ) = VERT(JJ,L0,ITs) - CENTR(JJ,L0,ITs)
        V2(JJ) = VERT(JJ,L1,ITs) - CENTR(JJ,L0,ITs)
        V3(JJ) = VERT(JJ,L1,ITs) - CENTR(JJ,L1,ITs)
        V4(JJ) = VERT(JJ,L2,ITs) - CENTR(JJ,L1,ITs)
 2000 continue
C     Trova la derivata dei vettori posizione
      DO 2010 JJ = 1, 3
        DV1(JJ) = DP(L0,JJ)
        DV2(JJ) = DP(L1,JJ)
        DV3(JJ) = DP(L1,JJ)
        DV4(JJ) = DP(L2,JJ)
 2010 continue
C
C     Corregge la derivata dei vettori posizione se il lato considerato
C     appartiene alla sfera che si muove (distinguendo se ITS appartiene
C     o no alla sfera che si muove)
      INDEX = 0
      IF(NS1.EQ.NSFE.AND.NS2.NE.NSFE) INDEX = 1
      IF(NS1.NE.NSFE.AND.NS2.EQ.NSFE) INDEX = 1
      IF(INDEX.EQ.1) THEN
        T12(1) = Sphere(1,NS2) - Sphere(1,NS1)
        T12(2) = Sphere(2,NS2) - Sphere(2,NS1)
        T12(3) = Sphere(3,NS2) - Sphere(3,NS1)
        DIST2 = T12(1)*T12(1) + T12(2)*T12(2) + T12(3)*T12(3)
        If(IOpt.eq.0) Then
          DO 2020 JJ = 1, 3
            DV1(JJ) = DV1(JJ) + (Sphere(4,NS1)**2 - Sphere(4,NS2)**2) *
     %                T12(IC) * T12(JJ) / (DIST2*DIST2)
            DV2(JJ) = DV2(JJ) + (Sphere(4,NS1)**2 - Sphere(4,NS2)**2) *
     %                T12(IC) * T12(JJ) / (DIST2*DIST2)
 2020 continue
          FACT = (Sphere(4,NS1)**2-Sphere(4,NS2)**2+DIST2)
     $         / (2.D0*DIST2)
          DV1(IC) = DV1(IC) - FACT
          DV2(IC) = DV2(IC) - FACT
        Else If(IOpt.eq.1) Then
          DO 2030 JJ = 1, 3
            DV1(JJ) = DV1(JJ) + Sphere(4,NS2) * T12(JJ) / DIST2
            DV2(JJ) = DV2(JJ) + Sphere(4,NS2) * T12(JJ) / DIST2
 2030 continue
        Else
          write(6,'("Illegal IOpt in DerBet.")')
          Call Abend()
        EndIf
      ENDIF
C
      INDEX = 0
      IF(NS1.EQ.NSFE.AND.NS3.NE.NSFE) INDEX = 1
      IF(NS1.NE.NSFE.AND.NS3.EQ.NSFE) INDEX = 1
      IF(INDEX.EQ.1) THEN
        T12(1) = Sphere(1,NS3) - Sphere(1,NS1)
        T12(2) = Sphere(2,NS3) - Sphere(2,NS1)
        T12(3) = Sphere(3,NS3) - Sphere(3,NS1)
        DIST2 = T12(1)*T12(1) + T12(2)*T12(2) + T12(3)*T12(3)
        If(IOpt.eq.0) Then
          DO 2040 JJ = 1, 3
            DV3(JJ) = DV3(JJ) + (Sphere(4,NS1)**2 - Sphere(4,NS3)**2) *
     %                T12(IC) * T12(JJ) / (DIST2*DIST2)
            DV4(JJ) = DV4(JJ) + (Sphere(4,NS1)**2 - Sphere(4,NS3)**2) *
     %                T12(IC) * T12(JJ) / (DIST2*DIST2)
 2040 continue
          FACT = (Sphere(4,NS1)**2-Sphere(4,NS3)**2+DIST2)
     $         / (2.D0*DIST2)
          DV3(IC) = DV3(IC) - FACT
          DV4(IC) = DV4(IC) - FACT
        Else If(IOpt.eq.1) Then
          DO 2050 JJ = 1, 3
            DV3(JJ) = DV3(JJ) + Sphere(4,NS3) * T12(JJ) / DIST2
            DV4(JJ) = DV4(JJ) + Sphere(4,NS3) * T12(JJ) / DIST2
 2050 continue
        Else
          write(6,'("Illegal IOpt in DerBet.")')
          Call Abend()
        EndIf
      ENDIF
C
C     Trova il vettore tangente al lato L0-L1: T0 = V2 x (V2 x V1)
      DO 2060 JJ = 1, 3
        P1(JJ) = V2(JJ)
        P2(JJ) = V1(JJ)
 2060 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2070 JJ = 1, 3
        P2(JJ) = P3(JJ)
 2070 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2080 JJ = 1, 3
        T0(JJ) = P3(JJ)
 2080 continue
      dnxxt0 = DNORM3
C     Trova il vettore tangente al lato L1-L2: T1 = V3 x (V3 x V4)
      DO 2090 JJ = 1, 3
        P1(JJ) = V3(JJ)
        P2(JJ) = V4(JJ)
 2090 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2100 JJ = 1, 3
        P2(JJ) = P3(JJ)
 2100 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2110 JJ = 1, 3
        T1(JJ) = P3(JJ)
 2110 continue
      dnxxt1 = DNORM3
C
C     Trova l'angolo Beta(L1)
      PROD = T0(1)*T1(1) + T0(2)*T1(2) + T0(3)*T1(3)
      BETA = PI - ACOS( PROD / (dnxxt0*dnxxt1) )
      COSBETA = cos(BETA)
      SENBETA = Sin(BETA)
C
C     Trova la derivata della tangente T0 :
C     DT0 = DV2 x (V2 x V1) + V2 x (DV2 x V1) + V2 x (V2 x DV1)
C
      DO 2120 JJ = 1, 3
        DT0(JJ) = 0.D0
 2120 continue
C
      DO 2130 JJ = 1, 3
        P1(JJ) = V2(JJ)
        P2(JJ) = V1(JJ)
 2130 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2140 JJ = 1, 3
        P1(JJ) = DV2(JJ)
        P2(JJ) = P3(JJ)
 2140 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2150 JJ = 1, 3
        DT0(JJ) = DT0(JJ) + P3(JJ)
 2150 continue
C
      DO 2160 JJ = 1, 3
        P1(JJ) = DV2(JJ)
        P2(JJ) = V1(JJ)
 2160 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2170 JJ = 1, 3
        P1(JJ) = V2(JJ)
        P2(JJ) = P3(JJ)
 2170 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2180 JJ = 1, 3
        DT0(JJ) = DT0(JJ) + P3(JJ)
 2180 continue
C
      DO 2190 JJ = 1, 3
        P1(JJ) = V2(JJ)
        P2(JJ) = DV1(JJ)
 2190 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2200 JJ = 1, 3
        P1(JJ) = V2(JJ)
        P2(JJ) = P3(JJ)
 2200 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2210 JJ = 1, 3
        DT0(JJ) = DT0(JJ) + P3(JJ)
 2210 continue
C
C     Trova la derivata della tangente T1 :
C     DT1 = DV3 x (V3 x V4) + V3 x (DV3 x V4) + V3 x (V3 x DV4)
C
      DO 2220 JJ = 1, 3
        DT1(JJ) = 0.D0
 2220 continue
C
      DO 2230 JJ = 1, 3
        P1(JJ) = V3(JJ)
        P2(JJ) = V4(JJ)
 2230 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2240 JJ = 1, 3
        P1(JJ) = DV3(JJ)
        P2(JJ) = P3(JJ)
 2240 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2250 JJ = 1, 3
        DT1(JJ) = DT1(JJ) + P3(JJ)
 2250 continue
C
      DO 2260 JJ = 1, 3
        P1(JJ) = DV3(JJ)
        P2(JJ) = V4(JJ)
 2260 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2270 JJ = 1, 3
        P1(JJ) = V3(JJ)
        P2(JJ) = P3(JJ)
 2270 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2280 JJ = 1, 3
        DT1(JJ) = DT1(JJ) + P3(JJ)
 2280 continue
C
      DO 2290 JJ = 1, 3
        P1(JJ) = V3(JJ)
        P2(JJ) = DV4(JJ)
 2290 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2300 JJ = 1, 3
        P1(JJ) = V3(JJ)
        P2(JJ) = P3(JJ)
 2300 continue
      CALL VECP(P1,P2,P3,DNORM3)
      DO 2310 JJ = 1, 3
        DT1(JJ) = DT1(JJ) + P3(JJ)
 2310 continue
C
C     Infine calcola la derivata dell'angolo Beta(L1)
      DO 2320 JJ = 1, 3
        P1(JJ) = T1(JJ) + COSBETA * dnxxt1 * T0(JJ) / dnxxt0
        P2(JJ) = T0(JJ) + COSBETA * dnxxt0 * T1(JJ) / dnxxt1
 2320 continue
      DO 2330 JJ = 1, 3
        DBETA = 0.D0
 2330 continue
      DO 2340 JJ = 1, 3
        DBETA = DBETA + ( P1(JJ) * DT0(JJ) + P2(JJ) * DT1(JJ) )
 2340 continue
      DBETA = DBETA / ( SENBETA * dnxxt0 * dnxxt1 )
      DA = Sphere(4,NS1) * Sphere(4,NS1) * DBETA
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nTs)
      END
*
      Function iRToInt(r)
      implicit real*8 (a-h, o-z)
c
c     floating to integer conversion.
c
      iRToInt = int(r)
      return
      end
