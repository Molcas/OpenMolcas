!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Deriva(IPrint,ToAng,NAt,nTs,NEsf,NEsfP,RSolv,Tessera,Vert,Centr,Sphere,ISphe,IntSph,NOrd,NVert,NewSph,DerTes,DerPunt, &
                  DerRad,DerCentr)

implicit real*8(a-h,o-z)
parameter(MxVert=20)
integer AlGe(63), Casca(MxVert)
dimension DISTKI(3), DERP(3), DERPT(3)
dimension Tessera(4,*), Vert(3,MxVert,*), Centr(3,MxVert,*)
dimension Sphere(4,*), ISphe(*)
dimension IntSph(MxVert,*), NewSph(2,NEsf), NOrd(*), NVert(*)
dimension DerTes(nTs,NAt,3), DerPunt(nTs,NAt,3,3)
dimension DerRad(NEsf,NAt,3), DerCentr(NEsf,NAt,3,3)
save Zero, One
data ZERO,One/0.0d0,1.0d0/

! Compute the derivatives of area and of representative point position
! for each tessera

! Le derivate contengono termini dovuti direttamente allo spostamento
! del centro della sfera NSJ, e termini "mediati" dagli
! spostamenti del centro e dal cambiamento del raggio delle sfere
! "aggiunte" (create da PEDRA, oltre a quelle originarie).
! Chiamiamo ITS la tessera, I la sfera a cui ITS appartiene, J la sfera
! che si muove. Le derivate che cerchiamo sono dS(ITS)/dx(J) e
! dQ(ITS)/dx(J) (indicando con x(J) una qualsiasi coordinata di J.
! Per prima cosa calcoliamo le derivate dell'area e del punto
! rappresentativo delle tessere DIRETTAMENTE tagliate da J, poi di quelle
! che appartengono a J e sono tagliate da altre sfere.
! Poi consideriamo le "sfere aggiunte" che vengono modificate dal
! movimento di J: per ogni sfera aggiunta si esaminano tutte le
! tessere e caso per caso si calcolano le derivate.
!
! In ogni caso le derivate totali dS/dx(J) e dQ/dx(J) si scrivono come una
! combinazione delle seguenti derivate parziali:
!      dS/dR(I)  ;  dS/dz(I)  ;  dS/dR(K)  ;  dS/dz(K)
!      dQ/dR(I)  ;  dQ/dz(I)  ;  dQ/dR(K)  ;  dQ/dz(K)
! dove, ancora una volta, z indica una qualunque coordinata e K un'altra
! sfera (la stessa J o una sfera aggiunta modificata dal movimento
! di J) che taglia S.
!
! NESFJ e' la sfera che sta attorno all'atomo NSJ: se NSJ non ha
! nessuna sfera, NESFJ = 0

! Initiate the derivative arrays
NDeg = 3*NAt
call FZero(DerTes(1,1,1),nTs*NDeg)
call FZero(DerPunt(1,1,1,1),3*nTs*NDeg)
call FZero(DerRad(1,1,1),NEsf*NDeg)
call FZero(DerCentr(1,1,1,1),3*NEsf*NDeg)
! The geometric quantities are expected to be in Angstrom
call DScal_(4*nTs,ToAng,Tessera(1,1),1)
call DScal_(nTs,ToAng,Tessera(4,1),4)
call DScal_(4*NEsf,ToAng,Sphere(1,1),1)
call DScal_(3*MxVert*nTs,ToAng,Vert(1,1,1),1)
call DScal_(3*MxVert*nTs,ToAng,Centr(1,1,1),1)

do NSJ=1,NAt
  do ICoord=1,3

    NESFJ = 0
    do I=1,NESFP
      if (NSJ == NORD(I)) NESFJ = I
    end do

    ! Memorizza in DERRAD(NS,NSJ,ICOORD) la derivata del raggio di ogni
    ! sfera NS e in DERCENTR(NS,NSJ,ICOORD,3) le derivate delle
    ! coordinate del centro di NS rispetto alla coord. ICOORD della
    ! sfera NSJ.
    !
    ! Se NS e' una sfera originaria queste derivate sono 0, tranne
    ! DERCENTR(NESFJ,NSJ,ICOORD,ICOORD)=1 :

    do NSFE=1,NESFP
      DERRAD(NSFE,NSJ,ICOORD) = ZERO
      do JJ=1,3
        DERCENTR(NSFE,NSJ,ICOORD,JJ) = ZERO
      end do
    end do
    if (NESFJ /= 0) DERCENTR(NESFJ,NSJ,ICOORD,ICOORD) = 1.0d0

    ! Se non c'e' nessuna sfera attorno all'atomo NSJ ...
    if (NESFJ == 0) goto 1000

    ! 1) Effetti diretti.
    ! Loop sulle tessere.
    do ITS=1,NTS
      DERTS = ZERO
      do JJ=1,3
        DERPT(JJ) = ZERO
      end do
      NV = NVERT(ITS)
      NSI = ISPHE(ITS)
      if (NSI == NESFJ) goto 200

      ! Derivate nel caso I non = J
      !
      !           dS/dx(J), dQ/dx(J)
      ! ITS ha un lato su J ?
      ICONT = 0
      do N=1,NV
        if (INTSPH(N,ITs) == NESFJ) ICONT = ICONT+1
      end do
      if (ICONT >= 1) then
        call dSd(0,ITS,ICOORD,NESFJ,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
        DERTS = DERTS+DERS
        do JJ=1,3
          DERPT(JJ) = DERPT(JJ)+DERP(JJ)
        end do
      end if
      goto 150
200   continue
      ! Derivate nel caso I = J

      ! Loop sulle sfere K che intersecano I
      do NSK=1,NESF
        if (NSK == NSI) goto 300
        ! ITS ha un lato su K ?
        ICONT = 0
        do N=1,NV
          if (INTSPH(N,ITs) == NSK) ICONT = ICONT+1
        end do
        if (ICONT == 0) goto 300
        call dSd(0,ITS,ICOORD,NSK,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
        DERTS = DERTS-DERS
        do JJ=1,3
          DERPT(JJ) = DERPT(JJ)-DERP(JJ)
        end do
300     continue
      end do
150   continue
      DERTES(ITS,NSJ,ICOORD) = DERTS
      do JJ=1,3
        DERPUNT(ITS,NSJ,ICOORD,JJ) = DERPT(JJ)
      end do
    end do

    ! 2) Effetti indiretti.
    ! Loop sulle sfere aggiunte
    do NSA=NESFP+1,NESF
      do II=1,63
        ALGE(II) = 0
      end do
      ! Costruiamo l'"albero genealogico" della sfera NSA
      ALGE(1) = NSA
      ALGE(2) = abs(NEWSPH(1,NSA))
      ALGE(3) = abs(NEWSPH(2,NSA))
      LIVEL = 3
      NUMBER = 2
510   continue
      NSUB = 1
      do II=LIVEL-NUMBER+1,LIVEL
        if (ALGE(II) > NESFP) then
          ALGE(LIVEL+NSUB) = abs(NEWSPH(1,ALGE(II)))
          ALGE(LIVEL+NSUB+1) = abs(NEWSPH(2,ALGE(II)))
        end if
        NSUB = NSUB+2
      end do
      NUMBER = NUMBER*2
      LIVEL = LIVEL+NUMBER
      if (NUMBER < 32) goto 510
      ! Si accerta che nell'ultimo livello ci siano solo sfere originarie
      do II=34,63
        if (ALGE(II) > NESFP) then
          write(6,10) NSA,ALGE(II)
        end if
      end do
      ! Quando un elemento di ALGE e' = NESFJ, costruisce la corrispondente
      ! "cascata" di sfere aggiunte che collega NESFJ a NSA
      do LIVEL=2,6
        MIN = 2**(LIVEL-1)
        MAX = 2**(LIVEL)-1
        do II=MIN,MAX
          if (ALGE(II) /= NESFJ) goto 700
          do K=1,MxVert
            CASCA(K) = 0
          end do
          CASCA(1) = NESFJ
          INDEX = II
          K = 2
          do LL=LIVEL,2,-1
            FACT = (INDEX-2**(LL-1))/2.d0
            INDEX = 2**(LL-2)+iRToInt(FACT)
            CASCA(K) = ALGE(INDEX)
            K = K+1
          end do
          ! Contiamo gli elementi diversi da 0 in CASCA
          ICONT = 0
          do K=1,MxVert
            if (CASCA(K) /= 0) ICONT = ICONT+1
          end do

          ! Costruiamo le derivate composte del raggio e delle coordinate di
          ! NSA (ultimo elemento di CASCA)
          ! rispetto alla coordinata ICOORD di NESFJ (primo elemento di CASCA)
          NS1 = CASCA(1)
          NS2 = CASCA(2)
          call dRdC(NS2,ICOORD,NS1,DR,RSolv,Sphere,NewSph)
          call dCdC(1,NS2,ICOORD,NS1,DX,Sphere,NewSph)
          call dCdC(2,NS2,ICOORD,NS1,DY,Sphere,NewSph)
          call dCdC(3,NS2,ICOORD,NS1,DZ,Sphere,NewSph)
          do I=3,ICONT
            DDR = ZERO
            DDX = ZERO
            DDY = ZERO
            DDZ = ZERO
            NS1 = CASCA(I-1)
            NS2 = CASCA(I)

            ! Derivata del raggio dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dRdR(NS2,NS1,DER,RSolv,Sphere,NewSph)
            DDR = DER*DR
            call dRdC(NS2,1,NS1,DER,RSolv,Sphere,NewSph)
            DDR = DDR+DER*DX
            call dRdC(NS2,2,NS1,DER,RSolv,Sphere,NewSph)
            DDR = DDR+DER*DY
            call dRdC(NS2,3,NS1,DER,RSolv,Sphere,NewSph)
            DDR = DDR+DER*DZ

            ! Derivata della coord. X dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dCdR(1,NS2,NS1,DER,Sphere,NewSph)
            DDX = DER*DR
            call dCdC(1,NS2,1,NS1,DER,Sphere,NewSph)
            DDX = DDX+DER*DX
            call dCdC(1,NS2,2,NS1,DER,Sphere,NewSph)
            DDX = DDX+DER*DY
            call dCdC(1,NS2,3,NS1,DER,Sphere,NewSph)
            DDX = DDX+DER*DZ

            ! Derivata della coord. Y dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dCdR(2,NS2,NS1,DER,Sphere,NewSph)
            DDY = DER*DR
            call dCdC(2,NS2,1,NS1,DER,Sphere,NewSph)
            DDY = DDY+DER*DX
            call dCdC(2,NS2,2,NS1,DER,Sphere,NewSph)
            DDY = DDY+DER*DY
            call dCdC(2,NS2,3,NS1,DER,Sphere,NewSph)
            DDY = DDY+DER*DZ

            ! Derivata della coord. Z dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dCdR(3,NS2,NS1,DER,Sphere,NewSph)
            DDZ = DER*DR
            call dCdC(3,NS2,1,NS1,DER,Sphere,NewSph)
            DDZ = DDZ+DER*DX
            call dCdC(3,NS2,2,NS1,DER,Sphere,NewSph)
            DDZ = DDZ+DER*DY
            call dCdC(3,NS2,3,NS1,DER,Sphere,NewSph)
            DDZ = DDZ+DER*DZ

            DR = DDR
            DX = DDX
            DY = DDY
            DZ = DDZ
          end do

          ! Se NS e' una sfera aggiunta, memorizza le derivate del raggio
          ! e delle coordinate del centro:
          DERRAD(NSA,NSJ,ICOORD) = DR
          DERCENTR(NSA,NSJ,ICOORD,1) = DX
          DERCENTR(NSA,NSJ,ICOORD,2) = DY
          DERCENTR(NSA,NSJ,ICOORD,3) = DZ

          ! Per ogni sfera aggiunta connessa a NESFJ, fa passare tutte le
          ! tessere, calcolando, se e' il caso, le derivate
          do ITS=1,NTS
            DERTS = ZERO
            do JJ=1,3
              DERPT(JJ) = ZERO
            end do
            NV = NVERT(ITS)
            NSI = ISPHE(ITS)
            ! Derivate delle tessere appartenenti a NSA
            if (NSI == NSA) then
              ! Derivate relative al raggio di NSI: 1) scaling
              DERTS = 2.d0*Tessera(4,ITS)*DR/Sphere(4,NSI)
              DERPT(1) = (Tessera(1,ITS)-Sphere(1,NSI))*DR/Sphere(4,NSI)
              DERPT(2) = (Tessera(2,ITS)-Sphere(2,NSI))*DR/Sphere(4,NSI)
              DERPT(3) = (Tessera(3,ITS)-Sphere(3,NSI))*DR/Sphere(4,NSI)
              ! Loop sulle altre sfere K che tagliano ITS
              do NSK=1,NESF
                if (NSK == NSI) goto 950
                ICONT = 0
                do N=1,NV
                  if (INTSPH(N,ITs) == NSK) ICONT = ICONT+1
                end do
                if (ICONT == 0) goto 950
                ! Derivate relative al raggio di NSI: 2) raggio delle altre sfere
                call dSd(1,ITS,IJunk,NSK,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
                DERTS = DERTS-DERS*Sphere(4,NSK)*DR/Sphere(4,NSI)
                do JJ=1,3
                  DERPT(JJ) = DERPT(JJ)-DERP(JJ)*Sphere(4,NSK)*DR/Sphere(4,NSI)
                end do
                ! Derivate relative al raggio di NSI: 3) coord. delle altre sfere
                DISTKI(1) = Sphere(1,NSK)-Sphere(1,NSI)
                DISTKI(2) = Sphere(2,NSK)-Sphere(2,NSI)
                DISTKI(3) = Sphere(3,NSK)-Sphere(3,NSI)
                FAC = DR/Sphere(4,NSI)
                do JJ=1,3
                  call dSd(0,ITS,JJ,NSK,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
                  if (DISTKI(JJ) /= ZERO) then
                    DERTS = DERTS-FAC*DERS*DISTKI(JJ)
                    do LL=1,3
                      DERPT(LL) = DERPT(LL)-FAC*DERP(LL)*DISTKI(JJ)
                    end do
                  end if
                end do
                ! Derivate relative alle coordinate di NSI
                call dSd(0,ITS,1,NSK,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
                DERTS = DERTS-DX*DERS
                do JJ=1,3
                  DERPT(JJ) = DERPT(JJ)-DX*DERP(JJ)
                end do
                call dSd(0,ITS,2,NSK,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
                DERTS = DERTS-DY*DERS
                do JJ=1,3
                  DERPT(JJ) = DERPT(JJ)-DY*DERP(JJ)
                end do
                call dSd(0,ITS,3,NSK,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
                DERTS = DERTS-DZ*DERS
                do JJ=1,3
                  DERPT(JJ) = DERPT(JJ)-DZ*DERP(JJ)
                end do
950             continue
              end do
              DERTES(ITS,NSJ,ICOORD) = DERTES(ITS,NSJ,ICOORD)+DERTS
              do JJ=1,3
                DERPUNT(ITS,NSJ,ICOORD,JJ) = DERPUNT(ITS,NSJ,ICOORD,JJ)+DERPT(JJ)
              end do
            else
              ! Derivate delle tessere tagliate da NSA
              ICONT = 0
              do N=1,NV
                if (INTSPH(N,ITs) == NSA) ICONT = ICONT+1
              end do
              if (ICONT == 0) goto 900
              ! Derivate relative al raggio di NSA
              call dSd(1,ITS,IJunk,NSA,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
              DERTS = DERS*DR
              do JJ=1,3
                DERPT(JJ) = DERP(JJ)*DR
              end do
              ! Derivate relative alle coordinate di NSA
              call dSd(0,ITS,1,NSA,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
              DERTS = DERTS+DX*DERS
              do JJ=1,3
                DERPT(JJ) = DERPT(JJ)+DX*DERP(JJ)
              end do
              call dSd(0,ITS,2,NSA,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
              DERTS = DERTS+DY*DERS
              do JJ=1,3
                DERPT(JJ) = DERPT(JJ)+DY*DERP(JJ)
              end do
              call dSd(0,ITS,3,NSA,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)
              DERTS = DERTS+DZ*DERS
              do JJ=1,3
                DERPT(JJ) = DERPT(JJ)+DZ*DERP(JJ)
              end do

              DERTES(ITS,NSJ,ICOORD) = DERTES(ITS,NSJ,ICOORD)+DERTS
              do JJ=1,3
                DERPUNT(ITS,NSJ,ICOORD,JJ) = DERPUNT(ITS,NSJ,ICOORD,JJ)+DERPT(JJ)
              end do
            end if
900         continue
          end do
700       continue
        end do
      end do
    end do

    ! Controlla che nessuna derivata sia eccessiva (per imprecisioni numeriche)
    do ITS=1,NTS
      if (abs(DERTES(ITS,NSJ,ICOORD)) > 10.d0) then
        if (IPrint >= 2) then
          write(6,*) 'DERIVATA GEOMETRICA ECCESSIVA TRASCURATA'
          write(6,*) 'TESSERA, SFERA, COORDINATA',ITS,NESFJ,ICOORD
        end if
        DERTES(ITS,NSJ,ICOORD) = 0.d0
      end if
    end do

    ! End of loop on atoms and coordinates

  end do
1000 continue
end do
! Put the geometric quantities in Bohr again
ToBohr = One/ToAng
call DScal_(4*nTs,ToBohr,Tessera(1,1),1)
call DScal_(nTs,ToBohr,Tessera(4,1),4)
call DScal_(4*NEsf,ToBohr,Sphere(1,1),1)
call DScal_(30*nTs,ToBohr,Vert(1,1,1),1)
call DScal_(30*nTs,ToBohr,Centr(1,1,1),1)

return

10 format(/,'SFERA ',I3,/,'AL LIVELLO 5 LA SFERA',I3," E' AGGIUNTA")

end subroutine Deriva
!====
subroutine DCDC(JJ,NSI,ICOORD,NESFJ,DC,Sphere,NewSph)

implicit real*8(A-H,O-Z)
dimension Sphere(4,*), NewSph(2,*)
dimension COORDJ(3), COORDK(3)

! Trova la derivata della coordinata JJ del centro della sfera
! NSI rispetto alla coordinata ICOORD di NSJ, che interseca NSI.
!
! La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e NSK
!
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1,NSI) < 0) .or. (NEWSPH(2,NSI) < 0)) goto 100
K = NEWSPH(1,NSI)
if (K == NESFJ) K = NEWSPH(2,NSI)
COORDJ(1) = Sphere(1,NESFJ)
COORDJ(2) = Sphere(2,NESFJ)
COORDJ(3) = Sphere(3,NESFJ)
COORDK(1) = Sphere(1,K)
COORDK(2) = Sphere(2,K)
COORDK(3) = Sphere(3,K)
D2 = (Sphere(1,NESFJ)-Sphere(1,K))**2+(Sphere(2,NESFJ)-Sphere(2,K))**2+(Sphere(3,NESFJ)-Sphere(3,K))**2
D = sqrt(D2)
DC = (Sphere(4,NESFJ)-Sphere(4,K))*(COORDJ(ICOORD)-COORDK(ICOORD))*(COORDJ(JJ)-COORDK(JJ))/(2.d0*D**3)
if (JJ == ICOORD) DC = DC+0.5d0-(Sphere(4,NESFJ)-Sphere(4,K))/(2.d0*D)
goto 200

100 continue
NSK = NEWSPH(1,NSI)
if (abs(NSK) == NESFJ) NSK = NEWSPH(2,NSI)
COORDJ(1) = Sphere(1,NESFJ)
COORDJ(2) = Sphere(2,NESFJ)
COORDJ(3) = Sphere(3,NESFJ)
COORDK(1) = Sphere(1,abs(NSK))
COORDK(2) = Sphere(2,abs(NSK))
COORDK(3) = Sphere(3,abs(NSK))
D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
D = sqrt(D2)
if (NSK > 0) then
  DC = Sphere(4,NESFJ)*(COORDJ(JJ)-COORDK(JJ))*(COORDJ(ICOORD)-COORDK(ICOORD))/D**3
  if (ICOORD == JJ) DC = DC+1.d0-Sphere(4,NESFJ)/D
else
  DC = -Sphere(4,abs(NSK))*(COORDK(JJ)-COORDJ(JJ))*(COORDK(ICOORD)-COORDJ(ICOORD))/D**3
  if (ICOORD == JJ) DC = DC+Sphere(4,abs(NSK))/D
end if

200 continue

return

end subroutine DCDC
!====
subroutine dCdR(JJ,NSI,NESFJ,DC,Sphere,NewSph)

implicit real*8(A-H,O-Z)
dimension Sphere(4,*), NewSph(2,*)
dimension COORDJ(3), COORDK(3)
data ZERO/0.0d0/

! Trova la derivata della coordinata JJ del centro della sfera
! NSI rispetto al raggio dellla sfera NSJ.
!
! La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e NSK
!
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1,NSI) < 0) .or. (NEWSPH(2,NSI) < 0)) goto 100
NSK = NEWSPH(1,NSI)
if (NSK == NESFJ) NSK = NEWSPH(2,NSI)
COORDJ(1) = Sphere(1,NESFJ)
COORDJ(2) = Sphere(2,NESFJ)
COORDJ(3) = Sphere(3,NESFJ)
COORDK(1) = Sphere(1,NSK)
COORDK(2) = Sphere(2,NSK)
COORDK(3) = Sphere(3,NSK)
D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2+(Sphere(2,NESFJ)-Sphere(2,NSK))**2+(Sphere(3,NESFJ)-Sphere(3,NSK))**2
D = sqrt(D2)
DC = -(COORDJ(JJ)-COORDK(JJ))/(2.d0*D)
goto 200

100 continue
NSK = NEWSPH(1,NSI)
if (abs(NSK) == NESFJ) NSK = NEWSPH(2,NSI)
DC = ZERO
if (NSK < ZERO) goto 200
COORDJ(1) = Sphere(1,NESFJ)
COORDJ(2) = Sphere(2,NESFJ)
COORDJ(3) = Sphere(3,NESFJ)
COORDK(1) = Sphere(1,NSK)
COORDK(2) = Sphere(2,NSK)
COORDK(3) = Sphere(3,NSK)
D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2+(Sphere(2,NESFJ)-Sphere(2,NSK))**2+(Sphere(3,NESFJ)-Sphere(3,NSK))**2
D = sqrt(D2)
DC = -(COORDJ(JJ)-COORDK(JJ))/D

200 continue

return

end subroutine dCdR
!====
subroutine dRdC(NSI,ICOORD,NESFJ,DR,RSolv,Sphere,NewSph)

implicit real*8(A-H,O-Z)
dimension Sphere(4,*), NewSph(2,*)
dimension COORDJ(3), COORDK(3)

! Trova la derivata del raggio della sfera NSI rispetto alla
! coordinata ICOORD (1=X, 2=Y, 3=Z) della sfera NSJ, che interseca
! NSI.
!
! La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e K
!
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1,NSI) < 0) .or. (NEWSPH(2,NSI) < 0)) goto 100
K = NEWSPH(1,NSI)
if (K == NESFJ) K = NEWSPH(2,NSI)
COORDJ(1) = Sphere(1,NESFJ)
COORDJ(2) = Sphere(2,NESFJ)
COORDJ(3) = Sphere(3,NESFJ)
COORDK(1) = Sphere(1,K)
COORDK(2) = Sphere(2,K)
COORDK(3) = Sphere(3,K)
D2 = (Sphere(1,NESFJ)-Sphere(1,K))**2+(Sphere(2,NESFJ)-Sphere(2,K))**2+(Sphere(3,NESFJ)-Sphere(3,K))**2
D = sqrt(D2)
B = 0.5d0*(D+Sphere(4,NESFJ)-Sphere(4,K))
RS = RSOLV
A = ((Sphere(4,NESFJ)+RS)**2+D2-(Sphere(4,K)+RS)**2)/D
DR = (2.d0*A*B-2.d0*B*D-A*D)*(COORDJ(ICOORD)-COORDK(ICOORD))/(4.d0*D2*(Sphere(4,NSI)+RS))
goto 200

100 continue
NSK = NEWSPH(1,NSI)
if (abs(NSK) == NESFJ) NSK = NEWSPH(2,NSI)
COORDJ(1) = Sphere(1,NESFJ)
COORDJ(2) = Sphere(2,NESFJ)
COORDJ(3) = Sphere(3,NESFJ)
COORDK(1) = Sphere(1,abs(NSK))
COORDK(2) = Sphere(2,abs(NSK))
COORDK(3) = Sphere(3,abs(NSK))
RI = Sphere(4,NSI)+RSOLV
RJ = Sphere(4,NESFJ)+RSOLV
RK = Sphere(4,abs(NSK))+RSOLV
DIFF = COORDJ(ICOORD)-COORDK(ICOORD)
D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
D = sqrt(D2)
FAC = Sphere(4,NESFJ)*(RJ*RJ-D*D-RK*RK)
if (NSK < 0) FAC = Sphere(4,abs(NSK))*(RK*RK-D*D-RJ*RJ)
DR = DIFF*FAC/(2.d0*D**3*RI)

200 continue

return

end subroutine dRdC
!====
subroutine dRdR(NSI,NESFJ,DR,RSolv,Sphere,NewSph)

implicit real*8(A-H,O-Z)
dimension Sphere(4,*), NewSph(2,*)

! Trova la derivata del raggio della sfera NSI rispetto al raggio
! della sfera NSJ.
!
! La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e NSK
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1,NSI) < 0) .or. (NEWSPH(2,NSI) < 0)) goto 100
NSK = NEWSPH(1,NSI)
if (NSK == NESFJ) NSK = NEWSPH(2,NSI)
RS = RSOLV
RJ = Sphere(4,NESFJ)+RS
RK = Sphere(4,NSK)+RS
RI = Sphere(4,NSI)+RS
D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2+(Sphere(2,NESFJ)-Sphere(2,NSK))**2+(Sphere(3,NESFJ)-Sphere(3,NSK))**2
D = sqrt(D2)
DR = (-3.d0*RJ*RJ+RK*RK+2.d0*RJ*RK+3.d0*D*RJ-D*RK)/(4.d0*D*RI)
goto 200

100 continue
NSK = NEWSPH(1,NSI)
if (abs(NSK) == NESFJ) NSK = NEWSPH(2,NSI)

if (NSK > 0) then
  RS = RSOLV
  RJ = Sphere(4,NESFJ)+RS
  RK = Sphere(4,NSK)+RS
  RI = Sphere(4,NSI)+RS
  D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2+(Sphere(2,NESFJ)-Sphere(2,NSK))**2+(Sphere(3,NESFJ)-Sphere(3,NSK))**2
  D = sqrt(D2)
  DR = (2.d0*D*RJ+2.d0*D*Sphere(4,NESFJ)-2.d0*RJ*Sphere(4,NESFJ)+D*D-RJ*RJ-RK*RK)/(2.d0*D*RI)
else
  RS = RSOLV
  RJ = Sphere(4,NESFJ)+RS
  RI = Sphere(4,NSI)+RS
  D2 = (Sphere(1,NESFJ)-Sphere(1,abs(NSK)))**2+(Sphere(2,NESFJ)-Sphere(2,abs(NSK)))**2+(Sphere(3,NESFJ)-Sphere(3,abs(NSK)))**2
  D = sqrt(D2)
  DR = (Sphere(4,abs(NSK))*RJ)/(D*RI)
end if
200 continue

return

end subroutine dRdR
!====
subroutine dSd(IOpt,ITS,ICOORD,NESFJ,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)

implicit real*8(A-H,O-Z)
parameter(MxVert=20)
dimension Sphere(4,*), ISphe(*), NewSph(2,*)
dimension IntSph(MxVert,*), NVert(*)
dimension Tessera(4,*), Vert(3,MxVert,*), Centr(3,MxVert,*)
dimension DP(MxVert,3), DERP(3), SUMDER(3), SUMVERT(3), PT(3)
data ZERO/0.0d0/

! Trova le derivate dell'area e del punto rappresentativo della
! tessera ITS rispetto
!
!   IOpt = 0 : alla coordinata ICOORD(1=X, 2=Y, 3=Z) dell'atomo NSJ.
!              (ex dsdx)
!   IOpt = 1 : al raggio della sfera NSJ. (ex dsdr)
!
! NESFJ e' la sfera attorno all'atomo NSJ.
! ITS appartiene alla sfera NSI (diversa da NESFJ, e intersecata
! da NESFJ).

DERTS = ZERO
NSI = ISPHE(ITS)
NV = NVERT(ITS)
do L=1,NV
  DP(L,1) = ZERO
  DP(L,2) = ZERO
  DP(L,3) = ZERO
end do
do L=1,NV
  if (INTSPH(L,ITs) /= NESFJ) goto 100
  ! L1 e' il lato di ITS che sta sulla sfera NESFJ
  L1 = L
  L0 = L-1
  if (L1 == 1) L0 = NV
  L2 = L1+1
  if (L1 == NV) L2 = 1
  L3 = L2+1
  if (L2 == NV) L3 = 1
  ! Trova la derivata del vertice L1: DP(L1)
  call DVer(IOpt,ICOORD,ITS,L0,L1,L2,DP(L1,1),DP(L1,2),DP(L1,3),Vert,Centr,nTs,Sphere,IntSph)
  ! Trova la derivata del vertice L2: DP(L2)
  call DVer(IOpt,ICOORD,ITS,L1,-L2,L3,DP(L2,1),DP(L2,2),DP(L2,3),Vert,Centr,nTs,Sphere,IntSph)

  ! 1) Derivata dell'area della tessera.

  ! Calcola i contributi dovuti ai lati L0-L1, L1-L2, L2-L3
  if (INTSPH(L0,ITs) /= NSI) then
    call DerPhi(IOpt,ICOORD,NESFJ,ITS,L0,L1,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
    DERTS = DERTS+DA
  end if
  call DerPhi(IOpt,ICOORD,NESFJ,ITS,L1,L2,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
  DERTS = DERTS+DA
  if (INTSPH(L2,ITs) /= ISPHE(ITS)) then
    call DerPhi(IOpt,ICOORD,NESFJ,ITS,L2,L3,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
    DERTS = DERTS+DA
  end if

  ! Calcola i contributi dovuti ai vertici L1 e L2
  call DerBet(IOpt,ICOORD,NESFJ,ITS,L0,L1,L2,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
  DERTS = DERTS-DA
  call DerBet(IOpt,ICOORD,NESFJ,ITS,L1,L2,L3,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
  DERTS = DERTS-DA
100 continue
end do
DERS = DERTS

! 2) Derivata del punto rappresentativo.

! Trova il punto rappresentativo riferito al centro della sfera
PT(1) = Tessera(1,ITS)-Sphere(1,NSI)
PT(2) = Tessera(2,ITS)-Sphere(2,NSI)
PT(3) = Tessera(3,ITS)-Sphere(3,NSI)
! Trova il modulo della somma dei vertici della tessera
do JJ=1,3
  SUMVERT(JJ) = ZERO
end do
do L=1,NV
  SUMVERT(1) = SUMVERT(1)+(VERT(1,L,ITs)-Sphere(1,NSI))
  SUMVERT(2) = SUMVERT(2)+(VERT(2,L,ITs)-Sphere(2,NSI))
  SUMVERT(3) = SUMVERT(3)+(VERT(3,L,ITs)-Sphere(3,NSI))
end do
SUM = ZERO
do JJ=1,3
  SUM = SUM+SUMVERT(JJ)*SUMVERT(JJ)
end do
SUM = sqrt(SUM)
! Trova la somma delle derivate dei vertici
do JJ=1,3
  SUMDER(JJ) = DP(L1,JJ)+DP(L2,JJ)
end do

PROD = ZERO
do JJ=1,3
  PROD = PROD+PT(JJ)*SUMDER(JJ)
end do
do JJ=1,3
  DERP(JJ) = (SUMDER(JJ)*Sphere(4,NSI)/SUM)-(PT(JJ)*PROD/(Sphere(4,NSI)*SUM))
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(NewSph)

end subroutine dSd
!====
subroutine DVer(IOpt,IC,ITS,L0,L,L2,DX,DY,DZ,Vert,Centr,nTs,Sphere,IntSph)

implicit real*8(A-H,O-Z)
parameter(MxVert=20)
dimension Sphere(4,*), IntSph(MxVert,*)
dimension Vert(3,MxVert,*), Centr(3,MxVert,*)
dimension P(3), P1(3), P2(3), P3(3)

! Trova la derivata della posizione del vertice L della tessera ITS
!
! IOpt = 0 : rispetto alla coordinata IC della sfera che, intersecando
!            ISPHE(ITS), forma il lato in esame. (ex derver)
! IOpt = 1 : rispetto al raggio della sfera aggiunta NSJ che,
!            intersecando la tessera ITS, crea il lato considerato.
!            (ex derver1, IC not referenced)
!
! Se stiamo considerando il primo vertice del lato : L > 0,
! il lato che si muove e' L-L2, e il vettore derivata e' diretto
! lungo la tangente al lato precedente (L0-L).
! Se stiamo considerando il secondo vertice del lato : L < 0,
! il lato che si muove e' L0-L, e il vettore derivata e' tangente
! al lato seguente (L-L2).

ZERO = 0.d0
if (L > 0) then
  L1 = L
  NSJ = INTSPH(L1,ITs)
else
  L1 = -L
  NSJ = INTSPH(L0,ITs)
end if

! Il vettore P indica la posizione del vertice rispetto al centro
! della sfera NSJ
P(1) = VERT(1,L1,ITs)-Sphere(1,NSJ)
P(2) = VERT(2,L1,ITs)-Sphere(2,NSJ)
P(3) = VERT(3,L1,ITs)-Sphere(3,NSJ)

! Trova il versore tangente al lato L0-L1 se stiamo considerando il
! primo vertice, L2-L1 se consideriamo il secondo
if (L > 0) then
  do JJ=1,3
    P1(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L0,ITs)
    P2(JJ) = VERT(JJ,L0,ITs)-CENTR(JJ,L0,ITs)
  end do
else
  do JJ=1,3
    P1(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L1,ITs)
    P2(JJ) = VERT(JJ,L2,ITs)-CENTR(JJ,L1,ITs)
  end do
end if
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P3(JJ) = P3(JJ)/DNORM3
end do
! Trova la derivata
PROD = P(1)*P3(1)+P(2)*P3(2)+P(3)*P3(3)
if (IOpt == 0) then
  if ((PROD == ZERO) .and. (P(IC) /= ZERO)) then
    write(6,'(a)') 'Stop in DVer.'
    call Abend()
  end if
  if (PROD == ZERO) PROD = 1.d0
  FACT = P(IC)/PROD
else if (IOpt == 1) then
  if (PROD == ZERO) then
    write(6,'(a)') 'Stop in DVer.'
    call Abend()
  end if
  FACT = Sphere(4,NSJ)/PROD
else
  FACT = 0.0d0
  write(6,'(a)') 'Illegal IOpt in DVer.'
  call Abend()
end if
DX = FACT*P3(1)
DY = FACT*P3(2)
DZ = FACT*P3(3)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nTs)

end subroutine DVer
!====
subroutine DerPhi(IOpt,IC,NESFJ,ITS,L1,L2,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)

implicit real*8(A-H,O-Z)
parameter(MxVert=20)
dimension Sphere(4,*), ISphe(*), IntSph(MxVert,*)
dimension Vert(3,MxVert,*), Centr(3,MxVert,*)
dimension DP(MxVert,3), V1(3), V2(3), T12(3)
dimension VEC1(3), VEC2(3), VEC3(3), VEC4(3)
save Zero, One, Small1
data Zero/0.0d0/,One/1.0d0/,Small1/1.d-6/

! Find the derivative of: Phi(L1) Cos[Theta(L1)]
!
! IOpt = 0 : refers to the L1 side of tessera ITS, knowing the derivatives
!            of the position of the vertices by which is delimited
!
! IOpt = 1 : refers to the L1 side of tessera ITS wrt the radius of the
!            NSJ sphere knowing the derivatives of the positions of the
!            vertices by which is delimited
!
! NS1 is the sphere to which the tessera belongs, NS2 is the sphere that
! creates the side L1 by intersecting NS1

Small = 1.0d-12
NS1 = ISPHE(ITS)
NS2 = INTSPH(L1,ITs)

! Finds the coordinates of the vertices wrt the center of the circle on
! which is the arc L1 and the radius of the circle

do JJ=1,3
  V1(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L1,ITs)
  V2(JJ) = VERT(JJ,L2,ITs)-CENTR(JJ,L1,ITs)
end do
RC2 = V1(1)**2+V1(2)**2+V1(3)**2
PROD = V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
COSPHI = PROD/RC2

! If the tessera is intersected very near its vertex the numerical
! approximation could lead to PROD>RC2 and COSPHI>1. Avoid this by:

Delta = 1.d-12
if (abs(CosPhi) > One) CosPhi = sign(One-Delta,CosPhi)
SENPHI = sqrt(One-COSPHI*COSPHI)

! Compute the derivative of Phi(L1)

do JJ=1,3
  VEC1(JJ) = V1(JJ)-COSPHI*V2(JJ)
  VEC2(JJ) = DP(L2,JJ)
  VEC3(JJ) = V2(JJ)-COSPHI*V1(JJ)
  VEC4(JJ) = DP(L1,JJ)
end do

! If the side is just that created by the moving sphere some components
! are corrected
! IOpt = 0: due to the derivative of the distance between sphere centers
! IOpt = 1: due to the derivative of the radius of sphere NESFJ

if (NS2 == NESFJ) then
  T12(1) = Sphere(1,NESFJ)-Sphere(1,NS1)
  T12(2) = Sphere(2,NESFJ)-Sphere(2,NS1)
  T12(3) = Sphere(3,NESFJ)-Sphere(3,NS1)
  DIST2 = T12(1)*T12(1)+T12(2)*T12(2)+T12(3)*T12(3)
  if (IOpt == 0) then
    FACT = (Sphere(4,NS1)**2-Sphere(4,NESFJ)**2+DIST2)/(2.d0*DIST2)
    VEC2(IC) = VEC2(IC)-FACT
    VEC4(IC) = VEC4(IC)-FACT
  else if (IOpt == 1) then
    do JJ=1,3
      VEC2(JJ) = VEC2(JJ)+Sphere(4,NESFJ)*T12(JJ)/DIST2
      VEC4(JJ) = VEC4(JJ)+Sphere(4,NESFJ)*T12(JJ)/DIST2
    end do
  else
    write(6,'(a)') 'Illegal IOpt in DerPhi.'
    call Abend()
  end if
end if

DPHI = 0.d0
do JJ=1,3
  DPHI = DPHI-(VEC1(JJ)*VEC2(JJ)+VEC3(JJ)*VEC4(JJ))
end do
if (abs(SenPhi) < Small) then
  if (abs(DPhi) > Small1) then
    write(6,'(a)') 'SenPhi small but not DPhi in DerPhi.'
    call Abend()
  end if
  DPhi = Zero
else
  DPHI = DPHI/(RC2*SENPHI)
end if

! Compute the cosine of the polar angle

DNORM1 = 0.d0
DNORM2 = 0.d0
V1(1) = VERT(1,L1,ITs)-Sphere(1,NS1)
V1(2) = VERT(2,L1,ITs)-Sphere(2,NS1)
V1(3) = VERT(3,L1,ITs)-Sphere(3,NS1)
V2(1) = Sphere(1,NS2)-Sphere(1,NS1)
V2(2) = Sphere(2,NS2)-Sphere(2,NS1)
V2(3) = Sphere(3,NS2)-Sphere(3,NS1)
do JJ=1,3
  DNORM1 = DNORM1+V1(JJ)*V1(JJ)
  DNORM2 = DNORM2+V2(JJ)*V2(JJ)
end do
DNORM1 = sqrt(DNORM1)
DNORM2 = sqrt(DNORM2)
COSTH = (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))/(DNORM1*DNORM2)

! If the side is not that formed by the moving sphere the derivative of
! the polar angle vanishes

DCOS = 0.d0
if (NS2 /= NESFJ) goto 100

! Otherwise:

do JJ=1,3
  DCOS = DCOS+V2(JJ)*DP(L1,JJ)
end do
if (IOpt == 0) DCOS = DCOS+V1(IC)-Sphere(4,NS1)*COSTH*V2(IC)/DNORM2
DCOS = DCOS/(Sphere(4,NS1)*DNORM2)
100 continue
DA = COSTH*DPHI+acos(COSPHI)*DCOS
DA = Sphere(4,NS1)*Sphere(4,NS1)*DA

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nTs)

end subroutine DerPhi
!====
subroutine DerBet(IOpt,IC,NSFE,ITS,L0,L1,L2,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)

implicit real*8(A-H,O-Z)
parameter(MxVert=20)
dimension Sphere(4,*), ISphe(*), IntSph(MxVert,*)
dimension Vert(3,MxVert,*), Centr(3,MxVert,*)
dimension DP(MxVert,3)
dimension V1(3), V2(3), V3(3), V4(3), DV1(3), DV2(3), DV3(3), DV4(3)
dimension P1(3), P2(3), P3(3), T0(3), T1(3), T12(3), DT0(3), DT1(3)

! Trova la derivata dell'angolo esterno del vertice L1 della tessera
! ITS, formato dai lati L0-L1 e L1-L2,
!
! IOpt = 0 : rispetto alla coordinata IC (ex derbeta)
! IOpt = 1 : rispetto al raggio (ex derbeta1)
!
! della sfera NSFE, conoscendo la derivata della posizione del vertice P1.

PI = 4.0d0*atan(1.0d0)
NS1 = ISPHE(ITS)
NS2 = INTSPH(L0,ITs)
NS3 = INTSPH(L1,ITs)

! Trova i vettori posizione dei vertici L0, L1 rispetto al centro
! dell'arco L0, e dei vertici L1, L2 rispetto al centro dell'arco L1
do JJ=1,3
  V1(JJ) = VERT(JJ,L0,ITs)-CENTR(JJ,L0,ITs)
  V2(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L0,ITs)
  V3(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L1,ITs)
  V4(JJ) = VERT(JJ,L2,ITs)-CENTR(JJ,L1,ITs)
end do
! Trova la derivata dei vettori posizione
do JJ=1,3
  DV1(JJ) = DP(L0,JJ)
  DV2(JJ) = DP(L1,JJ)
  DV3(JJ) = DP(L1,JJ)
  DV4(JJ) = DP(L2,JJ)
end do

! Corregge la derivata dei vettori posizione se il lato considerato
! appartiene alla sfera che si muove (distinguendo se ITS appartiene
! o no alla sfera che si muove)
INDEX = 0
if ((NS1 == NSFE) .and. (NS2 /= NSFE)) INDEX = 1
if ((NS1 /= NSFE) .and. (NS2 == NSFE)) INDEX = 1
if (INDEX == 1) then
  T12(1) = Sphere(1,NS2)-Sphere(1,NS1)
  T12(2) = Sphere(2,NS2)-Sphere(2,NS1)
  T12(3) = Sphere(3,NS2)-Sphere(3,NS1)
  DIST2 = T12(1)*T12(1)+T12(2)*T12(2)+T12(3)*T12(3)
  if (IOpt == 0) then
    do JJ=1,3
      DV1(JJ) = DV1(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS2)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
      DV2(JJ) = DV2(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS2)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
    end do
    FACT = (Sphere(4,NS1)**2-Sphere(4,NS2)**2+DIST2)/(2.d0*DIST2)
    DV1(IC) = DV1(IC)-FACT
    DV2(IC) = DV2(IC)-FACT
  else if (IOpt == 1) then
    do JJ=1,3
      DV1(JJ) = DV1(JJ)+Sphere(4,NS2)*T12(JJ)/DIST2
      DV2(JJ) = DV2(JJ)+Sphere(4,NS2)*T12(JJ)/DIST2
    end do
  else
    write(6,'(a)') 'Illegal IOpt in DerBet.'
    call Abend()
  end if
end if

INDEX = 0
if ((NS1 == NSFE) .and. (NS3 /= NSFE)) INDEX = 1
if ((NS1 /= NSFE) .and. (NS3 == NSFE)) INDEX = 1
if (INDEX == 1) then
  T12(1) = Sphere(1,NS3)-Sphere(1,NS1)
  T12(2) = Sphere(2,NS3)-Sphere(2,NS1)
  T12(3) = Sphere(3,NS3)-Sphere(3,NS1)
  DIST2 = T12(1)*T12(1)+T12(2)*T12(2)+T12(3)*T12(3)
  if (IOpt == 0) then
    do JJ=1,3
      DV3(JJ) = DV3(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS3)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
      DV4(JJ) = DV4(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS3)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
    end do
    FACT = (Sphere(4,NS1)**2-Sphere(4,NS3)**2+DIST2)/(2.d0*DIST2)
    DV3(IC) = DV3(IC)-FACT
    DV4(IC) = DV4(IC)-FACT
  else if (IOpt == 1) then
    do JJ=1,3
      DV3(JJ) = DV3(JJ)+Sphere(4,NS3)*T12(JJ)/DIST2
      DV4(JJ) = DV4(JJ)+Sphere(4,NS3)*T12(JJ)/DIST2
    end do
  else
    write(6,'(a)') 'Illegal IOpt in DerBet.'
    call Abend()
  end if
end if

! Trova il vettore tangente al lato L0-L1: T0 = V2 x (V2 x V1)
do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = V1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  T0(JJ) = P3(JJ)
end do
dnxxt0 = DNORM3
! Trova il vettore tangente al lato L1-L2: T1 = V3 x (V3 x V4)
do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = V4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  T1(JJ) = P3(JJ)
end do
dnxxt1 = DNORM3

! Trova l'angolo Beta(L1)
PROD = T0(1)*T1(1)+T0(2)*T1(2)+T0(3)*T1(3)
BETA = PI-acos(PROD/(dnxxt0*dnxxt1))
COSBETA = cos(BETA)
SENBETA = sin(BETA)

! Trova la derivata della tangente T0 :
! DT0 = DV2 x (V2 x V1) + V2 x (DV2 x V1) + V2 x (V2 x DV1)

do JJ=1,3
  DT0(JJ) = 0.d0
end do

do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = V1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = DV2(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT0(JJ) = DT0(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = DV2(JJ)
  P2(JJ) = V1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT0(JJ) = DT0(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = DV1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT0(JJ) = DT0(JJ)+P3(JJ)
end do

! Trova la derivata della tangente T1 :
! DT1 = DV3 x (V3 x V4) + V3 x (DV3 x V4) + V3 x (V3 x DV4)

do JJ=1,3
  DT1(JJ) = 0.d0
end do

do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = V4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = DV3(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT1(JJ) = DT1(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = DV3(JJ)
  P2(JJ) = V4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT1(JJ) = DT1(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = DV4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT1(JJ) = DT1(JJ)+P3(JJ)
end do

! Infine calcola la derivata dell'angolo Beta(L1)
do JJ=1,3
  P1(JJ) = T1(JJ)+COSBETA*dnxxt1*T0(JJ)/dnxxt0
  P2(JJ) = T0(JJ)+COSBETA*dnxxt0*T1(JJ)/dnxxt1
end do
do JJ=1,3
  DBETA = 0.d0
end do
do JJ=1,3
  DBETA = DBETA+(P1(JJ)*DT0(JJ)+P2(JJ)*DT1(JJ))
end do
DBETA = DBETA/(SENBETA*dnxxt0*dnxxt1)
DA = Sphere(4,NS1)*Sphere(4,NS1)*DBETA

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nTs)

end subroutine DerBet
!====
function iRToInt(r)

implicit real*8(a-h,o-z)

! floating to integer conversion.

iRToInt = int(r)

return

end function iRToInt
