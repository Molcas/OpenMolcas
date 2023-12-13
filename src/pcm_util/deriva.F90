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

subroutine Deriva(IPrint,NAt,nTs,NEsf,NEsfP,RSolv,Tessera,Vert,Centr,Sphere,ISphe,IntSph,NOrd,NVert,NewSph,DerTes,DerPunt,DerRad, &
                  DerCentr)

use PCM_Arrays, only: MxVert
use Constants, only: Zero, One, Two, Ten, Half, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IPrint, NAt, nTs, NEsf, NEsfP, ISphe(nTs), IntSph(MxVert,nTs), NOrd(NEsfP), NVert(nTs), &
                                 NewSph(2,NEsf)
real(kind=wp), intent(in) :: RSolv
real(kind=wp), intent(inout) :: Tessera(4,nTs), Vert(3,MxVert,nTs), Centr(3,MxVert,nTs), Sphere(4,NEsf)
real(kind=wp), intent(out) :: DerTes(nTs,NAt,3), DerPunt(nTs,NAt,3,3), DerRad(NEsf,NAt,3), DerCentr(NEsf,NAt,3,3)
integer(kind=iwp) :: AlGe(63), Casca(MxVert), I, ICONT, ICoord, idx, II, IJunk, IMAX, IMIN, ITS, JJ, K, LIVEL, LL, N, NESFJ, NS1, &
                     NS2, NSA, NSI, NSJ, NSK, NSUB, NUM, NV
real(kind=wp) :: DDR, DDX, DDY, DDZ, DER, DERP(3), DERPT(3), DERS, DERTS, DISTKI(3), DR, DX, DY, DZ, FAC, FACT, ToBohr

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

! Initialize the derivative arrays
DerTes(:,:,:) = Zero
DerPunt(:,:,:,:) = Zero
DerRad(:,:,:) = Zero
DerCentr(:,:,:,:) = Zero
! The geometric quantities are expected to be in Angstrom
Tessera(1:3,:) = Angstrom*Tessera(1:3,:)
Tessera(4,:) = Angstrom*Angstrom*Tessera(4,:)
Sphere(:,:) = Angstrom*Sphere(:,:)
Vert(:,:,:) = Angstrom*Vert(:,:,:)
Centr(:,:,:) = Angstrom*Centr(:,:,:)

outer: do NSJ=1,NAt
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

    DERRAD(1:NESFP,NSJ,ICOORD) = ZERO
    DERCENTR(1:NESFP,NSJ,ICOORD,:) = ZERO
    if (NESFJ /= 0) DERCENTR(NESFJ,NSJ,ICOORD,ICOORD) = One

    ! Se non c'e' nessuna sfera attorno all'atomo NSJ ...
    if (NESFJ == 0) cycle outer

    ! 1) Effetti diretti.
    ! Loop sulle tessere.
    do ITS=1,NTS
      DERTS = ZERO
      DERPT(:) = ZERO
      NV = NVERT(ITS)
      NSI = ISPHE(ITS)
      if (NSI == NESFJ) then
        ! Derivate nel caso I = J

        ! Loop sulle sfere K che intersecano I
        do NSK=1,NESF
          if (NSK == NSI) cycle
          ! ITS ha un lato su K ?
          ICONT = 0
          do N=1,NV
            if (INTSPH(N,ITs) == NSK) ICONT = ICONT+1
          end do
          if (ICONT == 0) cycle
          call dSd(0,ICOORD,NSK,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
          DERTS = DERTS-DERS
          DERPT(:) = DERPT(:)-DERP(:)
        end do
      else
        ! Derivate nel caso I non = J
        !
        !           dS/dx(J), dQ/dx(J)
        ! ITS ha un lato su J ?
        ICONT = 0
        do N=1,NV
          if (INTSPH(N,ITs) == NESFJ) ICONT = ICONT+1
        end do
        if (ICONT >= 1) then
          call dSd(0,ICOORD,NESFJ,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
          DERTS = DERTS+DERS
          DERPT(:) = DERPT(:)+DERP(:)
        end if
      end if
      DERTES(ITS,NSJ,ICOORD) = DERTS
      DERPUNT(ITS,NSJ,ICOORD,:) = DERPT(:)
    end do

    ! 2) Effetti indiretti.
    ! Loop sulle sfere aggiunte
    do NSA=NESFP+1,NESF
      ALGE(:) = 0
      ! Costruiamo l'"albero genealogico" della sfera NSA
      ALGE(1) = NSA
      ALGE(2) = abs(NEWSPH(1,NSA))
      ALGE(3) = abs(NEWSPH(2,NSA))
      LIVEL = 3
      NUM = 2
      do while (NUM < 32)
        NSUB = 1
        do II=LIVEL-NUM+1,LIVEL
          if (ALGE(II) > NESFP) then
            ALGE(LIVEL+NSUB) = abs(NEWSPH(1,ALGE(II)))
            ALGE(LIVEL+NSUB+1) = abs(NEWSPH(2,ALGE(II)))
          end if
          NSUB = NSUB+2
        end do
        NUM = NUM*2
        LIVEL = LIVEL+NUM
      end do
      ! Si accerta che nell'ultimo livello ci siano solo sfere originarie
      do II=34,63
        if (ALGE(II) > NESFP) then
          write(u6,10) NSA,ALGE(II)
        end if
      end do
      ! Quando un elemento di ALGE e' = NESFJ, costruisce la corrispondente
      ! "cascata" di sfere aggiunte che collega NESFJ a NSA
      do LIVEL=2,6
        IMIN = 2**(LIVEL-1)
        IMAX = 2**(LIVEL)-1
        do II=IMIN,IMAX
          if (ALGE(II) /= NESFJ) cycle
          CASCA(:) = 0
          CASCA(1) = NESFJ
          idx = II
          K = 2
          do LL=LIVEL,2,-1
            FACT = (idx-2**(LL-1))*Half
            idx = 2**(LL-2)+int(FACT)
            CASCA(K) = ALGE(idx)
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
          call dRdC(NS2,ICOORD,NS1,DR,RSolv,Sphere,NewSph(:,NS2))
          call dCdC(1,ICOORD,NS1,DX,Sphere,NewSph(:,NS2))
          call dCdC(2,ICOORD,NS1,DY,Sphere,NewSph(:,NS2))
          call dCdC(3,ICOORD,NS1,DZ,Sphere,NewSph(:,NS2))
          do I=3,ICONT
            DDR = ZERO
            DDX = ZERO
            DDY = ZERO
            DDZ = ZERO
            NS1 = CASCA(I-1)
            NS2 = CASCA(I)

            ! Derivata del raggio dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dRdR(NS2,NS1,DER,RSolv,Sphere,NewSph(:,NS2))
            DDR = DER*DR
            call dRdC(NS2,1,NS1,DER,RSolv,Sphere,NewSph(:,NS2))
            DDR = DDR+DER*DX
            call dRdC(NS2,2,NS1,DER,RSolv,Sphere,NewSph(:,NS2))
            DDR = DDR+DER*DY
            call dRdC(NS2,3,NS1,DER,RSolv,Sphere,NewSph(:,NS2))
            DDR = DDR+DER*DZ

            ! Derivata della coord. X dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dCdR(1,NS1,DER,Sphere,NewSph(:,NS2))
            DDX = DER*DR
            call dCdC(1,1,NS1,DER,Sphere,NewSph(:,NS2))
            DDX = DDX+DER*DX
            call dCdC(1,2,NS1,DER,Sphere,NewSph(:,NS2))
            DDX = DDX+DER*DY
            call dCdC(1,3,NS1,DER,Sphere,NewSph(:,NS2))
            DDX = DDX+DER*DZ

            ! Derivata della coord. Y dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dCdR(2,NS1,DER,Sphere,NewSph(:,NS2))
            DDY = DER*DR
            call dCdC(2,1,NS1,DER,Sphere,NewSph(:,NS2))
            DDY = DDY+DER*DX
            call dCdC(2,2,NS1,DER,Sphere,NewSph(:,NS2))
            DDY = DDY+DER*DY
            call dCdC(2,3,NS1,DER,Sphere,NewSph(:,NS2))
            DDY = DDY+DER*DZ

            ! Derivata della coord. Z dell'elemento I di CASCA rispetto
            ! alla coord. ICOORD dell'elemento 1 di CASCA
            call dCdR(3,NS1,DER,Sphere,NewSph(:,NS2))
            DDZ = DER*DR
            call dCdC(3,1,NS1,DER,Sphere,NewSph(:,NS2))
            DDZ = DDZ+DER*DX
            call dCdC(3,2,NS1,DER,Sphere,NewSph(:,NS2))
            DDZ = DDZ+DER*DY
            call dCdC(3,3,NS1,DER,Sphere,NewSph(:,NS2))
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
            DERPT(:) = ZERO
            NV = NVERT(ITS)
            NSI = ISPHE(ITS)
            ! Derivate delle tessere appartenenti a NSA
            if (NSI == NSA) then
              ! Derivate relative al raggio di NSI: 1) scaling
              DERTS = Two*Tessera(4,ITS)*DR/Sphere(4,NSI)
              DERPT(1) = (Tessera(1,ITS)-Sphere(1,NSI))*DR/Sphere(4,NSI)
              DERPT(2) = (Tessera(2,ITS)-Sphere(2,NSI))*DR/Sphere(4,NSI)
              DERPT(3) = (Tessera(3,ITS)-Sphere(3,NSI))*DR/Sphere(4,NSI)
              ! Loop sulle altre sfere K che tagliano ITS
              do NSK=1,NESF
                if (NSK == NSI) cycle
                ICONT = 0
                do N=1,NV
                  if (INTSPH(N,ITs) == NSK) ICONT = ICONT+1
                end do
                if (ICONT == 0) cycle
                ! Derivate relative al raggio di NSI: 2) raggio delle altre sfere
                call dSd(1,IJunk,NSK,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
                DERTS = DERTS-DERS*Sphere(4,NSK)*DR/Sphere(4,NSI)
                DERPT(:) = DERPT(:)-DERP(:)*Sphere(4,NSK)*DR/Sphere(4,NSI)
                ! Derivate relative al raggio di NSI: 3) coord. delle altre sfere
                DISTKI(1) = Sphere(1,NSK)-Sphere(1,NSI)
                DISTKI(2) = Sphere(2,NSK)-Sphere(2,NSI)
                DISTKI(3) = Sphere(3,NSK)-Sphere(3,NSI)
                FAC = DR/Sphere(4,NSI)
                do JJ=1,3
                  call dSd(0,JJ,NSK,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
                  if (DISTKI(JJ) /= ZERO) then
                    DERTS = DERTS-FAC*DERS*DISTKI(JJ)
                    DERPT(:) = DERPT(:)-FAC*DERP(:)*DISTKI(JJ)
                  end if
                end do
                ! Derivate relative alle coordinate di NSI
                call dSd(0,1,NSK,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
                DERTS = DERTS-DX*DERS
                DERPT(:) = DERPT(:)-DX*DERP(:)
                call dSd(0,2,NSK,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
                DERTS = DERTS-DY*DERS
                DERPT(:) = DERPT(:)-DY*DERP(:)
                call dSd(0,3,NSK,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
                DERTS = DERTS-DZ*DERS
                DERPT(:) = DERPT(:)-DZ*DERP(:)
              end do
              DERTES(ITS,NSJ,ICOORD) = DERTES(ITS,NSJ,ICOORD)+DERTS
              DERPUNT(ITS,NSJ,ICOORD,:) = DERPUNT(ITS,NSJ,ICOORD,:)+DERPT(:)
            else
              ! Derivate delle tessere tagliate da NSA
              ICONT = 0
              do N=1,NV
                if (INTSPH(N,ITs) == NSA) ICONT = ICONT+1
              end do
              if (ICONT == 0) cycle
              ! Derivate relative al raggio di NSA
              call dSd(1,IJunk,NSA,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
              DERTS = DERS*DR
              DERPT(:) = DERP(:)*DR
              ! Derivate relative alle coordinate di NSA
              call dSd(0,1,NSA,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
              DERTS = DERTS+DX*DERS
              DERPT(:) = DERPT(:)+DX*DERP(:)
              call dSd(0,2,NSA,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
              DERTS = DERTS+DY*DERS
              DERPT(:) = DERPT(:)+DY*DERP(:)
              call dSd(0,3,NSA,DERS,DERP,Tessera(:,ITS),Vert(:,:,ITS),Centr(:,:,ITS),Sphere,NSI,IntSph(:,ITS),NV)
              DERTS = DERTS+DZ*DERS
              DERPT(:) = DERPT(:)+DZ*DERP(:)

              DERTES(ITS,NSJ,ICOORD) = DERTES(ITS,NSJ,ICOORD)+DERTS
              DERPUNT(ITS,NSJ,ICOORD,:) = DERPUNT(ITS,NSJ,ICOORD,:)+DERPT(:)
            end if
          end do
        end do
      end do
    end do

    ! Controlla che nessuna derivata sia eccessiva (per imprecisioni numeriche)
    do ITS=1,NTS
      if (abs(DERTES(ITS,NSJ,ICOORD)) > Ten) then
        if (IPrint >= 2) then
          write(u6,*) 'DERIVATA GEOMETRICA ECCESSIVA TRASCURATA'
          write(u6,*) 'TESSERA, SFERA, COORDINATA',ITS,NESFJ,ICOORD
        end if
        DERTES(ITS,NSJ,ICOORD) = Zero
      end if
    end do

    ! End of loop on atoms and coordinates

  end do
end do outer
! Put the geometric quantities in Bohr again
ToBohr = One/Angstrom
Tessera(1:3,:) = ToBohr*Tessera(1:3,:)
Tessera(4,:) = ToBohr*ToBohr*Tessera(4,:)
Sphere(:,:) = ToBohr*Sphere(:,:)
Vert(:,:,:) = ToBohr*Vert(:,:,:)
Centr(:,:,:) = ToBohr*Centr(:,:,:)

return

10 format(/,'SFERA ',I3,/,'AL LIVELLO 5 LA SFERA',I3," E' AGGIUNTA")

end subroutine Deriva
