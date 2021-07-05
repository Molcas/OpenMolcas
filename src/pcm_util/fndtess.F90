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

subroutine FndTess(iPrint,ToAng,LcNAtm,Xs,Ys,Zs,Rs,pNs,nn)

use PCM_arrays, only: PCMSph, PCMTess, Vert, Centr, SSph, PCMDM, PCM_N, PCMiSph, NVert, IntSph, NewSph

implicit real*8(A-H,O-Z)
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "status.fh"
real*8 Xs(nn), Ys(nn), Zs(nn), Rs(nn)
integer pNs(nn)
real*8, allocatable :: Xt(:), Yt(:), Zt(:), At(:), pVert(:), pCentr(:), pSSph(:), CV(:)
integer, allocatable :: pNewS(:), pIntS(:), pNVert(:), pISph(:), JTR(:)

! Definition of solute cavity and computation of vertices,
! representative points and surfaces of the tesserae by the
! Gauss-Bonnet Theorem.

! Allocate space for X, Y, Z, Area, ISPHE (index of sphere to which
! tessera belongs); then allocate space for IntSph (indices of spheres
! cutting the tessera), NewSph (indices of spheres creating new smoothing
! spheres), SSph (surface of each sphere exposed to the solvent),
! Finally, allocate temporary space for NVert (number of vertices for
! any tessera), Vert (coordinates of vertices), Centr (center of arcs).

call mma_allocate(Xt,MxTs,Label='Xt')
call mma_allocate(Yt,MxTs,Label='Yt')
call mma_allocate(Zt,MxTs,Label='Zt')
call mma_allocate(At,MxTs,Label='At')
call mma_allocate(pVert,3*MxVert*MxTs,Label='pVert')
call mma_allocate(pCentr,3*MxVert*MxTs,Label='pCentr')
call mma_allocate(pSSph,MxSph,Label='pSSph')
call mma_allocate(CV,3000,Label='CV')

call mma_allocate(JTR,3*MxTs,Label='JTR')
call mma_allocate(pISph,MxTs,Label='pIShp')
call mma_allocate(pNVert,MxTs,Label='pNVert')
call mma_allocate(pIntS,MxVert*MxTs,Label='pIntS')
call mma_allocate(pNewS,2*MxSph,Label='pNewS')

Omega = RSlPar(2)
Ret = RSlPar(3)
Fro = RSlPar(4)
TsAre = RSlPar(7)
RSolv = RSlPar(19)
ITsNum = ISlPar(11)

call FndTess_(iPrint,ToAng,LcNAtm,MxSph,MxTs,Xs,Ys,Zs,Rs,Ret,Omega,Fro,RSolv,NSinit,NS,ITsNum,TsAre,nTs,Xt,Yt,Zt,At,pISph,pNVert, &
              pVert,pCentr,pIntS,pNewS,pSSph,JTR,CV)
!                                                                      *
!***********************************************************************
!     Re-allocate with actual dimensioning                             *
if (RctFld_Status /= Active) then

  ! Allocate PCM arrays
  call mma_allocate(PCMSph,4,NS,Label='PCMSph')
  call mma_allocate(PCMTess,4,nTs,Label='PCMTess')
  call mma_allocate(Vert,3,MxVert,nTs,Label='Vert')
  call mma_allocate(Centr,3,MxVert,nTs,Label='Centr')
  call mma_allocate(SSph,NS,Label='SSph')
  call mma_allocate(PCMDM,nTs,nTs,Label='PCMDM')
  call mma_allocate(PCM_N,NS,Label='PCM_N')
  call mma_allocate(PCMiSph,nTs,Label='PCMiSph')
  call mma_allocate(NVert,nTs,Label='NVert')
  call mma_allocate(IntSph,MxVert,nTs,Label='IntSph')
  call mma_allocate(NewSph,2,NS,Label='NewSph')

  nPCM_info_r = 4*NS+4*nTs+3*MxVert*nTs+3*MxVert*nTs+NS+nTs**2
  nPCM_info_i = NS+nTs+nTs+MxVert*nTs+2*NS+nTs**2
  nPCM_info = nPCM_info_r+nPCM_info_i

end if

! PCMSph
do iS=1,NS
  PCMSph(1,iS) = Xs(iS)
  PCMSph(2,iS) = Ys(iS)
  PCMSph(3,iS) = Zs(iS)
  PCMSph(4,iS) = Rs(iS)
end do

! PCMTess
do iTs=1,nTs
  PCMTess(1,iTs) = Xt(iTs)
  PCMTess(2,iTs) = Yt(iTs)
  PCMTess(3,iTs) = Zt(iTs)
  PCMTess(4,iTs) = At(iTs)
end do

! Vert
call dcopy_(3*MxVert*nTs,pVert,1,Vert,1)

! Centr
call dcopy_(3*MxVert*nTs,pCentr,1,Centr,1)

! SSph
call dcopy_(nS,pSSph,1,SSph,1)

! nOrd
call ICopy(nS,pNs,1,PCM_N,1)

! ISph
call ICopy(nTs,pISph,1,PCMiSph,1)

! NVert
call ICopy(nTs,pNVert,1,NVert,1)

! IntSph
call ICopy(MxVert*nTs,pIntS,1,IntSph,1)

! NewSph
call ICopy(2*NS,pNewS,1,NewSph,1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate temporary memory

call mma_deallocate(pNewS)
call mma_deallocate(pIntS)
call mma_deallocate(pNVert)
call mma_deallocate(pISph)
call mma_deallocate(JTR)
call mma_deallocate(CV)
call mma_deallocate(pSSPh)
call mma_deallocate(pCentr)
call mma_deallocate(pVert)
call mma_deallocate(At)
call mma_deallocate(Zt)
call mma_deallocate(Yt)
call mma_deallocate(Xt)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine FndTess
!====
subroutine FndTess_(IPrint,ToAng,NAT,MxSph,MxTs,XE,YE,ZE,RE,RET,Omega,FRO,RSolv,NEsfP,NEsf,ITsNum,TSAre,NTS,XCTs,YCTs,ZCTs,AS, &
                    ISphe,Nvert,Vert,Centr,IntSph,NewSph,SSfe,JTR,CV)
! Definition of solute cavity and computation of vertices,
! representative points and surfaces of the tesserae by the
! Gauss-Bonnet Theorem.

implicit real*8(A-H,O-Z)
parameter(MxVert=20)
dimension XCTs(*), YCTs(*), ZCTs(*), AS(*)
dimension ISphe(*)
dimension XE(*), YE(*), ZE(*), RE(*)
dimension SSFE(*), IntSph(MxVert,*), NVert(*)
dimension NewSph(2,*)
dimension Vert(3,MxVert,*), Centr(3,MxVert,*)
dimension JTR(3,*), CV(3,*)
dimension PP(3), PTS(3,MxVert), CCC(3,MxVert)
save Zero, First
data ZERO,FIRST/0.0d0,0.0174533d0/

! PEDRA works with Angstroms
do ISFE=1,NESFP
  XE(ISFE) = XE(ISFE)*ToAng
  YE(ISFE) = YE(ISFE)*ToAng
  ZE(ISFE) = ZE(ISFE)*ToAng
end do
NESF = NESFP
if (IPRINT == 2) write(6,800)

! creation of new spheres

do N=1,NESF
  NEWSPH(1,N) = 0
  NEWSPH(2,N) = 0
end do
ITYPC = 0
OMEGA = OMEGA*FIRST
SENOM = sin(OMEGA)
COSOM2 = (cos(OMEGA))**2
RTDD = RET+RSOLV
RTDD2 = RTDD*RTDD
NET = NESF
NN = 2
NE = NESF
NEV = NESF
GO TO 100
110 continue
NN = NE+1
NE = NET
if (NE > MxSph) then
  write(6,1111)
  call Abend()
end if
100 continue
do I=NN,NE
  NES = I-1
  do J=1,NES
    RIJ2 = (XE(I)-XE(J))**2+(YE(I)-YE(J))**2+(ZE(I)-ZE(J))**2
    RIJ = sqrt(RIJ2)
    RJD = RE(J)+RSOLV
    TEST1 = RE(I)+RJD+RSOLV
    if (RIJ >= TEST1) GO TO 130
    REG = max(RE(I),RE(J))
    REP = min(RE(I),RE(J))
    REG2 = REG*REG
    REP2 = REP*REP
    TEST2 = REP*SENOM+sqrt(REG2-REP2*COSOM2)
    if (RIJ <= TEST2) GO TO 130
    REGD2 = (REG+RSOLV)*(REG+RSOLV)
    TEST3 = (REGD2+REG2-RTDD2)/REG
    if (RIJ >= TEST3) GO TO 130
    do K=1,NEV
      if ((K == J) .or. (K == I)) GO TO 140
      RJK2 = (XE(J)-XE(K))**2+(YE(J)-YE(K))**2+(ZE(J)-ZE(K))**2
      if (RJK2 >= RIJ2) GO TO 140
      RIK2 = (XE(I)-XE(K))**2+(YE(I)-YE(K))**2+(ZE(I)-ZE(K))**2
      if (RIK2 >= RIJ2) GO TO 140
      RJK = sqrt(RJK2)
      RIK = sqrt(RIK2)
      SP = (RIJ+RJK+RIK)/2.0d0
      HH = 4*(SP*(SP-RIJ)*(SP-RIK)*(SP-RJK))/RIJ2
      REO = RE(K)*FRO
      if (K >= NE) REO = 0.0002d0
      REO2 = REO*REO
      if (HH < REO2) GO TO 130
140   continue
    end do
    REPD2 = (REP+RSOLV)**2
    TEST8 = sqrt(REPD2-RTDD2)+sqrt(REGD2-RTDD2)
    if (RIJ <= TEST8) GO TO 150
    REND2 = REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
    if (REND2 <= RTDD2) GO TO 130
    REN = sqrt(REND2)-RSOLV
    FC = REG/(RIJ-REG)
    TEST7 = REG-RE(I)
    KG = I
    KP = J
    if (TEST7 <= 0.000000001d0) GO TO 160
    KG = J
    KP = I
160 continue
    FC1 = FC+1.0d0
    XEN = (XE(KG)+FC*XE(KP))/FC1
    YEN = (YE(KG)+FC*YE(KP))/FC1
    ZEN = (ZE(KG)+FC*ZE(KP))/FC1
    ITYPC = 1
    GO TO 170
150 continue
    R2GN = RIJ-REP+REG
    RGN = R2GN/2.0d0
    FC = R2GN/(RIJ+REP-REG)
    FC1 = FC+1.0d0
    TEST7 = REG-RE(I)
    KG = I
    KP = J
    if (TEST7 <= 0.000000001d0) GO TO 180
    KG = J
    KP = I
180 continue
    XEN = (XE(KG)+FC*XE(KP))/FC1
    YEN = (YE(KG)+FC*YE(KP))/FC1
    ZEN = (ZE(KG)+FC*ZE(KP))/FC1
    REN = sqrt(REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ))-RSOLV
170 continue
    NET = NET+1
    XE(NET) = XEN
    YE(NET) = YEN
    ZE(NET) = ZEN
    RE(NET) = REN

    ! Nella matrice NEWSPH(2,NESF) sono memorizzati i numeri delle
    ! sfere "generatrici" della nuova sfera NET: se la nuova sfera e'
    ! del tipo A o B entrambi i numeri sono positivi, se e' di tipo
    ! C il numero della sfera "principale" e' negativo
    ! (per la definizione del tipo si veda JCC 11, 1047 (1990))

    if (ITYPC == 0) then
      NEWSPH(1,NET) = KG
      NEWSPH(2,NET) = KP
    elseif (ITYPC == 1) then
      NEWSPH(1,NET) = -KG
      NEWSPH(2,NET) = KP
    end if

130 continue
  end do
  NEV = NET
end do
if (NET /= NE) GO TO 110
NESF = NET

! Division of the surface into tesserae

VCav = ZERO
Scav = ZERO

! Controlla se ciascuna tessera e' scoperta o va tagliata

NN1 = 0
do NSFE=1,NESF
  XEN = XE(NSFE)
  YEN = YE(NSFE)
  ZEN = ZE(NSFE)
  REN = RE(NSFE)
  if ((ITsNum == 0) .and. (TsAre == 0.d0)) then
    IPtype = 2
    IPFlag = 0
    ITsNum = 60
  elseif ((ITsNum > 0) .and. (TsAre == 0.d0)) then
    IPFlag = 0
  elseif (TsAre > 0.d0) then
    IPFlag = 1
  end if
  call PolyGen(MxTs,IPtype,IPflag,TsAre,ITsNum,XEN,YEN,ZEN,REN,ITsEff,CV,JTR)
  do ITS=1,ITsEff
    N1 = JTR(1,ITS)
    N2 = JTR(2,ITS)
    N3 = JTR(3,ITS)
    PTS(1,1) = CV(1,N1)
    PTS(2,1) = CV(2,N1)
    PTS(3,1) = CV(3,N1)
    PTS(1,2) = CV(1,N2)
    PTS(2,2) = CV(2,N2)
    PTS(3,2) = CV(3,N2)
    PTS(1,3) = CV(1,N3)
    PTS(2,3) = CV(2,N3)
    PTS(3,3) = CV(3,N3)
    NV = 3
    do JJ=1,3
      PP(JJ) = ZERO
    end do

    ! Per ciascuna tessera, trova la porzione scoperta e ne
    ! calcola l'area con il teorema di Gauss-Bonnet; il punto
    ! rappresentativo e' definito come media dei vertici della porzione
    ! scoperta di tessera e passato in PP.
    ! I vertici di ciascuna tessera sono conservati in
    ! VERT(3,MxVert,MxTs), il numero di vertici di ciascuna tessera e'
    ! in NVERT(MxTs), e i centri dei cerchi di ciascun lato sono in
    ! CENTR(3,MxVert,MxTs). In INTSPH(MxVert,MxTs) sono registrate le
    ! sfere a cui appartengono i lati delle tessere.

    call TESSERA(iPrint,MxTs,Nesf,NSFE,NV,XE,YE,ZE,RE,IntSph,PTS,CCC,PP,AREA)
    if (AREA == 0.d0) goto 310
    NN1 = NN1+1
    NN = min(NN1,MxTs)
    XCTS(NN) = PP(1)
    YCTS(NN) = PP(2)
    ZCTS(NN) = PP(3)
    AS(NN) = AREA

    ISPHE(NN) = NSFE
    NVERT(NN) = NV
    do IV=1,NV
      do JJ=1,3
        VERT(JJ,IV,NN) = PTS(JJ,IV)
        CENTR(JJ,IV,NN) = CCC(JJ,IV)
      end do
    end do
    do IV=1,NV
      INTSPH(IV,NN) = INTSPH(IV,MxTs)
    end do
310 continue
  end do
end do
NTS = NN

! Verifica se due tessere sono troppo vicine
TEST = 0.02d0
TEST2 = TEST*TEST
do I=1,NTS-1
  if (AS(I) == ZERO) goto 400
  XI = XCTS(I)
  YI = YCTS(I)
  ZI = ZCTS(I)
  II = I+1
  do J=II,NTS
    if (ISPHE(I) == ISPHE(J)) goto 410
    if (AS(J) == ZERO) goto 410
    XJ = XCTS(J)
    YJ = YCTS(J)
    ZJ = ZCTS(J)
    RIJ = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
    if (RIJ > TEST2) goto 410

    ! La routine originaria sostituiva le due tessere troppo vicine con una
    ! sola tessera. Nel caso Gauss-Bonnet, anche i vertici delle tessere
    ! e i centri degli archi vengono memorizzati ed e' impossibile sostituirli
    ! nello stesso modo: percio' l'area della tessera piu' piccola verra'
    ! trascurata per evitare problemi nella autopolarizzazione.
    if (IPRINT == 2) write(6,1000) I,J,TEST2
    if (AS(I) < AS(J)) AS(I) = ZERO
    if (AS(I) >= AS(J)) AS(J) = ZERO
410 continue
  end do
400 continue
end do

! E' preferibile eliminare del tutto le tessere che per
! qualche motivo hanno AREA = 0, ridefinendo tutti gli
! indici: l'errore numerico cosi' introdotto e' in genere
! trascurabile e si evitano problemi di convergenza

! Define here the number of tesserae in electrostatic calculations
! to avoid problems in successive calcn.
ITS = 0
450 continue
ITS = ITS+1
if (AS(ITS) < 1.D-10) then
  do I=ITS,NTS-1
    AS(I) = AS(I+1)
    XCTS(I) = XCTS(I+1)
    YCTS(I) = YCTS(I+1)
    ZCTS(I) = ZCTS(I+1)
    ISPHE(I) = ISPHE(I+1)
    NVERT(I) = NVERT(I+1)
    do IV=1,MxVert
      INTSPH(IV,I) = INTSPH(IV,I+1)
      do IC=1,3
        VERT(IC,IV,I) = VERT(IC,IV,I+1)
        CENTR(IC,IV,I) = CENTR(IC,IV,I+1)
      end do
    end do
  end do
  NTS = NTS-1
  ITS = ITS-1
end if
if (ITS < NTS) goto 450
!***********************************************************************
! Calcola il volume della cavita' con la formula (t. di Gauss):
!            V=SOMMAsulleTESSERE{A r*n}/3
! dove r e' la distanza del punto rappresentativo dall'origine,
! n e' il versore normale alla tessera, A l'area della tessera,
! e * indica il prodotto scalare.
!***********************************************************************
VCav = ZERO
do ITS=1,NTS
  NSFE = ISPHE(ITS)
  ! Trova il versore normale
  XN = (XCTS(ITS)-XE(NSFE))/RE(NSFE)
  YN = (YCTS(ITS)-YE(NSFE))/RE(NSFE)
  ZN = (ZCTS(ITS)-ZE(NSFE))/RE(NSFE)
  ! Trova il prodotto scalare
  PROD = XCTS(ITS)*XN+YCTS(ITS)*YN+ZCTS(ITS)*ZN
  VCav = VCav+AS(ITS)*PROD/3.d0
end do
!***********************************************************************
! Stampa la geometria della cavita'
Scav = ZERO
do I=1,NESF
  SSFE(I) = ZERO
end do
do I=1,NTS
  K = ISPHE(I)
  SSFE(K) = SSFE(K)+AS(I)
end do
OMEGA = OMEGA/FIRST
if (IPRINT == 2) write(6,1100) OMEGA,RSOLV,RET,FRO,NESF
do I=1,NESF
  if (IPRINT == 2) write(6,1200) I,XE(I),YE(I),ZE(I),RE(I),SSFE(I)
  Scav = Scav+SSFE(I)
end do
if (IPRINT == 2) write(6,1300) NTS,Scav,VCav

! Trasforma i risultati in bohr
do I=1,NESF
  RE(I) = RE(I)/ToAng
  XE(I) = XE(I)/ToAng
  YE(I) = YE(I)/ToAng
  ZE(I) = ZE(I)/ToAng
end do
do I=1,NTS
  do J=1,NVERT(I)
    do L=1,3
      VERT(L,J,I) = VERT(L,J,I)/ToAng
      CENTR(L,J,I) = CENTR(L,J,I)/ToAng
    end do
  end do
end do
do I=1,NTS
  AS(I) = AS(I)/(ToAng*ToAng)
  XCTS(I) = XCTS(I)/ToAng
  YCTS(I) = YCTS(I)/ToAng
  ZCTS(I) = ZCTS(I)/ToAng
end do
if (IPRINT == 3) then
  write(6,1500)
  write(6,1600)
  write(6,1700) (I,ISPHE(I),AS(I),XCTS(I),YCTS(I),ZCTS(I),I=1,NTS)
end if
if (NN1 > MxTs) then
  write(6,1240) NN1,MxTs
  write(6,1112)
  call Abend()
end if

return

800 format(/,'**** POLARISABLE CONTINUUM MODEL - UNIVERSITIES OF NAPLES AND PISA *****')
1000 format(/,'ATTENZIONE: I CENTRI DELLE TESSERE ',I4,',',I4,' DISTANO MENO DI',F8.6,' A',/)
1100 format(/,'** CHARACTERISTICS OF THE CAVITY **',//,'  GEOMETRICAL PARAMETERS: OMEGA=',F8.3,' RSOLV=',F8.3,' RET=',F8.3, &
            ' FRO=',F8.3,//,'  TOTAL NUMBER OF SPHERES',I5,//,'  CENTERS AND RADII :',//, &
            '  SPHERE       CENTER  (X,Y,Z) (A)     RADIUS (A)','       AREA(A*A)')
1111 format(' Cavity: too many spheres; increase Omega')
1112 format('Too many tesserae.')
!1200 format(I5,4F15.9,F18.9)
1200 format(I5,4F10.3,F18.3)
1240 format(' NN1=',I10,' but MxTs=',I10,'.')
1300 format(/,' TOTAL NUMBER OF TESSERAE ',I8,//,' SURFACE AREA',F14.8,'(A**2)',8X,'CAVITY VOLUME',F14.8,' (A**3)')
!1400 format(/,'SFERA',I4,' COORD.',I4,'  DERIVATA   ',F10.6,/)
1500 format('1 *** SUDDIVISIONE DELLA SUPERFICIE  ***')
1600 format(' TESSERA  SFERA   AREA   X Y Z CENTRO TESSERA  X Y Z PUNTO NORMALE')
1700 format(2I4,7F12.7)
!1800 format(/,'**** END OF CAVITY DEFINITION ****',/)
! Avoid unused argument warnings
if (.false.) call Unused_integer(NAT)

end subroutine FndTess_
