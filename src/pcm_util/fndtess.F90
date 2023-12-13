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

subroutine FndTess(iPrint,Xs,Ys,Zs,Rs,pNs,m)

use PCM_arrays, only: Centr, IntSph, MxSph, MxTs, MxVert, NewSph, NVert, PCM_N, PCMDM, PCMiSph, PCMSph, PCMTess, PCMTess, SSph, Vert
use rctfld_module, only: iSLPar, nPCM_Info, nS, nSInit, nTS, rSLPar, rSolv
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Half, Angstrom, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, m, pNs(m)
real(kind=wp), intent(inout) :: Xs(m), Ys(m), Zs(m), Rs(m)
integer(kind=iwp) :: I, II, IPFlag, IPtype, iTs, iTsNum, ITsEff, ITYPC, J, K, KG, KP, N1, N2, N3, NE, NES, NET, NEV, NN, NN1, &
                     nPCM_info_i, nPCM_info_r, NSFE, NV
real(kind=wp) :: AREA, COSOM2, FC, FC1, Fro, HH, Omega, PP(3), PROD, R2GN, REG, REG2, REGD2, REN, REND2, REO, REO2, REP, REP2, &
                 REPD2, Ret, RGN, RIJ, RIJ2, RIK, RIK2, RJD, RJK, RJK2, RTDD, RTDD2, Scav, SENOM, SP, TEST, TEST1, TEST2, TEST3, &
                 TEST7, TEST8, TsAre, VCav, XEN, XI, XJ, XN, YEN, YI, YJ, YN, ZEN, ZI, ZJ, ZN
logical(kind=iwp) :: FIRST
integer(kind=iwp), allocatable :: JTR(:,:), pIntS(:,:), pISph(:), pNewS(:,:), pNVert(:)
real(kind=wp), allocatable :: At(:), CCC(:,:), CV(:,:), pCentr(:,:,:), pSSph(:), PTS(:,:), pVert(:,:,:), Xt(:), Yt(:), Zt(:)

! Definition of solute cavity and computation of vertices,
! representative points and surfaces of the tesserae by the
! Gauss-Bonnet Theorem.

! Allocate space for X, Y, Z, Area, pISph (index of sphere to which
! tessera belongs); then allocate space for IntSph (indices of spheres
! cutting the tessera), pNewS (indices of spheres creating new smoothing
! spheres), SSph (surface of each sphere exposed to the solvent),
! Finally, allocate temporary space for NVert (number of vertices for
! any tessera), Vert (coordinates of vertices), Centr (center of arcs).

call mma_allocate(Xt,MxTs,Label='Xt')
call mma_allocate(Yt,MxTs,Label='Yt')
call mma_allocate(Zt,MxTs,Label='Zt')
call mma_allocate(At,MxTs,Label='At')
call mma_allocate(pVert,3,MxVert,MxTs,Label='pVert')
call mma_allocate(pCentr,3,MxVert,MxTs,Label='pCentr')
call mma_allocate(pSSph,MxSph,Label='pSSph')
call mma_allocate(CV,3,MxSph,Label='CV')

call mma_allocate(JTR,3,MxTs,Label='JTR')
call mma_allocate(pISph,MxTs,Label='pIShp')
call mma_allocate(pNVert,MxTs,Label='pNVert')
call mma_allocate(pIntS,MxVert,MxTs,Label='pIntS')
call mma_allocate(pNewS,2,MxSph,Label='pNewS')

Omega = RSlPar(2)
Ret = RSlPar(3)
Fro = RSlPar(4)
TsAre = RSlPar(7)
RSolv = RSlPar(19)
ITsNum = ISlPar(11)

! PEDRA works with Angstroms
NS = NSinit
Xs(1:NS) = Xs(1:NS)*Angstrom
Ys(1:NS) = Ys(1:NS)*Angstrom
Zs(1:NS) = Zs(1:NS)*Angstrom
if (IPRINT == 2) write(u6,800)

! creation of new spheres

pNewS(:,1:NS) = 0
ITYPC = 0
OMEGA = OMEGA*DEG2RAD
SENOM = sin(OMEGA)
COSOM2 = (cos(OMEGA))**2
RTDD = RET+RSOLV
RTDD2 = RTDD*RTDD
NET = NS
NN = 2
NE = NET-1 ! just to get the loop started
NEV = NS
FIRST = .true.
do while (NET /= NE)
  if (FIRST) then
    NE = NS
    FIRST = .false.
  else
    NN = NE+1
    NE = NET
    if (NE > MxSph) then
      write(u6,1111)
      call Abend()
    end if
  end if
  do I=NN,NE
    NES = I-1
    middle: do J=1,NES
      RIJ2 = (Xs(I)-Xs(J))**2+(Ys(I)-Ys(J))**2+(Zs(I)-Zs(J))**2
      RIJ = sqrt(RIJ2)
      RJD = Rs(J)+RSOLV
      TEST1 = Rs(I)+RJD+RSOLV
      if (RIJ >= TEST1) cycle middle
      REG = max(Rs(I),Rs(J))
      REP = min(Rs(I),Rs(J))
      REG2 = REG*REG
      REP2 = REP*REP
      TEST2 = REP*SENOM+sqrt(REG2-REP2*COSOM2)
      if (RIJ <= TEST2) cycle middle
      REGD2 = (REG+RSOLV)*(REG+RSOLV)
      TEST3 = (REGD2+REG2-RTDD2)/REG
      if (RIJ >= TEST3) cycle middle
      do K=1,NEV
        if ((K == J) .or. (K == I)) cycle
        RJK2 = (Xs(J)-Xs(K))**2+(Ys(J)-Ys(K))**2+(Zs(J)-Zs(K))**2
        if (RJK2 >= RIJ2) cycle
        RIK2 = (Xs(I)-Xs(K))**2+(Ys(I)-Ys(K))**2+(Zs(I)-Zs(K))**2
        if (RIK2 >= RIJ2) cycle
        RJK = sqrt(RJK2)
        RIK = sqrt(RIK2)
        SP = (RIJ+RJK+RIK)*Half
        HH = 4*(SP*(SP-RIJ)*(SP-RIK)*(SP-RJK))/RIJ2
        REO = Rs(K)*FRO
        if (K >= NE) REO = 2.0e-4_wp
        REO2 = REO*REO
        if (HH < REO2) cycle middle
      end do
      REPD2 = (REP+RSOLV)**2
      TEST8 = sqrt(REPD2-RTDD2)+sqrt(REGD2-RTDD2)
      if (RIJ <= TEST8) then
        R2GN = RIJ-REP+REG
        RGN = R2GN*Half
        FC = R2GN/(RIJ+REP-REG)
        REN = sqrt(REGD2+RGN*(RGN-(REGD2+RIJ2-REPD2)/RIJ))-RSOLV
      else
        REND2 = REGD2+REG2-(REG/RIJ)*(REGD2+RIJ2-REPD2)
        if (REND2 <= RTDD2) cycle middle
        REN = sqrt(REND2)-RSOLV
        FC = REG/(RIJ-REG)
        ITYPC = 1
      end if
      FC1 = FC+One
      TEST7 = REG-Rs(I)
      if (TEST7 <= 1.0e-9_wp) then
        KG = I
        KP = J
      else
        KG = J
        KP = I
      end if
      XEN = (Xs(KG)+FC*Xs(KP))/FC1
      YEN = (Ys(KG)+FC*Ys(KP))/FC1
      ZEN = (Zs(KG)+FC*Zs(KP))/FC1
      NET = NET+1
      Xs(NET) = XEN
      Ys(NET) = YEN
      Zs(NET) = ZEN
      Rs(NET) = REN

      ! Nella matrice pNewS(2,NS) sono memorizzati i numeri delle
      ! sfere "generatrici" della nuova sfera NET: se la nuova sfera e'
      ! del tipo A o B entrambi i numeri sono positivi, se e' di tipo
      ! C il numero della sfera "principale" e' negativo
      ! (per la definizione del tipo si veda JCC 11, 1047 (1990))

      if (ITYPC == 0) then
        pNewS(1,NET) = KG
        pNewS(2,NET) = KP
      else if (ITYPC == 1) then
        pNewS(1,NET) = -KG
        pNewS(2,NET) = KP
      end if

    end do middle
    NEV = NET
  end do
  if (NET == NE) exit
end do
NS = NET

! Division of the surface into tesserae

VCav = ZERO
Scav = ZERO

! Controlla se ciascuna tessera e' scoperta o va tagliata

call mma_allocate(CCC,3,MxVert,Label='CCC')
call mma_allocate(PTS,3,MxVert,Label='PTS')

NN1 = 0
do NSFE=1,NS
  XEN = Xs(NSFE)
  YEN = Ys(NSFE)
  ZEN = Zs(NSFE)
  REN = Rs(NSFE)
  if ((ITsNum == 0) .and. (TsAre == Zero)) then
    IPtype = 2
    IPFlag = 0
    ITsNum = 60
  else if ((ITsNum > 0) .and. (TsAre == Zero)) then
    IPFlag = 0
  else if (TsAre > Zero) then
    IPFlag = 1
  end if
  call PolyGen(MxTs,MxSph,IPtype,IPflag,TsAre,ITsNum,XEN,YEN,ZEN,REN,ITsEff,CV,JTR)
  do ITS=1,ITsEff
    N1 = JTR(1,ITS)
    N2 = JTR(2,ITS)
    N3 = JTR(3,ITS)
    PTS(:,1) = CV(:,N1)
    PTS(:,2) = CV(:,N2)
    PTS(:,3) = CV(:,N3)
    NV = 3
    PP(:) = ZERO

    ! Per ciascuna tessera, trova la porzione scoperta e ne
    ! calcola l'area con il teorema di Gauss-Bonnet; il punto
    ! rappresentativo e' definito come media dei vertici della porzione
    ! scoperta di tessera e passato in PP.
    ! I vertici di ciascuna tessera sono conservati in
    ! pVERT(3,MxVert,MxTs), il numero di vertici di ciascuna tessera e'
    ! in pNVERT(MxTs), e i centri dei cerchi di ciascun lato sono in
    ! pCENTR(3,MxVert,MxTs). In pIntS(MxVert,MxTs) sono registrate le
    ! sfere a cui appartengono i lati delle tessere.

    call TESSERA(iPrint,NS,NSFE,NV,Xs,Ys,Zs,Rs,pIntS(:,MxTs),PTS,CCC,PP,AREA)
    if (AREA == Zero) cycle
    NN1 = NN1+1
    NN = min(NN1,MxTs)
    Xt(NN) = PP(1)
    Yt(NN) = PP(2)
    Zt(NN) = PP(3)
    At(NN) = AREA

    pISph(NN) = NSFE
    pNVERT(NN) = NV
    pVERT(:,1:NV,NN) = PTS(:,1:NV)
    pCENTR(:,1:NV,NN) = CCC(:,1:NV)
    pIntS(1:NV,NN) = pIntS(1:NV,MxTs)
  end do
end do
NTS = NN

call mma_deallocate(CCC)
call mma_deallocate(PTS)

! Verifica se due tessere sono troppo vicine
TEST = 0.02_wp
TEST2 = TEST*TEST
do I=1,NTS-1
  if (At(I) == ZERO) cycle
  XI = Xt(I)
  YI = Yt(I)
  ZI = Zt(I)
  II = I+1
  do J=II,NTS
    if (pISph(I) == pISph(J)) cycle
    if (At(J) == ZERO) cycle
    XJ = Xt(J)
    YJ = Yt(J)
    ZJ = Zt(J)
    RIJ = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
    if (RIJ > TEST2) cycle

    ! La routine originaria sostituiva le due tessere troppo vicine con una
    ! sola tessera. Nel caso Gauss-Bonnet, anche i vertici delle tessere
    ! e i centri degli archi vengono memorizzati ed e' impossibile sostituirli
    ! nello stesso modo: percio' l'area della tessera piu' piccola verra'
    ! trascurata per evitare problemi nella autopolarizzazione.
    if (IPRINT == 2) write(u6,1000) I,J,TEST2
    if (At(I) < At(J)) At(I) = ZERO
    if (At(I) >= At(J)) At(J) = ZERO
  end do
end do

! E' preferibile eliminare del tutto le tessere che per
! qualche motivo hanno AREA = 0, ridefinendo tutti gli
! indici: l'errore numerico cosi' introdotto e' in genere
! trascurabile e si evitano problemi di convergenza

! Define here the number of tesserae in electrostatic calculations
! to avoid problems in successive calcn.
ITS = 0
do while (ITS < NTS)
  ITS = ITS+1
  if (At(ITS) < 1.0e-10_wp) then
    do I=ITS,NTS-1
      At(I) = At(I+1)
      Xt(I) = Xt(I+1)
      Yt(I) = Yt(I+1)
      Zt(I) = Zt(I+1)
      pISph(I) = pISph(I+1)
      pNVERT(I) = pNVERT(I+1)
      pVERT(:,:,I) = pVERT(:,:,I+1)
      pCENTR(:,:,I) = pCENTR(:,:,I+1)
      pIntS(:,I) = pIntS(:,I+1)
    end do
    NTS = NTS-1
    ITS = ITS-1
  end if
end do
!***********************************************************************
! Calcola il volume della cavita' con la formula (t. di Gauss):
!            V=SOMMAsulleTESSERE{A r*n}/3
! dove r e' la distanza del punto rappresentativo dall'origine,
! n e' il versore normale alla tessera, A l'area della tessera,
! e * indica il prodotto scalare.
!***********************************************************************
VCav = ZERO
do ITS=1,NTS
  NSFE = pISph(ITS)
  ! Trova il versore normale
  XN = (Xt(ITS)-Xs(NSFE))/Rs(NSFE)
  YN = (Yt(ITS)-Ys(NSFE))/Rs(NSFE)
  ZN = (Zt(ITS)-Zs(NSFE))/Rs(NSFE)
  ! Trova il prodotto scalare
  PROD = Xt(ITS)*XN+Yt(ITS)*YN+Zt(ITS)*ZN
  VCav = VCav+At(ITS)*PROD/Three
end do
!***********************************************************************
! Stampa la geometria della cavita'
Scav = ZERO
pSSph(1:NS) = ZERO
do I=1,NTS
  K = pISph(I)
  pSSph(K) = pSSph(K)+At(I)
end do
OMEGA = OMEGA/DEG2RAD
if (IPRINT == 2) write(u6,1100) OMEGA,RSOLV,RET,FRO,NS
do I=1,NS
  if (IPRINT == 2) write(u6,1200) I,Xs(I),Ys(I),Zs(I),Rs(I),pSSph(I)
  Scav = Scav+pSSph(I)
end do
if (IPRINT == 2) write(u6,1300) NTS,Scav,VCav

! Trasforma i risultati in bohr
Rs(1:NS) = Rs(1:NS)/Angstrom
Xs(1:NS) = Xs(1:NS)/Angstrom
Ys(1:NS) = Ys(1:NS)/Angstrom
Zs(1:NS) = Zs(1:NS)/Angstrom
do I=1,NTS
  pVERT(:,1:pNVERT(I),I) = pVERT(:,1:pNVERT(I),I)/Angstrom
  pCENTR(:,1:pNVERT(I),I) = pCENTR(:,1:pNVERT(I),I)/Angstrom
end do
At(1:NTS) = At(1:NTS)/(Angstrom*Angstrom)
Xt(1:NTS) = Xt(1:NTS)/Angstrom
Yt(1:NTS) = Yt(1:NTS)/Angstrom
Zt(1:NTS) = Zt(1:NTS)/Angstrom
if (IPRINT == 3) then
  write(u6,1500)
  write(u6,1600)
  write(u6,1700) (I,pISph(I),At(I),Xt(I),Yt(I),Zt(I),I=1,NTS)
end if
if (NN1 > MxTs) then
  write(u6,1240) NN1,MxTs
  write(u6,1112)
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Re-allocate with actual dimensioning
if (.not. allocated(PCMSph)) then

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
PCMSph(1,:) = Xs(1:NS)
PCMSph(2,:) = Ys(1:NS)
PCMSph(3,:) = Zs(1:NS)
PCMSph(4,:) = Rs(1:NS)

! PCMTess
PCMTess(1,:) = Xt(1:nTs)
PCMTess(2,:) = Yt(1:nTs)
PCMTess(3,:) = Zt(1:nTs)
PCMTess(4,:) = At(1:nTs)

! Vert
Vert(:,:,:) = pVert(:,:,1:nTs)

! Centr
Centr(:,:,:) = pCentr(:,:,1:nTs)

! SSph
SSph(:) = pSSph(1:NS)

! nOrd
PCM_N(:) = pNs(1:NS)

! ISph
PCMiSph(:) = pISph(1:nTs)

! NVert
NVert(:) = pNVert(1:nTs)

! IntSph
IntSph(:,:) = pIntS(:,1:nTs)

! NewSph
NewSph(:,:) = pNewS(:,1:NS)
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

end subroutine FndTess
