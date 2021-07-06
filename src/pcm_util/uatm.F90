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

subroutine UATM(IOut,ICharg,NAt,NSfe,ToAng,Re,Alpha,C,IAn,NOrd,Chg,iPrint)
! New settings of radii for electrostatic cavity
! for HF/6-31(+)G* and ICOMP=4
! explicit values for C,N,O,F,S,Cl,Br,I, otherwise modified UFF radii

use Constants, only: Zero, Two, Four, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut, ICharg, NAt, IAn(*), iPrint
integer(kind=iwp), intent(inout) :: NSfe, NOrd(*)
real(kind=wp), intent(in) :: ToAng, Alpha, C(3,NAt)
real(kind=wp), intent(inout) :: Re(*), Chg(NAt)
integer(kind=iwp), parameter :: ISAX = 1000, MxBond = 12
logical(kind=iwp) :: AlBond, OKCHG, OKUAH
character(len=2) :: AtSymb
integer(kind=iwp) :: i1, IAddH, IAt, IBond(MxBond,ISAX), IBtype(MxBond,ISAX), IEffbn, IHNum(ISAX), IQTot, IRX, &
                     ITrBnd(MxBond,ISAX), ItrBtp(MxBond,ISAX), IUSE, iverif, J, JE, JHAT, JJ, NBond(ISAX), NH, NHI, NTotH, &
                     NTrBnd(ISAX), NUAH
real(kind=wp) :: area, DAl(ISAX), DAr, DDHYb, DDX, DH, DHyb(ISAX), DId, DRQ, DX(ISAX), FrQ, GX, PBO(MxBond,ISAX), QTot, RX, surf
integer(kind=iwp), parameter :: NOKUAH(9) = [6,7,8,9,15,16,17,35,53]
real(kind=wp), parameter :: Coeff(0:5) = [1.0_wp,0.9_wp,0.6_wp,0.3_wp,0.1_wp,0.0_wp], DRQM = 0.3_wp, &
                            R0(0:7) = [0.00_wp,1.00_wp,1.50_wp,1.98_wp,2.08_wp,2.35_wp,2.35_wp,2.35_wp], &
                            Gamma(0:7) = [0.00_wp,0.00_wp,0.09_wp,0.13_wp,0.15_wp,0.15_wp,0.15_wp,0.15_wp]
character, parameter :: BCH(3) = ['s','d','t'], HH1(10) = [' ','H','H','H','H','H','H','H','H','H'], &
                        HH2(10) = [' ',' ','2','3','4','5','6','7','8','9'], HY1(4) = ['*','s','s','s'], &
                        HY2(4) = [' ','p','p','p'], HY3(4) = [' ',' ','2','3']
integer(kind=iwp), external :: IRowAt
real(kind=wp), external :: AtNear, HybNew, Pauling

! The number of explicitly parametrized atoms is NUAH and their atomic
! numbers are in NOKUAH

NUAH = size(NOKUAH)

! find bonds

OkChg = .false.
AlBond = .false.
call FndBnd(IOut,0,AlBond,ToAng,MxBond,NAt,IAn,C,NBond,IBond,IBtype,PBO,Re)
NTotH = 0
do IAt=1,NAt
  if (IAn(IAt) == 1) NTotH = NTotH+1
  IHNum(IAt) = 0
  NTrBnd(IAt) = 0
  do JJ=1,NBond(IAt)
    JHAT = IBond(JJ,IAT)
    if (IAn(JHAT) == 1) then
      IHNum(IAt) = IHNum(IAt)+1
    else
      NTrBnd(IAt) = NTrBnd(IAt)+1
      ITrBnd(NTrBnd(IAt),IAt) = JHAT
      ITrBtp(NTrBnd(IAt),IAt) = IBType(JJ,IAt)
    end if
  end do
end do
area = Zero
do i1=1,nsfe
  iat = nord(i1)
  IEffbn = min(5,NTrBnd(IAt))
  surf = Four*Pi*re(I1)*re(I1)*coeff(IEffBn)
  area = area+surf
end do

if (iPrint > 5) then
  write(iOut,*)
  write(iOut,*)
  write(iOut,'(6X,A)') 'Polarized Continuum Model Cavity'
  write(iOut,'(6X,A)') '================================'
  write(iOut,*)
  write(iOut,'(6X,A)') ' Nord Group  Hybr  Charge Alpha Radius            Bonded to'
end if

! Assign Charge and Hybridization to atoms

QTot = Zero
do IAt=1,NAt
  OKUAH = .false.
  do iverif=1,NUAH
    if (IAN(IAT) == NOKUAH(IVerif)) OKUAH = .true.
  end do
  DHYB(IAt) = HybNew(OKUAH,OKCHG,IAt,IAn,NBond,IBond,IBtype,PBO,Chg(IAt))
  QTot = QTot+Chg(IAt)
end do

! Verify and possibly correct charges

if (QTot >= 1.0e-5_wp) then
  IQTot = int(QTot+0.1_wp)
!else if (QTot < 1.0e-5_wp)
else
  IQtot = int(QTot-0.1_wp)
end if
if (IQTot /= ICharg) then
  FrQ = real(ICharg,kind=wp)/real(NAt-NTotH,kind=wp)
  do IAt=1,NAt
    !write(IOut,'(6X," Atom",I3," Old Q=",F5.2," New Q=",F5.2)') IAt,Chg(IAt),FrQ
    if (IAn(IAt) /= 1) Chg(IAt) = FrQ
  end do
end if
! Determine Re, Alpha
NSfe = 0
do IAt=1,NAt
  DDHYb = Zero
  DRQ = Zero
  DDX = Zero
  IRX = IRowAt(IAn(IAt))
  RX = R0(IRX)
  GX = gamma(IRX)
  NHI = 1
  DH = Zero
  DAr = Zero
  DId = Zero
  DX(IAt) = Zero
  OKUAH = .false.
  do iverif=1,NUAH
    if (IAN(IAT) == NOKUAH(IVerif)) then
      OKUAH = .true.
      DX(IAt) = AtNear(IAt,IAn,NBond,IBond,DH,DId,DAr,Chg)
    end if
  end do
  DAl(IAt) = DAr-DId

  ! retain only H atoms not bonded to OKUAH atoms

  IUSE = 1
  if (IAn(IAt) == 1) IUSE = 0

  if (IUSE == 1) then
    NSfe = NSfe+1
    I1 = NSfe
    NOrd(I1) = IAt
    OKUAH = .false.
    Re(I1) = Pauling(IAn(IAt))
    IAddH = min(3,IHNum(IAt))
    Re(I1) = Re(I1)+real(IAddH,kind=wp)*Gx
    do iverif=1,NUAH
      if (IAN(IAT) == NOKUAH(IVerif)) OKUAH = .true.
    end do
    if (.not. okuah) goto 10
    if (Chg(IAt) < -1.0e-2_wp) then
      DRQ = CHg(IAt)*DrQM
      if (IAn(IAt) == 7) DrQ = Chg(IAT)*0.2_wp
    end if
    if (Chg(IAt) > 0.4_wp) then
      if (IAn(IAt) == 8) DRQ = -Chg(IAt)*0.26_wp
      if (IAn(IAt) == 15) DRQ = -Chg(IAt)*0.45_wp
      if (IAn(IAt) == 16) DRQ = -Chg(IAt)*0.55_wp
    end if
    NHI = int(DHyb(IAt))+1
    if (NHI < 4) then
      DDHyb = (Four-DHyb(IAt))*GX
      if (IAn(IAt) /= 6) DDHyb = DDHyb/Two
    end if
    DDX = GX*(DX(IAT)+DAl(IAt))
    Re(I1) = RX+DDX+DDHyb+DRQ
    ! protect for too small C atoms
    if ((IAn(IAT) == 6) .and. (Re(I1) < 1.5_wp)) Re(I1) = 1.5_wp
10  continue
    NH = IHNum(IAt)+1
    JE = min(4,NTrBnd(IAt))
    if (iPrint > 5) write(IOut,'(6X,1X,I3,2X,A2,2A1,3X,3A1,2X,F5.2,3X,F4.2,2X,F5.3,3X,6(1X,A2,A3,"[",A1,"]"))') &
      IAt,AtSymb(IAn(IAt)),HH1(NH),HH2(NH),HY1(NHI),HY2(NHI),HY3(NHI),Chg(IAt),Alpha,Re(I1),(AtSymb(IAn(ITrBnd(J,IAt))),'   ', &
      BCH(ITrBtp(J,IAt)),J=1,JE)
    if (NTrBnd(IAt) <= 4) goto 20
    JE = min(8,NTrBnd(IAt))
    if (iPrint > 5) write(IOut,'(6X,40X,4(1X,A2,A3,"[",A1,"]"))') (AtSymb(IAn(ITRBnd(JJ,IAt))),'   ',BCH(ITrBtp(JJ,IAt)),JJ=5,JE)

    if (NTrBnd(IAt) <= 8) goto 20
    JE = NTrBnd(IAt)
    if (iPrint > 5) write(IOut,'(6X,40X,4(1X,A2,A3,"[",A1,"]"))') (AtSymb(IAn(ITRBnd(JJ,IAt))),'   ',BCH(ITrBtp(JJ,IAt)),JJ=9,JE)
20  continue
  end if
end do
if (iPrint > 5) then
  write(IOut,'(6X,1X,78("-"))')
  write(IOut,*)
end if

return

end subroutine UATM
