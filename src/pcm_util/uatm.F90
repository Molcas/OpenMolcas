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

implicit real*8(A-H,O-Z)
logical AlBond, OKUAH, OKCHG
!character AtSymb*2,AppNum*3,BCH*3,HH1*10,HH2*10,HY1*4,HY2*4,HY3*4
character AtSymb*2, BCH*3, HH1*10, HH2*10, HY1*4, HY2*4, HY3*4
parameter(MxBond=12)
parameter(ISAX=1000)
dimension RE(*), C(3,NAt), IAn(*), NOrd(*), Chg(NAt)
dimension NBond(ISAX), IBond(MxBond,ISAX), IBtype(MxBond,ISAX)
dimension NTrBnd(ISAX), ITrBnd(MxBond,ISAX)
dimension ItrBtp(MxBond,ISAX)
dimension IHNum(ISAX), PBO(MxBond,ISAX)
dimension DHyb(ISAX), R0(0:7), gamma(0:7), NOKUAH(20), Coeff(0:5)
dimension DX(ISAX), DAl(ISAX)
!save R0,Gamma,DRQM,DAlDon,DAlAcc,NOKUAH,Coeff
save R0, Gamma, DRQM, NOKUAH, Coeff
data R0/0.0D+00,1.00D+00,1.50D+00,1.98D+00,2.08D+00,3*2.35D+00/
data Gamma/2*0.00D+00,9.00D-02,1.30D-01,4*1.50D-01/
!data DAlDon,DAlAcc,DRQM/2*1.0D-01,3.0D-01/
data DRQM/3.0D-01/
data NOKUAH/6,7,8,9,15,16,17,35,53,11*0/
data Coeff/1.0d+0,0.9d+0,0.6d+0,0.3d+0,0.1d+0,0.0d+0/

BCH = 'sdt'
HH1 = ' HHHHHHHHH'
HH2 = '  23456789'
HY1 = '*sss'
HY2 = ' ppp'
HY3 = '  23'
Pi = 4.0d0*atan(1.0d0)

! The number of explicitly parametrized atoms is NUAH and their atomic
! numbers are in NOKUAH

NUAH = 9

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
area = 0.0d0
do i1=1,nsfe
  iat = nord(i1)
  IEffbn = min(5,NTrBnd(IAt))
  surf = 4.0d0*Pi*re(I1)*re(I1)*coeff(IEffBn)
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

QTot = 0.0D+00
do IAt=1,NAt
  OKUAH = .false.
  do iverif=1,NUAH
    if (IAN(IAT) == NOKUAH(IVerif)) OKUAH = .true.
  end do
  DHYB(IAt) = HybNew(OKUAH,OKCHG,IAt,IAn,NBond,IBond,IBtype,PBO,Chg(IAt))
  QTot = QTot+Chg(IAt)
end do

! Verify and possibly correct charges

if (QTot >= 1.0D-05) then
  IQTot = int(QTot+1.0D-01)
!else if (QTot < 1.0D-05)
else
  IQtot = int(QTot-1.0D-01)
end if
if (IQTot /= ICharg) then
  FrQ = dble(ICharg)/dble(NAt-NTotH)
  do IAt=1,NAt
    !write(IOut,'(6X," Atom",I3," Old Q=",F5.2," New Q=",F5.2)') IAt,Chg(IAt),FrQ
    if (IAn(IAt) /= 1) Chg(IAt) = FrQ
  end do
end if
! Determine Re, Alpha
NSfe = 0
do IAt=1,NAt
  DDHYb = 0.0D+00
  DRQ = 0.0D+00
  DDX = 0.0D+00
  IRX = IRowAt(IAn(IAt))
  RX = R0(IRX)
  GX = gamma(IRX)
  NH = 0
  NHI = 0
  DH = 0.0D+00
  DAr = 0.0D+00
  DId = 0.0D+00
  DX(IAt) = 0.0D+00
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
    Re(I1) = Re(I1)+dble(IAddH)*Gx
    do iverif=1,NUAH
      if (IAN(IAT) == NOKUAH(IVerif)) OKUAH = .true.
    end do
    if (.not. okuah) goto 10
    if (Chg(IAt) < -1.0D-02) then
      DRQ = CHg(IAt)*DrQM
      if (IAn(IAt) == 7) DrQ = Chg(IAT)*2.0D-01
    end if
    if (Chg(IAt) > 4.0D-01) then
      if (IAn(IAt) == 8) DRQ = -Chg(IAt)*2.6D-01
      if (IAn(IAt) == 15) DRQ = -Chg(IAt)*4.5D-01
      if (IAn(IAt) == 16) DRQ = -Chg(IAt)*5.5D-01
    end if
    NHI = int(DHyb(IAt))
    NH = IHNum(IAt)
    if (NHI < 3) then
      DDHyb = (4.0D+00-DHyb(IAt))*GX
      if (IAn(IAt) /= 6) DDHyb = DDHyb/2.0D+00
    end if
    DDX = GX*(DX(IAT)+DAl(IAt))
    Re(I1) = RX+DDX+DDHyb+DRQ
    ! protect for too small C atoms
    if ((IAn(IAT) == 6) .and. (Re(I1) < 1.5D+00)) Re(I1) = 1.5D+00
10  continue
    NH = IHNum(IAt)
    JE = min(4,NTrBnd(IAt))
    if (iPrint > 5) write(IOut,'(6X,1X,I3,2X,A2,2A1,3X,3A1,2X,F5.2,3X,F4.2,2X,F5.3,3X,6(1X,A2,A3,"[",A1,"]"))') &
      IAt,AtSymb(IAn(IAt)),HH1(NH+1:NH+1),HH2(NH+1:NH+1),HY1(NHI+1:NHI+1),HY2(NHI+1:NHI+1),HY3(NHI+1:NHI+1),Chg(IAt),Alpha,Re(I1), &
      (AtSymb(IAn(ITrBnd(J,IAt))),'   ',BCH(ITrBtp(J,IAt):ITrBtp(J,IAt)),J=1,JE)
    if (NTrBnd(IAt) <= 4) goto 20
    JE = min(8,NTrBnd(IAt))
    if (iPrint > 5) write(IOut,'(6X,40X,4(1X,A2,A3,"[",A1,"]"))') (AtSymb(IAn(ITRBnd(JJ,IAt))),'   ', &
                                                                  BCH(ITrBtp(JJ,IAt):ITrBtp(JJ,IAt)),JJ=5,JE)

    if (NTrBnd(IAt) <= 8) goto 20
    JE = NTrBnd(IAt)
    if (iPrint > 5) write(IOut,'(6X,40X,4(1X,A2,A3,"[",A1,"]"))') (AtSymb(IAn(ITRBnd(JJ,IAt))),'   ', &
                                                                  BCH(ITrBtp(JJ,IAt):ITrBtp(JJ,IAt)),JJ=9,JE)
20  continue
  end if
end do
if (iPrint > 5) then
  write(IOut,'(6X,1X,78("-"))')
  write(IOut,*)
end if

return

end subroutine UATM
