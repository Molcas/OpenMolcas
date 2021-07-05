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

function AtNear(IAt,IAn,NBond,IBond,DH,DId,DAr,Chg)

implicit real*8(A-H,O-Z)
!logical AddH
parameter(MxBond=12)
dimension IAn(*), NBond(*), IBond(MxBond,*), Chg(*)

AtNear = 0.0D+00
NearAt = 0
NH = 0
NNit = 0
NF = 0
NX = 0
NId = 0
NCSP3 = 0
NCSP2 = 0
NCSP = 0
DX = 0.0D+00
IAnI = IAn(IAt)
NbI = NBond(IAt)
do j=1,NbI
  Jat = ibond(j,IAt)
  IAnJ = Ian(JAt)
  NbJ = NBond(JAt)
  ChgJ = Chg(JAt)
  ! Hydrogen
  if (IAnJ == 1) NH = NH+1
  ! Carbon
  if (IAnJ == 6) then
    if (nbj == 4) NCSP3 = NCSP3+1
    if (nbj == 3) NCSP2 = NCSP2+1
    if (nbj == 2) NCSP = NCSP+1
  end if
  ! Nitrogen
  if (IAnJ == 7) then
    if (NBJ < 4) then
      if (IAnI == IAnJ) NId = NId+1
      NX = NX+1
      NNit = NNit+1
    else
      NX = NX+1
    end if
  end if
  ! Oxygen
  if (IAnJ == 8) then
    NOX = 1
    if ((IAnI == IAnJ) .and. (abs(ChgJ) < 1.0D-01)) NId = NId+1
    if (CHGJ > 4.0D-01) then
      NOX = 0
      DX = 5.0D-01
    end if
    if (CHGJ < -4.0D-01) NOX = 3
    NX = NX+NOX
  end if
  ! Fluorine
  if (IAnJ == 9) then
    if (IAnI == IAnJ) NId = NId+1
    NF = NF+1
  end if
  ! Positive Phosphorus
  if ((IAnJ == 15) .and. (NBJ == 4)) DX = DX-1.5d0
  ! Positive Sulfur
  if ((IAnJ == 16) .and. (NBJ == 3)) DX = DX-2.7d0
end do
NearAt = NH
NAr = 0
if ((IAn(IAt) == 6) .and. (abs(Chg(IAt)) < 4.0D-01)) NearAt = NearAt+NH
if ((IAn(IAt) == 6) .and. (NBond(IAt) == 3) .and. (NCSP2 >= 1)) then
  if ((NH+NCSP3) /= 2) NAR = NAlpAr(IAt,IAn,NBond,IBond,Chg)
end if
if ((IAn(IAt) == 6) .and. (NBond(IAt) == 4)) NearAt = NearAt-NX-NCSP2-NNit+NF+NCAlph(IAt,NH,NCSP3,IAn,NBond,IBond,Chg)
if (NH > 3) NearAt = NearAt-1
AtNear = dble(NearAt)-DX
DH = dble(NH)
DAr = 5.0D-01*dble(NAr)
DId = 5.0D-01*dble(NId)

return

end function AtNear
!====
function Hybnew(OKUAH,OKCHG,IAt,IAn,NBond,IBond,IBType,PBO,Chg)

implicit real*8(A-H,O-Z)
logical PIAT, OKUAH, OKChg
parameter(MxBond=12)
dimension NBond(*), IBond(MxBond,*), IBtype(MxBond,*), IAn(*)
dimension PBO(MxBond,*)

NumatI = IAn(IAt)
HybNew = 3.0D+00
!Chg = 0.0D+00
IQQ = 0
if (.not. OKUAH) then
  HybNew = 0.0D+00
  return
end if

! Hydrogen

if (IAn(IAt) == 1) then
  HybNew = 0.0D+00
  if ((NBond(IAT) == 0) .and. (.not. OKCHG)) Chg = 1.0D+00
end if

! Carbon

if (IAN(IAt) == 6) then
  NBTot = 0
  PBtot = 0.0D+00
  do jj=1,NBond(IAt)
    NBtot = NBTot+IBType(jj,IAt)
    PBTot = PBTot+PBO(jj,IAt)
  end do
  if ((NBond(IAt) == 3) .and. (NBTot >= 4)) HybNew = 2.0D+00
  if ((NBond(IAt) == 3) .and. (PBtot > 3.7d+0)) Hybnew = 2.0d+00
  if (NBond(IAt) == 2) HybNew = 1.0D+00
end if

! Nitrogen, Phosphorus, Arsenic, Antimony

if (IColAt(NumAtI) == 5) then
  NBIA = NBond(IAt)
  goto(10,20,30,40) NBIA
10 continue
  HybNew = 1.0D+00
  ICC = IBond(1,IAt)
  ICN = IAn(IBond(1,IAt))
  if ((ICN == 6) .and. (NBond(ICC) == 1) .and. (.not. OKCHG)) Chg = -1.0D+00
  goto 50
20 continue
  if (PIAT(IAt,IAn,NBond,IBond)) then
    HybNew = 2.0D+00
  else
    if (.not. OKCHG) Chg = -1.0D+00
  end if
  goto 50
30 continue
  if (PIAT(IAt,IAn,NBond,IBond)) HybNew = 2.0D+00
  goto 50
40 continue
  NBII = 0
  do ii=1,NBIA
    NumII = IAn(IBond(II,IAt))
    if ((NumII == 6) .or. (NumII == 1)) NBII = NBII+1
  end do
  if ((NBII > 3) .and. (.not. OKCHG)) Chg = 1.0D+00
50 continue
end if

! Oxygen, Sulfur, Selenium, Tellurium

if (IColAt(NumAtI) == 6) then

  ! tricoordinated (positively charged only if bonded to C and/or H only)

  if (NBond(IAt) == 3) then
    NBII = 0
    do ii=1,3
      NumII = IAn(IBond(II,IAt))
      if ((NumII == 1) .or. (NumII == 6)) NBII = NBII+1
    end do
    if (NBII == 3) then
      if (.not. OKCHG) Chg = 1.0D+00
      HybNew = 3.0D+00
    end if
  end if

  ! bicoordinated

  if (NBond(IAt) == 2) then
    if (.not. OKCHG) Chg = 0.0D+00
    HybNew = 3.0D+00

    ! test for protonated aldehydes, ketones and similar for S

    do jj=1,2
      JAt = IBond(JJ,IAt)
      NumJ = IAn(JAt)
      NBJ = NBond(JAt)
      if ((NumJ == 6) .and. (NBJ == 3)) then
        IQQ = 0
        do kk=1,3
          KAt = IBond(kk,JAt)
          NumK = IAn(KAt)
          NBK = NBond(KAt)
          if ((NumK == 6) .and. (NBK == 4)) IQQ = IQQ+1
        end do
      end if
    end do
    if (IQQ >= 2) then
      if (.not. OKCHG) Chg = 1.0D+00
      HybNew = 2.0D+00
    end if
  end if
  if (NBond(IAt) == 1) then
    HybNew = 2.0D+00
    Jat = IBond(1,IAt)
    icc = IAn(JAt)
    if (icc == 1) then
      HybNew = 3.0D+00
      if (.not. OKCHG) Chg = -1.0D+00
    end if
    if (icc == 6) then
      if (NBond(JAt) == 4) then
        if (.not. OKCHG) CHg = -1.0D+00
        HybNew = 3.0D+00
      end if
      ISP2 = 0
      JSP2 = 0
      ICARB = 0
      if (NBond(JAt) == 3) JSP2 = 1
      do K=1,NBond(JAt)
        KAT = IBond(K,JAt)
        if ((IAN(KAT) == 6) .and. (NBond(KAT) == 3)) ISP2 = ISP2+1
        if ((IAN(KAT) == 8) .and. (NBond(KAt) == 1)) ICARB = ICARB+1
      end do
      if (ISP2 > 1) then
        if (.not. OKCHG) Chg = -1.0D+00
        HybNew = 3.0D+00
      end if
      if ((JSp2 == 1) .and. (ICARB == 2)) then
        if (.not. OKCHG) Chg = -5.0D-01
        HybNew = 3.0D+00
      end if
    end if
  end if
end if

! Halogens

if ((IColAt(NumAtI) == 7) .and. (NBond(IAt) == 0)) then
  if (.not. OKCHG) Chg = -1.0D+00
end if

!write(IOut,'(I3,2F3.1)') IAn(IAt),HybNew,Chg
return

end function Hybnew
!====
integer function NCAlph(IAt,NHI,NCSP3I,IAn,NBond,IBond,Chg)

implicit real*8(A-H,O-Z)
parameter(MxBond=12)
dimension NBond(*), IBond(MxBond,*), IAn(*), chg(*)

NHetI = 4-NHI-NCSP3I
NCP1 = 0
NCM1 = 0
!write(IOut,'("Atom",I2," NHet=",I1)') IAt,NHetI
do jj=1,4
  jat = IBond(jj,iat)
  ianj = ian(jat)
  nbj = nbond(jat)
  nhj = 0
  ncsp3j = 0
  nhetj = 0
  iplj = 0
  if ((ianj == 6) .and. (nbj == 4)) then
    do kk=1,4
      kat = IBond(kk,JAt)
      iank = ian(KAt)
      nbk = nbond(kat)
      if (iank == 1) NHJ = NHJ+1
      if ((iank == 6) .and. (nbk == 4)) NCSP3J = NCSP3J+1
      if (chg(kat) > 4.0d-01) iplj = 1
    end do
    NHETJ = 4-NHJ-NCSP3J
    if ((NHetI >= 0) .and. (NHetJ == 0)) NCP1 = NCP1+1
    if ((NHetI == 0) .and. (NHetJ > 0) .and. (iplj == 0)) NCM1 = NCM1+1
  end if
  !write(IOut,*) JAt,NHetJ,NCP1,NCM1
end do
NCAlph = NCP1-NCM1

return

end function NCAlph
!====
integer function NAlPar(IAt,IAn,NBond,IBond,Chg)

implicit real*8(A-H,O-Z)
parameter(MxBond=12)
dimension NBond(*), IBond(MxBond,*), IAn(*), Chg(*)

NAlPar = -1
NArI = 0
do jj=1,3
  NCSP2J = 0
  NHetJ = 0
  jat = IBond(jj,iat)
  ianj = ian(jat)
  nbj = nbond(jat)
  if ((ianj == 7) .and. (nbj >= 3)) NCSP2J = 2
  if ((ianj == 6) .and. (nbj == 3)) then
    do kk=1,3
      kat = IBond(kk,JAt)
      iank = ian(KAt)
      nbk = nbond(Kat)
      chgk = chg(KAt)
      if (chgk < 4.0D-01) then
        if ((iank == 6) .and. (nbk == 3)) NCSP2J = NCSP2J+1
        if ((iank == 8) .or. (iank == 9)) NHetJ = NHetJ+1
        if ((iank == 17) .or. (iank == 35) .or. (iank == 53)) NHetJ = NHetJ+1
        if (iank == 7) then
          if (nbk <= 2) NHetJ = NHetJ+1
          if (nbk >= 3) NCSP2J = NCSP2J+1
        end if
      end if
    end do
  end if
  if ((NCSP2J >= 2) .and. (NHetJ == 0)) NarI = NArI+1
end do
if (NArI >= 2) NAlPar = 1

return

end function NAlPar
!====
function IPBO(ToAng,IA,JA,RIJ,BondOr)
! Return the order of the bond between atoms of atomic numbers IA
! and JA, separated by RIJ, or 0 if there is no bond.

implicit real*8(A-H,O-Z)

! Generate connectivity based on bond distances alone.  The criteria
! is whether the distances is no more than 30% longer than the sum of
! covalent radii of involved atoms. For the moment all bond types are
! determined using Pauling bond orders.
! Note that RIJ is multiplied by ToAng, i.e. if it is already in Ang.
! ToAng must be set = 1

IPBO = 0
R1IJ = RIJ*ToAng
R0IJ = RCov97(IA,JA)
!if (R1IJ > R0IJ*1.3d0) return
BondOr = exp((R0IJ-R1IJ)/0.3d0)
if (BondOr < 2.0D-01) return
IBondO = int(Sngl(BondOr+0.5d0))
IBondO = max(IBondO,1)
IBondO = min(IBondO,3)
IPBO = IBondO

return

end function IPBO
!====
integer function IColAt(NumAt)

dimension ICol(0:108)
save ICol
data ICol/0, &
          1,2, &
          1,2,3,4,5,6,7,8, &
          1,2,3,4,5,6,7,8, &
          1,2, &
            9,9,9,9,9,9,9,9,9,9, &
              3,4,5,6,7,8, &
          1,2, &
            9,9,9,9,9,9,9,9,9,9, &
              3,4,5,6,7,8, &
          1,2, &
            52*9/

IcolAt = ICol(NumAt)

return

end function IColAt
!====
logical function PiAt(IAt,IAn,NBond,IBond)

parameter(MxBond=12)
dimension IAn(*), NBond(*), IBond(MxBond,*)

PiAt = .false.
ISPI = -1
do J=1,NBond(IAT)
  JAt = IBond(J,IAt)
  IAnJ = IAn(JAt)
  ITpJ = IColAt(IAnJ)
  NBnJ = NBond(JAT)
  NPiJ = 0
  do K=1,NBnJ
    KAt = IBond(K,JAt)
    if ((IAn(KAt) == 6) .and. (NBond(KAt) == 3)) NPiJ = NPiJ+1
  end do
  if ((IAnJ == 6) .and. (NBnJ == 3)) then
    if (NPiJ >= 2) then
      ISPI = ISPI+2
    else
      ISPI = ISPI+1
    end if
  end if
  if ((ITpJ == 5) .and. (NBnJ == 2)) ISPI = ISPI+1
  if ((ITpJ == 5) .and. (NPiJ >= 2)) ISPI = ISPI+1
end do
if (ISPI > 0) PiAt = .true.

return

end function PiAt
!====
integer function IRowAt(NumAt)

dimension IRow(0:108)
save IRow
data IROW/0,2*1,8*2,8*3,18*4,18*5,32*6,22*7/

IRowAt = IRow(NumAt)

return

end function IRowAt
!====
function AtSymb(I)
! ATSYMB(I) is the Atomic Symbol corresponding to the Atomic number I
! if I=0  Atsymb=Bq
! if I=-1 Atsymb=X

character*2 AtSymb
#include "periodic_table.fh"

if (I > 0) then
  AtSymb = PTab(i)
  return
end if
if (I == -1) AtSymb = ' X'
if (I == 0) AtSymb = 'Bq'

return

end function AtSymb
