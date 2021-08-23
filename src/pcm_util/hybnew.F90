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

function Hybnew(OKUAH,OKCHG,MxBond,IAt,IAn,NBond,IBond,IBType,PBO,Chg)

use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Hybnew
logical(kind=iwp), intent(in) :: OKUAH, OKCHG
integer(kind=iwp), intent(in) :: MxBond, IAt, IAN(*), NBond(*), IBond(MxBond,*), IBtype(MxBond,*)
real(kind=wp), intent(in) :: PBO(MxBond,*)
real(kind=wp), intent(inout) :: Chg
integer(kind=iwp) :: ICARB, ICC, ICN, ii, IQQ, ISP2, JAt, jj, JSP2, K, KAt, kk, NBIA, NBII, NBJ, NBK, NBTot, NumatI, NumII, NumJ, &
                     NumK
real(kind=wp) :: PBtot
logical(kind=iwp) :: PIAT
integer(kind=iwp), external :: IColAt

NumatI = IAn(IAt)
HybNew = Three
!Chg = Zero
IQQ = 0
if (.not. OKUAH) then
  HybNew = Zero
  return
end if

! Hydrogen

if (IAn(IAt) == 1) then
  HybNew = Zero
  if ((NBond(IAT) == 0) .and. (.not. OKCHG)) Chg = One
end if

! Carbon

if (IAN(IAt) == 6) then
  NBTot = 0
  PBtot = Zero
  do jj=1,NBond(IAt)
    NBtot = NBTot+IBType(jj,IAt)
    PBTot = PBTot+PBO(jj,IAt)
  end do
  if ((NBond(IAt) == 3) .and. (NBTot >= 4)) HybNew = Two
  if ((NBond(IAt) == 3) .and. (PBtot > 3.7_wp)) Hybnew = Two
  if (NBond(IAt) == 2) HybNew = One
end if

! Nitrogen, Phosphorus, Arsenic, Antimony

if (IColAt(NumAtI) == 5) then
  NBIA = NBond(IAt)
  select case (NBIA)
    case (1)
      HybNew = One
      ICC = IBond(1,IAt)
      ICN = IAn(IBond(1,IAt))
      if ((ICN == 6) .and. (NBond(ICC) == 1) .and. (.not. OKCHG)) Chg = -One
    case (2)
      if (PIAT(MxBond,IAt,IAn,NBond,IBond)) then
        HybNew = Two
      else
        if (.not. OKCHG) Chg = -One
      end if
    case (3)
      if (PIAT(MxBond,IAt,IAn,NBond,IBond)) HybNew = Two
    case (4)
      NBII = 0
      do ii=1,NBIA
        NumII = IAn(IBond(II,IAt))
        if ((NumII == 6) .or. (NumII == 1)) NBII = NBII+1
      end do
      if ((NBII > 3) .and. (.not. OKCHG)) Chg = One
  end select
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
      if (.not. OKCHG) Chg = One
      HybNew = Three
    end if
  end if

  ! bicoordinated

  if (NBond(IAt) == 2) then
    if (.not. OKCHG) Chg = Zero
    HybNew = Three

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
      if (.not. OKCHG) Chg = One
      HybNew = Two
    end if
  end if
  if (NBond(IAt) == 1) then
    HybNew = Two
    Jat = IBond(1,IAt)
    icc = IAn(JAt)
    if (icc == 1) then
      HybNew = Three
      if (.not. OKCHG) Chg = -One
    end if
    if (icc == 6) then
      if (NBond(JAt) == 4) then
        if (.not. OKCHG) CHg = -One
        HybNew = Three
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
        if (.not. OKCHG) Chg = -One
        HybNew = Three
      end if
      if ((JSp2 == 1) .and. (ICARB == 2)) then
        if (.not. OKCHG) Chg = -Half
        HybNew = Three
      end if
    end if
  end if
end if

! Halogens

if ((IColAt(NumAtI) == 7) .and. (NBond(IAt) == 0)) then
  if (.not. OKCHG) Chg = -One
end if

!write(IOut,'(I3,2F3.1)') IAn(IAt),HybNew,Chg

return

end function Hybnew
