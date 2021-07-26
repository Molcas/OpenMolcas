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

function AtNear(MxBond,IAt,IAn,NBond,IBond,DH,DId,DAr,Chg)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: AtNear
integer(kind=iwp), intent(in) :: MxBond, IAt, IAn(*), NBond(*), IBond(MxBond,*)
real(kind=wp), intent(out) :: DH, DId, DAr
real(kind=wp), intent(in) :: Chg(*)
integer(kind=iwp) :: IAnI, IAnJ, j, Jat, NAr, NbI, NbJ, NCSP, NCSP2, NCSP3, NearAt, NF, NH, NId, NNit, NOX, NX
real(kind=wp) :: ChgJ, DX
integer(kind=iwp), external :: NAlPAr, NCAlph

AtNear = Zero
NearAt = 0
NH = 0
NNit = 0
NF = 0
NX = 0
NId = 0
NCSP3 = 0
NCSP2 = 0
NCSP = 0
DX = Zero
IAnI = IAn(IAt)
NbI = NBond(IAt)
do j=1,NbI
  JAt = IBond(j,IAt)
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
    if ((IAnI == IAnJ) .and. (abs(ChgJ) < 0.1_wp)) NId = NId+1
    if (CHGJ > 0.4_wp) then
      NOX = 0
      DX = Half
    end if
    if (CHGJ < -0.4_wp) NOX = 3
    NX = NX+NOX
  end if
  ! Fluorine
  if (IAnJ == 9) then
    if (IAnI == IAnJ) NId = NId+1
    NF = NF+1
  end if
  ! Positive Phosphorus
  if ((IAnJ == 15) .and. (NBJ == 4)) DX = DX-1.5_wp
  ! Positive Sulfur
  if ((IAnJ == 16) .and. (NBJ == 3)) DX = DX-2.7_wp
end do
NearAt = NH
NAr = 0
if ((IAn(IAt) == 6) .and. (abs(Chg(IAt)) < 0.4_wp)) NearAt = NearAt+NH
if ((IAn(IAt) == 6) .and. (NBond(IAt) == 3) .and. (NCSP2 >= 1)) then
  if ((NH+NCSP3) /= 2) NAR = NAlpAr(MxBond,IAt,IAn,NBond,IBond,Chg)
end if
if ((IAn(IAt) == 6) .and. (NBond(IAt) == 4)) NearAt = NearAt-NX-NCSP2-NNit+NF+NCAlph(MxBond,IAt,NH,NCSP3,IAn,NBond,IBond,Chg)
if (NH > 3) NearAt = NearAt-1
AtNear = real(NearAt,kind=wp)-DX
DH = real(NH,kind=wp)
DAr = Half*real(NAr,kind=wp)
DId = Half*real(NId,kind=wp)

return

end function AtNear
