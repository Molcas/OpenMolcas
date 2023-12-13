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

subroutine PickPoints(nPick,Pick,DPick,nEPP,EPCo,Coo,dLimmo,BS)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nEPP
integer(kind=iwp), intent(out) :: nPick, Pick(nEPP)
real(kind=wp), intent(out) :: DPick(nEPP)
real(kind=wp), intent(in) :: EPCo(3,nEPP), Coo(3), dLimmo(2), BS
integer(kind=iwp) :: iP
real(kind=wp) :: Distad, xtwo, ytwo, ztwo

nPick = 0
do iP=1,nEPP
  xtwo = (EPCo(1,iP)-Coo(1))**2
  ytwo = (EPCo(2,iP)-Coo(2))**2
  ztwo = (EPCo(3,iP)-Coo(3))**2
  Distad = sqrt(xtwo+ytwo+ztwo)
  if ((Distad < dLimmo(2)*BS) .and. (Distad > dLimmo(1)*BS)) then
    nPick = nPick+1
    Pick(nPick) = iP
    DPick(nPick) = Distad
  end if
end do

return

end subroutine PickPoints
