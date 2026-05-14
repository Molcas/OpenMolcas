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

subroutine AtmLst(Cart,nAtom,Coor,mAtom)

use Symmetry_Info, only: nIrrep, iOper
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom, mAtom
real(kind=wp), intent(in) :: Cart(3,nAtom)
real(kind=wp), intent(out) :: Coor(3,mAtom)
integer(kind=iwp) :: iAtom, iEnd, ig, iGo, iSt
real(kind=wp) :: r(3), x, y, z
logical(kind=iwp) :: New

!call RecPrt(' In AtmLst:Cart',' ',Cart,3,nAtom)

! Loop over list of symmetry unique centers

iSt = 1
do iAtom=1,nAtom
  iEnd = iSt
  Coor(:,iSt) = Cart(:,iAtom)

  ! Loop over the operators of the point group

  do ig=1,nIrrep-1
    r(1) = merge(-One,One,btest(iOper(ig),0))
    r(2) = merge(-One,One,btest(iOper(ig),1))
    r(3) = merge(-One,One,btest(iOper(ig),2))
    x = r(1)*Cart(1,iAtom)
    y = r(2)*Cart(2,iAtom)
    z = r(3)*Cart(3,iAtom)

    New = .true.
    do iGo=iSt,iEnd
      if (New .and. (x == Coor(1,iGo)) .and. (y == Coor(2,iGo)) .and. (z == Coor(3,iGo))) New = .false.
    end do
    if (New) then
      iEnd = iEnd+1
      Coor(1,iEnd) = x
      Coor(2,iEnd) = y
      Coor(3,iEnd) = z
    end if
  end do      ! End loop over operators
  iSt = iEnd+1
end do        ! End loop over centers

!call RecPrt(' In AtmLst: Coor',' ',Coor,3,mAtom)

return

end subroutine AtmLst
