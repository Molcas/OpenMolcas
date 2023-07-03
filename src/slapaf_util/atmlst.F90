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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Cart(3,nAtom), Coor(3,mAtom), r(3)
logical New

!call RecPrt(' In AtmLst:Cart',' ',Cart,3,nAtom)

! Loop over list of symmetry unique centers

iSt = 1
do iAtom=1,nAtom
  iEnd = iSt
  call dcopy_(3,Cart(1,iAtom),1,Coor(1,iSt),1)

  ! Loop over the operators of the point group

  do ig=1,nIrrep-1
    r(1) = One
    if (iand(iOper(ig),1) /= 0) r(1) = -One
    r(2) = One
    if (iand(iOper(ig),2) /= 0) r(2) = -One
    r(3) = One
    if (iand(iOper(ig),4) /= 0) r(3) = -One
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
