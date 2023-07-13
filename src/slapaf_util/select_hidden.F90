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

subroutine Select_Hidden(mTtAtm,nHidden,Coord,HiddenCoord,iHiddenAN,nKept,iPL)
! Select among the hidden atoms the ones to be kept

use Slapaf_Info, only: rHidden
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mTtAtm, nHidden, iPL
real(kind=wp), intent(in) :: Coord(3,mTtAtm), HiddenCoord(3,nHidden)
integer(kind=iwp), intent(inout) :: iHiddenAN(nHidden), nKept
integer(kind=iwp) :: iAN, iAtom, iHid
real(kind=wp) :: Dist, dMax, X, Y, Z

! Criteria: dMin < distance < rHidden

dMax = rHidden ! bohr
do iHid=1,nHidden
  X = HiddenCoord(1,iHid)
  Y = HiddenCoord(2,iHid)
  Z = HiddenCoord(3,iHid)
  iAN = iHiddenAN(iHid)
  iAtom = 0
  do
    iAtom = iAtom+1
    Dist = sqrt((X-Coord(1,iAtom))**2+(Y-Coord(2,iAtom))**2+(Z-Coord(3,iAtom))**2)
    if (Dist <= dMax) then
      iHiddenAN(iHid) = -iAN
      nKept = nKept+1
    end if
    if ((iAtom >= mTtAtm) .or. (iHiddenAN(iHid) > 0)) exit
  end do
end do

! The end

if ((iPL > 3) .and. (nKept > 0)) write(u6,'(A,i3,A)') ' Select_Hidden: ',nKept,' hidden atoms are kept'

return

end subroutine Select_Hidden
