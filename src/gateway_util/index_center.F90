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

function Index_Center(iCnt,iR,cIndex,iAtoms,nAtoms)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Index_Center
integer(kind=iwp), intent(in) :: iCnt, iR, nAtoms
integer(kind=iwp), intent(inout) :: cIndex(2,nAtoms), iAtoms
integer(kind=iwp) :: i

outer: do
  Index_Center = 0
  do i=1,iAtoms
    if ((cIndex(1,i) == iCnt) .and. (cIndex(2,i) == iR)) then
      Index_Center = i
      exit outer
    end if
  end do

  iAtoms = iAtoms+1
  cIndex(1,iAtoms) = iCnt
  cIndex(2,iAtoms) = iR
end do outer

return

end function Index_Center
