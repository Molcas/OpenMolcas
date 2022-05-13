!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
subroutine ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)

use nq_pdft, only: lft, lGGA
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mAO, mGrid, nMOs
real(kind=wp), intent(out) :: TabSO2(nMOs,mAO,mGrid)
real(kind=wp), intent(in) :: TabSO(mAO,mGrid,nMOs)
integer(kind=iwp) :: iEnd, iGrid, jAO

! TabSO : mAO*mGrid x nMOs
! TabSO2: nMOs x mAO*nGrid

! loop over first and optionally second derivatives of the SOs
! this defines the length of nAO to 3 or 9.
if (lft .and. lGGA) then
  iEnd = 9
else
  iEnd = 3
end if

do iGrid=1,mGrid
  do jAO=1,iEnd
    TabSO2(:,jAO,iGrid) = TabSO(jAO+1,iGrid,:)
  end do
end do

return

end subroutine ConvertTabSO
