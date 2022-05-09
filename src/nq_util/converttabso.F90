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
integer(kind=iwp) :: mAO, mGrid, nMOs
real(kind=wp) :: TabSO2(nMOs,mAO*mGrid), TabSO(mAO,mGrid,nMOs)
integer(kind=iwp) :: iAO, iEnd, iGrid, iOff, iSt, jAO, nAOGrid

nAOGrid = mAO*mGrid
! TabSO : mAO*mGrid x nMOs
! TabSO2: nMOs x mAO*nGrid

! loop over first and optionally second derivatives of the SOs
! this defines the length of nAO to 3 or 9.
iSt = 1
if (lft .and. lGGA) then
  iEnd = 9
else
  iEnd = 3
end if

do iGrid=1,mGrid

  do jAO=iSt,iEnd

    iOff = (iGrid-1)*mAO+jAO

    iAO = jAO+1
    call DCopy_(nMOs,TabSO(iAO,iGrid,1),nAOGrid,TabSO2(:,iOff),1)
  end do
end do

return

end subroutine ConvertTabSO
