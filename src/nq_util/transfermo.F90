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
subroutine TransferMO(MOas,TabMO,mAO,mGrid,nMOs,iAO)
! Purpose:
! Transferring MO information to MOas to be used in dgemm.
! It records the MO values on each grid point, too.
! But the difference from TransActMO is that the first and
! the second elements are the values of the first MO at grid
! point 1 and grid point 2.

use nq_Info, only: mIrrep, mOrb, nOrbt, OffBasFro
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mAO, mGrid, nMOs, iAO
real(kind=wp), intent(out) :: MOas(mGrid,nOrbt)
real(kind=wp), intent(in) :: TabMO(mAO,mGrid,nMOs)
integer(kind=iwp) :: iIrrep, IOff1, iOff2

IOff2 = 0
do iIrrep=0,mIrrep-1
  IOff1 = OffBasFro(iIrrep)
  MOas(:,IOff2+1:IOff2+mOrb(iIrrep)) = TabMO(iAO,:,IOff1+1:IOff1+mOrb(iIrrep))
  IOff2 = IOff2+mOrb(iIrrep)
end do

return

end subroutine TransferMO
