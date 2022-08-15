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
subroutine TransActMO(MOs,TabMO,mAO,mGrid,nMOs)
! Purpose:
! Trasnferring active orbitals to the MOs array.
! It records the MO values on each grid point.
! The first and the second elements are the MO values
! of the first and the second active MO at grid point 1.

use nq_Info, only: IOff_Ash, IOff_BasAct, mIrrep, nAsh, NASHT
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mAO, mGrid, nMOs
real(kind=wp), intent(out) :: MOs(NASHT,mGrid)
real(kind=wp), intent(in) :: TabMO(mAO,mGrid,nMOs)
integer(kind=iwp) :: iGrid, iIrrep, IOff1, IOff2

do iGrid=1,mGrid
  do iIrrep=0,mIrrep-1
    IOff1 = IOff_Ash(iIrrep)
    IOff2 = IOff_BasAct(iIrrep)
    MOs(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = TabMO(1,iGrid,IOff2+1:IOff2+nAsh(iIrrep))
  end do
end do

return

end subroutine TransActMO
