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
integer(kind=iwp) :: mAO, mGrid, nMOs
real(kind=wp) :: MOs(mGrid*NASHT), TabMO(mAO,mGrid,nMOs)
! Input: mAO mGrid nMOs TabMO
! Output: MOs
integer(kind=iwp) :: iGrid, iIrrep, IOff1, iOff2, iOff3, nGridPi

nGridPi = mAO*mGrid
do iGrid=1,mGrid
  IOff1 = (iGrid-1)*NASHT
  do iIrrep=0,mIrrep-1
    IOff2 = IOff_Ash(iIrrep)+1
    IOff3 = IOff_BasAct(iIrrep)+1
    call DCopy_(nAsh(iIrrep),TabMO(1,iGrid,IOff3),nGridPi,MOs(IOff1+IOff2),1)
  end do
end do

return

end subroutine TransActMO
