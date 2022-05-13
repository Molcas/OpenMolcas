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
subroutine TransActMO2(MOs,MOas,mGrid)
! Purpose:
! obtaining an active MO array with a structure of MOs in TransActMO
! from an MO array with a structure of that in TransferMO

use nq_Info, only: iOff_Ash, mIrrep, nAsh, NASHT, nIsh, nOrbt, OffOrb
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mGrid
real(kind=wp), intent(out) :: MOs(NASHT,mGrid)
real(kind=wp), intent(in) :: MOas(mGrid,nOrbt)
integer(kind=iwp) :: iGrid, iIrrep, IOff1, IOff2

do iGrid=1,mGrid
  do iIrrep=0,mIrrep-1
    IOff2 = iOff_Ash(iIrrep)
    IOff1 = OffOrb(iIrrep)+nIsh(iIrrep)
    MOs(IOff2+1:IOff2+nAsh(iIrrep),iGrid) = MOas(iGrid,IOff1+1:IOff1+nAsh(iIrrep))
  end do
end do

return

end subroutine TransActMO2
