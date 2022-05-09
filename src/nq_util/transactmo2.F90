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
integer(kind=iwp) :: mGrid
real(kind=wp) :: MOs(mGrid*NASHT), MOas(mGrid*nOrbt)
! Input: mGrid MOas
! Output: MOs
integer(kind=iwp) :: iGrid, iIrrep, IOff1, iOff2, iOff3

do iGrid=1,mGrid
  IOff3 = (iGrid-1)*nAsht
  do iIrrep=0,mIrrep-1
    IOff2 = IOff3+iOff_Ash(iIrrep)+1
    IOff1 = (OffOrb(iIrrep)+nIsh(iIrrep))*mGrid+iGrid
    call DCopy_(nAsh(iIrrep),MOas(iOff1),mGrid,MOs(IOff2),1)
  end do
end do

return

end subroutine TransActMO2
