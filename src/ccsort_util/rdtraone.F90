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

subroutine RdTraOne()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     read the header of the transformed one-electron integral file    *
!                                                                      *
!***********************************************************************

use ccsort_global, only: nBasX, nDelX, nFroX, nSymX
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iDisk, LuTraOne, nOrbX(8), TocTraOne(64)
real(kind=wp) :: EcorX
character(len=LenIn8), allocatable :: BsLbl(:)

LuTraOne = 3

call DaName(LuTraOne,'TRAONE')

!ulf
iDisk = 0

call mma_allocate(BsLbl,MxOrb,label='BsLbl')
call WR_MOTRA_Info(LuTraOne,2,iDisk,TocTraOne,64,EcorX,nSymX,nBasX,nOrbX,nFroX,nDelX,8,BsLbl,LenIn8*MxOrb)
call mma_deallocate(BsLbl)

call Daclos(LuTraOne)

return

end subroutine RdTraOne
