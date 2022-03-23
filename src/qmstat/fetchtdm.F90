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

subroutine FetchTDM(nB,nS,iBigT,TDMchar)

use qmstat_global, only: nState
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nB, nS, iBigT
character(len=6) :: TDMchar
#include "WrkSpc.fh"
integer(kind=iwp) :: iDisk, indx, iS1, iS2, kaunter, Lu, nSize
integer(kind=iwp), allocatable :: iTocBig(:)
integer(kind=iwp), external :: IsFreeUnit

iDisk = 0
kaunter = 0
nSize = nTri_Elem(nB)
indx = 0
Lu = 72
Lu = IsFreeUnit(Lu)
call DaName(Lu,TDMchar)
call mma_allocate(iTocBig,nTri_Elem(nState),label='iTocBig')
call iDaFile(Lu,2,iTocBig,nTri_Elem(nState),iDisk)
do iS1=1,nS
  do iS2=1,iS1
    kaunter = kaunter+1
    iDisk = iTocBig(kaunter)
    call dDaFile(Lu,2,Work(iBigT+indx),nSize,iDisk)
    indx = indx+nSize
  end do
end do
call mma_deallocate(iTocBig)
call DaClos(Lu)

return

end subroutine FetchTDM
