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

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nB, nS, iBigT
character(len=6) :: TDMchar
#include "maxi.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iDisk, indx, iS1, iS2, iTocBig(MxStOT), kaunter, Lu, nSize !IFG
integer(kind=iwp), external :: IsFreeUnit

iDisk = 0
kaunter = 0
nSize = nB*(nB+1)/2
indx = 0
Lu = 72
Lu = IsFreeUnit(Lu)
call DaName(Lu,TDMchar)
call iDaFile(Lu,2,iTocBig,MxStOT,iDisk)
do iS1=1,nS
  do iS2=1,iS1
    kaunter = kaunter+1
    iDisk = iTocBig(kaunter)
    call dDaFile(Lu,2,Work(iBigT+indx),nSize,iDisk)
    indx = indx+nSize
  end do
end do
call DaClos(Lu)

return

end subroutine FetchTDM
