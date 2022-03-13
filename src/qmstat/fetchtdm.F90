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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "WrkSpc.fh"
dimension iTocBig(MxStOT)
character TDMchar*6

iDisk = 0
kaunter = 0
nSize = nB*(nB+1)/2
index = 0
Lu = 72
Lu = IsFreeUnit(Lu)
call DaName(Lu,TDMchar)
call iDaFile(Lu,2,iTocBig,MxStOT,iDisk)
do iS1=1,nS
  do iS2=1,iS1
    kaunter = kaunter+1
    iDisk = iTocBig(kaunter)
    call dDaFile(Lu,2,Work(iBigT+index),nSize,iDisk)
    index = index+nSize
  end do
end do
call DaClos(Lu)

return

end subroutine FetchTDM
