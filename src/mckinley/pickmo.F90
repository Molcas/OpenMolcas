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

subroutine PickMO(COUT,nOut,nAcO,icmp,iBasi,iBasn,jBasj,jBasn,kBask,kBasn,lBasl,lBasn,iaoii)

use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use pso_stuff
use Symmetry_Info, only: nIrrep

implicit real*8(a-h,o-z)
#include "etwas.fh"
#include "real.fh"
real*8 COUT(nOut)
integer iCmp(4), iBas(4), nBs(4)
integer iAOii(4)

iBas(1) = iBasi
iBas(2) = jBasj
iBas(3) = kBask
iBas(4) = lBasl
nBs(1) = iBasn
nBs(2) = jBasn
nBs(3) = kBasn
nBs(4) = lBasn
ip2 = 1

do iCnt=3,4
  ipC = 0
  do iIrrep=0,nIrrep-1
    iOrb = nIsh(iIrrep)
    do iAsh=1,nAsh(iIrrep)
      jj = iCmp(iCnt)
      do i1=1,jj
        iSO = iAOtSO(iAOii(iCnt)+i1,iIrrep)+iBas(iCnt)-1
        if (iSO > 0) then
          ip1 = ipC+(iOrb+iAsh-1)*nBas(iIrrep)+iSO
          call dcopy_(nBs(iCnt),CMO(ip1,1),1,COUT(ip2),1)
        else
          call dcopy_(nBs(iCnt),[0.0d0],0,COUT(ip2),1)
        end if
        ip2 = ip2+nBs(iCnt)
      end do
    end do
    ipc = ipc+nBas(iIrrep)**2
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nAcO)

end subroutine PickMO
