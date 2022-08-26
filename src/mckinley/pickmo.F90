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

subroutine PickMO(COUT,nOut,icmp,iBasi,iBasn,jBasj,jBasn,kBask,kBasn,lBasl,lBasn,iaoii)

use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use pso_stuff, only: CMO
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nOut, iCmp(4), iBasi, iBasn, jBasj, jBasn, kBask, kBasn, lBasl, lBasn, iAOii(4)
real(kind=wp) :: COUT(nOut)
#include "etwas.fh"
integer(kind=iwp) :: i1, iAsh, iBas(4), iCnt, iIrrep, iOrb, ip1, ip2, ipC, iSO, jj, nBs(4)

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
          call dcopy_(nBs(iCnt),[Zero],0,COUT(ip2),1)
        end if
        ip2 = ip2+nBs(iCnt)
      end do
    end do
    ipc = ipc+nBas(iIrrep)**2
  end do
end do

return

end subroutine PickMO
