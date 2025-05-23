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

subroutine PickMO(COUT,nOut,nSD,iSD4)

use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use pso_stuff, only: CMO
use Symmetry_Info, only: nIrrep
use Etwas, only: nAsh, nIsh
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nOut, nSD, iSD4(0:nSD,4)
real(kind=wp), intent(_OUT_) :: COUT(nOut)
integer(kind=iwp) :: i1, iAOii(4), iAsh, iBas(4), iCmp(4), iCnt, iIrrep, iOrb, ip1, ip2, ipC, iSO, jj, nBs(4)

iCmp(:) = iSD4(2,:)
iBas(:) = iSD4(8,:)
nBs(:) = iSD4(19,:)
iAOii(:) = iSD4(7,:)
ip2 = 1

do iCnt=3,4
  ipC = 0
  do iIrrep=0,nIrrep-1
    iOrb = nIsh(iIrrep)
    do iAsh=1,nAsh(iIrrep)
      jj = iCmp(iCnt)
      do i1=1,jj
        iSO = iAOtSO(iAOii(iCnt)+i1,iIrrep)+iBas(iCnt)
        if (iSO > 0) then
          ip1 = ipC+(iOrb+iAsh-1)*nBas(iIrrep)+iSO
          COUT(ip2:ip2+nBs(iCnt)-1) = CMO(ip1:ip1+nBs(iCnt)-1,1)
        else
          COUT(ip2:ip2+nBs(iCnt)-1) = Zero
        end if
        ip2 = ip2+nBs(iCnt)
      end do
    end do
    ipc = ipc+nBas(iIrrep)**2
  end do
end do

return

end subroutine PickMO
