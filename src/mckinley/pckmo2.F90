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

subroutine PckMO2(COUT,icmpi,iBasi,jcmpj,jBasj,iAOi,jAOj)

use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use pso_stuff, only: CMO
use Symmetry_Info, only: nIrrep
use Etwas, only: nAsh, nIsh
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: COUT(*)
integer(kind=iwp), intent(in) :: icmpi, iBasi, jcmpj, jBasj, iAOi, jAOj
integer(kind=iwp) :: i1, iaoii(4), iAsh, iCmp(4), iCnt, iIrrep, iOrb, ip1, ip2, ipC, iSO, jj, nBs(4)

nBs(1) = iBasi
nBs(2) = jBasj
iAOii(1) = iAOi
iAOii(2) = jAOj
icmp(1) = icmpi
icmp(2) = jcmpj
ip2 = 1

do iCnt=1,2
  ipC = 0
  do iIrrep=0,nIrrep-1
    iOrb = nIsh(iIrrep)
    do iAsh=1,nAsh(iIrrep)
      jj = iCmp(iCnt)
      do i1=1,jj
        iSO = iAOtSO(iAOii(iCnt)+i1,iIrrep)
        if (iSO > 0) then
          ip1 = ipC+(iOrb+iAsh-1)*nBas(iIrrep)+iso
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

end subroutine PckMO2
