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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine DIN(Dens)

use Basis_Info, only: nBas
use pso_stuff, only: CMO, nDens
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Dens(nDens)
#include "etwas.fh"
integer(kind=iwp) :: iBas, iIrr, ip, ip1, ip2, ipD, jBas, nTemp2
real(kind=wp), allocatable :: Temp2(:)

!                                                                      *
!***********************************************************************
!                                                                      *
nTemp2 = 0
do iIrr=0,nIrrep-1
  nTemp2 = max(nTemp2,nBas(iIrr))
end do

call mma_allocate(Temp2,nTemp2**2,Label='Temp2')

ip = 1
ipD = 0
do iIrr=0,nIrrep-1

  if (nBas(iIrr) == 0) cycle

  call DGEMM_('N','T',nBas(iIrr),nBas(iIrr),nIsh(iIrr),One,CMO(ip,1),nBas(iIrr),CMO(ip,1),nBas(iIrr),Zero,Temp2,nBas(iIrr))
  do iBas=1,nBas(iIrr)
    do jBas=1,iBas-1
      ip1 = (iBas-1)*nBas(iIrr)+jBas
      ip2 = iBas*(iBas-1)/2+jBas
      Dens(ipD+ip2) = Temp2(ip1)*Four
    end do
    ip1 = (iBas-1)*nBas(iIrr)+iBas
    ip2 = iBas*(iBas+1)/2
    Dens(ipD+ip2) = Temp2(ip1)*Two
  end do
  ip = ip+nBas(iIrr)**2
  ipd = ipD+nBas(iIrr)*(nBas(iIrr)+1)/2

end do

call mma_deallocate(Temp2)

return

end subroutine DIN
