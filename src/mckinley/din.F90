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

use Index_Functions, only: nTri_Elem
use Basis_Info, only: nBas
use pso_stuff, only: CMO, nDens
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Dens(nDens)
#include "etwas.fh"
integer(kind=iwp) :: iBas, iIrr, ip, ip1, ip2, ipD, nTemp2
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
  ip1 = 0
  ip2 = ipD
  do iBas=1,nBas(iIrr)
    Dens(ip2+1:ip2+iBas-1) = Four*Temp2(ip1+1:ip1+iBas-1)
    Dens(ip2+iBas) = Two*Temp2(ip1+iBas)
    ip1 = ip1+nBas(iIrr)
    ip2 = ip2+iBas
  end do
  ip = ip+nBas(iIrr)**2
  ipD = ipD+nTri_Elem(nBas(iIrr))

end do

call mma_deallocate(Temp2)

return

end subroutine DIN
