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

subroutine DAN(Dens)

use Basis_Info, only: nBas
use pso_stuff
use Symmetry_Info, only: nIrrep

implicit real*8(a-h,o-z)
#include "etwas.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8 Dens(nDens)
integer na(0:7), ipcm(0:7)
real*8, allocatable :: Temp1(:), Temp2(:), Temp3(:)
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

ipD = 0
nnA = 0
ndenssq = 0
ipCC = 1
do i=0,nIrrep-1
  nDenssq = ndenssq+nBas(i)**2
  nA(i) = nnA
  ipcm(i) = ipCC
  nnA = nnA+nAsh(i)
  ipCC = ipCC+nBas(i)**2
end do

call mma_allocate(Temp1,nDensSQ,Label='Temp1')
call mma_allocate(Temp2,nDensSQ,Label='Temp2')
call mma_allocate(Temp3,nDensSQ,Label='Temp3')

do iS=0,nIrrep-1
  Temp1(:) = Zero
  if (nBas(is) > 0) then
    do iB=1,nAsh(iS)
      iiB = nA(iS)+iB
      do jB=1,nAsh(iS)
        jjB = nA(iS)+jB
        ijB = iTri(iiB,jjB)
        ip1 = nBas(iS)*(nISh(iS)+iB-1)+nIsh(is)+jb
        Temp1(ip1) = G1(ijB,1)
      end do
    end do

    call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),1.0d0,CMO(ipCM(iS),1),nBas(is),Temp1,nBas(is),0.0d0,Temp3,nBas(is))
    call DGEMM_('N','T',nBas(is),nBas(is),nBas(is),1.0d0,Temp3,nBas(is),CMO(ipCM(is),1),nBas(is),0.0d0,Temp2,nBas(is))

    do iBas=1,nBas(iS)
      do jBas=1,iBas
        ip1 = (iBas-1)*nBas(iS)+jBas
        ip2 = iTri(iBas,jBas)
        Fact = 2.0d0
        if (iBas == jBas) Fact = 1.0d0
        Dens(ipD+ip2) = Temp2(ip1)*Fact
      end do
    end do
    ipD = ipD+nBas(is)*(nBas(is)+1)/2
  end if
end do

call mma_deallocate(Temp3)
call mma_deallocate(Temp2)
call mma_deallocate(Temp1)

return

end subroutine DAN
