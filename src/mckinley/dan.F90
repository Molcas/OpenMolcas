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

use Index_Functions, only: iTri, nTri_Elem
use Basis_Info, only: nBas
use pso_stuff, only: CMO, G1, nDens
use Symmetry_Info, only: nIrrep
use Etwas, only: nAsh, nIsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Dens(nDens)
integer(kind=iwp) :: i, iB, iBas, iiB, ijB, ip1, ip2, ipCC, ipcm(0:7), ipD, iS, jB, jjB, na(0:7), ndenssq, nnA
real(kind=wp), allocatable :: Temp1(:), Temp2(:), Temp3(:)

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

    call DGEMM_('N','N',nBas(is),nBas(is),nBas(is),One,CMO(ipCM(iS),1),nBas(is),Temp1,nBas(is),Zero,Temp3,nBas(is))
    call DGEMM_('N','T',nBas(is),nBas(is),nBas(is),One,Temp3,nBas(is),CMO(ipCM(is),1),nBas(is),Zero,Temp2,nBas(is))

    ip1 = 0
    ip2 = ipD
    do iBas=1,nBas(iS)
      Dens(ip2+1:ip2+iBas-1) = Two*Temp2(ip1+1:ip1+iBas-1)
      Dens(ip2+iBas) = Temp2(ip1+iBas)
      ip2 = ip2+iBas
      ip1 = ip1+nBas(iS)
    end do
    ipD = ipD+nTri_Elem(nBas(is))
  end if
end do

call mma_deallocate(Temp3)
call mma_deallocate(Temp2)
call mma_deallocate(Temp1)

return

end subroutine DAN
