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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine GetDmatAO(DMO,DAO,nDMO,nDAO)
! Purpose: calculate the active 1RDM in AO basis given that in MO basis

use Index_Functions, only: iTri
use MCLR_Data, only: CMO, ipMat, nA, nDens
use input_mclr, only: nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDMO, nDAO
real(kind=wp), intent(in) :: DMO(nDMO)
real(kind=wp), intent(out) :: DAO(nDAO)
integer(kind=iwp) :: i, iA, iAA, ij, iS, j, jA, jAA, nbas_tot, nLCMO
real(kind=wp), allocatable :: D1(:), NatCMO(:), OCCU(:)

!                                                                      *
!***********************************************************************
!                                                                      *
nLCMO = sum(nbas(1:nSym)**2)
nbas_tot = sum(nbas(1:nSym))

call mma_allocate(D1,nLCMO)
D1(:) = Zero
! First, converting DMO into D1
! similar to computing D_K from G1q, as done in out_pt2.f
!********************************************************
do iS=1,nSym
  do iA=1,nash(is)
    do jA=1,nash(is)
      i = iA+nish(is)
      j = jA+nish(is)
      iAA = iA+na(is)
      jAA = jA+na(is)
      D1(ipmat(is,is)+i-1+(j-1)*nbas(is)) = DMO(iTri(iAA,jAA))
    end do
  end do
end do

!********************************************************
call mma_allocate(OCCU,nbas_tot,Label='OCCU')
call mma_allocate(NatCMO,nDens,Label='NatCMO')

call NatOrb_MCLR(D1,CMO,NatCMO,OCCU)
call dmat_MCLR(NatCMO,OCCU,DAO)
ij = 0
do iS=1,nSym
  do i=1,nbas(is)
    DAO(1:i-1) = Half*DAO(1:i-1)
    ij = ij+i+1
  end do
end do
call mma_deallocate(D1)
call mma_deallocate(OCCU)
call mma_deallocate(NatCMO)

end subroutine GetDmatAO
