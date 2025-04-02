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
use MCLR_Data, only: CMO
use MCLR_Data, only: ipMat, nA, nDens2
use input_mclr, only: nSym, nAsh, nBas, nIsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half

implicit none
#include "SysDef.fh"
! Input
integer nDMO, nDAO
real*8, dimension(nDMO) :: DMO
! Output
real*8, dimension(nDAO) :: DAO
! Auxiliaries
real*8, dimension(:), allocatable :: D1, OCCU, NatCMO
integer nLCMO, iS, i, j, iAA, jAA, nbas_tot, ij, iA, jA

!                                                                      *
!***********************************************************************
!                                                                      *
nLCMO = 0
nbas_tot = 0
do iS=1,nSym
  nbas_tot = nbas_tot+nbas(is)
  nLCMO = nLCMO+nBas(is)**2
end do

call mma_allocate(D1,nLCMO)
call FZero(D1,nLCMO)
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
call mma_allocate(NatCMO,ndens2,Label='NatCMO')

call NatOrb_MCLR(D1,CMO,NatCMO,OCCU)
call dmat_MCLR(NatCMO,OCCU,DAO)
ij = 0
do iS=1,nSym
  do i=1,nbas(is)
    do j=1,i-1
      ij = ij+1
      DAO(ij) = Half*DAO(ij)
    end do
    ij = ij+1
  end do
end do
call mma_deallocate(D1)
call mma_deallocate(OCCU)
call mma_deallocate(NatCMO)

end subroutine GetDmatAO
