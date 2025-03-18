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
      Subroutine GetDmatAO(DMO,DAO,nDMO,nDAO)
      use Arrays, only: CMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Half
      use MCLR_Data, only: ipMat, nA, nDens2
      use input_mclr, only: nSym,nAsh,nBas,nIsh
      Implicit None
#include "SysDef.fh"
!*****Purpose: calculate the active 1RDM in AO basis given that in MO
!*****         basis
!*****Input
      INTEGER nDMO,nDAO
      Real*8,DIMENSION(nDMO)::DMO
!*****Output
      Real*8,DIMENSION(nDAO)::DAO
!*****Auxiliaries
      Real*8,DIMENSION(:),Allocatable::D1,OCCU,NatCMO
      INTEGER nLCMO,iS,i,j,iAA,jAA,nbas_tot,ij,iA,jA, iTri
!                                                                      *
!***********************************************************************
!                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!                                                                      *
!***********************************************************************
!                                                                      *
      nLCMO=0
      nbas_tot=0
      DO iS=1,nSym
       nbas_tot=nbas_tot+nbas(is)
       nLCMO=nLCMO+nBas(is)**2
      END DO

      Call mma_allocate(D1,nLCMO)
      Call FZero(D1,nLCMO)
!*****First, converting DMO into D1
!*****similar to computing D_K from G1q, as done in out_pt2.f
!********************************************************
      DO iS=1,nSym
       Do iA=1,nash(is)
        do jA=1,nash(is)
         i=iA+nish(is)
         j=jA+nish(is)
         iAA=iA+na(is)
         jAA=jA+na(is)
         D1(ipmat(is,is)+i-1+(j-1)*nbas(is))=DMO(itri(iAA,jAA))
        end do
       End Do
      END DO

!********************************************************
      Call mma_allocate(OCCU,nbas_tot,Label='OCCU')
      Call mma_allocate(NatCMO,ndens2,Label='NatCMO')

      Call NatOrb_MCLR(D1,CMO,NatCMO,OCCU)
      Call dmat_MCLR(NatCMO,OCCU,DAO)
      ij=0
      DO iS=1,nSym
       Do i=1,nbas(is)
        do j=1,i-1
         ij=ij+1
         DAO(ij)=Half*DAO(ij)
        end do
        ij=ij+1
       End Do
      END DO
      Call mma_deallocate(D1)
      Call mma_deallocate(OCCU)
      Call mma_deallocate(NatCMO)

      End Subroutine GetDmatAO
