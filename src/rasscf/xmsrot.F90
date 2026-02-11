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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
      Subroutine XMSRot(CMO,FI,FA)
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on May. 21, 2020, created this file.               *
! ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
      use rasscf_global, only: LROOTS, NAC
      use general_data, only: NTOT1,NTOT2

      Implicit None

!*****Input
      Real*8,DIMENSION(NTOT1):: FI,FA
      Real*8,Dimension(NTOT2)::CMO
!*****Auxiliary quantities
      Real*8,DIMENSION(:,:),Allocatable::FckO
!*****FckO:  Fock matrix for MO
      Real*8,DIMENSION(:,:),Allocatable::FckS,EigVec
!*****FckS:  Fock matrix for states
      Real*8,DIMENSION(:,:,:),Allocatable::GDMat
!*****GDMat: density matrix or transition density matrix

!     Allocating Memory
      CALL mma_allocate(GDMat,lRoots*(lRoots+1)/2,NAC,NAC)
      CALL mma_allocate(FckO,NAC,NAC)
      CALL mma_allocate(FckS,lRoots,lRoots)
      CALL mma_allocate(EigVec,lRoots,lRoots)
!
      CALL CalcFckO(CMO,FI,FA,FckO)

      CALL GetGDMat(GDMAt)
!
      CALL CalcFckS(FckO,GDMat,FckS)
!
      CALL CalcEigVec(FckS,lRoots,EigVec)
!
      call printmat('ROT_VEC','XMS-PDFT',eigvec,lroots,lroots,7,8,'N')
!     Deallocating Memory
      CALL mma_deallocate(GDMat)
      CALL mma_deallocate(FckO)
      CALL mma_deallocate(FckS)
      CALL mma_deallocate(EigVec)

      End Subroutine XMSRot
