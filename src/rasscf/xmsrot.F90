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
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on May. 21, 2020, created this file.               *
! ****************************************************************

subroutine XMSRot(CMO,FI,FA)

use Index_Functions, only: nTri_Elem
use rasscf_global, only: LROOTS, NAC
use general_data, only: NTOT1, NTOT2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
real(kind=wp) :: CMO(NTOT2), FI(NTOT1), FA(NTOT1)
real(kind=wp), allocatable :: EigVec(:,:), FckO(:,:), FckS(:,:), GDMat(:,:,:)

!FckO:  Fock matrix for MO
!FckS:  Fock matrix for states
!GDMat: density matrix or transition density matrix

! Allocating Memory
call mma_allocate(GDMat,nTri_Elem(lRoots),NAC,NAC)
call mma_allocate(FckO,NAC,NAC)
call mma_allocate(FckS,lRoots,lRoots)
call mma_allocate(EigVec,lRoots,lRoots)

call CalcFckO(CMO,FI,FA,FckO)

call GetGDMat(GDMAt)

call CalcFckS(FckO,GDMat,FckS)

call CalcEigVec(FckS,lRoots,EigVec)

call printmat('ROT_VEC','XMS-PDFT',eigvec,lroots,lroots,7,8,'N')
! Deallocating Memory
call mma_deallocate(GDMat)
call mma_deallocate(FckO)
call mma_deallocate(FckS)
call mma_deallocate(EigVec)

end subroutine XMSRot
