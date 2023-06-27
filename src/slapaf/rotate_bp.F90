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
! Copyright (C) 2023, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Rotate_BP
!
!> @brief
!>   Rotate \p BP to match \p Ref
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Rotate the vectors in a branching plane (\p BP) such that the plane
!> matches the reference (\p Ref). Do so in an "optimal" manner, meaning
!> that the vectors change as little as possible. This is achieved by
!> the SVD. The output value \p Phi is the "total" cosine for the angle
!> between the planes.
!>
!> @param[in,out] BP    Branching plane to transform
!> @param[in]     Ref   Reference branching plane
!> @param[in]     nDim  Dimension of the vectors
!> @param[in]     nVec  Number of vectors defining the planes (usually 2)
!> @param[out]    Phi   Product of principal cosines (total cosine)
!***********************************************************************

subroutine Rotate_BP(BP,Ref,nDim,nVec,Phi)

use linalg_mod, only: Gram_Schmidt, mult
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nVec
real(kind=wp), intent(inout) :: BP(nDim,nVec)
real(kind=wp), intent(in) :: Ref(nDim,nVec)
real(kind=wp), intent(out) :: Phi
integer(kind=iwp) :: m
real(kind=wp), allocatable :: Amat(:,:), Aux(:,:), OB1(:,:), OB2(:,:), Sing(:), Umat(:,:), Vmat(:,:)

! Get orthormal bases for the two planes

call mma_allocate(OB1,nDim,2,Label='OB1')
call mma_allocate(OB2,nDim,2,Label='OB2')
call Gram_Schmidt(BP,nVec,OB1,m)
if (m /= nVec) call WarningMessage(2,'Rotate_BP: Linear dependence in BP')
call Gram_Schmidt(Ref,nVec,OB2,m)
if (m /= nVec) call WarningMessage(2,'Rotate_BP: Linear dependence in Ref')

! Find principal angles and vectors

call mma_allocate(Amat,nVec,nVec,Label='Amat')
call mma_allocate(Umat,nVec,nVec,Label='Umat')
call mma_allocate(Vmat,nVec,nVec,Label='Vmat')
call mma_allocate(Sing,nVec,Label='Sing')
call mult(OB1,OB2,Amat,.true.,.false.)
call full_svd(nVec,nVec,Amat,Umat,Vmat,Sing)
Sing(:) = min(Sing,One)
Phi = product(Sing)

! Transform bases to principal vectors

call mma_allocate(Aux,nDim,nVec,Label='Aux')
call mult(OB1,Umat,Aux,.false.,.false.)
OB1(:,:) = Aux
call mult(OB2,Vmat,Aux,.false.,.true.)
OB2(:,:) = Aux
call mma_deallocate(Aux)

! Transform BP by replacing the OB1 basis with OB2
! (note that this is not a general transform, it is only valid for vectors within the BP)

call mult(OB1,BP,Amat,.true.,.false.)
call mult(OB2,Amat,BP,.false.,.false.)

call mma_deallocate(OB1)
call mma_deallocate(OB2)
call mma_deallocate(Amat)
call mma_deallocate(Umat)
call mma_deallocate(Vmat)
call mma_deallocate(Sing)

end subroutine Rotate_BP
