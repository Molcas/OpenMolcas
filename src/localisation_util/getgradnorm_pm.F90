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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetGradNorm_PM(nAtoms,nOrb2Loc,PA,GradNorm)
! Thomas Bondo Pedersen, December 2005.
!
! Purpose: compute the gradient norm of the Pipek-Mezey functional.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Four, Eight
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: GradNorm
real(kind=wp), allocatable :: Rmat(:,:)
integer(kind=iwp) :: iAtom, i,j
real(kind=wp) :: Fun, Rjj

call mma_Allocate(RMat,nOrb2Loc,nOrb2Loc,Label='RMat')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GradientNorm - needed for all optimization methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the R matrix used to calculate the gradient norm
RMat(:,:) = Zero
do iAtom=1,nAtoms
  do j=1,nOrb2Loc
    Rjj = PA(j,j,iAtom)
    do i=1,nOrb2Loc
      Rmat(i,j) = Rmat(i,j)+PA(i,j,iAtom)*Rjj
    end do
  end do
end do

!The Gradient Norm (always positive!) is used as threshold criterium
GradNorm = Zero
do i=1,nOrb2Loc-1
  do j=i+1,nOrb2Loc
    GradNorm = GradNorm + (Rmat(i,j) - Rmat(j,i))**2 !square to avoid negative norm
  end do
end do
GradNorm = Four*sqrt(GradNorm) !sqrt to complement previous **2

call mma_Deallocate(RMat)

if (Debug) then
    write(u6,*) ' '
    write(u6,*) 'In GetGrad_PM'
    write(u6,*) '-------------'
    write(u6,*) ' '
    write(u6,*) 'GradNorm = ',gradnorm
    ! Trace of Rmat is the value of the PM functional
    Fun=Zero
    do i=1,nOrb2Loc
        Fun = Fun+Rmat(i,i)
    end do
    call RecPrt('RMat',' ',RMat(:,:), nOrb2Loc, nOrb2Loc)
    write(u6,*) 'PM_functional = Tr(R) = ',Fun
end if



end subroutine GetGradNorm_PM
