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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Transpose_TDM
!
!> @brief
!>   Transpose a transition density matrix in place.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Transpose a transition density matrix, stored in symmetry blocks,
!> replacing the original matrix. The matrices contain only the symmetry
!> blocks that match the total symmetry of the transition.
!>
!> @param[in,out] TDM       Transition density matrix
!> @param[in]     Symmetry  Symmetry of the transition
!***********************************************************************

subroutine Transpose_TDM(TDM,Symmetry)

use Symmetry_Info, only: Mul, nIrrep
use rassi_data, only: NBASF
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: TDM(*)
integer(kind=iwp), intent(in) :: Symmetry
integer(kind=iwp) :: i, iBlock(0:8), iSym1, iSym2, j, nTot
real(kind=wp), allocatable :: Tmp(:)

! Compute the location of all the stored symmetry blocks
nTot = 0
iBlock(0) = 0
do iSym1=1,nIrrep
  nTot = nTot+nBasF(iSym1)*nBasF(Mul(Symmetry,iSym1))
  iBlock(iSym1) = nTot
end do
! Make a copy so we can transpose in place
call mma_Allocate(Tmp,nTot,Label='Tmp')
Tmp(:) = TDM(1:nTot)
! Transpose symmetry block (a,b) onto symmetry block (b,a)
do iSym1=1,nIrrep
  iSym2 = Mul(Symmetry,iSym1)
  do i=1,nBasF(iSym2)
    do j=1,nBasF(iSym1)
      TDM(iBlock(iSym2-1)+(j-1)*nBasF(iSym2)+i) = Tmp(iBlock(iSym1-1)+(i-1)*nBasF(iSym1)+j)
    end do
  end do
end do
call mma_deAllocate(Tmp)

end subroutine Transpose_TDM
