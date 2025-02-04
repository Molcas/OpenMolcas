!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Mk_iSO2Ind(iSO2Sh,iSO2Ind,nSO,nShell)

use Basis_Info, only: nBas_Aux
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSO, iSO2Sh(nSO), nShell
integer(kind=iwp), intent(out) :: iSO2Ind(nSO)
integer(kind=iwp) :: iB, iIrrep, iSh, iSO
integer(kind=iwp), allocatable :: nTemp(:)

call mma_allocate(nTemp,nShell,Label='nTemp')

iSO = 0
! Loop over the irreps
do iIrrep=0,nIrrep-1

  nTemp(:) = 0
  ! Loop over the canonical index of the auxiliary functions, starting from 1 for the first irrep.
  do iB=1,nBas_Aux(iIrrep)
    iSO = iSO+1
    iSh = iSO2Sh(iSO)  ! Pick up the corresponding shell index
    nTemp(iSh) = nTemp(iSh)+1  ! Increment the count of basis functions for shell ish
    ! This table translates the global index of a basis function, iSO, to the local index of the shell to which it belongs.
    iSO2Ind(iSO) = nTemp(iSh)
  end do

end do

call mma_deallocate(nTemp)

end subroutine Mk_iSO2Ind
