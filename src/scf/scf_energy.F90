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

subroutine SCF_Energy(FstItr,E1_,E2_,EV)

use Interfaces_SCF, only: PMat_SCF
use InfSCF, only: nD
use stdalloc, only: mma_allocate, mma_deallocate
use MxDM, only: MxIter

implicit none
logical FstItr
real*8 E1_, E2_, EV
real*8, allocatable :: XCf(:,:)
integer nXCf

! Allocate memory for coefficients for minimized densities.

nXCf = MxIter
call mma_allocate(XCf,nXCf,nD)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute
! (1) density matrix,
! (2) two-electron part of Fock matrix, and
! (3) variational energies.
!                                                                      *
!***********************************************************************
!                                                                      *
! (1) 1-particle density

call DMat(XCf,nXCf,nD)
!                                                                      *
!***********************************************************************
!                                                                      *
! (2) 2-electron part of the Fock matrix, and the potential of the
!     external field (none linear or bi-linear).

!     Affects data in position iPsLst

call PMat_SCF(FstItr,XCf,nXCf,nD)

call mma_deallocate(XCf)
!                                                                      *
!***********************************************************************
!                                                                      *
! (3) the energy.

call EneClc(E1_,E2_,EV)

return

end subroutine SCF_Energy
