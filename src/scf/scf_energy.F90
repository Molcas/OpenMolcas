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
      Subroutine SCF_Energy(FstItr,E1_,E2_,EV)
      Use Interfaces_SCF, Only: PMat_SCF
      use InfSCF, Only: nD
      use stdalloc, only: mma_allocate, mma_deallocate
      use MxDM, only: MxIter
      Implicit None
      Logical FstItr
      Real*8 E1_, E2_, EV

      Real*8, Allocatable :: XCf(:,:)
      Integer nXCf
!
!     Allocate memory for coefficients for minimized densities.
!
      nXCf = MxIter
      Call mma_allocate(XCf,nXCf,nD)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute
!     (1) density matrix,
!     (2) two-electron part of Fock matrix, and
!     (3) variational energies.
!                                                                      *
!***********************************************************************
!                                                                      *
!     (1) 1-particle density
!
      Call DMat(XCf,nXCf,nD)
!                                                                      *
!***********************************************************************
!                                                                      *
!     (2) 2-electron part of the Fock matrix, and the potential of the
!         external field (none linear or bi-linear).
!
!         Affects data in position iPsLst
!
      Call PMat_SCF(FstItr,XCf,nXCf,nD)
!
      Call mma_deallocate(XCf)
!                                                                      *
!***********************************************************************
!                                                                      *
!     (3) the energy.
!
      Call EneClc(E1_,E2_,EV)
!
      Return
      End Subroutine SCF_Energy
