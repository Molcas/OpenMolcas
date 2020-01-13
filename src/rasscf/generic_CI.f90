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
! Copyright (C) 2020, Oskar Weser                                      *
!***********************************************************************


module generic_CI
    implicit none
    private
    public :: CI_solver_factory

    type, abstract :: CI_runner
    end type

    type, abstract :: CI_runner_args
    end type

contains


!>  @brief
!>    Create a CI-solver from building blocks.
!>
!>  @author Oskar Weser
!>
!>  @paramin[in] CMO MO coefficients
!>  @paramin[in] DIAF DIAGONAL of Fock matrix useful for NECI
!>  @paramin[in] D1I_MO Inactive 1-dens matrix
!>  @paramin[in] TUVX Active 2-el integrals
!>  @paramin[inout] F_In Fock matrix from inactive density
!>  @paramin[inout] D1S_MO Average spin 1-dens matrix
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
!>  @paramin[in] fake_run  If true the CI step is not performed, but
!>    the RDMs are read from previous run.
      subroutine CI_solver_factory(
     &          CMO, DIAF, D1I_AO, D1A_AO, TUVX, F_IN,
     &          D1S_MO, DMAT, PSMAT, PAMAT)
      use general_data, only : iSpin, ntot, ntot1, ntot2, nAsh, nBas
      use rasscf_data, only : iter, lRoots, nRoots, S, KSDFT, EMY,
     &    rotmax, Ener, Nac, nAcPar, nAcpr2

      use gugx_data, only : IfCAS
      use gas_data, only : ngssh, iDoGas, nGAS, iGSOCCX

      end subroutine

end module
