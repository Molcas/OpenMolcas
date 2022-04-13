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

#include "macros.fh"

!> This module defines an abstract class for CI-solvers.
!> If you inherit from CI_solver_t and override the deferred methods,
!> your initialization and cleanup will be automatically called.
module generic_CI
    use general_data, only : ntot, ntot1, ntot2
    use rasscf_data, only : nAcPar, nAcpr2, nroots
    use definitions, only: wp
    implicit none
    private
    public :: CI_solver_t

    type, abstract :: CI_solver_t
    contains
      procedure(CI_run_t), deferred :: run
      procedure(CI_cleanup_t), deferred :: cleanup
    end type

    abstract interface
!>  @brief
!>      Interface for procedure pointers to CI-solvers
!>
!>  @author Oskar Weser
!>
!>  @param[in] actual_iter The actual iteration number starting at 0.
!>      This means 0 is 1A, 1 is 1B, 2 is 2 and so on.
!>  @param[in] iroot specified roots for SA-CASSCF, e.g. 1,3,9,...
!>  @param[in] weight weights specified for roots for SA-CASSCF
!>  @param[in] CMO MO coefficients
!>  @param[in] DIAF DIAGONAL of Fock matrix useful for NECI
!>  @param[in] D1I_MO Inactive 1-dens matrix
!>  @param[in] TUVX Active 2-el integrals
!>  @param[in,out] F_In Fock matrix from inactive density
!>  @param[in,out] D1S_MO Average spin 1-dens matrix
!>  @param[out] DMAT Average 1 body density matrix
!>  @param[out] PSMAT Average symm. 2-dens matrix
!>  @param[out] PAMAT Average antisymm. 2-dens matrix

        subroutine CI_run_t(this, actual_iter, ifinal, iroot, weight, &
                            CMO, DIAF, D1I_AO, D1A_AO, TUVX, F_IN, &
                            D1S_MO, DMAT, PSMAT, PAMAT)
            import :: CI_solver_t, ntot, ntot1, ntot2, nAcPar, nAcpr2, nroots,&
                      wp

            class(CI_solver_t), intent(in) :: this
            integer, intent(in) :: actual_iter, iroot(nroots), ifinal
            real(wp), intent(in) :: weight(nroots), &
                                  CMO(nTot2), DIAF(nTot), D1I_AO(nTot2), &
                                  D1A_AO(nTot2), TUVX(nAcpr2)
            real(wp), intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
            real(wp), intent(out) :: DMAT(nAcpar), PSMAT(nAcpr2), PAMAT(nAcpr2)
        end subroutine

!>  @brief
!>      Interface to cleanup routine for CI-solvers
!>
!>  @details
!>  It would be nice to have automatic finalization.
!>  Unfortunately we have to support old compilers.
!>
!>  This method will be called at the end of rasscf.
!>
!>  @author Oskar Weser
        subroutine CI_cleanup_t(this)
            import :: CI_solver_t
            class(CI_solver_t), intent(inout) :: this
        end subroutine

    end interface

    contains

end module generic_CI
