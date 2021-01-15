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
!> I you inherit from CI_solver_t and override the deferred methods,
!> your initialization and cleanup will be automatically called.
module generic_CI
    use definitions, only: wp
    use general_data, only : ntot, ntot1, ntot2
    use rasscf_data, only : nAcPar, nAcpr2
    implicit none
    private
    public :: CI_solver_t, unused

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
!>  @paramin[in] actual_iter The actual iteration number starting at 0.
!>      This means 0 is 1A, 1 is 1B, 2 is 2 and so on.
!>  @paramin[in] CMO MO coefficients
!>  @paramin[in] DIAF DIAGONAL of Fock matrix useful for NECI
!>  @paramin[in] D1I_MO Inactive 1-dens matrix
!>  @paramin[in] TUVX Active 2-el integrals
!>  @paramin[inout] F_In Fock matrix from inactive density
!>  @paramin[inout] D1S_MO Average spin 1-dens matrix
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
        subroutine CI_run_t(this, actual_iter, CMO, DIAF, D1I_AO, D1A_AO, TUVX, &
                            F_IN, D1S_MO, DMAT, PSMAT, PAMAT)
            import :: CI_solver_t, ntot, ntot1, ntot2, nAcPar, nAcpr2

            class(CI_solver_t), intent(in) :: this
            integer, intent(in) :: actual_iter
            real*8, intent(in) :: CMO(nTot2), DIAF(nTot), D1I_AO(nTot2), &
                                  D1A_AO(nTot2), TUVX(nAcpr2)
            real*8, intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
            real*8, intent(out) :: DMAT(nAcpar), PSMAT(nAcpr2), PAMAT(nAcpr2)
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

    subroutine unused(CI_solver)
        class(CI_solver_t), intent(in) :: CI_solver
        unused_var(CI_solver)

        block
            real(wp), parameter :: PI = 4._wp * atan(1._wp)
            write(*, *) my_cos([0._wp, 2._wp * PI])
        end block

        contains
            real(wp) impure elemental function my_cos(x)
                real(wp), intent(in) :: x
                my_cos = cos(x)
            end function
    end subroutine

end module generic_CI
