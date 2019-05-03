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
! Copyright (C) 2014, Giovanni Li Manni                                *
!               2019, Oskar Weser                                      *
!***********************************************************************
module fcidump_transformations
  use stdalloc, only : mma_allocate, mma_deallocate

  use general_data, only : nActEl, nAsh, ntot, ntot1, ntot2, nBas, nSym
  use rasscf_data, only : nAcPar, core_energy => Emy
  use index_symmetry, only : one_el_idx_flatten
  implicit none
  private
  public :: get_orbital_E, fold_Fock
contains

!>  @brief
!>    Get the orbital energies.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  If it is the first iteration (iter == 1), the one electron
!>  energies are read from the InpOrb.
!>  Otherwise orbital_energies is a copy of DIAF.
!>
!>  @param[in] iter
!>  @param[in] DIAF
!>  @param[out] orbital_energies
  subroutine get_orbital_E(iter, DIAF, orbital_energies)
    implicit none
    integer, intent(in) :: iter
    real(kind=8), intent(in) :: DIAF(:)
    real(kind=8), intent(inout) :: orbital_energies(:)

    if (iter == 1) then
      call read_orbital_energies(nSym, nBas, orbital_energies)
    else
      orbital_energies(:) = DIAF(:)
    end if
  contains
    subroutine read_orbital_energies(nSym, nBas, orbital_energies)
      implicit none
      integer, intent(in) :: nSym, nBas(:)
      real(kind=8), intent(inout) :: orbital_energies(:)
      real(kind=8) :: Dummy(1)
      integer :: LuInpOrb = 10, iDummy(1), err
      character(*), parameter ::  FnInpOrb = 'INPORB'
      character(80) :: VecTit
      logical :: okay
      call f_Inquire(FnInpOrb, okay)
      if (okay) then
        call RdVec(FnInpOrb,LuInpOrb,'E',nSym,nBas,nBas, &
          Dummy, Dummy, orbital_energies, iDummy, &
          VecTit, 0, err)
      else
        Write (6,*) 'RdCMO: Error finding MO file'
        call QTrace()
        call Abend()
      end if
    end subroutine read_orbital_energies
  end subroutine get_orbital_E


!>  @brief
!>    Fold the Fock-Matrix
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The fock_table gets filled with the Fock matrix elements
!>  whose absolute value is larger than cutoff.
!>  The values are given by
!>  \f[ < i | F | j > \f]
!>  The index is given by i and j.
!>
!>  @param[in] CMO The MO coefficients.
!>  @param[in] D1I_MO The inactive one-body density matrix in MO-space
!>  @param[inout] F_In
!>  @param[out] folded_Fock
!>    \f[\sum_{\sigma\rho} {In}^D_{\sigma\rho} (g_{\mu\nu\sigma\rho})  \f]
! TODO(Oskar): write LaTeX
  subroutine fold_Fock(CMO, D1A_MO, F_In, folded_Fock)
    implicit none
    real(8), intent(in) :: CMO(nTot2), D1A_MO(nTot2)
    real(8), intent(inout) :: F_In(nTot1)
    real(8), intent(out) :: folded_Fock(:)
    integer :: i, n
    real(8) :: core_E_per_act_el
    real(8), allocatable :: &
! active one-body density matrix in AO-space
        D1A_AO(:),&
! inactive one-body density matrix in AO-space
        D1I_AO(:)

    call mma_allocate(D1A_AO, nTot2)
    call get_D1A_RASSCF(CMO, D1A_MO, D1A_AO)
    call mma_allocate(D1I_AO, ntot2)
    call get_D1A_RASSCF(CMO, D1I_AO)
! SGFCIN has side effects and EMY/core_energy is set in this routine.
! Besides F_In will contain the one electron contribution afterwards.
    call SGFCIN(CMO, folded_fock, F_In, D1I_AO, D1A_AO, D1A_AO)
    call mma_deallocate(D1A_AO)
    call mma_deallocate(D1I_AO)

! Remove the one electron contribution in the diagonal elements.
    if (nActEl /= 0) then
      core_E_per_act_el = core_energy / dble(nActEl)
    else
      core_E_per_act_el = 0.0d0
    end if

    do i = 1, sum(nAsh)
      n = one_el_idx_flatten(i, i)
      folded_Fock(n) = folded_Fock(n) - core_E_per_act_el
    end do
  end subroutine fold_Fock

end module fcidump_transformations
