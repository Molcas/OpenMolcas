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
  use rasscf_data, only : nAcPar, core_energy => Emy, nAc
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
!>  If it is the 1A or 1B iteration (0 <= actual_iter <= 1), the one electron
!>  energies are read from the InpOrb.
!>  Otherwise orbital_energies is a copy of DIAF.
!>
!>  @param[in] iter
!>  @param[in] DIAF
!>  @param[out] orbital_energies
  subroutine get_orbital_E(actual_iter, DIAF, orbital_energies)
    integer, intent(in) :: actual_iter
    real*8, intent(in) :: DIAF(:)
    real*8, intent(out) :: orbital_energies(:)

    orbital_energies = 0.d0
    if (0 <= actual_iter .or. actual_iter <= 1) then
      call read_orbital_energies(nSym, nBas, orbital_energies)
    else
      orbital_energies(:) = DIAF(:)
    end if
  contains
    subroutine read_orbital_energies(nSym, nBas, orbital_energies)
      integer, intent(in) :: nSym, nBas(:)
      real*8, intent(inout) :: orbital_energies(:)
      real*8 :: Dummy(1)
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
!>    Generate inactive Fock-matrix in active MO space
!>
!>  @author Oskar Weser
!>
!>  @details
!>  Generate the Fock-matrix for the frozen and inactive orbitals.
!>  in the basis of the active MOs as obtained from ::SGFCIN.
!>  Has the sideeffect of setting ::EMY to the core energy.
!>
!>  @param[in] CMO The MO coefficients.
!>  @param[in] D1I_AO The inactive one-body density matrix in AO-space
!>    \f[D^{\text{AO}, I} = 2 C (C^I)^\dagger \f]
!>    See ::get_D1I_rasscf.
!>  @param[in] D1A_AO The active one-body density matrix in AO-space
!>    \f[ D^{\text{AO}, A} = C^A D^A (C^A)^\dagger \f]
!>    See ::get_D1A_rasscf.
!>  @param[in] D1S_MO The active spin density matrix in MO-space
!>    \f[ D^A_\alpha - D^A_\beta \f]
!>  @param[inout] FI The inactive Fock matrix in AO-space
!>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu})\f]
!>    In output FI contains also the core energy added to
!>    the diagonal elements.
!>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu}) + \frac{E^{(0)}}{n_el} \delta_{\mu\nu} \f]
!>  @param[out] folded_Fock The inactive Fock matrix
!>    in the basis of the active MOs as obtained from ::SGFCIN.
  subroutine fold_Fock(CMO, D1I_AO, D1A_AO, D1S_MO, F_In, folded_Fock)
    real*8, intent(in) :: CMO(nTot2), D1A_AO(nTot2), D1I_AO(nTot2), D1S_MO(nAcPar)
    real*8, intent(inout) :: F_In(nTot1)
    real*8, intent(out) :: folded_Fock(nAcPar)
    integer :: i, n
    real*8 :: core_E_per_act_el,&
! one-body spin density matrix in AO-space
        D1S_AO(nTot2),&
! blocked one-body spin density matrix in MO-space
        D1S_MO_blocked(nAcPar)

    D1S_MO_blocked(:nAcPar) = D1S_MO(:nAcPar)
    IF (nAsh(1) /= nAc) call DBlock(D1S_MO_blocked)
    call get_D1A_RASSCF(CMO, D1S_MO_blocked, D1S_AO)
! SGFCIN has side effects and EMY/core_energy is set in this routine.
! Besides F_In will contain the one electron contribution afterwards.
    call SGFCIN(CMO, folded_fock, F_In, D1I_AO, D1A_AO, D1S_AO)

! Remove the core energy in the diagonal elements.
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
