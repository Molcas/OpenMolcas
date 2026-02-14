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

use general_data, only: nActEl, nAsh, nBas, nSym, ntot1, ntot2
use rasscf_global, only: nAc, nAcPar, Emy
use index_symmetry, only: one_el_idx_flatten
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

public :: fold_Fock, get_orbital_E

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
!>  @param[in] actual_iter
!>  @param[in] DIAF
!>  @param[out] orbital_energies
subroutine get_orbital_E(actual_iter,DIAF,orbital_energies)

  integer(kind=iwp), intent(in) :: actual_iter
  real(kind=wp), intent(in) :: DIAF(:)
  real(kind=wp), intent(out) :: orbital_energies(:)

  orbital_energies = Zero
  if ((0 <= actual_iter) .or. (actual_iter <= 1)) then
    call read_orbital_energies(nSym,nBas,orbital_energies)
  else
    orbital_energies(:) = DIAF(:)
  end if

  contains

  subroutine read_orbital_energies(nSym,nBas,orbital_energies)

    use Definitions, only: u6

    integer(kind=iwp), intent(in) :: nSym, nBas(:)
    real(kind=wp), intent(inout) :: orbital_energies(:)
    integer(kind=iwp) :: err, iDummy(1)
    real(kind=wp) :: Dummy(1)
    character(len=80) :: VecTit
    logical(kind=iwp) :: okay
    integer(kind=iwp), parameter :: LuInpOrb = 10
    character(len=*), parameter :: FnInpOrb = 'INPORB'

    call f_Inquire(FnInpOrb,okay)
    if (okay) then
      call RdVec(FnInpOrb,LuInpOrb,'E',nSym,nBas,nBas,Dummy,Dummy,orbital_energies,iDummy,VecTit,0,err)
    else
      write(u6,*) 'RdCMO: Error finding MO file'
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
!>  Has the sideeffect of setting \p EMY to the core energy.
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
!>  @param[in,out] F_In The inactive Fock matrix in AO-space
!>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu})\f]
!>    In output FI contains also the core energy added to
!>    the diagonal elements.
!>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} -
!>      \frac{1}{2} g_{\mu\sigma\rho\nu}) + \frac{E^{(0)}}{n_el} \delta_{\mu\nu} \f]
!>  @param[out] folded_Fock The inactive Fock matrix
!>    in the basis of the active MOs as obtained from ::SGFCIN.
subroutine fold_Fock(CMO,D1I_AO,D1A_AO,D1S_MO,F_In,folded_Fock)

  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), intent(in) :: CMO(nTot2), D1A_AO(nTot2), D1I_AO(nTot2), D1S_MO(nAcPar)
  real(kind=wp), intent(inout) :: F_In(nTot1)
  real(kind=wp), intent(out) :: folded_Fock(nAcPar)
  integer(kind=iwp) :: i, n
  real(kind=wp) :: core_E_per_act_el
  real(kind=wp), allocatable :: D1S_AO(:), D1S_MO_blocked(:)

  ! D1S_AO: one-body spin density matrix in AO-space
  ! D1S_MO: blocked one-body spin density matrix in MO-space
  call mma_allocate(D1S_AO,nTot2,Label='D1S_AO')
  call mma_allocate(D1S_MO_blocked,nAcPar,Label='D1S_MO_blocked')

  D1S_MO_blocked(:nAcPar) = D1S_MO(:nAcPar)
  if (nAsh(1) /= nAc) call DBlock(D1S_MO_blocked)
  call get_D1A_RASSCF(CMO,D1S_MO_blocked,D1S_AO)
  ! SGFCIN has side effects and EMY/core_energy is set in this routine.
  ! Besides F_In will contain the one electron contribution afterwards.
  call SGFCIN(CMO,folded_fock,F_In,D1I_AO,D1A_AO,D1S_AO)

  call mma_deallocate(D1S_AO)
  call mma_deallocate(D1S_MO_blocked)

  ! Remove the core energy in the diagonal elements.
  if (nActEl /= 0) then
    core_E_per_act_el = Emy/real(nActEl,kind=wp)
  else
    core_E_per_act_el = Zero
  end if

  do i=1,sum(nAsh)
    n = one_el_idx_flatten(i,i)
    folded_Fock(n) = folded_Fock(n)-core_E_per_act_el
  end do

end subroutine fold_Fock

end module fcidump_transformations
