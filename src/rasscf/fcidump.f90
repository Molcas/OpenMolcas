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
module fcidump
  use stdalloc, only : mma_allocate, mma_deallocate

  use general_data, only : nSym, nAsh, ntot, ntot1, ntot2
  use rasscf_data, only : core_energy => EMY, nAcPr2, nAcPar
  use gas_data, only : ngssh

  use fcidump_transformations, only : get_orbital_E, fold_Fock
  use fcidump_reorder, only : get_P_GAS, get_P_inp, ReOrInp, ReOrFlag
  use fcidump_dump, only : make_fcidumps

  implicit none
  logical :: DumpOnly = .false.
  private
  public :: transform_and_dump, DumpOnly
  save

contains

  subroutine transform_and_dump(iter, CMO, DIAF, D1I_MO, TUVX, F_IN)
    implicit none
    integer, intent(in) :: iter
    real(kind=8), intent(in) :: CMO(nTot2), DIAF(nTot), TUVX(nAcpr2),&
        D1I_MO(nTot2)
    real(kind=8), intent(inout) :: F_In(nTot1)

    real(kind=8), allocatable ::  folded_Fock(:), orbital_E(:)
    integer :: permutation(sum(nAsh(:nSym)))

    select case (ReOrFlag)
      case (2:)
        permutation = get_P_inp(ReOrInp)
        call mma_deallocate(ReOrInp)
      case (-1)
        permutation = get_P_GAS(nGSSH)
    end select

    call mma_allocate(orbital_E, size(DIAF))
    call get_orbital_E(iter, DIAF, orbital_E)

    call mma_allocate(folded_Fock, nAcPar)
! fold_fock works inplace on F_In
    call fold_Fock(CMO, F_In, D1I_MO, folded_Fock)

    if (ReOrFlag /= 0) then
      call make_fcidumps(orbital_E, folded_Fock, TUVX, core_energy, permutation)
    else
      call make_fcidumps(orbital_E, folded_Fock, TUVX, core_energy)
    end if

    call mma_deallocate(orbital_E)
    call mma_deallocate(folded_Fock)
  end subroutine transform_and_dump

end module fcidump
