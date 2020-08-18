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
! Copyright (C) 2019, Oskar Weser                                      *
!***********************************************************************



!>  @brief
!>    Master module for fcidump.
module fcidump
  use rasscf_data, only : nacpar
  use general_data, only : nAsh, nTot, nTot1, nTot2
  use fcidump_tables, only : OrbitalTable, FockTable, TwoElIntTable,&
    mma_allocate, mma_deallocate, fill_orbitals, fill_fock, fill_2ElInt
  use fcidump_transformations, only : get_orbital_E, fold_Fock
  use fcidump_reorder, only : reorder
  use fcidump_dump, only : dump_ascii, dump_hdf5
  implicit none
  private
  public :: make_fcidumps, transform, DumpOnly, cleanup
  logical :: DumpOnly = .false.
  save
contains

  subroutine make_fcidumps(ascii_path, h5_path, orbital_energies, folded_Fock,&
                           TUVX, core_energy, permutation)
    use general_data, only : nSym, nAsh
    character(*), intent(in) :: ascii_path, h5_path
    real*8, intent(in) :: orbital_energies(:), folded_Fock(:), TUVX(:), core_energy
    integer, intent(in), optional :: permutation(:)
    type(OrbitalTable) :: orbital_table
    type(FockTable) :: fock_table
    type(TwoElIntTable) :: two_el_table
    integer :: orbsym(sum(nAsh(:nSym))), n, j

    call mma_allocate(fock_table, nacpar)
    call mma_allocate(two_el_table, size(TUVX))
    call mma_allocate(orbital_table, sum(nAsh))

    call fill_orbitals(orbital_table, orbital_energies)
    call fill_fock(fock_table, folded_Fock)
    call fill_2ElInt(two_el_table, TUVX)

    n = 1
    do j = 1, nSym
      orbsym(n : n + nAsh(j) - 1) = j
      n = n + nAsh(j)
    end do

    if (present(permutation)) then
      call reorder(orbital_table, fock_table, two_el_table, orbsym, permutation)
    end if

    call dump_ascii(ascii_path, core_energy, orbital_table, fock_table, &
                    & two_el_table, orbsym)
    call dump_hdf5(h5_path, core_energy, orbital_table, fock_table, &
                    & two_el_table, orbsym)

    call mma_deallocate(fock_table)
    call mma_deallocate(two_el_table)
    call mma_deallocate(orbital_table)
  end subroutine make_fcidumps

  subroutine transform(actual_iter, CMO, DIAF, D1I_AO, D1A_AO, D1S_MO, &
                       F_IN, orbital_E, folded_Fock)
    integer, intent(in) :: actual_iter
    real*8, intent(in) :: DIAF(nTot),&
      CMO(nTot2),&
      D1I_AO(nTot2),&
      D1A_AO(nTot2),&
      D1S_MO(nAcPar)
    real*8, intent(inout) :: F_IN(nTot1)
    real*8, intent(out) :: orbital_E(nTot), folded_Fock(nAcPar)

    call get_orbital_E(actual_iter, DIAF, orbital_E)
    call fold_Fock(CMO, D1I_AO, D1A_AO, D1S_MO, F_In, folded_Fock)
  end subroutine transform

  subroutine cleanup()
    use fcidump_reorder, only : fcidump_reorder_cleanup => cleanup
    call fcidump_reorder_cleanup()
  end subroutine
end module fcidump
