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
  use fcidump_tables
  use fcidump_transformations
  use fcidump_reorder
  use fcidump_dump
  implicit none
  private
  public :: make_fcidumps
  save
contains

  subroutine make_fcidumps(iter, nacpar, nAsh, TUVX, DIAF, &
        CMO, F_IN, D1I_MO, core_energy, permutation)
    implicit none
    integer, intent(in) :: iter, nacpar, nAsh(:)
    real(8), intent(in) :: TUVX(:), DIAF(:), core_energy, &
      CMO(:), F_IN(:), D1I_MO(:)
    integer, intent(in), optional :: permutation(:)
    real(8), allocatable :: orbital_energies(:), folded_Fock(:)
    type(OrbitalTable) :: orbital_table
    type(FockTable) :: fock_table
    type(TwoElIntTable) :: two_el_table

    call mma_allocate(fock_table, nacpar)
    call mma_allocate(two_el_table, size(TUVX))
    call mma_allocate(orbital_table, sum(nAsh))

    call mma_allocate(orbital_energies, size(DIAF))
    call get_orbital_E(iter, DIAF, orbital_energies)
    call fill_orbitals(orbital_table, orbital_energies)
    call mma_deallocate(orbital_energies)

    call mma_allocate(folded_Fock, nAcPar)
    call get_folded_Fock(CMO, F_In, D1I_MO, folded_Fock)
    call fill_fock(fock_table, folded_Fock)
    call mma_deallocate(folded_Fock)

    call fill_2ElInt(two_el_table, TUVX)

    if (present(permutation)) then
      call reorder(orbital_table, fock_table, two_el_table, permutation)
    end if

    call dump_ascii(core_energy, orbital_table, fock_table, two_el_table)
    call dump_hdf5(core_energy, orbital_table, fock_table, two_el_table)

    call mma_deallocate(fock_table)
    call mma_deallocate(two_el_table)
    call mma_deallocate(orbital_table)
  end subroutine make_fcidumps
end module fcidump
