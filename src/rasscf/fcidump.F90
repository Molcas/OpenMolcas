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
!               2026, Nike Dattani                                     *
!               2026, Jaafar Mehrez                                    *
!***********************************************************************

!> @brief
!>   Master module for fcidump.
module fcidump

use rasscf_global, only: nacpar
use general_data, only: nTot, nTot1, nTot2
use fcidump_tables, only: fill_2ElInt, fill_fock, fill_orbitals, FockTable, mma_allocate, mma_deallocate, OrbitalTable, &
                          TwoElIntTable
use fcidump_transformations, only: fold_Fock, get_orbital_E
use fcidump_reorder, only: reorder
use fcidump_dump, only: dump_ascii, dump_fort55, dump_hdf5
use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: DumpOnly = .false.
integer(kind=iwp) :: DmpMode = 0 ! 0=FCIDUMP only, 1=fort55 only

public :: cleanup, DumpOnly, DmpMode, make_fcidumps, transform

contains

subroutine make_fcidumps(ascii_path,h5_path,orbital_energies,folded_Fock,TUVX,core_energy,permutation,fort55_path)

  use general_data, only: nSym, nAsh

  character(len=*), intent(in) :: ascii_path, h5_path
  real(kind=wp), intent(in) :: orbital_energies(:), folded_Fock(:), TUVX(:), core_energy
  integer(kind=iwp), intent(in), optional :: permutation(:)
  character(len=*), intent(in), optional :: fort55_path
  integer(kind=iwp) :: i, j, n
  integer(kind=iwp), allocatable :: energy_perm(:), inv_perm(:)
  type(OrbitalTable) :: orbital_table
  type(FockTable) :: fock_table
  type(TwoElIntTable) :: two_el_table
  integer(kind=iwp), allocatable :: orbsym(:)

  call mma_allocate(orbsym,sum(nAsh(:nSym)))
  call mma_allocate(fock_table,nacpar)
  call mma_allocate(two_el_table,size(TUVX))
  call mma_allocate(orbital_table,sum(nAsh))

  call fill_orbitals(orbital_table,orbital_energies)
  call fill_fock(fock_table,folded_Fock)
  call fill_2ElInt(two_el_table,TUVX)

  n = 1
  do j=1,nSym
    orbsym(n:n+nAsh(j)-1) = j
    n = n+nAsh(j)
  end do

  if (present(permutation)) call reorder(orbital_table,fock_table,two_el_table,orbsym,permutation)

  if (present(fort55_path) .and. DmpMode == 1) then
    call mma_allocate(energy_perm,sum(nAsh(:nSym)),Label='energy_perm')
    call mma_allocate(inv_perm,sum(nAsh(:nSym)),Label='inv_perm')
    call energy_sort_permutation(orbital_energies,energy_perm)
    call reorder(orbital_table,fock_table,two_el_table,orbsym,energy_perm)
    call dump_fort55(fort55_path,core_energy,orbital_table,fock_table,two_el_table,orbsym)
    do i=1,size(energy_perm)
      inv_perm(energy_perm(i)) = i
    end do
    call reorder(orbital_table,fock_table,two_el_table,orbsym,inv_perm)
    call mma_deallocate(energy_perm)
    call mma_deallocate(inv_perm)
  end if

  if (DmpMode == 0) then
    call dump_ascii(ascii_path,core_energy,orbital_table,fock_table,two_el_table,orbsym)
    call dump_hdf5(h5_path,core_energy,orbital_table,fock_table,two_el_table,orbsym)
  end if

  call mma_deallocate(orbsym)
  call mma_deallocate(fock_table)
  call mma_deallocate(two_el_table)
  call mma_deallocate(orbital_table)

end subroutine make_fcidumps

!> @brief
!>  Compute a permutation that sorts active orbitals by energy.
!>
!> @author Jaafar Mehrez
!>
!> @details
!>   MRCC expect orbitals in energy order. This routine extracts
!>   the active orbital energies and returns a permutation array
!>   permutation(i) = new position (1-based) of original orbital i.
subroutine energy_sort_permutation(orbital_energies,permutation)

  use general_data, only: nSym, nAsh, nBas, nFro, nIsh
  use Definitions, only: wp, iwp

  real(kind=wp), intent(in) :: orbital_energies(:)
  integer(kind=iwp), intent(out) :: permutation(:)

  integer(kind=iwp) :: nmo, i, j, isym, offset, count
  integer(kind=iwp), allocatable :: idx_map(:)
  real(kind=wp), allocatable :: active_energies(:)

  nmo = sum(nAsh(:nSym))
  call mma_allocate(active_energies,nmo,Label='active_energies')
  call mma_allocate(idx_map,nmo,Label='idx_map')

  count = 1
  offset = 1
  do isym=1,nSym
    do i=1,nAsh(isym)
      active_energies(count) = orbital_energies(offset + nFro(isym) + nIsh(isym) + i - 1)
      idx_map(count) = count
      count = count + 1
    end do
    offset = offset + nBas(isym)
  end do

  do i=1,nmo-1
    do j=i+1,nmo
      if (active_energies(idx_map(i)) > active_energies(idx_map(j)) + 1.0e-10_wp) then
        offset = idx_map(i)
        idx_map(i) = idx_map(j)
        idx_map(j) = offset
      else if (abs(active_energies(idx_map(i)) - active_energies(idx_map(j))) <= 1.0e-10_wp) then
        if (idx_map(i) > idx_map(j)) then
          offset = idx_map(i)
          idx_map(i) = idx_map(j)
          idx_map(j) = offset
        end if
      end if
    end do
  end do

  do i=1,nmo
    permutation(idx_map(i)) = i
  end do

  call mma_deallocate(active_energies)
  call mma_deallocate(idx_map)

end subroutine energy_sort_permutation

subroutine transform(actual_iter,CMO,DIAF,D1I_AO,D1A_AO,D1S_MO,F_IN,orbital_E,folded_Fock)

  integer(kind=iwp), intent(in) :: actual_iter
  real(kind=wp), intent(in) :: CMO(nTot2), DIAF(nTot), D1I_AO(nTot2), D1A_AO(nTot2), D1S_MO(nAcPar)
  real(kind=wp), intent(inout) :: F_IN(nTot1)
  real(kind=wp), intent(out) :: orbital_E(nTot), folded_Fock(nAcPar)

  call get_orbital_E(actual_iter,DIAF,orbital_E)
  call fold_Fock(CMO,D1I_AO,D1A_AO,D1S_MO,F_In,folded_Fock)

end subroutine transform

subroutine cleanup()

  use fcidump_reorder, only: fcidump_reorder_cleanup => cleanup

  call fcidump_reorder_cleanup()

end subroutine cleanup

end module fcidump
