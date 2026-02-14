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

!> @brief
!>   Master module for fcidump.
module fcidump

use rasscf_global, only: nacpar
use general_data, only: nTot, nTot1, nTot2
use fcidump_tables, only: fill_2ElInt, fill_fock, fill_orbitals, FockTable, mma_allocate, mma_deallocate, OrbitalTable, &
                          TwoElIntTable
use fcidump_transformations, only: fold_Fock, get_orbital_E
use fcidump_reorder, only: reorder
use fcidump_dump, only: dump_ascii, dump_hdf5
use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: DumpOnly = .false.

public :: cleanup, DumpOnly, make_fcidumps, transform

contains

subroutine make_fcidumps(ascii_path,h5_path,orbital_energies,folded_Fock,TUVX,core_energy,permutation)

  use general_data, only: nSym, nAsh

  character(len=*), intent(in) :: ascii_path, h5_path
  real(kind=wp), intent(in) :: orbital_energies(:), folded_Fock(:), TUVX(:), core_energy
  integer(kind=iwp), intent(in), optional :: permutation(:)
  integer(kind=iwp) :: j, n
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

  call dump_ascii(ascii_path,core_energy,orbital_table,fock_table,two_el_table,orbsym)
  call dump_hdf5(h5_path,core_energy,orbital_table,fock_table,two_el_table,orbsym)

  call mma_deallocate(orbsym)
  call mma_deallocate(fock_table)
  call mma_deallocate(two_el_table)
  call mma_deallocate(orbital_table)

end subroutine make_fcidumps

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
