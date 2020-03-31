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
module fcidump_tables
  use stdalloc, only : mma_allocate, mma_deallocate
  use index_symmetry, only : one_el_idx, two_el_idx
  implicit none
  private
  public :: FockTable, TwoElIntTable, OrbitalTable, mma_allocate, &
    mma_deallocate, length, print, fill_orbitals, fill_fock, fill_2ElInt, &
    cutoff_default, unused
  save

  type :: FockTable
    sequence
    real*8, allocatable, dimension(:) :: values ! <i | F | j >
    integer, allocatable, dimension(:, :) :: index ! i, j
    real*8 :: cutoff
    integer :: length
  end type FockTable

  type :: TwoElIntTable
    sequence
    real*8, allocatable, dimension(:) :: values ! <ij| 1/r_{12} |kl>
    integer, allocatable, dimension(:, :) :: index ! i, j, k, l
    real*8 :: cutoff
    integer :: length
  end type TwoElIntTable

  type :: OrbitalTable
    sequence
    real*8, allocatable, dimension(:) :: values ! <i| F |i>
    integer, allocatable, dimension(:) :: index ! i
  end type OrbitalTable

  real*8, parameter :: cutoff_default = 1.0d-11

  interface mma_allocate
    module procedure FockTable_allocate, TwoElIntTable_allocate, &
      OrbitalTable_allocate
  end interface

  interface mma_deallocate
    module procedure FockTable_deallocate, TwoElIntTable_deallocate, &
      OrbitalTable_deallocate
  end interface

  interface length
    module procedure FockTable_length, TwoElIntTable_length, &
      OrbitalTable_length
  end interface

  interface print
    module procedure FockTable_print, TwoElIntTable_print, OrbitalTable_print
  end interface

  interface unused
    module procedure FockTable_unused, TwoElIntTable_unused, OrbitalTable_unused
  end interface

contains

  subroutine OrbitalTable_allocate(orbital_table, n)
    implicit none
    integer, intent(in) :: n
    type(OrbitalTable), intent(inout) :: orbital_table
    call mma_allocate(orbital_table%values, n)
    call mma_allocate(orbital_table%index, n)
  end subroutine OrbitalTable_allocate

  subroutine OrbitalTable_deallocate(orbital_table)
    implicit none
    type(OrbitalTable), intent(inout) :: orbital_table
    call mma_deallocate(orbital_table%values)
    call mma_deallocate(orbital_table%index)
  end subroutine OrbitalTable_deallocate

!>  @brief
!>    Fill in orbital energies
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The orbitals table gets filled with the orbital energies from DIAF.
!>  If it is the first iteration (iter == 1) then the one electron
!>  energies are read from the InpOrb.
!>
!>  @param[in,out] orbitals Core
!>  @param[in] DIAF
!>  @param[in] iter
  subroutine fill_orbitals(table, orbital_energies)
    use general_data, only : nBas, nSym, nAsh, nFro, nIsh
    implicit none
    type(OrbitalTable), intent(inout) :: table
    real*8, intent(in) :: orbital_energies(:)
    integer :: i, n, iSym, iOff

    iOff = 0
    n = 1
    do iSym = 1, nSym
      if (nAsh(iSym) > 0) then
        do i = 1, nAsh(isym)
          table%index(n) = n
          table%values(n) = orbital_energies(ioff + nFro(iSym) + nIsh(iSym) + i)
          n = n + 1
        enddo
      end if
      iOff   = iOff + nBas(iSym)
    end do
  end subroutine fill_orbitals

  pure integer function OrbitalTable_length(table)
    type(OrbitalTable), intent(in) :: table
    OrbitalTable_length = size(table%values)
  end function OrbitalTable_length

  subroutine OrbitalTable_print(table)
    type(OrbitalTable), intent(in) :: table
    integer :: i
    do i = 1, length(table)
      write(6, '(E15.7, I7)') table%values(i), table%index(i)
    end do
  end subroutine OrbitalTable_print

  subroutine OrbitalTable_unused(table)
    type(OrbitalTable), intent(in) :: table
    integer :: n
    if (.false.) n = length(table)
  end subroutine OrbitalTable_unused

  subroutine FockTable_allocate(fock_table, n)
    implicit none
    integer, intent(in) :: n
    type(FockTable), intent(inout) :: fock_table
    call mma_allocate(fock_table%values, n)
    call mma_allocate(fock_table%index, 2, n)
  end subroutine FockTable_allocate

  subroutine FockTable_deallocate(fock_table)
    implicit none
    type(FockTable), intent(inout) :: fock_table
    call mma_deallocate(fock_table%values)
    call mma_deallocate(fock_table%index)
  end subroutine FockTable_deallocate

!>  @brief
!>    Fill in fock matrix elements.
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
!>  @param[in,out] fock_table
!>  @param[in] CMO The occupation number vector in MO-space.
!>  @param[in] F_In
!>    \f[\sum_{\sigma\rho} {In}^D_{\sigma\rho} (g_{\mu\nu\sigma\rho})  \f]
!>  @param[in] D1I_MO The inactive one-body density matrix in MO-space
!>  @param[in] cutoff Optional parameter that is set by default to
!>    fciqmc_tables::cutoff_default.
  subroutine fill_fock(fock_table, Fock, cutoff)
    use general_data, only : nActEl, nAsh, ntot, ntot1, ntot2
    use rasscf_data, only : nAcPar
    implicit none
    real*8, intent(in) :: Fock(:)
    type(FockTable), intent(inout) :: fock_table
    real*8, optional, intent(in) :: cutoff

    integer :: i, n
    real*8 :: cutoff_

    cutoff_ = merge(cutoff, cutoff_default, present(cutoff))

    n = 0
    do i = 1, size(Fock)
      if (abs(Fock(i)) >= cutoff_) then
        n = n + 1
        call one_el_idx(i, fock_table%index(:, n))
        fock_table%values(n) = Fock(i)
      end if
    end do
    fock_table%length = n
    fock_table%cutoff = cutoff_
  end subroutine fill_fock

  pure integer function FockTable_length(table)
    implicit none
    type(FockTable), intent(in) :: table
    FockTable_length = table%length
  end function FockTable_length

  subroutine FockTable_print(table)
    type(FockTable), intent(in) :: table
    integer :: i, j
    do j = 1, length(table)
      write(6, '(E15.7, I7, I7)') table%values(j), (table%index(i, j), i=1, 2)
    end do
  end subroutine FockTable_print

  subroutine FockTable_unused(table)
    type(FockTable), intent(in) :: table
    integer :: n
    if (.false.) n = length(table)
  end subroutine FockTable_unused

  subroutine TwoElIntTable_allocate(table, n)
    implicit none
    integer, intent(in) :: n
    type(TwoElIntTable), intent(inout) :: table
    call mma_allocate(table%values, n)
    call mma_allocate(table%index, 4, n)
  end subroutine TwoElIntTable_allocate

  subroutine TwoElIntTable_deallocate(table)
    implicit none
    type(TwoElIntTable), intent(inout) :: table
    call mma_deallocate(table%values)
    call mma_deallocate(table%index)
  end subroutine TwoElIntTable_deallocate

!>  @brief
!>    Fill in two electron integrals
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The two_el_table gets filled with those two electron integrals
!>  whose absolute value is larger than cutoff.
!>  The values are given by
!>  \f[ < i, j | \frac{1}{r_{1,2}} | k, l > \f]
!>  The index is given by i, j, k, and l.
!>
!>  @param[in,out] two_el_table
!>  @param[in] TUVX
!>  @param[in] cutoff Optional parameter that is set by default to
!>    fciqmc_tables::cutoff_default.
  subroutine fill_2ElInt(two_el_table, TUVX, cutoff)
    implicit none
    real*8, intent(in) :: TUVX(:)
    type(TwoElIntTable), intent(inout) :: two_el_table
    integer :: i, n
    real*8, optional :: cutoff
    real*8 :: cutoff_

    integer, parameter :: max_test = 20
    integer :: l_twoel_test

    cutoff_ = merge(cutoff, cutoff_default, present(cutoff))

    n = 0
    do i = 1, size(TUVX)
      if (abs(TUVX(i)) >= cutoff_) then
        n = n + 1
        call two_el_idx(i, two_el_table%index(:, n))
        two_el_table%values(n) = TUVX(i)
      end if
    end do
    two_el_table%length = n
    two_el_table%cutoff = cutoff_

! ========== For testing purposes FROM HERE =============
    l_twoel_test = min(max_test, length(two_el_table))
    call Add_Info('TwoEl Integral element Input', &
      TUVX(:l_twoel_test), l_twoel_test, 8)
! ========== For testing purposes TO HERE ===============
  end subroutine fill_2ElInt

  pure integer function TwoElIntTable_length(table)
    type(TwoElIntTable), intent(in) :: table
    TwoElIntTable_length = table%length
  end function TwoElIntTable_length

  subroutine TwoElIntTable_print(table)
    type(TwoElIntTable), intent(in) :: table
    integer :: i, j
    do j = 1, length(table)
      write(6, '(E15.7, I7, I7, I7, I7)') &
        table%values(j), (table%index(i, j), i=1, 4)
    end do
  end subroutine TwoElIntTable_print

  subroutine TwoElIntTable_unused(table)
    type(TwoElIntTable), intent(in) :: table
    integer :: n
    if (.false.) n = length(table)
  end subroutine TwoElIntTable_unused
end module fcidump_tables
