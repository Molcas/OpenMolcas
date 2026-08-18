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

use index_symmetry, only: one_el_idx, two_el_idx
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
private

type :: FockTable
  real(kind=wp), allocatable :: values(:)    ! <i| F | j>
  integer(kind=iwp), allocatable :: idx(:,:) ! i, j
  real(kind=wp) :: cutoff
  integer(kind=iwp) :: length
end type FockTable

type :: TwoElIntTable
  real(kind=wp), allocatable :: values(:)    ! <ij| 1/r_{12} |kl>
  integer(kind=iwp), allocatable :: idx(:,:) ! i, j, k, l
  real(kind=wp) :: cutoff
  integer(kind=iwp) :: length
end type TwoElIntTable

type :: OrbitalTable
  real(kind=wp), allocatable :: values(:)  ! <i| F |i>
  integer(kind=iwp), allocatable :: idx(:) ! i
end type OrbitalTable

real(kind=wp), parameter :: cutoff_default = 1.0e-11_wp

interface mma_allocate
  module procedure :: FockTable_allocate, TwoElIntTable_allocate, OrbitalTable_allocate
end interface mma_allocate

interface mma_deallocate
  module procedure :: FockTable_deallocate, TwoElIntTable_deallocate, OrbitalTable_deallocate
end interface mma_deallocate

interface length
  module procedure :: FockTable_length, TwoElIntTable_length, OrbitalTable_length
end interface length

interface print_table
  module procedure :: FockTable_print, TwoElIntTable_print, OrbitalTable_print
end interface print_table

public :: cutoff_default, fill_2ElInt, fill_fock, fill_orbitals, FockTable, length, mma_allocate, mma_deallocate, OrbitalTable, &
          print_table, TwoElIntTable

contains

subroutine OrbitalTable_allocate(orbital_table,n)

  type(OrbitalTable), intent(inout) :: orbital_table
  integer(kind=iwp), intent(in) :: n

  call mma_allocate(orbital_table%values,n)
  call mma_allocate(orbital_table%idx,n)

end subroutine OrbitalTable_allocate

subroutine OrbitalTable_deallocate(orbital_table)

  type(OrbitalTable), intent(inout) :: orbital_table

  call mma_deallocate(orbital_table%values)
  call mma_deallocate(orbital_table%idx)

end subroutine OrbitalTable_deallocate

!>  @brief
!>    Fill in orbital energies
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The orbitals table gets filled with the orbital energies from orbital_energies.
!>
!>  @param[in,out] table
!>  @param[in] orbital_energies
subroutine fill_orbitals(table,orbital_energies)

  use general_data, only: nAsh, nBas, nFro, nIsh, nSym

  type(OrbitalTable), intent(inout) :: table
  real(kind=wp), intent(in) :: orbital_energies(:)
  integer(kind=iwp) :: i, iOff, iSym, n

  iOff = 0
  n = 1
  do iSym=1,nSym
    if (nAsh(iSym) > 0) then
      do i=1,nAsh(isym)
        table%idx(n) = n
        table%values(n) = orbital_energies(ioff+nFro(iSym)+nIsh(iSym)+i)
        n = n+1
      end do
    end if
    iOff = iOff+nBas(iSym)
  end do

end subroutine fill_orbitals

pure function OrbitalTable_length(table)

  integer(kind=iwp) :: OrbitalTable_length
  type(OrbitalTable), intent(in) :: table

  OrbitalTable_length = size(table%values)

end function OrbitalTable_length

subroutine OrbitalTable_print(table)

  type(OrbitalTable), intent(in) :: table
  integer(kind=iwp) :: i

  do i=1,length(table)
    write(u6,'(ES15.7, I7)') table%values(i),table%idx(i)
  end do

end subroutine OrbitalTable_print

subroutine FockTable_allocate(fock_table,n)

  type(FockTable), intent(inout) :: fock_table
  integer(kind=iwp), intent(in) :: n

  call mma_allocate(fock_table%values,n)
  call mma_allocate(fock_table%idx,2,n)

end subroutine FockTable_allocate

subroutine FockTable_deallocate(fock_table)

  type(FockTable), intent(inout) :: fock_table

  call mma_deallocate(fock_table%values)
  call mma_deallocate(fock_table%idx)

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
!>  @param[in] Fock
!>    \f[\sum_{\sigma\rho} {In}^D_{\sigma\rho} (g_{\mu\nu\sigma\rho})  \f]
!>  @param[in] cutoff Optional parameter that is set by default to
!>    fciqmc_tables::cutoff_default.
subroutine fill_fock(fock_table,Fock,cutoff)

  type(FockTable), intent(inout) :: fock_table
  real(kind=wp), intent(in) :: Fock(:)
  real(kind=wp), optional, intent(in) :: cutoff
  integer(kind=iwp) :: i, n
  real(kind=wp) :: cutoff_

  cutoff_ = cutoff_default
  if (present(cutoff)) cutoff_ = cutoff

  n = 0
  do i=1,size(Fock)
    if (abs(Fock(i)) >= cutoff_) then
      n = n+1
      call one_el_idx(i,fock_table%idx(:,n))
      fock_table%values(n) = Fock(i)
    end if
  end do
  fock_table%length = n
  fock_table%cutoff = cutoff_

end subroutine fill_fock

pure function FockTable_length(table)

  integer(kind=iwp) :: FockTable_length
  type(FockTable), intent(in) :: table

  FockTable_length = table%length

end function FockTable_length

subroutine FockTable_print(table)

  type(FockTable), intent(in) :: table
  integer(kind=iwp) :: i, j

  do j=1,length(table)
    write(u6,'(ES15.7, I7, I7)') table%values(j),(table%idx(i,j),i=1,2)
  end do

end subroutine FockTable_print

subroutine TwoElIntTable_allocate(table,n)

  type(TwoElIntTable), intent(inout) :: table
  integer(kind=iwp), intent(in) :: n

  call mma_allocate(table%values,n)
  call mma_allocate(table%idx,4,n)

end subroutine TwoElIntTable_allocate

subroutine TwoElIntTable_deallocate(table)

  type(TwoElIntTable), intent(inout) :: table

  call mma_deallocate(table%values)
  call mma_deallocate(table%idx)

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
subroutine fill_2ElInt(two_el_table,TUVX,cutoff)

  type(TwoElIntTable), intent(inout) :: two_el_table
  real(kind=wp), intent(in) :: TUVX(:)
  real(kind=wp), optional :: cutoff
  integer(kind=iwp) :: i, l_twoel_test, n
  real(kind=wp) :: cutoff_
  integer(kind=iwp), parameter :: max_test = 20

  cutoff_ = cutoff_default
  if (present(cutoff)) cutoff_ = cutoff

  n = 0
  do i=1,size(TUVX)
    if (abs(TUVX(i)) >= cutoff_) then
      n = n+1
      call two_el_idx(i,two_el_table%idx(:,n))
      two_el_table%values(n) = TUVX(i)
    end if
  end do
  two_el_table%length = n
  two_el_table%cutoff = cutoff_

  ! ========== For testing purposes FROM HERE =============
  l_twoel_test = min(max_test,length(two_el_table))
  call Add_Info('TwoEl Integral element Input',TUVX(:l_twoel_test),l_twoel_test,8)
  ! ========== For testing purposes TO HERE ===============

end subroutine fill_2ElInt

pure function TwoElIntTable_length(table)

  integer(kind=iwp) :: TwoElIntTable_length
  type(TwoElIntTable), intent(in) :: table

  TwoElIntTable_length = table%length

end function TwoElIntTable_length

subroutine TwoElIntTable_print(table)

  type(TwoElIntTable), intent(in) :: table
  integer(kind=iwp) :: i, j

  do j=1,length(table)
    write(u6,'(ES15.7, I7, I7, I7, I7)') table%values(j),(table%idx(i,j),i=1,4)
  end do

end subroutine TwoElIntTable_print

end module fcidump_tables
