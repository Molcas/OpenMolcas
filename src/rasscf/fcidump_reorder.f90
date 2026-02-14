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

module fcidump_reorder

use fcidump_tables, only: FockTable, length, OrbitalTable, TwoElIntTable
use sorting_funcs, only: leq_i
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
private

! n==0: Don't reorder.
! n>=2: User defined permutation with n non-fixed point elements.
! n==-1: Use GAS sorting scheme.
integer(kind=iwp) :: ReOrFlag = 0
integer(kind=iwp), allocatable :: ReOrInp(:)

public :: cleanup, get_P_GAS, get_P_inp, reorder, ReOrFlag, ReOrInp

interface reorder
  module procedure :: FockTable_reorder, TwoElIntTable_reorder, OrbitalTable_reorder, ALL_reorder
end interface reorder

contains

subroutine OrbitalTable_reorder(orbitals,P)

  type(OrbitalTable), intent(inout) :: orbitals
  integer(kind=iwp), intent(in) :: P(:)
  integer(kind=iwp) :: i

  do i=1,length(orbitals)
    orbitals%idx(i) = P(orbitals%idx(i))
  end do

end subroutine OrbitalTable_reorder

subroutine FockTable_reorder(fock,P)

  type(FockTable), intent(inout) :: fock
  integer(kind=iwp), intent(in) :: P(:)
  integer(kind=iwp) :: i, j

  do j=1,length(fock)
    do i=1,2
      fock%idx(i,j) = P(fock%idx(i,j))
    end do
  end do

end subroutine FockTable_reorder

subroutine TwoElIntTable_reorder(two_el_table,P)

  type(TwoElIntTable), intent(inout) :: two_el_table
  integer(kind=iwp), intent(in) :: P(:)
  integer(kind=iwp) :: i, j

  do j=1,length(two_el_table)
    do i=1,4
      two_el_table%idx(i,j) = P(two_el_table%idx(i,j))
    end do
  end do
  do j=1,length(two_el_table)
    do i=1,4
      two_el_table%idx(i,j) = P(two_el_table%idx(i,j))
    end do
  end do

end subroutine TwoElIntTable_reorder

function get_P_GAS(ngssh) result(P)

  use sorting, only: argsort
  use general_data, only: nSym
  use gas_data, only: nGAS

  integer(kind=iwp), intent(in) :: ngssh(:,:)
  integer(kind=iwp) :: P(sum(ngssh))
  integer(kind=iwp) :: i, iGAS, iSym
  integer(kind=iwp), allocatable :: X(:)

  call mma_allocate(X,sum(ngssh))
  X(:) = [(((iGAS,i=1,ngssh(iGAS,iSym)),iGAS=1,nGAS),iSym=1,nSym)]
  P(:) = argsort(X,leq_i)
  call mma_deallocate(X)

end function get_P_GAS

function get_P_inp(ReOrInp) result(P)

  use sorting, only: sort
  use general_data, only: nAsh

  integer(kind=iwp) :: P(sum(nAsh))
  integer(kind=iwp), intent(in) :: ReOrInp(:)
  integer(kind=iwp) :: i
  integer(kind=iwp), allocatable :: change_idx(:)

  call mma_allocate(change_idx,size(ReOrInp))
  P(:) = [(i,i=1,size(P))]
  change_idx(:) = ReOrInp
  call sort(change_idx,leq_i)
  P(change_idx) = ReOrInp
  call mma_deallocate(change_idx)

end function get_P_inp

subroutine ALL_reorder(orbitals,fock,two_el_table,orbsym,P)

  type(OrbitalTable), intent(inout) :: orbitals
  type(FockTable), intent(inout) :: fock
  type(TwoElIntTable), intent(inout) :: two_el_table
  integer(kind=iwp), intent(inout) :: orbsym(:)
  integer(kind=iwp), intent(in) :: P(:)

  call reorder(orbitals,P)
  call reorder(fock,P)
  call reorder(two_el_table,P)
  orbsym(:) = orbsym(P)

end subroutine ALL_reorder

subroutine cleanup()

  use fcidump_tables, only: mma_deallocate

  call mma_deallocate(ReOrInp,safe='*')

end subroutine cleanup

end module fcidump_reorder
