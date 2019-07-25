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
  use stdalloc, only : mma_allocate, mma_deallocate
  use fcidump_tables, only :  FockTable, TwoElIntTable, OrbitalTable,&
    mma_allocate, mma_deallocate, length
  implicit none
  private
  public :: reorder, get_P_GAS, get_P_inp, ReOrFlag, ReOrInp, cleanup
  save
! n==0: Don't reorder.
! n>=2: User defined permutation with n non-fixed point elements.
! n==-1: Use GAS sorting scheme.
  integer :: ReOrFlag = 0
  integer, allocatable :: ReOrInp(:)

  interface reorder
    module procedure FockTable_reorder, TwoElIntTable_reorder, &
      OrbitalTable_reorder, ALL_reorder
  end interface

contains

  subroutine OrbitalTable_reorder(orbitals, P)
    implicit none
    type(OrbitalTable), intent(inout) :: orbitals
    integer, intent(in) :: P(:)
    integer :: i
    do i = 1, length(orbitals)
      orbitals%index(i) = P(orbitals%index(i))
    end do
  end subroutine

  subroutine FockTable_reorder(fock, P)
    implicit none
    type(FockTable), intent(inout) :: fock
    integer, intent(in) :: P(:)
    integer :: i, j
    do j = 1, length(fock)
      do i = 1, 2
        fock%index(i, j) = P(fock%index(i, j))
      end do
    end do
  end subroutine

  subroutine TwoElIntTable_reorder(two_el_table, P)
    implicit none
    type(TwoElIntTable), intent(inout) :: two_el_table
    integer, intent(in) :: P(:)
    integer :: i, j
    do j = 1, length(two_el_table)
      do i = 1, 4
        two_el_table%index(i, j) = P(two_el_table%index(i, j))
      end do
    end do
  end subroutine

  function get_P_GAS(ngssh) result(P)
    use sorting, only : argsort
    implicit none
    integer, intent(in) :: ngssh(:, :)
    integer :: P(sum(ngssh)), X(sum(ngssh))
    integer :: iGAS, iSym, iOrb, bounds(2)
    bounds = shape(ngssh)
    iOrb = 1
    do iSym = 1, bounds(2)
      do iGAS = 1, bounds(1)
        X(iOrb : iOrb + ngssh(iGAS, iSym)) = iGAS
        iOrb = iOrb + ngssh(iGAS, iSym) + 1
      end do
    end do
    P = argsort(X)
  end function

  function get_P_inp(ReOrInp) result(P)
    use sorting, only : sort
    use general_data, only : nAsh
    implicit none
    integer, intent(in) :: ReOrInp(:)
    integer :: P(sum(nAsh)), change_idx(size(ReOrInp)), i
    P = [(i, i = 1, size(P))]
    change_idx = ReOrInp
    call sort(change_idx)
    P(change_idx) = ReOrInp
  end function

  subroutine ALL_reorder(orbitals, fock, two_el_table, P)
    implicit none
    type(OrbitalTable), intent(inout) :: orbitals
    type(FockTable), intent(inout) :: fock
    type(TwoElIntTable), intent(inout) :: two_el_table
    integer, intent(in) :: P(:)
    call reorder(orbitals, P)
    call reorder(fock, P)
    call reorder(two_el_table, P)
  end subroutine

  subroutine cleanup()
    implicit none
    call mma_deallocate(ReOrInp)
  end subroutine
end module fcidump_reorder
