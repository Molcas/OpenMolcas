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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************

#include "macros.fh"

module rhodyn_utils
! module contains some auxiliary routines

use linalg_mod, only: abort_
use definitions, only: wp, iwp

implicit none

public :: mult, removeColumn, removeLineAndColumn, dashes, transform

interface mult
  module procedure mult_2D, multZ_2D
end interface
interface removeLineAndColumn
  module procedure removeLineAndColumnZ, removeLineAndColumnR
end interface
interface removeColumn
  module procedure removeColumnZ, removeColumnR
end interface
interface transform
  module procedure transformZ, transformR
end interface

contains

subroutine dashes(length)
  ! print the line of dashes of length 'length'
  ! to the standard output fileunit 6

  integer(kind=iwp) :: i, l, fileunit
  integer(kind=iwp), intent(in), optional :: length

  if (present(length)) then
    l = length
  else
    l = 72
  end if
  fileunit = 6
  do i=1,l
    write(fileunit,'(A)',advance='no') '-'
  end do
  write(fileunit,*)

end subroutine dashes

subroutine mult_2D(a,b,c,transpA,transpB)

  real(kind=wp), intent(in) :: a(:,:), b(:,:)
  real(kind=wp), intent(out) :: c(:,:)
  logical, intent(in), optional :: transpA, transpB
  logical :: transpA_, transpB_
  integer(kind=iwp) :: m, n, k, k1, k2

  if (present(transpA)) then
    transpA_ = transpA
  else
    transpA_ = .false.
  end if
  if (present(transpB)) then
    transpB_ = transpB
  else
    transpB_ = .false.
  end if
  m = size(a,merge(2,1,transpA_))
  ASSERT(m == size(c,1))
  n = size(b,merge(1,2,transpB_))
  ASSERT(n == size(c,2))
  k1 = size(a,merge(1,2,transpA_))
  k2 = size(b,merge(2,1,transpB_))
  ASSERT(k1 == k2)
  k = k1
  call dgemm_(merge('T','N',transpA_),merge('T','N',transpB_),m,n,k,1.d0,a,size(a,1),b,size(b,1),0.d0,c,size(c,1))

end subroutine mult_2D

subroutine multZ_2D(a,b,c,transpA,transpB)

  complex(kind=wp), intent(in) :: a(:,:), b(:,:)
  complex(kind=wp), intent(out) :: c(:,:)
  logical, intent(in), optional :: transpA, transpB
  logical :: transpA_, transpB_
  integer(kind=iwp) :: m, n, k, k1, k2

  if (present(transpA)) then
    transpA_ = transpA
  else
    transpA_ = .false.
  end if
  if (present(transpB)) then
    transpB_ = transpB
  else
    transpB_ = .false.
  end if
  m = size(a,merge(2,1,transpA_))
  ASSERT(m == size(c,1))
  n = size(b,merge(1,2,transpB_))
  ASSERT(n == size(c,2))
  k1 = size(a,merge(1,2,transpA_))
  k2 = size(b,merge(2,1,transpB_))
  ASSERT(k1 == k2)
  k = k1
  call zgemm_(merge('C','N',transpA_),merge('C','N',transpB_),m,n,k,(1.0d0,0.0d0),a,size(a,1),b,size(b,1),(0.0d0,0.0d0),c,size(c,1))

end subroutine multZ_2D

subroutine removeLineAndColumnZ(a,remLCarray)

  use stdalloc, only: mma_allocate, mma_deallocate

  complex(kind=wp), dimension(:,:), allocatable, intent(inout) :: a
  integer(kind=iwp), dimension(:), intent(in) :: remLCarray
  ! temp variables
  complex(kind=wp), dimension(:), allocatable :: b
  logical, dimension(:,:), allocatable :: mask
  integer(kind=iwp) :: i, j, l, m, n, k

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remLCarray)
  call mma_allocate(mask,n,m,label='mask')
  ! mask creation
  l = 0
  mask = .false.
  do i=1,n
    do j=1,m
      if (any(remLCarray == i) .and. any(remLCarray == j)) then
        mask(i,j) = .true.
        l = l+1
      end if
    end do
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,k,k,label='a')
  ! copy back
  a(:,:) = reshape(b,(/k,k/))
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeLineAndColumnZ

subroutine removeLineAndColumnR(a,remLCarray)

  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: a
  integer(kind=iwp), dimension(:), intent(in) :: remLCarray
  ! temp variables
  real(kind=wp), dimension(:), allocatable :: b
  logical, dimension(:,:), allocatable :: mask
  integer(kind=iwp) :: i, j, l, m, n, k

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remLCarray)
  call mma_allocate(mask,n,m,label='a')
  ! mask creation
  l = 0
  mask = .false.
  do i=1,n
    do j=1,m
      if (any(remLCarray == i) .and. any(remLCarray == j)) then
        mask(i,j) = .true.
        l = l+1
      end if
    end do
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,k,k,label='a')
  ! copy back
  a(:,:) = reshape(b,(/k,k/))
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeLineAndColumnR

subroutine removeColumnZ(a,remCarray)

  use stdalloc, only: mma_allocate, mma_deallocate

  complex(kind=wp), dimension(:,:), allocatable, intent(inout) :: a
  integer(kind=iwp), dimension(:), intent(in) :: remCarray
  ! temp variables
  complex(kind=wp), dimension(:), allocatable :: b
  logical, dimension(:,:), allocatable :: mask
  integer(kind=iwp) :: j, l, m, n, k

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remCarray)
  call mma_allocate(mask,n,m,label='mask')
  ! mask creation
  l = 0
  mask = .false.
  do j=1,m
    if (any(remCarray == j)) then
      mask(:,j) = .true.
      l = l+n
    end if
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,n,k,label='a')
  ! copy back
  a(:,:) = reshape(b,(/n,k/))
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeColumnZ

subroutine removeColumnR(a,remCarray)

  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: a
  integer(kind=iwp), dimension(:), intent(in) :: remCarray
  ! temp variables
  real(kind=wp), dimension(:), allocatable :: b
  logical, dimension(:,:), allocatable :: mask
  integer(kind=iwp) :: j, l, m, n, k

  ! sizes
  n = size(a,1)
  m = size(a,2)
  k = size(remCarray)
  call mma_allocate(mask,n,m,label='mask')
  ! mask creation
  l = 0
  mask = .false.
  do j=1,m
    if (any(remCarray == j)) then
      mask(:,j) = .true.
      l = l+n
    end if
  end do
  call mma_allocate(b,l,label='b')
  ! copy remaining elements into temp b
  b(:) = pack(a,mask)
  call mma_deallocate(a)
  call mma_allocate(a,n,k,label='a')
  ! copy back
  a(:,:) = reshape(b,(/n,k/))
  call mma_deallocate(b)
  call mma_deallocate(mask)

end subroutine removeColumnR

subroutine transformZ(a,u,b,order)
  ! perform orthogonal transformation with transformation matrix
  ! order = True (default): direct transform : b = u^C * a * u
  ! order = False         : inverse transform: b = u   * a * u^C

  use stdalloc, only: mma_allocate, mma_deallocate

  complex(kind=wp), dimension(:,:), intent(in) :: a, u
  complex(kind=wp), dimension(:,:), intent(out) :: b
  complex(kind=wp), dimension(:,:), allocatable :: temp
  logical, optional, intent(in) :: order
  logical :: order_
  integer(kind=iwp) :: m, n

  if (present(order)) then
    order_ = order
  else
    order_ = .true.
  end if
  m = size(u,1)
  n = size(u,2)
  if (order_) then
    ASSERT(m == size(a,1))
    call mma_allocate(temp,n,m,label='temp')
    call multZ_2D(u,a,temp,.true.,.false.)
    call multZ_2D(temp,u,b,.false.,.false.)
  else
    ASSERT(n == size(a,1))
    call mma_allocate(temp,m,n,label='temp')
    call multZ_2D(u,a,temp,.false.,.false.)
    call multZ_2D(temp,u,b,.false.,.true.)
  end if
  call mma_deallocate(temp)

end subroutine transformZ

subroutine transformR(a,u,b,order)
  ! perform orthogonal transformation with transformation matrix
  ! order = True (default): direct transform : b = u^T * a * u
  ! order = False       : inverse transform: b = u   * a * u^T
  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp), dimension(:,:), intent(in) :: a, u
  real(kind=wp), dimension(:,:), intent(out) :: b
  real(kind=wp), dimension(:,:), allocatable :: temp
  logical, optional, intent(in) :: order
  logical :: order_
  integer(kind=iwp) :: m, n

  if (present(order)) then
    order_ = order
  else
    order_ = .true.
  end if
  m = size(u,1)
  n = size(u,2)
  if (order_) then
    ASSERT(m == size(a,1))
    call mma_allocate(temp,n,m,label='temp')
    call mult_2D(u,a,temp,.true.,.false.)
    call mult_2D(temp,u,b,.false.,.false.)
  else
    ASSERT(n == size(a,1))
    call mma_allocate(temp,m,n,label='temp')
    call mult_2D(u,a,temp,.false.,.false.)
    call mult_2D(temp,u,b,.false.,.true.)
  end if
  call mma_deallocate(temp)

end subroutine transformR

subroutine sortci(N1,A,WR,C,print_level)

  integer(kind=iwp) :: N1, INFO, LWORK
  integer(kind=iwp), intent(in) :: print_level
  real(kind=wp), dimension(N1,N1) :: A, C, diag, B
  real(kind=wp), dimension(N1) :: WR
  real(kind=wp), dimension(2*N1) :: WORK
  integer(kind=iwp) :: k, l

  B = A
  LWORK = 2*N1
  call dsyev_('V','U',N1,A,N1,WR,WORK,LWORK,INFO)
  if (INFO /= 0) then
    write(6,*) 'ERROR in sortci'
    call Abend()
  end if
  call dsyev_('V','U',N1,A,N1,WR,WORK,LWORK,INFO)
  C = A
  if (print_level > 3) then
    call transform(B,C,diag)
    call dashes()
    write(6,*) 'Printout the diagonalized matrix:'
    call dashes()
    do k=1,10
      write(6,*) (diag(k,l),l=1,10)
    end do
  end if

end subroutine sortci

end module rhodyn_utils
