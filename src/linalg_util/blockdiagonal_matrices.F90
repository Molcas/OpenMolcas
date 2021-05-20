!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICEnSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Oskar Weser                                      *
!***********************************************************************

module blockdiagonal_matrices

#include "intent.fh"

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

public :: t_blockdiagonal, new, delete, from_raw, from_symm_raw, to_raw, blocksizes

type :: t_blockdiagonal
  real(kind=wp), allocatable :: blck(:,:)
end type

interface new
  module procedure :: block_new
end interface

interface delete
  module procedure :: block_delete
end interface

! Private extensions to mma interfaces
interface cptr2loff
  module procedure :: block_cptr2loff
end interface
interface mma_allocate
  module procedure :: block_mma_allo_1D, block_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: block_mma_free_1D
end interface

contains

subroutine block_new(blocks,blocksizes)
  type(t_blockdiagonal), allocatable, intent(_OUT_) :: blocks(:)
  integer(kind=iwp), intent(in) :: blocksizes(:)
  integer(kind=iwp) :: i, L
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A)
      import :: wp
      real(kind=wp), allocatable :: A(:,:)
    end subroutine c_null_alloc
  end interface
# endif

  if (allocated(blocks)) call block_delete(blocks)
  call mma_allocate(blocks,size(blocksizes),label='blocks')
  do i=1,size(blocks)
#   ifdef _GARBLE_
    ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
    call c_null_alloc(blocks(i)%blck)
#   endif
    L = blocksizes(i)
    call mma_allocate(blocks(i)%blck,L,L,label='Block')
  end do

# include "macros.fh"
  unused_proc(mma_allocate(blocks,[0,0]))

end subroutine block_new

subroutine block_delete(blocks)
  type(t_blockdiagonal), allocatable, intent(_OUT_) :: blocks(:)
  integer(kind=iwp) :: i

  do i=1,size(blocks)
    call mma_deallocate(blocks(i)%blck)
  end do
  call mma_deallocate(blocks)
end subroutine block_delete

subroutine from_raw(S_buffer,S)
  real(kind=wp), intent(in) :: S_buffer(:)
  type(t_blockdiagonal), intent(_OUT_) :: S(:) ! not allocatable, but has allocatable components
  integer(kind=iwp) :: i_block, offset, col, block_size, idx_block

  idx_block = 1
  do i_block=1,size(S)
    block_size = size(S(i_block)%blck,1)
    do col=1,block_size
      offset = idx_block+(col-1)*block_size
      S(i_block)%blck(:,col) = S_buffer(offset:offset-1+block_size)
    end do
    idx_block = idx_block+(block_size**2)
  end do
end subroutine from_raw

subroutine to_raw(S,S_buffer)
  type(t_blockdiagonal), intent(in) :: S(:)
  real(kind=wp), intent(out) :: S_buffer(:)
  integer(kind=iwp) :: i_block, offset, col, block_size, idx_block

  idx_block = 1
  do i_block=1,size(S)
    block_size = size(S(i_block)%blck,1)
    do col=1,block_size
      offset = idx_block+(col-1)*block_size
      S_buffer(offset:offset-1+block_size) = S(i_block)%blck(:,col)
    end do
    idx_block = idx_block+(block_size**2)
  end do
end subroutine to_raw

subroutine from_symm_raw(S_buffer,S)
  real(kind=wp), intent(in) :: S_buffer(:)
  type(t_blockdiagonal), intent(_OUT_) :: S(:) ! not allocatable, but has allocatable components
  integer(kind=iwp) :: i_block, block_size, idx_block

  idx_block = 1
  do i_block=1,size(S)
    block_size = size(S(i_block)%blck,1)
    if (block_size > 0) then
      call square(S_buffer(idx_block:),S(i_block)%blck,1,block_size,block_size)
    end if
    idx_block = idx_block+(block_size**2+block_size)/2
  end do
end subroutine from_symm_raw

pure function blocksizes(A) result(res)
  type(t_blockdiagonal), intent(in) :: A(:)
  integer(kind=iwp) :: res(size(A))
  integer(kind=iwp) :: i

  res(:) = [(size(A(i)%blck,1),i=1,size(A))]
end function blocksizes

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define block_cptr2loff, block_mma_allo_1D, block_mma_allo_1D_lim, block_mma_free_1D
#define _TYPE_ type(t_blockdiagonal)
#  define _FUNC_NAME_ block_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ block_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'blk_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module blockdiagonal_matrices
