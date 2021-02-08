************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICEnSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Oskar Weser                                      *
************************************************************************
#include "intent.h"
      module blockdiagonal_matrices
        use stdalloc, only : mma_allocate, mma_deallocate
        implicit none
        private
        public :: t_blockdiagonal, new, delete, from_raw,
     &    from_symm_raw, to_raw, blocksizes
        save

        type :: t_blockdiagonal
          real*8, allocatable :: block(:, :)
        end type

        interface new
          module procedure block_new
        end interface

        interface delete
          module procedure block_delete
        end interface

      contains

        subroutine block_new(blocks, blocksizes)
          type(t_blockdiagonal), allocatable, intent(out) :: blocks(:)
          integer, intent(in) :: blocksizes(:)
          integer :: i, L, err

          if (allocated(blocks)) deallocate(blocks)
          allocate(blocks(size(blocksizes)), stat=err)
          if (err /= 0) call abort_('Allocation failed in '//
     &        'blockdiagonal_matrices::new')
          do i = 1, size(blocks)
            L = blocksizes(i)
            call mma_allocate(blocks(i)%block, L, L, label='Block')
          end do
        end subroutine

        subroutine block_delete(blocks)
          type(t_blockdiagonal), allocatable :: blocks(:)
          integer :: i

          do i = 1, size(blocks)
            call mma_deallocate(blocks(i)%block)
          end do
          deallocate(blocks)
        end subroutine

        subroutine from_raw(S_buffer, S)
          implicit none
          real*8, intent(in) :: S_buffer(:)
          type(t_blockdiagonal), intent(_OUT_) :: S(:)
          integer :: i_block, offset, col, block_size, idx_block

          idx_block = 1
          do i_block = 1, size(S)
            block_size = size(S(i_block)%block, 1)
            do col = 1, block_size
              offset = idx_block + (col - 1) * block_size
              S(i_block)%block(:, col) =
     &          S_buffer(offset : offset - 1 + block_size)
            end do
            idx_block = idx_block + (block_size**2)
          end do
        end subroutine

        subroutine to_raw(S, S_buffer)
          implicit none
          type(t_blockdiagonal), intent(in) :: S(:)
          real*8, intent(out) :: S_buffer(:)
          integer :: i_block, offset, col, block_size, idx_block

          idx_block = 1
          do i_block = 1, size(S)
            block_size = size(S(i_block)%block, 1)
            do col = 1, block_size
              offset = idx_block + (col - 1) * block_size
              S_buffer(offset : offset - 1 + block_size) =
     &            S(i_block)%block(:, col)
            end do
            idx_block = idx_block + (block_size**2)
          end do
        end subroutine

        subroutine from_symm_raw(S_buffer, S)
          implicit none
          real*8, intent(in) :: S_buffer(:)
          type(t_blockdiagonal), intent(_OUT_) :: S(:)
          integer :: i_block, block_size, idx_block

          idx_block = 1
          do i_block = 1, size(S)
            block_size = size(S(i_block)%block, 1)
            if (block_size > 0) then
              call square(S_buffer(idx_block:), S(i_block)%block,1,
     &                    block_size, block_size)
            end if
            idx_block = idx_block + (block_size**2 + block_size) / 2
          end do
        end subroutine

        pure function blocksizes(A) result(res)
          implicit none
          type(t_blockdiagonal), intent(in) :: A(:)
          integer :: res(size(A))

          integer :: i
          res(:) = [(size(A(i)%block, 1), i = 1, size(A))]
        end function

        subroutine abort_(message)
          character(len=*), intent(in) :: message
          call WarningMessage(2, message)
          call Abend()
        end subroutine
      end module blockdiagonal_matrices
