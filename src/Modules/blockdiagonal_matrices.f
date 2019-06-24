      module blockdiagonal_matrices
        use stdalloc, only : mma_allocate, mma_deallocate
        implicit none
        private
        public :: t_blockdiagonal, new, delete, fill_from_buffer,
     &    fill_from_symm_buffer, fill_to_buffer
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
          type(t_blockdiagonal), intent(out) :: blocks(:)
          integer, intent(in) :: blocksizes(:)
          integer :: i, L

          do i = 1, size(blocks)
            L = blocksizes(i)
            call mma_allocate(blocks(i)%block, L, L)
          end do
        end subroutine

        subroutine block_delete(blocks)
          type(t_blockdiagonal) :: blocks(:)
          integer :: i

          do i = 1, size(blocks)
            call mma_deallocate(blocks(i)%block)
          end do
        end subroutine

        subroutine fill_from_buffer(S_buffer, S)
          implicit none
          real*8, intent(in) :: S_buffer(:)
          type(t_blockdiagonal), intent(inout) :: S(:)
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

        subroutine fill_to_buffer(S, S_buffer)
          implicit none
          real*8, intent(out) :: S_buffer(:)
          type(t_blockdiagonal), intent(in) :: S(:)
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

        subroutine fill_from_symm_buffer(S_buffer, S)
          implicit none
          real*8, intent(in) :: S_buffer(:)
          type(t_blockdiagonal), intent(inout) :: S(:)
          integer :: i_block, block_size, idx_block

          idx_block = 1
          do i_block = 1, size(S)
            block_size = size(S(i_block)%block, 1)
            if (block_size > 0) then
              call square(S_buffer(idx_block), S(i_block)%block,1,
     &                    block_size, block_size)
            end if
            idx_block = idx_block + (block_size**2 + block_size) / 2
          end do
        end subroutine

      end module blockdiagonal_matrices
