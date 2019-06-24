      module blockdiagonal_matrices
        use stdalloc, only : mma_allocate, mma_deallocate
        implicit none
        private
        public :: t_blockdiagonal, new, delete
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



      end module blockdiagonal_matrices
