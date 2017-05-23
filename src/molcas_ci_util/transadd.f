************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine transadd(n,a,lda)
      ! Add the transpose to a square matrix, A := A + A^T

      implicit none

      ! dummy arguments
      integer, intent(in) :: n, lda
      real*8, intent(inout) :: a(lda,*)

      ! small buffer block
      integer, parameter :: nblksz = 8

      ! number of blocks, remaining elements
      integer :: nblk, nrem
      integer :: iblk, jblk, ista, jsta, i, j

      nblk = n / nblksz
      nrem = n - nblk * nblksz

      ! loop over diagonal blocks
      do iblk = 1,nblk
        ista = nblksz * (iblk-1) + 1
        ! add transposed upper to lower triangular part
        do i = 1,nblksz
          do j = 1,i
            a(ista+i-1,ista+j-1) = a(ista+i-1,ista+j-1) +
     &       a(ista+j-1,ista+i-1)
          end do
        end do
        ! copy lower to upper part
        do i = 1,nblksz
          do j = 1,i-1
            a(ista+j-1,ista+i-1) = a(ista+i-1,ista+j-1)
          end do
        end do
      end do
      ! the tail block
      ista = nblksz * nblk + 1
      ! add transposed upper to lower triangular part
      do i = 1,nrem
        do j = 1,i
          a(ista+i-1,ista+j-1) = a(ista+i-1,ista+j-1) +
     &     a(ista+j-1,ista+i-1)
        end do
      end do
      ! copy lower to upper part
      do i = 1,nrem
        do j = 1,i-1
          a(ista+j-1,ista+i-1) = a(ista+i-1,ista+j-1)
        end do
      end do

      ! loop over lower triangular blocks
      do iblk = 1,nblk
        ista = nblksz * (iblk-1) + 1
        do jblk = 1,iblk-1
          jsta = nblksz * (jblk-1) + 1
          ! add transposed upper to lower triangular block
          do j = 1,nblksz
            do i = 1,nblksz
              a(ista+i-1,jsta+j-1) = a(ista+i-1,jsta+j-1) +
     &         a(jsta+j-1,ista+i-1)
            end do
          end do
          ! copy transpose of lower to upper triangular block
          do j = 1,nblksz
            do i = 1,nblksz
              a(jsta+i-1,ista+j-1) = a(ista+j-1,jsta+i-1)
            end do
          end do
        end do
      end do
      ! bottom blocks
      ista = nblksz * nblk + 1
      do jblk = 1,iblk-1
        jsta = nblksz * (jblk-1) + 1
        ! add transposed upper to lower triangular block
        do j = 1,nblksz
          do i = 1,nrem
            a(ista+i-1,jsta+j-1) = a(ista+i-1,jsta+j-1) +
     &       a(jsta+j-1,ista+i-1)
          end do
        end do
        ! copy transpose of lower to upper triangular block
        do j = 1,nrem
          do i = 1,nblksz
            a(jsta+i-1,ista+j-1) = a(ista+j-1,jsta+i-1)
          end do
        end do
      end do

      end subroutine
