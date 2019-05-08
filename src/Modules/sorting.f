************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Oskar Weser                                      *
************************************************************************
      module sorting
      implicit none
      private
      public :: sort, argsort

! When to switch to naive sort, which scales with n^2,
! but requires less overhead.
      integer :: naive_sort_trsh = 30

      interface swap
        module procedure I_swap, R_swap
      end interface

      interface naive_sort
        module procedure I_naive_sort, R_naive_sort
      end interface

      interface quicksort
        module procedure R_quicksort, I_quicksort
      end interface

!>  @brief
!>    Sort the array a inplace.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is a generic procedure that accepts integer and
!>  real inputs.
!>  Uses the quicksort algorithm.
!>
!>  \code{.unparsed}
!>  A = [2, 1, 3]
!>  call sort(A)
!>  A .eq. [1, 2, 3]
!>  \endcode
!>
!>  @param[in] a 1D-Array to be sorted
      interface sort
        module procedure R_sort, I_sort
      end interface

!>  @brief
!>    Returns the indices that would sort an array.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is a generic procedure that accepts integer and
!>  real inputs.
!>  Uses the quicksort algorithm.
!>
!>  \code{.unparsed}
!>  argsort([10, 9, 8, 7, 6])
!>  ! Returns [5, 4, 3, 2, 1]
!>  \endcode
!>
!>  @param[in] a 1D-Array to be sorted
      interface argsort
        module procedure I_argsort, R_argsort
      end interface

      contains

      subroutine I_swap(a, b)
        integer, intent(inout) :: a, b
        integer :: t
        t = a; a = b; b = t
      end subroutine

      subroutine R_swap(a, b)
        real*8, intent(inout) :: a, b
        real*8 :: t
        t = a; a = b; b = t
      end subroutine

      subroutine I_naive_sort(V, increasing)
        implicit none
        integer, intent(inout) :: V(:)
        integer :: t
        logical, intent(in) :: increasing
        integer :: i, j

        if (increasing) then
          do i = 2, size(V)
            t = V(i)
            j = i
            do while (j > 1 .and. V(j - 1) > t)
              V(j) = V(j - 1)
              j = j - 1
            end do
            V(j) = t
          end do
        else
          do i = 2, size(V)
            t = V(i)
            j = i
            do while (j > 1 .and. V(j - 1) < t)
              V(j) = V(j - 1)
              j = j - 1
            end do
            V(j) = t
          end do
        end if
      end subroutine I_naive_sort

      subroutine R_naive_sort(V, increasing)
        implicit none
        real*8, intent(inout) :: V(:)
        real*8 :: t
        logical, intent(in) :: increasing
        integer :: i, j

        if (increasing) then
          do i = 2, size(V)
            t = V(i)
            j = i
            do while (j > 1 .and. V(j - 1) > t)
              V(j) = V(j - 1)
              j = j - 1
            end do
            V(j) = t
          end do
        else
          do i = 2, size(V)
            t = V(i)
            j = i
            do while (j > 1 .and. V(j - 1) < t)
              V(j) = V(j - 1)
              j = j - 1
            end do
            V(j) = t
          end do
        end if
      end subroutine R_naive_sort

      recursive subroutine I_quicksort(V, increasing)
      implicit none
      integer, intent(inout) :: V(:)
      integer :: pivot
      logical, intent(in) :: increasing
      integer :: first, last, i, j

      if (size(V) > naive_sort_trsh) then
        first = 1; last = size(V, 1)
        i = first; j = last

        pivot = V((first + last) / 2)

        do
          if (increasing) then
            do while (V(i) > pivot)
              i = i + 1
            end do
            do while (pivot > V(j))
              j = j - 1
            end do
          else
            do while (V(i) < pivot)
              i = i + 1
            end do
            do while (pivot < V(j))
              j = j - 1
            end do
          end if
          if (i >= j) exit
          call swap(V(i), V(j))
          i = i + 1
          j = j - 1
        end do

        if (first < i - 1) call I_quicksort(V(: i - 1), increasing)
        if (j + 1 < last) call I_quicksort(V(j + 1 :), increasing)
      else
        call naive_sort(V, increasing)
      end if
      end subroutine I_quicksort

      recursive subroutine R_quicksort(V, increasing)
      implicit none
      real*8, intent(inout) :: V(:)
      real*8 :: pivot
      logical, intent(in) :: increasing
      integer :: first, last, i, j

      if (size(V) > naive_sort_trsh) then
        first = 1; last = size(V, 1)
        i = first; j = last

        pivot = V((first + last) / 2)

        do
          if (increasing) then
            do while (V(i) > pivot)
              i = i + 1
            end do
            do while (pivot > V(j))
              j = j - 1
            end do
          else
            do while (V(i) < pivot)
              i = i + 1
            end do
            do while (pivot < V(j))
              j = j - 1
            end do
          end if
          if (i >= j) exit
          call swap(V(i), V(j))
          i = i + 1
          j = j - 1
        end do

        if (first < i - 1) call R_quicksort(V(: i - 1), increasing)
        if (j + 1 < last) call R_quicksort(V(j + 1 :), increasing)
      else
        call naive_sort(V, increasing)
      end if
      end subroutine R_quicksort


      subroutine I_sort(V, increasing)
        implicit none
        integer, intent(inout):: V(:)
        logical, optional :: increasing
        logical :: increasing_
        increasing_ = merge(increasing, .true., present(increasing))
        call quicksort(V, increasing_)
      end subroutine

      subroutine R_sort(V, increasing)
        implicit none
        real*8, intent(inout):: V(:)
        logical, optional :: increasing
        logical :: increasing_
        increasing_ = merge(increasing, .true., present(increasing))
        call quicksort(V, increasing_)
      end subroutine

      function I_argsort(A, increasing) result(idx)
        implicit none
        integer, intent(inout) :: A(:)
        logical, optional :: increasing
        logical :: increasing_
        integer :: idx(size(A)), i
        increasing_ = merge(increasing, .true., present(increasing))
        idx = [(i, i = 1, size(idx))]
        call own_quicksort(idx, increasing_)
        contains
          recursive subroutine own_quicksort(idx, increasing)
            implicit none
            integer, intent(inout) :: idx(:)
            integer :: pivot

            logical, intent(in) :: increasing
            integer :: i, j

            if (size(idx) > naive_sort_trsh) then
              i = 1; j = size(idx)
              pivot = A(idx((1 + size(idx)) / 2))

              do
                if (increasing) then
                  do while (A(idx(i)) < pivot)
                     i = i + 1
                  end do
                  do while (pivot < A(idx(j)))
                     j = j - 1
                  end do
                else
                  do while (A(idx(i)) > pivot)
                     i = i + 1
                  end do
                  do while (pivot > A(idx(j)))
                     j = j - 1
                  end do
                end if
                if (i >= j) exit
                call swap(idx(i), idx(j))
                i = i + 1
                j = j - 1
              end do

              if (2 < i) call quicksort(idx(: i  - 1), increasing)
              if (j + 1 < size(idx)) then
                call quicksort(idx(j + 1 : ), increasing)
              end if
            else
              call own_naive_sort(idx, increasing)
            end if
          end subroutine own_quicksort

          subroutine own_naive_sort(idx, increasing)
            implicit none
            integer, intent(inout) :: idx(:)
            logical, intent(in) :: increasing
            integer :: i, j, t

            if (increasing) then
              do i = 2, size(idx)
                t = idx(i)
                j = i
                do while (j > 1 .and. A(idx((j-1))) > A(t))
                  idx(j) = idx(j - 1)
                  j = j - 1
                end do
                idx(j) = t
              end do
            else
              do i = 2, size(idx)
                t = idx(i)
                j = i
                do while (j > 1 .and. A(idx((j-1))) < A(t))
                  idx(j) = idx(j - 1)
                  j = j - 1
                end do
                idx(j) = t
              end do
            end if
          end subroutine own_naive_sort
      end function I_argsort

      function R_argsort(A, increasing) result(idx)
        implicit none
        real*8, intent(inout) :: A(:)
        logical, optional :: increasing
        logical :: increasing_
        integer :: idx(size(A)), i
        increasing_ = merge(increasing, .true., present(increasing))
        idx = [(i, i = 1, size(idx))]
        call own_quicksort(idx, increasing_)
        contains
          recursive subroutine own_quicksort(idx, increasing)
            implicit none
            integer, intent(inout) :: idx(:)
            real*8 :: pivot

            logical, intent(in) :: increasing
            integer :: i, j

            if (size(idx) > naive_sort_trsh) then
              i = 1; j = size(idx)
              pivot = A(idx((1 + size(idx)) / 2))

              do
                if (increasing) then
                  do while (A(idx(i)) < pivot)
                     i = i + 1
                  end do
                  do while (pivot < A(idx(j)))
                     j = j - 1
                  end do
                else
                  do while (A(idx(i)) > pivot)
                     i = i + 1
                  end do
                  do while (pivot > A(idx(j)))
                     j = j - 1
                  end do
                end if
                if (i >= j) exit
                call swap(idx(i), idx(j))
                i = i + 1
                j = j - 1
              end do

              if (2 < i) call quicksort(idx(: i  - 1), increasing)
              if (j + 1 < size(idx)) then
                call quicksort(idx(j + 1 : ), increasing)
              end if
            else
              call own_naive_sort(idx, increasing)
            end if
          end subroutine own_quicksort

          subroutine own_naive_sort(idx, increasing)
            implicit none
            integer, intent(inout) :: idx(:)
            logical, intent(in) :: increasing
            integer :: i, j, t

            if (increasing) then
              do i = 2, size(idx)
                t = idx(i)
                j = i
                do while (j > 1 .and. A(idx(j - 1)) > A(t))
                  idx(j) = idx(j - 1)
                  j = j - 1
                end do
                idx(j) = t
              end do
            else
              do i = 2, size(idx)
                t = idx(i)
                j = i
                do while (j > 1 .and. A(idx(j - 1)) < A(t))
                  idx(j) = idx(j - 1)
                  j = j - 1
                end do
                idx(j) = t
              end do
            end if
          end subroutine own_naive_sort
      end function R_argsort
      end module
