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
        public ::
     &    sort, argsort,  tSortScheme, sort_schemes, naive_sort_trsh

        type :: tSortSchemeValues
          integer ::
     &      mergesort = 1,
     &      quicksort = 2
        end type
        type(tSortSchemeValues), parameter ::
! TODO: Dear fellow MOLCAS developer of the future:
! Please replace the following explicit constructur
! with the default constructor
!     instance = type()
! as soon as possible.
! As of July 2019 the Sun compiler requires explicit construction
! for parameter variables. (Which is wrong IMHO.)
     &    sort_schemes =
     &      tSortSchemeValues(mergesort = 1, quicksort = 2)

        type :: tSortScheme
          integer :: val = sort_schemes%mergesort
        end type

! When to switch to naive sort, which scales with n^2,
! but requires less overhead.
        integer :: naive_sort_trsh = 30

        interface
          logical pure function compare_int(a, b)
            integer, intent(in) :: a, b
          end function

          logical pure function compare_real(x, y)
            real*8, intent(in) :: x, y
          end function
        end interface

!>  @brief
!>    Returns the indices that would sort a vector.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The vector is sorted until f(v(i), v(i + 1)) is true for all i.
!>  If you want to sort the vector increasingly, just use
!>  \code{.unparsed}
!>  logical pure function le(i, j)
!>    integer, intent(in) :: i, j
!>    le = i <= j
!>  end function
!>  \endcode
!>  If you want to change the sorting algorithm, you can do so on the fly.
!>  The default is always a stable sorting algorithm.
!>  \code{.unparsed}
!>   call argsort(V, le, tSortScheme(sort_schemes%mergesort))
!>   call argsort(V, le, tSortScheme(sort_schemes%quicksort))
!>  \endcode
!>
!>  @param[in] V Integer or real vector to be sorted
!>  @param[in] compare A logical pure function of two integer or real arguments.
!>  @param[in] scheme The sorting algorithm to use.
        interface argsort
          module procedure I1D_argsort, R1D_argsort
        end interface

        contains


        function I1D_argsort(V, compare, scheme) result(idx)
          integer, intent(inout) :: V(:)
          procedure(compare_int) :: compare
          type(tSortScheme), intent(in), optional :: scheme
          type(tSortScheme)  :: scheme_
          integer :: idx(lbound(V, 1):ubound(V, 1)), i
          if (present(scheme)) then
            scheme_ = scheme
          else
! Use default constructor, whatever the default is.
            scheme_ = tSortScheme()
          end if

          idx = [(i, i = lbound(V, 1), ubound(V, 1))]

          call sort(idx, my_compare, scheme_)

          contains
            logical pure function my_compare(x, y)
              integer, intent(in) :: x, y
              my_compare = compare(V(x), V(y))
            end function
        end function I1D_argsort


        function R1D_argsort(V, compare, scheme) result(idx)
          real*8, intent(inout) :: V(:)
          procedure(compare_real) :: compare
          type(tSortScheme), intent(in), optional :: scheme
          type(tSortScheme)  :: scheme_
          integer :: idx(lbound(V, 1):ubound(V, 1)), i
          if (present(scheme)) then
            scheme_ = scheme
          else
! Use default constructor, whatever the default is.
            scheme_ = tSortScheme()
          end if

          idx = [(i, i = lbound(V, 1), ubound(V, 1))]

          call sort(idx, my_compare, scheme_)

          contains
            logical pure function my_compare(x, y)
              integer, intent(in) :: x, y
              my_compare = compare(V(x), V(y))
            end function
        end function R1D_argsort

!>  @brief
!>    Sort the array V inplace.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The array is sorted until f(v(i), v(i + 1)) is true for all i.
!>  If you want to sort an array increasingly, just use
!>  \code{.unparsed}
!>  logical pure function le(i, j)
!>    integer, intent(in) :: i, j
!>    le = i <= j
!>  end function
!>  \endcode
!>  The real power of this routine shines, when the array
!>  is treated as an index array of other matrices or tensors.
!>  If you want to sort a 2D-matrix according to the column sum,
!>  just use:
!>  \code{.unparsed}
!>   integer :: col_idx(lbound(matrix, 2):ubound(matrix, 2))
!>   col_idx = [(i, i = lbound(col_idx, 1), ubound(col_idx, 1))]
!>   call sort(col_idx, col_sum)
!>
!>   contains
!>
!>    logical pure function col_sum(i, j)
!>      integer, intent(in) :: i, j
!>      col_sum = sum(abs(matrix(:, i))) <= sum(abs(matrix(:, j)))
!>    end function
!>  \endcode
!>  If you want to change the sorting algorithm, you can do so on the fly.
!>  The default is always a stable sorting algorithm.
!>  \code{.unparsed}
!>   call sort(col_idx, col_sum, tSortScheme(sort_schemes%mergesort))
!>   call sort(col_idx, col_sum, tSortScheme(sort_schemes%quicksort))
!>  \endcode
!>
!>  @param[in] V Integer 1D-Array to be sorted
!>  @param[in] compare A logical pure function of two integer arguments.
!>  @param[in] scheme The sorting algorithm to use.
        subroutine sort(V, compare, scheme)
          integer, intent(inout) :: V(:)
          procedure(compare_int) :: compare
          type(tSortScheme), intent(in), optional :: scheme
          type(tSortScheme)  :: scheme_
          if (present(scheme)) then
            scheme_ = scheme
          else
! Use default constructor, whatever the default is.
            scheme_ = tSortScheme()
          end if

          select case (scheme_%val)
            case(sort_schemes%mergesort)
              call mergesort(V, compare)
            case(sort_schemes%quicksort)
              call quicksort(V, compare)
          end select
        end subroutine sort

        subroutine naive_sort(V, compare)
          integer, intent(inout) :: V(:)
          procedure(compare_int) :: compare

          integer :: i, j, t

          do i = lbound(V, 1) + 1, ubound(V, 1)
            t = V(i)
            j = i
            do while (j > 1 .and. .not. compare(V(j - 1), t))
              V(j) = V(j - 1)
              j = j - 1
            end do
            V(j) = t
          end do
        end subroutine naive_sort

        subroutine swap(a, b)
          integer, intent(inout) :: a, b
          integer :: t
          t = a; a = b; b = t
        end subroutine

!>  @brief
!>    Merge two sorted arrays into sorted array.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  A and B are in-parameters and C is an out-parameter, if
!>  non overlapping memory is used.
!>  Either A or B is allowed to overlap with the second half
!>  of the output array C.
!>
!>  @param[inout] A Sorted 1D-array to be merged.
!>  @param[inout] B Sorted 1D-array to be merged.
!>  @param[inout] C Merged and sorted 1D-array.
!>  @param[in] compare A logical pure function of two integer arguments.
        subroutine merge_(A, B, C, compare)
! The target attribute is there to prevent the compiler from
! assuming non overlapping memory.
          integer, target, intent(inout) :: A(:), B(:), C(:)
          procedure(compare_int) :: compare

          integer :: i, j, k

          if (size(A) + size(B) > size(C)) stop

          i = lbound(A, 1)
          j = lbound(B, 1)
          do k = lbound(C, 1), ubound(C, 1)
            if (i <= ubound(A, 1) .and. j <= ubound(B, 1)) then
              if (compare(A(i), B(j))) then
                C(k) = A(i)
                i = i + 1
              else
                C(k) = B(j)
                j = j + 1
              end if
            else if (i <= ubound(A, 1)) then
              C(k) = A(i)
              i = i + 1
            else if (j <= ubound(B, 1)) then
              C(k) = B(j)
              j = j + 1
            end if
          end do
        end subroutine merge_


        recursive subroutine mergesort(A, compare)
          integer, intent(inout) :: A(:)
          procedure(compare_int) :: compare

          integer, allocatable :: work(:)
          integer :: half
          half = (ubound(A, 1) - lbound(A, 1)) / 2 + 1
          allocate(work(half))
          call mergesort_work(A, compare, work)
          deallocate(work)
        end subroutine mergesort

        recursive subroutine mergesort_work(A, compare, work)
          integer, intent(inout) :: A(:)
          procedure(compare_int) :: compare
          integer, intent(inout) :: work(:)

          integer :: half
          half = (ubound(A, 1) - lbound(A, 1)) / 2 + 1
          if (size(A) < 2) then
            continue
          else if (size(A) == 2) then
            call naive_sort(A, compare)
          else
            call mergesort_work(A( : half), compare, work)
            call mergesort_work(A(half + 1 :), compare, work)
            if (.not. compare(A(half), A(half + 1))) then
              work(1 : half) = A(1 : half)
              call merge_(work(1 : half), A(half + 1:), A, compare)
            endif
          end if
        end subroutine mergesort_work


        recursive subroutine quicksort(idx, compare)
          integer, intent(inout) :: idx(:)
          procedure(compare_int) :: compare

          integer :: i, j, pivot

          if (size(idx) > naive_sort_trsh) then
            i = lbound(idx, 1)
            j = ubound(idx, 1)
            pivot = idx((j - i) / 2 + 1)

            do
              do while (.not. compare(idx(i), pivot))
                 i = i + 1
              end do
              do while (.not. compare(pivot, idx(j)))
                 j = j - 1
              end do
              if (i >= j) exit
              call swap(idx(i), idx(j))
              i = i + 1
              j = j - 1
            end do

            if (lbound(idx, 1) + 1 < i) then
              call quicksort(idx(: i  - 1), compare)
            end if
            if (j + 1 < ubound(idx, 1)) then
              call quicksort(idx(j + 1 : ), compare)
            end if
          else
            call naive_sort(idx, compare)
          end if
        end subroutine quicksort
      end module sorting
