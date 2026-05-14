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
! Copyright (C) 2020, Oskar Weser                                      *
!               2021,2023, Ignacio Fdez. Galvan                        *
!***********************************************************************

module test_sorting_mod
    use fruit
    use sorting, only: sort, argsort
    use isotopes, only: maxatomnum, ptab
    use definitions, only: wp
    implicit none
    private
    public :: test_sort_ints, test_sort_isotopes, test_sort_reals

    ! not using the ElementList from the isotopes module because
    ! it is protected and proper initialization would require
    ! access to the isotopes_data.txt file
    ! (putting it here and not in the subroutine because of the lex_alphabet_leq problem, see below)
    type Element
      character(len=2) :: symbol
    end type Element
    type(Element), allocatable :: ElementList(:)

contains

    subroutine test_sort_ints
        integer, parameter :: n_test = 10
        integer :: numbers(n_test), idx(n_test)

        numbers = [10, 13, 78, -1, 2, 5, 4, 0, 0, -3]
        idx(:) = argsort(numbers, leq)

        call assert_true(all(numbers(idx( : size(numbers) - 1)) <= numbers(idx(2 :))))
        call assert_equals(idx, [10, 4, 8, 9, 5, 7, 6, 1, 2, 3], n_test)

        call sort(numbers, leq)

        call assert_true(all(numbers( : size(numbers) - 1) <= numbers(2 :)))
        call assert_equals(numbers, [-3, -1, 0, 0, 2, 4, 5, 10, 13, 78], n_test)

        contains

        ! This could also use sorting_funcs::leq_i
        logical pure function leq(i, j)
            integer, intent(in) :: i, j
            leq = i <= j
        end function
    end subroutine

    subroutine test_sort_reals
        integer, parameter :: n_test = 10
        real(kind=wp) :: numbers(n_test)
        integer :: idx(n_test)

        numbers = [0.471_wp, 0.117_wp, -0.357_wp, 0.318_wp, -0.696_wp, -0.426_wp, 0.854_wp, 0.343_wp, 0.697_wp, -0.570_wp]
        idx = argsort(numbers, leq)

        call assert_true(all(numbers(idx( : size(numbers) - 1)) <= numbers(idx(2 :))))
        call assert_equals(idx, [5, 10, 6, 3, 2, 4, 8, 1, 9, 7], n_test)

        contains

        ! This could also use sorting_funcs::leq_r
        logical pure function leq(i, j)
            real(kind=wp), intent(in) :: i, j
            leq = i <= j
        end function
    end subroutine

    subroutine test_sort_isotopes
        integer, allocatable :: idx(:)
        integer :: i

        allocate(ElementList(MaxAtomNum))
        do i=1,MaxAtomNum
          ElementList(i)%symbol = adjustl(PTab(i))
        end do

        allocate(idx(lbound(elementlist, 1) : ubound(elementlist, 1)))

        idx(:) = [(i, i = lbound(idx, 1), ubound(idx, 1))]

        call sort(idx, lex_alphabet_leq)

        print * , elementlist(idx(: 10))%symbol
        call assert_true(all(elementlist(idx(: 10))%symbol &
            == ['Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'B ', 'Ba']))

    end subroutine

    ! this should be internal to test_sort_isotopes,
    ! but the nvidia/pgi compiler chokes on it (segmentation fault)
    logical pure function lex_alphabet_leq(i, j)
        integer, intent(in) :: i, j
        lex_alphabet_leq = elementlist(i)%symbol <= elementlist(j)%symbol
    end function

end module

program test_sorting
    use fruit
    use test_sorting_mod

    implicit none
    integer :: failed_count, i, seed_size

    call random_seed(size=seed_size)
    call random_seed(put=[(i, i = 1, seed_size)])
    call init_fruit()
    call init_linalg()
    call inimem()

    call test_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) error stop

contains

    subroutine test_driver()
        call run_test_case(test_sort_ints, "test_sort_ints")
        call run_test_case(test_sort_reals, "test_sort_reals")
        call run_test_case(test_sort_isotopes, "test_sort_isotopes")
    end subroutine
end program test_sorting
