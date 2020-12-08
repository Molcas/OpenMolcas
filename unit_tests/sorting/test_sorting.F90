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
!***********************************************************************

#include "compiler_features.h"

module test_sorting_mod
    use fruit
    use sorting, only: sort
    use isotopes, only: initialize_isotopes, elementlist, Element_t
    use definitions, only: wp
    implicit none
    private
    public :: test_sort_isotopes, test_sort_reals

contains

    subroutine test_sort_reals
        integer, parameter :: n_test = 10
        integer :: numbers(n_test)

        numbers = [10, 13, 78, -1, 2, 5, 4, 0, 0, -3]
        call sort(numbers, leq)

        call assert_true(all(numbers( : size(numbers) - 1) <= numbers(2 :)))
        call assert_equals(numbers, [-3, -1, 0, 0, 2, 4, 5, 10, 13, 78], n_test)

        contains

        logical pure function leq(i, j)
            integer, intent(in) :: i, j
            leq = i <= j
        end function
    end subroutine

    subroutine test_sort_isotopes
        integer, allocatable :: idx(:)

        integer :: i

        call initialize_isotopes()

        allocate(idx(lbound(elementlist, 1) : ubound(elementlist, 1)))
        idx(:) = [(i, i = lbound(idx, 1), ubound(idx, 1))]

        call sort(idx, lex_alphabet_leq)

        call assert_true(all(elementlist(idx(: 10))%symbol &
            == ['Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'B ', 'Ba']))


        contains

            logical pure function lex_alphabet_leq(i, j)
                integer, intent(in) :: i, j
                lex_alphabet_leq = elementlist(i)%symbol <= elementlist(j)%symbol
            end function

    end subroutine

end module

program test_sorting
    use fruit
    use test_sorting_mod

    implicit none
    integer :: failed_count, i
    integer, parameter :: seed_size = 50
    integer, parameter :: seed(seed_size) = [(i, i = 1, seed_size)]

    call random_seed(put=seed)
    call init_fruit()
    call inimem()

    call test_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) error stop

contains

    subroutine test_driver()
        call run_test_case(test_sort_reals, "test_sort_reals")
        call run_test_case(test_sort_isotopes, "test_sort_isotopes")
    end subroutine
end program test_sorting
