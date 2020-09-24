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

module test_fortran_strings_mod
    use fruit
    use fortran_strings, only: StringWrapper_t, split, str
    implicit none
    private
    public :: test_split

contains

    subroutine test_split()
        character(*), parameter :: path = '/home/mustermann/file.txt'
        character(*), parameter :: names = 'Schroedinger Born Oppenheimer'

        type(StringWrapper_t), allocatable :: splitted(:)

        call split(path, '/', splitted)
        call assert_true(str(splitted(1)%str) == '')
        call assert_true(str(splitted(2)%str) == 'home')
        call assert_true(str(splitted(3)%str) == 'mustermann')
        call assert_true(str(splitted(4)%str) == 'file.txt')

        call split(names, ' ', splitted)
        call assert_true(str(splitted(1)%str) == 'Schroedinger')
        call assert_true(str(splitted(2)%str) == 'Born')
        call assert_true(str(splitted(3)%str) == 'Oppenheimer')
    end subroutine

end module test_fortran_strings_mod

program test_fortran_strings
    use fruit
    use test_fortran_strings_mod

    implicit none
    integer :: failed_count
    integer, parameter :: seed(20) = &
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]


    call init_fruit()
    call random_seed(put=seed)
    call inimem()

    call test_fortran_strings_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) error stop

contains

    subroutine test_fortran_strings_driver()
        call run_test_case(test_split, "test_split")
    end subroutine
end program test_fortran_strings
