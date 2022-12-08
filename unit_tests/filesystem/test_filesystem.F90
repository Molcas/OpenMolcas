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
! Copyright (C) 2021, Oskar Weser                                      *
!***********************************************************************

#include "compiler_features.h"

module test_filesystem_mod
    use fruit
    use filesystem, only: basename, inquire_, mkdir_, remove_, chdir_, getcwd_, copy_
    implicit none
    private
    public :: test_filesystem

contains

    subroutine test_filesystem()
        integer :: err
        character(*), parameter :: test_dir = 'test_dir'
        character(len=512) :: cwd, root

        call getcwd_(root)

        if (inquire_(test_dir)) then
            call remove_(test_dir, err)
            call assert_true(err == 0)
        end if
        call mkdir_(test_dir, err)
        call assert_true(err == 0)
        call assert_true(inquire_(test_dir))

        call chdir_(test_dir)
            call getcwd_(cwd)
            call assert_true(basename(cwd) == test_dir)

            block
                integer :: file_id, err
                character(*), parameter :: test_file = 'asdf', new_file = 'hallo Wrzlbrmft'
                call assert_false(inquire_(test_file))
                open(newunit=file_id, file=test_file)
                    write(file_id, '(A)') 'Hello World'
                close(file_id)
                call assert_true(inquire_(test_file))
                call copy_(test_file, new_file)
                call assert_true(inquire_(new_file))
                call remove_(new_file)
                call remove_(test_file, err)
                call assert_true(err == 0)
                call assert_false(inquire_(test_file))

                call copy_(test_file, new_file, err)
                call assert_true(err /= 0)
            end block
        call chdir_(root)

        call remove_(test_dir, err)
        call assert_true(err == 0)
        call assert_false(inquire_(test_dir))
    end subroutine

end module test_filesystem_mod

program test_filesystem_prog
    use fruit
    use test_filesystem_mod, only: test_filesystem

    implicit none
    integer :: failed_count, i, seed_size

    call random_seed(size=seed_size)
    call random_seed(put=[(i, i = 1, seed_size)])
    call init_fruit()
    call inimem()

    call run_test_case(test_filesystem, "test_filesystem")

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) error stop

end program test_filesystem_prog
