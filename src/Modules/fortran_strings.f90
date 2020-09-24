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
! Copyright (C) 2019, Oskar Weser                                      *
!***********************************************************************

#include "molcastypes.fh"

module fortran_strings
    use, intrinsic :: iso_c_binding, only: c_ptr, MOLCAS_C_INT, c_f_pointer
    implicit none
    save
    private
    public :: str, to_lower, to_upper, operator(.in.), split, &
        count_char, StringWrapper_t, Cptr_to_str

    ! This type exists to have an array of string pointers
    ! and to allow unequally sized strings.
    ! NOTE: Due to old compilers this had to be
    ! character(1), dimension(:), allocatable
    ! if possible change it to
    ! character(:), allocatable
    type :: StringWrapper_t
        character(1), allocatable :: str(:)
    end type


    !>  @brief
    !>    Convert to Fortran string
    !>
    !>  @author Oskar Weser
    !>
    !>  @details
    !>  It is a generic procedure that accepts Fortran integer,
    !>  Fortran real, and Fortran arrays with single character elements.
    !>
    !>  @param[in] A Fortran integer or real, or character array.
    interface str
        module procedure I_to_str, R_to_str, character_array_to_str
    end interface

    interface
        pure function strlen_c(c_string) bind(C, name='strlen_wrapper')
            import :: c_ptr, MOLCAS_C_INT
            type(c_ptr), intent(in) :: c_string
            integer(MOLCAS_C_INT) :: strlen_c
        end function
    end interface

    interface operator(.in.)
        module procedure substr_in_str
    end interface

    character(*), parameter :: &
        UPPERCASE_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', &
        lowercase_chars = 'abcdefghijklmnopqrstuvwxyz'

    contains

    pure function I_to_str(i) result(str)
        character(:), allocatable :: str
        integer, intent(in) :: i
        character(range(i) + 2) :: tmp
        write(tmp, '(I0)') I
        str = trim(tmp)
    end function

    pure function R_to_str(x) result(str)
        character(:), allocatable :: str
        real*8, intent(in) :: x
        character(range(x) + 2) :: tmp
        write(tmp, '(I0)') x
        str = trim(tmp)
    end function

    pure function character_array_to_str(array) result(res)
        character(1), intent(in) :: array(:)
        character(len=:), allocatable :: res
        integer :: i, L

        L = size(array)
        allocate(character(len=L) :: res)
        do i = 1, L
            res(i:i) = array(i)
        end do
    end function

    !> Convert C string pointer to Fortran string.
    function Cptr_to_str(c_string) result(res)
        type(c_ptr), intent(in) :: c_string
        character(len=:), allocatable :: res
        character, pointer, dimension(:) :: string
        integer :: i, L
        L = int(strlen_c(c_string))
        allocate(character(len=L) :: res)
        call c_f_pointer(c_string, string, [L])
        do i = 1, L
            res(i:i) = string(i)
        end do
    end function

    !> Changes a string to upper case
    pure function to_upper (in_str) result (string)
        character(*), intent(in) :: in_str
        character(len(in_str)) :: string
        integer :: ic, i, L

        L = len_trim(in_str)
        do i = 1, L
            ic = index(lowercase_chars, in_str(i:i))
            if (ic > 0) then
                string(i:i) = UPPERCASE_chars(ic:ic)
            else
                string(i:i) = in_str(i:i)
            end if
        end do
        string(L + 1: ) = ' '
    end function to_upper

    !> Changes a string to lower case
    pure function to_lower (in_str) result (string)
        character(*), intent(in) :: in_str
        character(len(in_str)) :: string
        integer :: ic, i, L

        L = len_trim(in_str)
        do i = 1, L
            ic = index(UPPERCASE_chars, in_str(i:i))
            if (ic > 0) then
                string(i:i) = lowercase_chars(ic:ic)
            else
                string(i:i) = in_str(i:i)
            end if
        end do
        string(L + 1: ) = ' '
    end function to_lower

    logical pure function substr_in_str(substring, string)
        character(*), intent(in) :: string, substring

        substr_in_str = index(string, substring) /= 0
    end function

    !> @brief
    !> Split a string at delimiter.
    pure subroutine split(str, delimiter, res)
        character(*), intent(in) :: str
        character(1), intent(in) :: delimiter
        type(StringWrapper_t), allocatable, intent(out) :: res(:)

        integer :: i, n, low

        allocate(res(count_char(str, delimiter) + 1))

        low = 1; n = 1
        do i = 1, len(str)
            if (str(i : i) == delimiter) then
                res(n)%str = str(low : i - 1)
                n = n + 1
                low = i + 1
            end if
        end do

        if (n == size(res)) then
            res(n)%str = str(low : )
        end if
    end subroutine

    !> @brief
    !> Count the occurence of a character in a string.
    pure function count_char(str, char) result(c)
        character(*), intent(in) :: str
        character(1), intent(in) :: char
        integer :: c
        integer :: i

        c = 0
        do i = 1, len(str)
            if (str(i : i) == char) c = c + 1
        end do
    end function
end module
