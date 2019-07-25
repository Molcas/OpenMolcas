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
***********************************************************************/

      module fortran_strings
        implicit none
        save
        private
        public :: str, to_lower, to_upper
#include "molcastypes.fh"

!>  @brief
!>    Convert to Fortran string
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is a generic procedure that accepts Fortran integer and
!>  real and C char* pointer arguments..
!>
!>  @param[in] A fortran integer or real, or a C char* pointer.
        interface str
          module procedure I_to_str, R_to_str, Cptr_to_str
        end interface

        interface
          function strlen_c(c_string) bind(C, name='strlen_wrapper')
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr) :: c_string
            integer(MOLCAS_C_INT) :: strlen_c
          end function
        end interface

        character(*), parameter ::
     &      UPPERCASE_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
     &      lowercase_chars = 'abcdefghijklmnopqrtsuvwxyz'

        contains

        function I_to_str(i) result(str)
          character(:), allocatable :: str
          integer, intent(in) :: i
          character(range(i) + 2) :: tmp
          write(tmp, '(I0)') I
          str = trim(tmp)
        end function

        function R_to_str(x) result(str)
          character(:), allocatable :: str
          real*8, intent(in) :: x
          character(range(x) + 2) :: tmp
          write(tmp, '(I0)') x
          str = trim(tmp)
        end function

!> Convert C string pointer to Fortran string.
        function Cptr_to_str(c_string) result(res)
          use, intrinsic :: iso_c_binding
          implicit none
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
          implicit none
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
          implicit none
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

      end module
