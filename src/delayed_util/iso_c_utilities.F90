!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module ISO_C_UTILITIES

use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, c_ptr, c_size_t

implicit none
private

character(kind=c_char), save, target :: dummy_string(1) = "?"

public :: C_F_STRING

contains

function C_F_STRING(CPTR) result(FPTR)
  ! Convert a null-terminated C string into a Fortran character array pointer

  interface ! strlen is a standard C function from <string.h>
    ! int strlen(char *string)
    function strlen(string) result(len) bind(C,name="strlen")
      import :: c_size_t, c_ptr
      type(c_ptr), value :: string ! A C pointer
      integer(kind=c_size_t) :: len
    end function
  end interface

  type(c_ptr), intent(in) :: CPTR ! The C address
  character(kind=c_char), pointer :: FPTR(:)

  if (c_associated(CPTR)) then
    call c_f_pointer(FPTR=FPTR,CPTR=CPTR,shape=[strlen(CPTR)])
  else
    ! To avoid segfaults, associate FPTR with a dummy target:
    FPTR => dummy_string
  end if

end function

end module ISO_C_UTILITIES
