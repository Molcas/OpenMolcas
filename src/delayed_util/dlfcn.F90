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

module DLFCN

use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_int, c_ptr
use ISO_C_UTILITIES, only: C_F_STRING

implicit none
private

! Valid modes for mode in DLOpen, obtained from the output of a C program:
integer, parameter :: RTLD_LAZY = 1, RTLD_GLOBAL = 256

! Struct for DLAddr
type, bind(C) :: DL_info
  type(c_ptr) :: dli_fname
  type(c_funptr) :: dli_fbase
  type(c_ptr) :: dli_sname
  type(c_funptr) :: dli_saddr
end type DL_info

public :: DL_info, DLAddr, DLClose, DLError, DLOpen, DLSym, RTLD_GLOBAL, RTLD_LAZY

interface ! All we need is interfaces for the prototypes in <dlfcn.h>

  function DLOpen(file,mode) result(handle) bind(C,name="dlopen")
    ! void *dlopen(const char *file, int mode);

    import :: c_char, c_int, c_ptr
    character(kind=c_char), dimension(*), intent(in) :: file
    ! C strings should be declared as character arrays
    integer(c_int), value :: mode
    type(c_ptr) :: handle

  end function DLOpen

  function DLSym(handle,name) result(funptr) bind(C,name="dlsym")
    ! void *dlsym(void *handle, const char *name);

    import :: c_char, c_funptr, c_ptr
    type(c_ptr), value :: handle
    character(kind=c_char), dimension(*), intent(in) :: name
    type(c_funptr) :: funptr ! A function pointer

  end function DLSym

  function DLClose(handle) result(status) bind(C,name="dlclose")
    ! int dlclose(void *handle);

    import :: c_int, c_ptr
    type(c_ptr), value :: handle
    integer(c_int) :: status

  end function DLClose

  function DLError() result(error) bind(C,name="dlerror")
    ! char *dlerror(void);

    import :: c_ptr
    type(c_ptr) :: error

  end function DLError

  ! dladdr is a Glibc extension, not POSIX
  function DLAddr(funptr,info) result(output) bind(C,name="dladdr")
    ! int dladdr(void *addr, Dl_info *info)

    import :: c_funptr, c_int, c_ptr
    type(c_funptr), value :: funptr ! A function pointer
    type(c_ptr), value :: info
    integer(c_int) :: output

  end function DLAddr

end interface

end module DLFCN
