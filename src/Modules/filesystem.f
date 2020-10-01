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

#include "compiler_features.h"
#include "molcastypes.fh"

      module filesystem
      use, intrinsic :: iso_c_binding, only: c_char, MOLCAS_C_INT,
     &      C_INT, c_ptr, C_NULL_CHAR
      use fortran_strings, only: split, StringWrapper_t, Cptr_to_str,
     &  str
      implicit none
      private
      public :: getcwd_, chdir_, symlink_, get_errno_, strerror_,
     &    mkdir_, remove_, real_path, basename
      interface
      subroutine getcwd_c(path, n, err) bind(C, name="getcwd_wrapper")
        import :: c_char, MOLCAS_C_INT
        implicit none
        character(len=1, kind=c_char), intent(out) :: path(*)
        integer(MOLCAS_C_INT), intent(in) :: n
        integer(MOLCAS_C_INT), intent(out) :: err
      end subroutine

      subroutine chdir_c(path, err) bind(C, name="chdir_wrapper")
        import :: c_char, MOLCAS_C_INT
        implicit none
        character(len=1, kind=c_char), intent(in) :: path(*)
        integer(MOLCAS_C_INT), intent(out) :: err
      end subroutine

      subroutine symlink_c(to, from, err)bind(C, name="symlink_wrapper")
        import :: c_char, MOLCAS_C_INT
        implicit none
        character(len=1, kind=c_char), intent(in) :: to(*), from(*)
        integer(MOLCAS_C_INT), intent(out) :: err
      end subroutine

      subroutine mkdir_c(path, mode, err) bind(C, name="mkdir_wrapper")
        import :: c_char, MOLCAS_C_INT
        implicit none
        character(len=1, kind=c_char), intent(in) :: path(*)
        integer(kind=MOLCAS_C_INT), intent(in) :: mode
        integer(kind=MOLCAS_C_INT), intent(out) :: err
      end subroutine mkdir_c

      function get_errno_c() bind(C, name="get_errno")
        import :: MOLCAS_C_INT
        implicit none
        integer(MOLCAS_C_INT) :: get_errno_c
      end function

#ifdef C_PTR_BINDING
      function strerror_c(errno) bind(C, name="strerror")
        import :: c_ptr, C_INT
        implicit none
        integer(C_INT), value, intent(in) :: errno
        type(c_ptr) :: strerror_c
       end function
#endif

      subroutine remove_c(path, err) bind(C, name="remove_wrapper")
        import :: c_char, MOLCAS_C_INT
        implicit none
        integer(kind=MOLCAS_C_INT), intent(out) :: err
        character(len=1, kind=c_char), intent(in) :: path(*)
      end subroutine remove_c

      end interface

      contains

      subroutine getcwd_(path, err)
        character(len=*), intent(out) :: path
        integer, intent(out), optional :: err
        integer(MOLCAS_C_INT) :: c_err
        call getcwd_c(path, len(path, MOLCAS_C_INT), c_err)
        if (present(err)) err = int(c_err)
      end subroutine

      subroutine chdir_(path, err)
        character(len=*), intent(in) :: path
        integer, intent(out), optional :: err
        integer(MOLCAS_C_INT) :: c_err
        call chdir_c(trim(path)//C_NULL_CHAR, c_err)
        if (present(err)) err = int(c_err)
      end subroutine

      subroutine symlink_(to, from, err)
        character(len=*), intent(in) :: to, from
        integer, intent(out), optional :: err
        integer(MOLCAS_C_INT) :: c_err
        call symlink_c(
     &    trim(to)//C_NULL_CHAR, trim(from)//C_NULL_CHAR, c_err)
        if (present(err)) err = int(c_err)
      end subroutine

      subroutine mkdir_(path, err)
        character(len=*), intent(in) :: path
        integer, optional, intent(out) :: err
        integer(MOLCAS_C_INT) :: loc_err
        call mkdir_c(trim(path)//C_NULL_CHAR,
     &        int(o'772', MOLCAS_C_INT), loc_err)
        if (present(err)) err = loc_err
      end subroutine

      integer function get_errno_()
        get_errno_ = int(get_errno_c())
      end function

!> Return Error String from Error number
      function strerror_(errnum) result(res)
        integer, intent(in) :: errnum
        character(len=:), allocatable :: res
#ifdef C_PTR_BINDING
        res = Cptr_to_str(strerror_c(int(errnum, C_INT)))
#else
        integer :: rc
        character(len=80) :: errstr
        integer, external :: aixerr
        errstr = ''
        rc = aixerr(errstr)
        res = trim(errstr)
#endif
      end function

      subroutine remove_(path, err)
        character(len=*) :: path
        integer, optional, intent(out) :: err
        integer(MOLCAS_C_INT) :: loc_err
        call remove_c(trim(path)//C_NULL_CHAR, loc_err)
        if (present(err)) err = int(loc_err)
      end subroutine

      function real_path(molcas_name) result(path)
        character(len=*), intent(in) :: molcas_name
        character(len=:), allocatable :: path
        character(len=1024) :: buffer
        integer :: L
        call prgmtranslate_master(molcas_name, buffer, L)
        path = buffer(:L)
      end function

      function basename(path) result(res)
        character(len=*), intent(in) :: path
        character(len=:), allocatable :: res
        type(StringWrapper_t), allocatable :: names(:)

        call split(path, '/', names)

        if (size(names(size(names))%str) /= 0) then
        ! Base is a normal file that is not a directory .../.../basename
            res = str(names(size(names))%str)
        else
        ! Base is itself a directory .../.../basename/
            res = str(names(size(names) - 1)%str)
        end if
      end function
      end module filesystem
