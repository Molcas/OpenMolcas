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

module filesystem

#include "compiler_features.h"

use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr, c_null_char
use fortran_strings, only: split, StringWrapper_t, Cptr_to_str, str
use linalg_mod, only: abort_
use Definitions, only: iwp, MOLCAS_C_INT

implicit none
private

public :: getcwd_, chdir_, symlink_, get_errno_, strerror_, mkdir_, &
  remove_, real_path, basename, inquire_, copy_

interface
  subroutine getcwd_c(path,n,err) bind(C,name='getcwd_wrapper')
    import :: c_char, MOLCAS_C_INT
    character(len=1,kind=c_char), intent(out) :: path(*)
    integer(kind=MOLCAS_C_INT), intent(in) :: n
    integer(kind=MOLCAS_C_INT), intent(out) :: err
  end subroutine getcwd_c

  subroutine chdir_c(path,err) bind(C,name='chdir_wrapper')
    import :: c_char, MOLCAS_C_INT
    character(len=1,kind=c_char), intent(in) :: path(*)
    integer(kind=MOLCAS_C_INT), intent(out) :: err
  end subroutine chdir_c

  subroutine symlink_c(to,from,err) bind(C,name='symlink_wrapper')
    import :: c_char, MOLCAS_C_INT
    character(len=1,kind=c_char), intent(in) :: to(*), from(*)
    integer(kind=MOLCAS_C_INT), intent(out) :: err
  end subroutine symlink_c

  subroutine mkdir_c(path,mode,err) bind(C,name='mkdir_wrapper')
    import :: c_char, MOLCAS_C_INT
    character(len=1,kind=c_char), intent(in) :: path(*)
    integer(kind=MOLCAS_C_INT), intent(in) :: mode
    integer(kind=MOLCAS_C_INT), intent(out) :: err
  end subroutine mkdir_c

  function get_errno_c() bind(C,name='get_errno')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: get_errno_c
  end function get_errno_c

# ifdef C_PTR_BINDING
  function strerror_c(errno) bind(C,name='strerror')
    import :: c_ptr, c_int
    type(c_ptr) :: strerror_c
    integer(kind=c_int), value, intent(in) :: errno
  end function strerror_c
# endif

  subroutine remove_c(path,err) bind(C,name='remove_wrapper')
    import :: c_char, MOLCAS_C_INT
    character(len=1,kind=c_char), intent(in) :: path(*)
    integer(kind=MOLCAS_C_INT), intent(out) :: err
  end subroutine remove_c

  subroutine copy_c(src,dst,err) bind(C,name='copy')
    import :: c_char, MOLCAS_C_INT
    character(len=1,kind=c_char), intent(in) :: src(*),dst(*)
    integer(kind=MOLCAS_C_INT), intent(out) :: err
  end subroutine copy_c

  function access_c(path) bind(C,name='access_wrapper')
    import :: c_char, MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT)  :: access_c
    character(len=1,kind=c_char), intent(in) :: path(*)
  end function access_c

end interface

contains

!> @brief
!> Return the current working directory.
subroutine getcwd_(path,err)
  character(len=*), intent(out) :: path
  integer(kind=iwp), intent(out), optional :: err
  integer(kind=MOLCAS_C_INT) :: c_err
  call getcwd_c(path,len(path,MOLCAS_C_INT),c_err)
  if (present(err)) err = int(c_err)
end subroutine getcwd_

!> @brief
!> Change the current directory to `path`.
subroutine chdir_(path,err)
  character(len=*), intent(in) :: path
  integer(kind=iwp), intent(out), optional :: err
  integer(kind=MOLCAS_C_INT) :: c_err
  call chdir_c(trim(path)//c_null_char,c_err)
  if (present(err)) err = int(c_err)
end subroutine chdir_

!> @brief
!> Create a weak symlink_ `from` -> `to`.
subroutine symlink_(to,from,err)
  character(len=*), intent(in) :: to, from
  integer(kind=iwp), intent(out), optional :: err
  integer(kind=MOLCAS_C_INT) :: c_err
  call symlink_c(trim(to)//c_null_char,trim(from)//c_null_char,c_err)
  if (present(err)) err = int(c_err)
end subroutine symlink_

!> @brief
!> Create a directory at `path`.
subroutine mkdir_(path,err)
  character(len=*), intent(in) :: path
  integer(kind=iwp), optional, intent(out) :: err
  integer(kind=MOLCAS_C_INT) :: loc_err
  call mkdir_c(trim(path)//c_null_char,int(o'772',MOLCAS_C_INT),loc_err)
  if (present(err)) err = loc_err
end subroutine mkdir_

!> @brief
!> Return the global error number.
function get_errno_()
  integer(kind=iwp) :: get_errno_
  get_errno_ = int(get_errno_c())
end function get_errno_

!> @brief
!> Return Error String from Error number
function strerror_(errnum) result(res)
  character(len=:), allocatable :: res
  integer(kind=iwp), intent(in) :: errnum
# ifndef C_PTR_BINDING
  integer(kind=iwp) :: rc
  character(len=80) :: errstr
  integer(kind=iwp), external :: aixerr
  errstr = ''
  rc = aixerr(errstr)
  res = trim(errstr)
# else
  res = Cptr_to_str(strerror_c(int(errnum,c_int)))
# endif
end function strerror_

!> @brief
!> Remove the file `path`.
subroutine remove_(path,err)
  character(len=*) :: path
  integer(kind=iwp), optional, intent(out) :: err
  integer(kind=MOLCAS_C_INT) :: loc_err
  call remove_c(trim(path)//c_null_char,loc_err)
  if (present(err)) err = int(loc_err)
end subroutine remove_

!> @brief
!> Translate internal Molcas filenames to system paths.
function real_path(molcas_name) result(path)
  character(len=:), allocatable :: path
  character(len=*), intent(in) :: molcas_name
  character(len=1024) :: buffer
  integer(kind=iwp) :: L
  call prgmtranslate_master(molcas_name,buffer,L)
  path = buffer(:L)
end function real_path

!> @brief
!> Return the basename of `path`.
function basename(path) result(res)
  character(len=:), allocatable :: res
  character(len=*), intent(in) :: path
  type(StringWrapper_t), allocatable :: names(:)

  call split(path,'/',names)

  if (size(names(size(names))%str) /= 0) then
    ! Base is a normal file that is not a directory .../.../basename
    res = str(names(size(names))%str)
  else
    ! Base is itself a directory .../.../basename/
    res = str(names(size(names)-1)%str)
  end if
end function basename

!> @brief
!> Return true if the file or directory `path` exists.
function inquire_(path)
  character(len=*), intent(in) :: path
  logical(kind=iwp) :: inquire_
  inquire_ = access_c(trim(path)//c_null_char) == 0
end function inquire_

!> @brief
!> Copy file from src to dst
!>
!> @details
!> This function is not guaranteed to be thread-safe.
!> As long as `dst` is different for two processes calling at the same
!> time, it is safe to use in process-based parallelisation.
subroutine copy_(src, dst, err)
  character(len=*), intent(in) :: src, dst
  integer(kind=iwp), intent(out), optional :: err
  integer(kind=iwp) :: err_
  call copy_c(trim(src)//c_null_char, trim(dst)//c_null_char, err_)
  if (present(err)) then
    err = err_
  else if (err_ /= 0) then
    call abort_('Error in copy')
  end if
end subroutine copy_

end module filesystem
