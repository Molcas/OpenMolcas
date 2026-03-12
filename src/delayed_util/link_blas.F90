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
! Copyright (C) 2015,2017, Ignacio Fdez. Galvan                        *
!***********************************************************************

module link_blas

use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_procpointer, c_funloc, c_funptr, c_int, c_loc, c_null_char, &
                                       c_null_funptr, c_null_ptr, c_ptr
use iso_c_utilities, only: C_F_string
use dlfcn, only: DL_Info, DLAddr, DLClose, DLError, DLOpen, DLSym, RTLD_global, RTLD_lazy
#include "f1.fh"
use Definitions, only: iwp, u6

implicit none
private

type(c_ptr), allocatable :: handles(:)

! Initializing procedure pointers is a F2008 feature, not supported by all compilers.
! When it is implemented, the pointers below should look like:
!
!   procedure(int_dasum), pointer :: lb_dasum=>int_dasum
!
! and then the exact placement of the initialization call (in start.f) is not critical

#include "f2.fh"

public :: lb_close, lb_initialize

contains

!#define _DEBUGPRINT_

!===============================================================================

subroutine lb_initialize(lib,prlev)

  character(len=*), intent(in) :: lib
  integer(kind=iwp), intent(in) :: prlev
  integer(kind=iwp) :: i
  logical(kind=iwp) :: loaded, try_load
  character(kind=c_char,len=1024) :: libname
  type(c_funptr) :: funptr = c_null_funptr
  type(DL_Info), target :: info
  character(len=1024), allocatable :: libs(:)

  try_load = (lib /= 'Internal')

  call lb_close()

  !*************************************************
  ! Try to load all the libraries specified by "lib"
  !*************************************************
  call split_string(lib,libs)

  allocate(handles(size(libs)))
  handles(:) = c_null_ptr

  loaded = .false.
  if (try_load) then
    do i=1,size(libs)
      libname = libs(i)
      if (len_trim(libs(i)) > 0) then
        handles(i) = dlopen(trim(libname)//c_null_char,int(ior(rtld_global,rtld_lazy),kind=c_int))
        loaded = (loaded .or. c_associated(handles(i)))
        if (prlev > 0) then
          if (c_associated(handles(i))) then
            write(u6,*) trim(libname),' loaded'
          else
            write(u6,*) c_f_string(dlerror())
          end if
        end if
      end if
    end do
  end if

  deallocate(libs)

  if (loaded) then
    !******************************************************************
    ! Associate all BLAS and LAPACK routines with the library functions
    !******************************************************************
#   include "f3.fh"

  else
    !**************************************
    ! Or use the fallback internal routines
    !**************************************
    if (prlev > 0) write(u6,*) 'Using internal BLAS+LAPACK'
#   include "f4.fh"
  end if

  if (prlev > 0) then
#   include "f5.fh"
  end if

end subroutine lb_initialize

!===============================================================================

subroutine lb_close()

  integer(kind=iwp) :: i, rc

  if (allocated(handles)) then
    do i=1,size(handles)
      if (c_associated(handles(i))) then
        rc = dlclose(handles(i))
        if (rc /= 0) write(u6,*) c_f_string(dlerror())
      end if
    end do
    deallocate(handles)
  end if

end subroutine lb_close

!===============================================================================

subroutine split_string(string,array)

  character(len=*) :: string
  integer(kind=iwp) :: i, lim, n, offset
  character(len=*), allocatable :: array(:)

  n = 1
  do i=1,len(string)
    if (string(i:i) == ':') n = n+1
  end do
  allocate(array(n))
  offset = 0
  do i=1,n-1
    lim = offset+index(string(offset+1:),':')
    array(i) = string(offset+1:lim-1)
    offset = lim
  end do
  array(n) = string(offset+1:)

end subroutine split_string

!===============================================================================

function link_func(funname)

  type(c_funptr) :: link_func
  character(kind=c_char,len=*), intent(in) :: funname
  integer(kind=iwp) :: i
  logical(kind=iwp) :: success

  !****************************************************
  ! Try to associate a function on all loaded libraries
  !****************************************************
  success = .false.
  link_func = c_null_funptr
  do i=1,size(handles)
    link_func = dlsym(handles(i),trim(funname)//'_'//c_null_char)
    if (c_associated(link_func)) then
      success = .true.
      exit
#   ifdef _DEBUGPRINT_
    else
      write(u6,*) c_f_string(dlerror())
#   endif
    end if
  end do
# ifdef _DEBUGPRINT_
  if (.not. success) write(u6,*) 'no ',trim(funname),' found'
# endif

end function link_func

end module link_blas
