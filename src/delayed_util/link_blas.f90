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
  use iso_c_binding
  use iso_c_utilities
  use dlfcn
#include "f1.fh"

  implicit none

  type(c_ptr), dimension(:), allocatable, private :: handles
!
! Initializing procedure pointers is a F2008 feature, not supported by all compilers.
! When it is implemented, the pointers below should look like:
!
!   procedure(int_dasum), pointer :: lb_dasum=>int_dasum
!
! and then the exact placement of the initialization call (in start.f) is not critical
!
#include "f2.fh"

contains

!define _DEBUGPRINT_

!===============================================================================

  subroutine lb_initialize(lib,prlev)
    character(len=*), intent(in) :: lib
    integer, intent(in) :: prlev
    character(len=1024), dimension(:), allocatable :: libs
    character(kind=c_char,len=1024) :: libname
    type(c_funptr) :: funptr=c_null_funptr
    type(DL_Info), target :: info
    integer :: i
    logical :: try_load,loaded

    try_load = (lib /= 'Internal')

    call lb_close()

!***************************************************
!   Try to load all the libraries specified by "lib"
!***************************************************
    call split_string(lib,libs)

    allocate(handles(size(libs)))
    handles(:)=c_null_ptr

    loaded=.false.
    if (try_load) then
      do i=1,size(libs)
        libname=libs(i)
        if (len_trim(libs(i)) > 0) then
          handles(i)=dlopen(trim(libname)//c_null_char, int(ior(rtld_global,rtld_lazy),kind=c_int))
          loaded=(loaded .or. c_associated(handles(i)))
          if (prlev > 0) then
            if (c_associated(handles(i))) then
              write (6,*) trim(libname),' loaded'
            else
              write(6,*) c_f_string(dlerror())
            end if
          end if
        end if
      end do
    end if

    deallocate(libs)

!********************************************************************
!   Associate all BLAS and LAPACK routines with the library functions
!********************************************************************
    if (loaded) then
#include "f3.fh"
!
!
!****************************************
!   Or use the fallback internal routines
!****************************************
    else
      if (prlev > 0) then
        write(6,*) 'Using internal BLAS+LAPACK'
      end if
#include "f4.fh"
    end if

    if (prlev > 0) then
#include "f5.fh"
    end if

  end subroutine lb_initialize

!===============================================================================

  subroutine lb_close()
    integer :: i,rc

    if (allocated(handles)) then
      do i=1,size(handles)
        if (c_associated(handles(i))) then
          rc = dlclose(handles(i))
          if (rc /= 0) then
            write(6,*) c_f_string(dlerror())
          end if
        end if
      end do
      deallocate(handles)
    end if

  end subroutine lb_close

!===============================================================================

  subroutine split_string(string,array)
    character(len=*) :: string
    character(len=*), dimension(:), allocatable :: array
    integer :: i,n,offset,lim

    n=1
    do i=1,len(string)
      if (string(i:i) == ':') n=n+1
    end do
    allocate(array(n))
    offset=0
    do i=1,n-1
      lim=offset+index(string(offset+1:),':')
      array(i)=string(offset+1:lim-1)
      offset=lim
    end do
    array(n)=string(offset+1:)

  end subroutine split_string

!===============================================================================

  function link_func(funname)
    character(kind=c_char,len=*) :: funname
    type(c_funptr) :: link_func
    integer :: i
    logical :: success

!******************************************************
!   Try to associate a function on all loaded libraries
!******************************************************
    success=.false.
    link_func=c_null_funptr
    do i=1,size(handles)
      link_func=dlsym(handles(i),trim(funname)//'_'//c_null_char)
      if (c_associated(link_func)) then
        success=.true.
        exit
#ifdef _DEBUGPRINT_
      else
        write(6,*) c_f_string(dlerror())
#endif
      end if
    end do
#ifdef _DEBUGPRINT_
    if (.not. success) then
      write(6,*) 'no ',trim(funname),' found'
    end if
#endif

  end function link_func

end module link_blas
