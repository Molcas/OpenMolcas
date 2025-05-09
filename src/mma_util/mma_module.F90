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

module mma_module

use, intrinsic :: iso_c_binding, only: c_char, c_ptr
use Definitions, only: wp, iwp, MOLCAS_C_INT, MOLCAS_C_REAL

implicit none
private

real(kind=wp) :: Work(1)
logical(kind=iwp) :: MemStat = .false.

! Fortran interfaces for C functions

interface

  function allocmem(ref,size_) bind(C,name='allocmem_')
    import :: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: allocmem
    real(kind=MOLCAS_C_REAL), intent(in) :: ref(*)
    integer(kind=MOLCAS_C_INT), intent(out) :: size_
  end function allocmem

  function c_getmem(name_,Op,dtyp,offset,len_) bind(C,name='c_getmem_')
    import :: c_char, MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_getmem
    character(kind=c_char), intent(in) :: name_(*), Op(*), dtyp(*)
    integer(kind=MOLCAS_C_INT), intent(inout) :: offset
    integer(kind=MOLCAS_C_INT), intent(in) :: len_
  end function c_getmem

  function cptr2woff(etyp,ptr) bind(C,name='cptr2woff_')
    import :: c_char, c_ptr, MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: cptr2woff
    character(kind=c_char), intent(in) :: etyp(*)
    type(c_ptr), value, intent(in) :: ptr
  end function

  function woff2cptr(etyp,offset) bind(C,name='woff2cptr_')
    import :: c_char, c_ptr, MOLCAS_C_INT
    type(c_ptr) :: woff2cptr
    character(kind=c_char), intent(in) :: etyp
    integer(kind=MOLCAS_C_INT), value, intent(in) :: offset
  end function woff2cptr

end interface

public :: allocmem, c_getmem, cptr2woff, MemStat, woff2cptr, Work

end module mma_module
