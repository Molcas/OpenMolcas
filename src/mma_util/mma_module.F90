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
use Definitions, only: wp, iwp, MOLCAS_C_INT

implicit none
private

logical(kind=iwp) :: MemStat = .false.

! Fortran interfaces for C functions

interface

  function allocmem(size_) bind(C,name='allocmem_')
    import :: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: allocmem
    integer(kind=MOLCAS_C_INT), intent(out) :: size_
  end function allocmem

  function c_getmem(name_,Op,dtyp,offset,len_) bind(C,name='c_getmem_')
    import :: c_char, MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_getmem
    character(kind=c_char), intent(in) :: name_(*), Op(*), dtyp(*)
    integer(kind=MOLCAS_C_INT), intent(inout) :: offset
    integer(kind=MOLCAS_C_INT), intent(in) :: len_
  end function c_getmem

  function cptr2woff(ptr) bind(C,name='cptr2woff_')
    import :: c_ptr, MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: cptr2woff
    type(c_ptr), value, intent(in) :: ptr
  end function cptr2woff

  function woff2cptr(offset) bind(C,name='woff2cptr_')
    import :: c_ptr, MOLCAS_C_INT
    type(c_ptr) :: woff2cptr
    integer(kind=MOLCAS_C_INT), value, intent(in) :: offset
  end function woff2cptr

end interface

public :: allocmem, c_getmem, cptr2woff, MemStat
#ifdef _GARBLE_
public :: garble

contains

subroutine garble(ipos,length,vartyp)

  use, intrinsic :: iso_c_binding, only: c_f_pointer
# ifdef _FPE_TRAP_
  use, intrinsic :: IEEE_Arithmetic, only: IEEE_signaling_NaN, IEEE_value
# endif
  use Definitions, only: byte, RtoB

  integer(kind=iwp), intent(in) :: ipos, length
  character(len=*), intent(in) :: vartyp
  integer(kind=iwp) :: ioff1, ioff2, foff1
  type(c_ptr) :: cptr
  integer(kind=iwp), pointer :: ibuf(:)
  integer(kind=byte), pointer :: i1buf(:)
  real(kind=wp), pointer :: rbuf(:)
  integer(kind=iwp), parameter :: igarbage = huge(igarbage)
  integer(kind=byte), parameter :: i1garbage = huge(i1garbage)
# ifndef _FPE_TRAP_
  real(kind=wp), parameter :: dgarbage = huge(dgarbage)
# else
  real(kind=wp) :: dgarbage

  dgarbage = IEEE_value(dgarbage,IEEE_signaling_NaN)
# endif

  ! Here we overwrite the underlying memory, using C pointers
  ! Do not try this at home!

  select case (vartyp)
    case ('REAL')
      cptr = woff2cptr(ipos-1)
      call c_f_pointer(cptr,rbuf,[length])
      rbuf(1:length) = dgarbage
      nullify(rbuf)
    case ('INTE')
      cptr = woff2cptr(ipos-1)
      call c_f_pointer(cptr,ibuf,[length])
      ibuf(1:length) = igarbage
      nullify(ibuf)
    case ('CHAR')
      ioff1 = (ipos-1)/RtoB+1
      ioff2 = mod(ipos-1,RtoB)+1
      foff1 = (ipos+length-2)/RtoB+1
      cptr = woff2cptr(ipos-1)
      call c_f_pointer(cptr,i1buf,[(foff1-ioff1+1)*RtoB])
      i1buf(ioff2:ioff2+length-1) = i1garbage
      nullify(i1buf)
    case default
      call Abend()
  end select

end subroutine garble

#endif

end module mma_module
