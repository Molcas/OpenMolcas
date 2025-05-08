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
! Copyright (C) 2012, Victor P. Vysotskiy                              *
!               2025, Ignacio Fdez. Galvan                             *
!***********************************************************************

#include "compiler_features.h"
#ifdef _GARBLE_

subroutine garble(ipos,length,vartyp)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr
use Definitions, only: wp, iwp, byte, RtoB

implicit none
integer(kind=iwp), intent(in) :: ipos, length
character(len=*), intent(in) :: vartyp
integer(kind=iwp) :: ioff1, ioff2, foff1
type(c_ptr) :: cptr
integer(kind=iwp), pointer :: ibuf(:)
integer(kind=byte), pointer :: i1buf(:)
real(kind=wp), pointer :: rbuf(:)
integer(kind=iwp), parameter :: igarbage = huge(igarbage)
integer(kind=byte), parameter :: i1garbage = huge(i1garbage)
real(kind=wp), parameter :: dgarbage = huge(dgarbage)
interface
  function woff2cptr(etyp,offset)
    import :: c_ptr, iwp
    type(c_ptr) :: woff2cptr
    character, intent(in) :: etyp
    integer(kind=iwp), value, intent(in) :: offset
  end function woff2cptr
end interface

! Here we overwrite the underlying memory, using C pointers
! Do not try this at home!

select case (vartyp)
  case ('REAL')
    cptr = woff2cptr('R',ipos-1)
    call c_f_pointer(cptr,rbuf,[length])
    rbuf(1:length) = dgarbage
    nullify(rbuf)
  case ('INTE')
    cptr = woff2cptr('I',ipos-1)
    call c_f_pointer(cptr,ibuf,[length])
    ibuf(1:length) = igarbage
    nullify(ibuf)
  case ('CHAR')
    ioff1 = (ipos-1)/RtoB+1
    ioff2 = mod(ipos-1,RtoB)+1
    foff1 = (ipos+length-2)/RtoB+1
    cptr = woff2cptr('C',ipos-1)
    call c_f_pointer(cptr,i1buf,[(foff1-ioff1+1)*RtoB])
    i1buf(ioff2:ioff2+length-1) = i1garbage
    nullify(i1buf)
end select

end subroutine garble

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(garble)

#endif
