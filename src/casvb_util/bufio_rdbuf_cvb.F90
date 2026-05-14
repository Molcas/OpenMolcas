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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine bufio_rdbuf_cvb()

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use casvb_global, only: file_id, ibuf, ibuffer, lbuf, nbuf, nword
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp) :: ioffset

if (nbuf < ibuf) then
  ibuffer(:) = 0
  return
end if
ioffset = (ibuf-1)*lbuf/RtoI
call bufio_rdbuf_cvb_internal(ibuffer)

return

! This is to allow type punning without an explicit interface
contains

subroutine bufio_rdbuf_cvb_internal(ibuffer)
  integer(kind=iwp), target :: ibuffer(*)
  real(kind=wp), pointer :: buffer(:)
  call c_f_pointer(c_loc(ibuffer(1)),buffer,[nword])
  call rdlow_cvb(buffer,nword,file_id,ioffset+1)
  nullify(buffer)
end subroutine bufio_rdbuf_cvb_internal

end subroutine bufio_rdbuf_cvb
