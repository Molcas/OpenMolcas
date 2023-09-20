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

subroutine bufio_wrbuf_cvb()

use casvb_global, only: file_id, ibuf, ibuffer, lbuf, nbuf, nword

implicit real*8(a-h,o-z)
#include "idbl_cvb.fh"

call bufio_wrbuf_internal(ibuffer)

return

! This is to allow type punning without an explicit interface
contains

subroutine bufio_wrbuf_internal(ibuffer)
  use iso_c_binding
  integer, target :: ibuffer(*)
  real*8, pointer :: buffer(:)
  if (ibuf == 0) return

  ioffset = (ibuf-1)*lbuf/idbl
  call c_f_pointer(c_loc(ibuffer(1)),buffer,[nword])
  call wrlow_cvb(buffer,nword,file_id,ioffset+1)
  nullify(buffer)
  if (ibuf > nbuf) then
    nbuf = ibuf
  end if
  return
end subroutine bufio_wrbuf_internal

end subroutine bufio_wrbuf_cvb
