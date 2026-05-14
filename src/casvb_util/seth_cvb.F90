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

subroutine seth_cvb(ivec,n)
! Buffered integer IO with integer offset (ncnt)

use casvb_global, only: ibuf, ibuffer, lbuf, ncnt
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, ivec(n)
integer(kind=iwp) :: i_max, i_min, ibuf_max, ibuf_min, ivec_max, ivec_offs, jbuf
logical(kind=iwp) :: full_buffer
logical(kind=iwp), parameter :: debug = .false.

if (n <= 0) return

ibuf_min = ncnt/lbuf+1
ibuf_max = (ncnt+n-1)/lbuf+1
ivec_offs = 1
do jbuf=ibuf_min,ibuf_max
  i_min = max(1,ncnt+1-(jbuf-1)*lbuf)
  i_max = min(lbuf,ncnt+n-(jbuf-1)*lbuf)
  full_buffer = (i_min == 1) .and. (i_max == lbuf)

  if (ibuf /= jbuf) then
    call bufio_wrbuf_cvb()
    call bufio_chbuf_cvb(jbuf)
    if (.not. full_buffer) call bufio_rdbuf_cvb()
  end if
  ivec_max = ivec_offs+i_max-i_min
  ibuffer(i_min:i_max) = ivec(ivec_offs:ivec_max)
  ivec_offs = ivec_max+1
end do

if (debug) then
  write(u6,*) ' wrbis :',n,ncnt
  write(u6,'(40i4)') ivec
end if
ncnt = ncnt+n

return

end subroutine seth_cvb
