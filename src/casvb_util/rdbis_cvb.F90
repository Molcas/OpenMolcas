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

subroutine rdbis_cvb(ivec,n,ioffset)
! Buffered integer IO with integer offset

use casvb_global, only: ibuf, ibuffer, lbuf
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(out) :: ivec(n)
integer(kind=iwp), intent(inout) :: ioffset
integer(kind=iwp) :: i_max, i_min, ibuf_max, ibuf_min, ivec_max, ivec_offs, jbuf
logical(kind=iwp), parameter :: debug = .false.

if (n <= 0) return

ibuf_min = ioffset/lbuf+1
ibuf_max = (ioffset+n-1)/lbuf+1
ivec_offs = 1
do jbuf=ibuf_min,ibuf_max
  i_min = max(1,ioffset+1-(jbuf-1)*lbuf)
  i_max = min(lbuf,ioffset+n-(jbuf-1)*lbuf)

  if (ibuf /= jbuf) then
    ! Following line only needed if reads and writes can be mixed
    call bufio_wrbuf_cvb()
    call bufio_chbuf_cvb(jbuf)
    call bufio_rdbuf_cvb()
  end if
  ivec_max = ivec_offs+i_max-i_min
  ivec(ivec_offs:ivec_max) = ibuffer(i_min:i_max)
  ivec_offs = ivec_max+1
end do

if (debug) then
  write(u6,*) ' rdbis :',n,ioffset
  write(u6,'(40i4)') ivec
end if
ioffset = ioffset+n

return

end subroutine rdbis_cvb
