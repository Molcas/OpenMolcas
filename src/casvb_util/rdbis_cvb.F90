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

implicit real*8(a-h,o-z)
#include "idbl_cvb.fh"
dimension ivec(n)
logical debug
data debug/.false./

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
  call imove_cvb(ibuffer(i_min),ivec(ivec_offs),i_max-i_min+1)
  ivec_offs = ivec_offs+i_max-i_min+1
end do

if (debug) then
  write(6,*) ' rdbis :',n,ioffset
  write(6,'(40i4)') ivec
end if
ioffset = ioffset+n

return

end subroutine rdbis_cvb
