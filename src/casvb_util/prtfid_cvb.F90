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

subroutine prtfid_cvb(chr,fileid)

implicit real*8(a-h,o-z)
#include "io_cvb.fh"
character*200 line
character*(*) chr

line = chr
call mkfn_cvb(fileid,ibf)
call appendchr_cvb(line,' file ',0)
call appendchr_cvb(line,filename(ibf),1)
call appendchr_cvb(line,'.',0)
write(6,'(a)') line(1:len_trim_cvb(line))

return

end subroutine prtfid_cvb
