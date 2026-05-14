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

use casvb_global, only: filename
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: chr
real(kind=wp), intent(in) :: fileid
integer(kind=iwp) :: ibf
character(len=200) :: line

line = chr
call mkfn_cvb(fileid,ibf)
call appendchr_cvb(line,' file ',0)
call appendchr_cvb(line,filename(ibf),1)
call appendchr_cvb(line,'.',0)
write(u6,'(a)') trim(line)

return

end subroutine prtfid_cvb
