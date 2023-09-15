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

subroutine rdline_init_cvb(variat)

implicit real*8(a-h,o-z)
! BLANKDELIM signifies whether blanks are used to delimit fields:
logical blankdelim, variat
#include "luinp_cvb.fh"
#include "rdline.fh"
parameter(nblank=2)
character*1 blanks(nblank)
save blankdelim
save blanks
data blankdelim/.true./
data blanks/' ',','/

if (variat) return
rewind(inp)
2100 read(inp,'(a)',end=2200) line
lenline = len_trim_cvb(line)
call strip_blanks_cvb(line,lenline,blanks,nblank,blankdelim)
call upper_case_cvb(line,lenline)
if (line(1:6) /= '&CASVB') goto 2100
!if (.not. (((.not. blankdelim) .and. (line(1:10) == '&CASVB&END')) .or. (blankdelim .and. (line(1:11) == '&CASVB &END')))) goto 2100

return

2200 write(6,*) ' WARNING: Initiation string not found in input file.'

return

end subroutine rdline_init_cvb
