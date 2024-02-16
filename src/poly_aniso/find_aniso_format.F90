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

subroutine find_aniso_format(old_aniso_format)

use Definitions, only: u5, u6

implicit none
logical :: old_aniso_format
character(len=280) :: line
integer :: LINENR

old_aniso_format = .false.

!=========== End of default settings====================================
rewind(u5)
50 read(u5,'(A280)',end=998) LINE
call NORMAL(LINE)
if (LINE(1:11) /= '&POLY_ANISO') Go To 50
LINENR = 0
100 read(u5,'(A280)',end=998) line
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') Go To 100
if (LINE == ' ') Go To 100
if (LINE(1:4) /= 'OLDA') Go To 100
if ((LINE(1:4) == 'END ') .or. (LINE(1:4) == '    ')) Go To 200

if (line(1:4) == 'OLDA') then

  old_aniso_format = .true.

  LINENR = LINENR+1
  Go To 100
end if

200 continue

Go To 190
!------ errors ------------------------------
998 continue
write(u6,*) ' READIN: Unexpected End of input file.'

190 continue
write(u6,*) 'find_aniso_format::  old_aniso_format=',old_aniso_format

return

end subroutine find_aniso_format
