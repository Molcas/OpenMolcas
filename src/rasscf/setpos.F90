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

subroutine SetPos(LUnit,KeyIn,Line,iRc)

use PrintLevel, only: TERSE
use output_ras, only: IPRLOC
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: LUNIT, iRC
character(len=*) KeyIn, Line
integer(kind=iwp) :: IPRLEV, KLen
character(len=16) :: Command, Key
#include "warnings.h"

! Read until, and including, a line beginning with a particular
! string in an ASCII file, assumed already opened, with unit
! number LUnit. That line is returned.
! Key lengths up to 16 bytes can be used, it is determined by
! the size of the input variable.

IPRLEV = IPRLOC(1)
iRc = _RC_ALL_IS_WELL_
KLen = min(16,len(KeyIn))
Key = ' '
Command = ' '
rewind(LUnit)

Key(1:KLen) = KeyIn(1:KLen)
call upcase(Key)
10 continue
read(LUnit,'(A)',end=9910,err=9920) Line
Command(1:KLen) = Line(1:KLen)
call upcase(Command)
if (Command /= Key) goto 10
return

!---  Error exits ----------------------
9910 continue
if (IPRLEV >= TERSE) then
  write(u6,*) ' SETPOS: Attempt to find an input line beginning'
  write(u6,*) ' with the keyword "',KeyIn,'" failed.'
end if
!call Quit(_RC_INPUT_ERROR_)
iRc = _RC_INPUT_ERROR_
return
9920 continue
if (IPRLEV >= TERSE) then
  write(u6,*) ' SETPOS: Attempt to find an input line beginning'
  write(u6,*) ' with the keyword "',KeyIn,'" failed.'
end if
!call Quit(_RC_INPUT_ERROR_)
iRc = _RC_INPUT_ERROR_

end subroutine SetPos
