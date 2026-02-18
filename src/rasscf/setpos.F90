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
integer(kind=iwp) :: IPRLEV, istatus, KLen
character(len=16) :: Command, Key
#include "warnings.h"

! Read until, and including, a line beginning with a particular
! string in an ASCII file, assumed already opened, with unit
! number LUnit. That line is returned.
! Key lengths up to 16 bytes can be used, it is determined by
! the size of the input variable.

IPRLEV = IPRLOC(1)
iRc = _RC_ALL_IS_WELL_
KLen = min(len(Key),len(KeyIn))
Key = ' '
Command = ' '
rewind(LUnit)

Key = KeyIn(1:KLen)
call upcase(Key)
do
  read(LUnit,'(A)',iostat=istatus) Line
  if (istatus /= 0) then
    if (IPRLEV >= TERSE) then
      if (istatus < 0) then
        write(u6,*) ' SETPOS: Attempt to find an input line beginning'
      else
        write(u6,*) ' SETPOS: Attempt to find an input line beginning'
      end if
      write(u6,*) ' with the keyword "',KeyIn,'" failed.'
    end if
    !call Quit(_RC_INPUT_ERROR_)
    iRc = _RC_INPUT_ERROR_
    exit
  end if
  Command(1:KLen) = Line(1:KLen)
  call upcase(Command)
  if (Command == Key) exit
end do

end subroutine SetPos
