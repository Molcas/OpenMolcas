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

subroutine FindErrorLine()

use getline_mod, only: iGetLine, MyUnit
use Definitions, only: iwp, u6

implicit none
character(len=180) :: line
integer(kind=iwp) :: isave, istatus, lunit

lunit = myunit
isave = igetline
rewind(lunit)
do
  read(lunit,'(a)',iostat=istatus) Line
  if (istatus > 0) exit
  call UpCase(Line)
  Line = adjustl(Line)
  if (Line(1:1) == '&') then
    line = line(2:)
    exit
  end if
end do
if (istatus <= 0) then
  igetline = 0
  write(u6,'(a,a,a)') ' >>>>> Input file for module ',line(1:index(line,' ')),' <<<<<'
  do
    read(lunit,'(A)',iostat=istatus) line
    if (istatus /= 0) exit
    igetline = igetline+1
    if (igetline == isave) then
      write(u6,*) '******   Error  *******'
      write(u6,'(a)') line
      write(u6,'(a)')
      call WarningMessage(2,'Error in FindErrorLine')
      call Quit_OnUserError()
    end if
    if (isave-igetline <= 50) write(u6,'(a)') line
  end do
  !  write(u6,'(a)') ' >>>>> Input error <<<<<'
  !  rewind(lunit)
  !  igetline = 0
  !end do
end if
call WarningMessage(1,'FindErrorLine: Error in input was not located;  Please, check it manually!')

return

end subroutine FindErrorLine
