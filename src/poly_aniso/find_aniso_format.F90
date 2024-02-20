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

use Definitions, only: iwp, u5, u6

implicit none
logical(kind=iwp), intent(out) :: old_aniso_format
character(len=280) :: line
integer(kind=iwp) :: istatus, LINENR

old_aniso_format = .false.

!=========== End of default settings====================================
rewind(u5)
do
  read(u5,'(A280)',iostat=istatus) LINE
  if (istatus < 0) then
    call Error()
    write(u6,*) 'find_aniso_format::  old_aniso_format=',old_aniso_format
    return
  end if
  call NORMAL(LINE)
  if (LINE(1:11) == '&POLY_ANISO') exit
end do
LINENR = 0
do
  read(u5,'(A280)',iostat=istatus) line
  if (istatus < 0) then
    call Error()
    write(u6,*) 'find_aniso_format::  old_aniso_format=',old_aniso_format
    exit
  end if
  LINENR = LINENR+1
  call NORMAL(LINE)
  if ((LINE(1:1) == '*') .or. (LINE == ' ')) cycle
  select case (LINE(1:4))
    case ('END ','    ')
      exit

    case ('OLDA')

      old_aniso_format = .true.

      LINENR = LINENR+1
  end select
end do

write(u6,*) 'find_aniso_format::  old_aniso_format=',old_aniso_format

return

contains

subroutine Error()

  write(u6,*) ' FIND_ANISO_FORMAT: Unexpected End of input file.'

end subroutine Error

end subroutine find_aniso_format
