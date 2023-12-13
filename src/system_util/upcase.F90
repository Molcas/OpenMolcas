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
! Copyright (C) Per Ake Malmqvist                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine uppercases a text string.                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-AAke Malmqvist                                          *
!          Lund University                                             *
!                                                                      *
!***********************************************************************

subroutine UpCase(string)

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: string
integer(kind=iwp), save :: ifset = 0, itab(0:255)
integer(kind=iwp) :: i, ii
character(len=*), parameter :: up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', &
                               lw = 'abcdefghijklmnopqrstuvwxyz'

if (ifset == 0) then
  ifset = 1
  do i=0,255
    itab(i) = i
  end do
  do ii=1,26
    i = ichar(lw(ii:ii))
    itab(i) = ichar(up(ii:ii))
  end do
end if

do ii=1,len(string)
  i = ichar(string(ii:ii))
  string(ii:ii) = char(itab(i))
end do

return

end subroutine UpCase
