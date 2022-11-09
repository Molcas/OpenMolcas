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

function nToken(Line)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nToken
character(len=*), intent(in) :: Line
integer(kind=iwp) :: iLine, nLine
logical(kind=iwp) :: On

nLine = len(Line)
nToken = 0
On = .true.

do iLine=1,nLine-1
  if (Line(iLine:iLine) == ' ') then
    On = .true.
  else
    if (On) nToken = nToken+1
    On = .false.
  end if
end do

end function nToken
