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

subroutine center_text(line)
!SVC: centers the text of a line

use Definitions, only: iwp

implicit none
character(len=*), intent(inout) :: line
character(len=len(line)) :: text
integer(kind=iwp) :: linewidth, textwidth, textoffset

text = adjustl(line)
linewidth = len(line)
textwidth = len_trim(text)
textoffset = (linewidth-textwidth)/2
if (textoffset > 0) then
  line = ' '
  line(textoffset+1:textoffset+textwidth) = text
end if

end subroutine center_text
