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

implicit none
character(len=*) :: line
character(len=100) :: text
integer :: linewidth, textwidth, textoffset

text = adjustl(line)
linewidth = len(line)
textwidth = len_trim(text)
textoffset = (linewidth-textwidth)/2
line = ' '
line(textoffset+1:textoffset+textwidth) = text

end subroutine center_text
