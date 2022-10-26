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

integer function nToken(Line)

character*(*) Line
logical On

nLine = len(Line)
nToken = 0
On = .true.

iLine = 1
99 continue
if (iLine == nLine) return
if (Line(iLine:iLine) == ' ') then
  On = .true.
else
  if (On) nToken = nToken+1
  On = .false.
end if
iLine = iLine+1
Go To 99

!write(6,*) 'Error in nToken!'
!call Abend()

end function nToken
