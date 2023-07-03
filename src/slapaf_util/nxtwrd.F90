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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine NxtWrd(Line,if,iE)
!***********************************************************************
!                                                                      *
! Object:                                                              *
!                                                                      *
! Called from: DefInt                                                  *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!***********************************************************************

implicit real*8(A-H,O-Z)
character*(*) Line

nChar = len(Line)
! Find first non-blank character
10 continue
if ((if == 0) .or. (if > nChar)) then
  call WarningMessage(2,'NxtWrd: (iF == 0) .or. (iF > nChar)')
  write(6,*) 'nChar=',nChar
  write(6,*) 'iF,iE=',if,iE
  call Abend()
end if
if (Line(if:if) == ' ') then
  if = if+1
  if (if >= nChar) then
    if = nChar
    iE = -1
    return
  end if
  Go To 10
end if
! Find the end of the present word
iE = if+1
20 continue
if (Line(iE:iE) /= ' ') then
  iE = iE+1
  if (iE > nChar) then
    iE = nChar
    return
  end if
  Go To 20
else
  iE = iE-1
end if

return

end subroutine NxtWrd
