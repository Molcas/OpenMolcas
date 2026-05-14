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

subroutine NxtWrd(Line,i_F,iE)
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

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Line
integer(kind=iwp), intent(inout) :: i_F, iE
integer(kind=iwp) :: nChar

nChar = len(Line)
! Find first non-blank character
do
  if ((i_F == 0) .or. (i_F > nChar)) then
    call WarningMessage(2,'NxtWrd: (i_F == 0) .or. (i_F > nChar)')
    write(u6,*) 'nChar=',nChar
    write(u6,*) 'i_F,iE=',i_F,iE
    call Abend()
  end if
  if (Line(i_F:i_F) /= ' ') exit
  i_F = i_F+1
  if (i_F >= nChar) then
    i_F = nChar
    iE = -1
    return
  end if
end do
! Find the end of the present word
iE = i_F+1
do
  if (Line(iE:iE) /= ' ') then
    iE = iE+1
    if (iE > nChar) then
      iE = nChar
      return
    end if
  else
    iE = iE-1
    exit
  end if
end do

return

end subroutine NxtWrd
