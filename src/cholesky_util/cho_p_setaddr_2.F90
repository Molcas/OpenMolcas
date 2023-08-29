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

subroutine Cho_P_SetAddr_2(InfRed,InfVec,MaxRed,MaxVec,N2,nSym,irc)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: MaxRed, MaxVec, N2, nSym
integer(kind=iwp), intent(inout) :: InfRed(MaxRed), InfVec(MaxVec,N2,nSym)
integer(kind=iwp), intent(out) :: irc

irc = 0

if (MaxRed > 0) then
  InfRed(1) = 0
else
  irc = 1
  return
end if

if ((MaxVec > 0) .and. (N2 >= 4)) then
  InfVec(1,3:4,1:nSym) = 0
else
  irc = 2
  return
end if

end subroutine Cho_P_SetAddr_2
