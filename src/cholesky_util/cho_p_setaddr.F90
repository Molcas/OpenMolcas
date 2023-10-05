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

subroutine Cho_P_SetAddr()
!
! Purpose: set initial disk adresses for local as well as global
!          reduced sets.

use Cholesky, only: Cho_Real_Par, InfRed, InfRed_G, InfVec, InfVec_G, LuPri, MaxRed, MaxVec, nSym, XnPass
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: irc
character(len=*), parameter :: SecNam = 'Cho_P_SetAddr'

if (Cho_Real_Par) then

  ! The variable XnPass must be zero: restart is not possible.
  ! ----------------------------------------------------------

  if (XnPass /= 0) call Cho_Quit('XnPass>0 error in '//SecNam,104)

  ! Global.
  ! -------

  call Cho_P_SetAddr_2(InfRed_G,InfVec_G,MaxRed,MaxVec,size(InfVec,2),nSym,irc)
  if (irc /= 0) then
    write(Lupri,*) SecNam,': Cho_P_SetAddr_2 returned ',irc
    call Cho_Quit('Error in '//SecNam,104)
  end if

end if

! Local.
! ------

call Cho_SetAddr(InfRed,InfVec,MaxRed,MaxVec,size(InfVec,2),nSym)

end subroutine Cho_P_SetAddr
