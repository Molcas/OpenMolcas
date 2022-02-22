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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine EvalERFun(ERFun,ERFunC,CMO,nOcc,nSym,Timing)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: ERFun
real(kind=wp), intent(_OUT_) :: ERFunC(*)
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp), intent(in) :: nSym, nOcc(nSym)
logical(kind=iwp), intent(in) :: Timing
integer(kind=iwp) :: irc
character(len=80) :: Txt
character(len=*), parameter :: SecNam = 'EvalERFun'

irc = 0
call Cho_Get_ER(irc,CMO,nOcc,ERFunC,ERFun,Timing)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'Cho_Get_ER returned',irc
  call SysAbendMsg(SecNam,'ER evaluation failed!',Txt)
end if

end subroutine EvalERFun
