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

implicit none
real*8 ERFun
real*8 ERFunC(*), CMO(*)
integer nSym
integer nOcc(nSym)
logical Timing
character*9 SecNam
parameter(SecNam='EvalERFun')
character*80 Txt
integer irc

irc = 0
call Cho_Get_ER(irc,CMO,nOcc,ERFunC,ERFun,Timing)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'Cho_Get_ER returned',irc
  call SysAbendMsg(SecNam,'ER evaluation failed!',Txt)
end if

end subroutine EvalERFun
