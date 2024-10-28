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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003-2005, Valera Veryazov                             *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine DmpLLs(iDskPt)

use LnkLst, only: DmpLst, Init_LLs, LLDelt, LLdGrd, LLGrad, LLx, LLy
use SCFFiles, only: LuDel, LuDgd, LuGrd, Lux, Luy
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iDskPt(5)

if (Init_LLs) then
  !call StatLLs()
  call DmpLst(LLGrad,LuGrd,iDskPt(1))
  call DmpLst(LLDgrd,LuDGd,iDskPt(2))
  call DmpLst(LLDelt,LuDel,iDskPt(3))
  call DmpLst(LLy,Lux,iDskPt(4))
  call DmpLst(LLx,Luy,iDskPt(5))
else
  write(u6,*) '****** W A R N I N G ! ******'
  write(u6,*) ' Linked list already killed!'
end if

return

end subroutine DmpLLs
