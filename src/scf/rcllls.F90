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

!#define _DEBUGPRINT_
subroutine RclLLs(iDskPt)

use InfSCF, only: MemRsv
use LnkLst, only: Init_LLs, LLDelt, LLdGrd, LLGrad, LLlGrd, LLx, LLy, RclLst
use SCFFiles, only: LuDel, LuDgd, LuGrd, LulGd, Lux, Luy
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: iDskPt(6)

call RclLst(LLGrad,LuGrd,iDskPt(1),MemRsv)
call RclLst(LLlGrd,LulGd,iDskPt(2),MemRsv)
call RclLst(LLDgrd,LuDGd,iDskPt(3),MemRsv)
call RclLst(LLDelt,LuDel,iDskPt(4),MemRsv)
call RclLst(LLy,Lux,iDskPt(5),MemRsv)
call RclLst(LLx,Luy,iDskPt(6),MemRsv)
Init_LLs = .true.
#ifdef _DEBUGPRINT_
call StatLLs()
#endif

return

end subroutine RclLLs
