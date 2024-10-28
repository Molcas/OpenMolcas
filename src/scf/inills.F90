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

subroutine IniLLs()
! initialize the diverse linked lists

use LnkLst, only: IniLst, Init_LLs, LLDelt, LLdGrd, LLGrad, LLlist, LLx, LLy
use InfSCF, only: MxOptm

implicit none

LLlist = 0
LLGrad = 0
call IniLst(LLGrad,20)
call IniLst(LLDgrd,20)
call IniLst(LLDelt,20)
call IniLst(LLy,20)
call IniLst(LLx,MxOptm)
Init_LLs = .true.

end subroutine IniLLs
