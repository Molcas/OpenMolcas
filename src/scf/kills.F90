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

subroutine KiLLs()
! dispose the diverse linked lists

use LnkLst, only: LLGrad, LLdGrd, LLDelt, LLy, LLx, Init_LLs, KilLst

implicit none

if (Init_LLs) then
  call KilLst(LLGrad)
  call KilLst(LLDgrd)
  call KilLst(LLDelt)
  call KilLst(LLy)
  call KilLst(LLx)
  Init_LLs = .false.
else
  write(6,*) '****** W A R N I N G ! ******'
  write(6,*) ' Linked list already killed!'
end if

return

end subroutine KiLLs
