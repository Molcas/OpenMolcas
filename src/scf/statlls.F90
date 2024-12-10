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

subroutine StatLLS()

use LnkLst, only: Init_LLs, LLDelt, LLdGrd, LLGrad, LLlGrd, LLx, LLy, StlLst
use Definitions, only: u6

implicit none

if (Init_LLs) then
  call StlLst(LLGrad)
  call StlLst(LLlGrd)
  call StlLst(LLDgrd)
  call StlLst(LLDelt)
  call StlLst(LLy)
  call StlLst(LLx)
else
  write(u6,*) '****** W A R N I N G ! ******'
  write(u6,*) ' Linked lists are not there!'
end if

return

end subroutine StatLLS
