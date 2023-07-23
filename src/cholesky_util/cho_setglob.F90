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

subroutine Cho_SetGlob()
!
! Purpose: define entries in choglob.fh

use Definitions, only: iwp

implicit none
#include "choglob.fh"
integer(kind=iwp) :: iSym, N
integer(kind=iwp), parameter :: iLarge = 999999

nnShl_G = 0
mmBstRT_G = 0

N = 8*nLoc_G
call iZero(iiBstR_G,N)
call iZero(nnBstR_G,N)
call iZero(nnBstRT_G,nLoc_G)
call iZero(NumCho_G,8)
NumChT_G = 0

do iSym=1,8
  LuCho_G(iSym) = -iLarge
end do
LuRed_G = -iLarge
LuRst_G = -iLarge

end subroutine Cho_SetGlob
