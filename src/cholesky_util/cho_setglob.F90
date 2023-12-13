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
! Purpose: define entries in Cholesky.

use Cholesky, only: iiBstR_G, LuCho_G, LuRed_G, LuRst_G, mmBstRT_G, nnBstR_G, nnBstRT_G, nnShl_G, NumCho_G, NumChT_G
use Definitions, only: iwp

implicit none
integer(kind=iwp), parameter :: iLarge = 999999

nnShl_G = 0
mmBstRT_G = 0

iiBstR_G(:,:) = 0
nnBstR_G(:,:) = 0
nnBstRT_G(:) = 0
NumCho_G(:) = 0
NumChT_G = 0

LuCho_G(:) = -iLarge
LuRed_G = -iLarge
LuRst_G = -iLarge

end subroutine Cho_SetGlob
