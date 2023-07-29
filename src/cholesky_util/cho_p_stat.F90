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

subroutine Cho_P_Stat()

use Cholesky, only: Cho_Real_Par, LuRed, LuRed_G, nSym, NumCho, NumCho_G, NumChT, NumChT_G
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iTmp, jTmp

if (Cho_Real_Par) then
  call Cho_P_IndxSwp()
  call iSwap(nSym,NumCho,1,NumCho_G,1)
  jTmp = NumChT
  NumChT = NumChT_G
  iTmp = LuRed
  LuRed = LuRed_G
  call Cho_Stat()
  LuRed = iTmp
  NumChT = jTmp
  call iSwap(nSym,NumCho,1,NumCho_G,1)
  call Cho_P_IndxSwp()
else
  call Cho_Stat()
end if

end subroutine Cho_P_Stat
