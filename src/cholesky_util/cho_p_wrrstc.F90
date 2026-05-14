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

subroutine Cho_P_WrRstC(iPass)
!
! Purpose: write global restart info to disk.

use Cholesky, only: Cho_Real_Par, LuRst, LuRst_G, nSym, NumCho, NumCho_G, TMISC
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iPass
integer(kind=iwp) :: iTmp
real(kind=wp) :: c1, c2, w1, w2

call CWTime(c1,w1)
if (Cho_Real_Par) then
  call Cho_P_IndxSwp()
  call iSwap(nSym,NumCho,1,NumCho_G,1)
  iTmp = LuRst
  LuRst = LuRst_G
  call Cho_WrRstC(iPass)
  LuRst = iTmp
  call iSwap(nSym,NumCho,1,NumCho_G,1)
  call Cho_P_IndxSwp()
else
  call Cho_WrRstC(iPass)
end if
call CWTime(c2,w2)
tMisc(1,3) = tMisc(1,3)+c2-c1
tMisc(2,3) = tMisc(2,3)+w2-w1

end subroutine Cho_P_WrRstC
