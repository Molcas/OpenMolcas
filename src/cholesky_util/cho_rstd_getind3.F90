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

subroutine Cho_RstD_GetInd3(iSP2F,l_iSP2F)

use Cholesky, only: LuRed, nnBstRT, nnShl, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: l_iSP2F
integer(kind=iwp), intent(out) :: iSP2F(l_iSP2F)
integer(kind=iwp) :: iAdr, iOpt

iOpt = 2
iAdr = nSym*nnShl+2*nnBstRT(1)
call iDAFile(LuRed,iOpt,iSP2F,l_iSP2F,iAdr)

end subroutine Cho_RstD_GetInd3
