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

subroutine Cho_P_QualSwp()

use Cholesky, only: iQuAB, iQuAB_L, nQual, nQual_L, nSym, pTemp
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: scr(8)

! Swap nQual array  local <-> global
scr(1:nSym) = nQual_L(1:nSym)
nQual_L(1:nSym) = nQual(1:nSym)
nQual(1:nSym) = scr(1:nSym)

! Swap pointer for iQuAB  local <-> global
pTemp => iQuAB
iQuAB => iQuAB_L
iQuAB_L => pTemp

end subroutine Cho_P_QualSwp
