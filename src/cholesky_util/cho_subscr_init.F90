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

subroutine Cho_SubScr_Init()
!
! Purpose: initialize screening in vector subtraction.

use Cholesky, only: DSPNm, DSubScr, nnBstR, nnShl, nSym
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iSym, l_DSubScr

l_DSubScr = nnBstR(1,1)
do iSym=2,nSym
  l_DSubScr = max(l_DSubScr,nnBstR(iSym,1))
end do
call mma_allocate(DSubScr,l_DSubScr,Label='DSubScr')

call mma_allocate(DSPNm,nnShl,Label='DSPNm')

end subroutine Cho_SubScr_Init
