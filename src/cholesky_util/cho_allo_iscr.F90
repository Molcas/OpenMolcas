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

subroutine Cho_Allo_iScr(DoDummy)
!
! Purpose: allocate iScr array for reading and reordering vectors.
!          If (DoDummy): make dummy (length 1) allocation.

use Cholesky, only: iScr, nnBstR, nSym
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: DoDummy
integer(kind=iwp) :: iSym, l_iScr

if (DoDummy) then
  l_iScr = 1
else
  l_iScr = nnBstR(1,1)
  do iSym=2,nSym
    l_iScr = max(l_iScr,nnBstR(iSym,1))
  end do
end if
call mma_allocate(iScr,l_iScr,Label='iScr')

end subroutine Cho_Allo_iScr
