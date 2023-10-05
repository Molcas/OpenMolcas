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

subroutine Cho_MaxDX(Diag,Dmax)
!
! Purpose: get max. diagonal elements in each sym. block,
!          qualified diagonals excluded.

use Cholesky, only: iiBstR, IndRed, iQuAB, nnBstR, nQual, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Diag(*)
real(kind=wp), intent(out) :: Dmax(nSym)
integer(kind=iwp) :: iab, iQ, iRab, jRab, jRab1, jRab2, jSym, MxQ
real(kind=wp), allocatable :: ExQ(:)

MxQ = nQual(1)
do jSym=2,nSym
  MxQ = max(MxQ,nQual(jSym))
end do
call mma_allocate(ExQ,MxQ,Label='ExQ')

do jSym=1,nSym

  Dmax(jSym) = Zero
  if (nQual(jSym) < 1) cycle  ! next symm

  do iQ=1,nQual(jSym)
    iab = IndRed(iQuAB(iQ,jSym),2) ! addr in 1st red set
    ExQ(iQ) = Diag(iab)
    Diag(iab) = Zero
  end do

  jRab1 = iiBstr(jSym,2)+1
  jRab2 = jRab1+nnBstR(jSym,2)-1
  do jRab=jRab1,jRab2
    iRab = IndRed(jRab,2) ! addr in 1st red set
    Dmax(jSym) = max(Dmax(jSym),Diag(iRab))
  end do

  ! Restore the qualified
  do iQ=1,nQual(jSym)
    iab = IndRed(iQuAB(iQ,jSym),2) ! addr in 1st red set
    Diag(iab) = ExQ(iQ)
  end do

end do

call mma_deallocate(ExQ)

end subroutine Cho_MaxDX
