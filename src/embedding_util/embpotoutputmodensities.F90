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
! Copyright (C) Thomas Dresselhaus                                     *
!***********************************************************************

! Actually the incoming densities are in an AO basis and treated as such...
subroutine embPotOutputMODensities(nAtoms,nSym,densInact,densAct,nBasPerSym,nBasTotSquare)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nSym, nBasPerSym(nSym), nBasTotSquare
real(kind=wp), intent(in) :: densInact(nBasTotSquare), densAct(nBasTotSquare)
integer(kind=iwp) :: i, iCnt, iCnt2, iCol, iDensDimPck, iRow
real(kind=wp), allocatable :: totDens(:), totDensP(:)

iDensDimPck = 0
do i=1,nSym
  iDensDimPck = iDensDimPck+((nBasPerSym(i)*(nBasPerSym(i)+1))/2)
end do

call mma_allocate(totDens,nBasTotSquare,label='TotD')
call mma_allocate(totDensP,iDensDimPck,label='TotDp')

call dcopy_(nBasTotSquare,densInact,1,totDens,1)
call daxpy_(nBasTotSquare,One,densAct,1,totDens,1)
!Pack density
iCnt = 0
iCnt2 = 0
do i=1,nSym
  do iRow=1,nBasPerSym(i)
    do iCol=1,iRow
      iCnt = iCnt+1
      if (iRow == iCol) then
        totDensP(iCnt) = totDens((iRow-1)*nBasPerSym(i)+iCol+iCnt2)
      else
        totDensP(iCnt) = Two*totDens((iRow-1)*nBasPerSym(i)+iCol+iCnt2)
      end if
    end do
  end do
  iCnt2 = iCnt2+nBasPerSym(i)*nBasPerSym(i)
end do

call embPotOutput(nAtoms,totDensP)

call mma_deallocate(totDens)
call mma_deallocate(totDensP)

end subroutine embPotOutputMODensities
