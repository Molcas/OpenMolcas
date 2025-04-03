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

subroutine ReLoad(A,idsym,NBAS1,NBAS2)

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nDens2
use input_mclr, only: nSym
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8 A(*)
integer idSym
integer nbas2(nsym), nbas1(nsym)
real*8, allocatable :: ATemp(:)
integer iS, jS, j

call mma_allocate(ATemp,ndens2,Label='ATemp')

do iS=1,nsym
  js = Mul(is,idsym)
  if (min(nbas1(is),nbas2(is)) < 1) cycle
  do j=0,min(nbas2(js),nbas1(js))-1
    call dcopy_(min(nbas1(is),nbas2(is)),A(ipMat(is,js)+j*nbas1(is)),1,ATemp(ipmat(is,js)+j*nbas2(is)),1)
  end do
end do
call dcopy_(ndens2,ATemp,1,A,1)
call mma_deallocate(ATemp)

end subroutine ReLoad
