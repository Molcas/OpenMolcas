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
use MCLR_Data, only: ipMat, nDens
use input_mclr, only: nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: A(*)
integer(kind=iwp), intent(in) :: idSym, nbas2(nsym), nbas1(nsym)
integer(kind=iwp) :: iS, j, jS, m, n1, n2
real(kind=wp), allocatable :: ATemp(:)

call mma_allocate(ATemp,nDens,Label='ATemp')

do iS=1,nsym
  js = Mul(is,idsym)
  m = min(nbas1(is),nbas2(is))
  do j=1,min(nbas2(js),nbas1(js))
    n1 = ipMat(is,js)+(j-1)*nbas1(is)
    n2 = ipMat(is,js)+(j-1)*nbas2(is)
    ATemp(n2:n2+m-1) = A(n1:n1+m-1)
  end do
end do
A(1:nDens) = ATemp(:)
call mma_deallocate(ATemp)

end subroutine ReLoad
