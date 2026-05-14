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

subroutine TTMul(A,B,C,nRowA,nColA,nRowB)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nRowA, nColA, nRowB
real(kind=wp), intent(in) :: A(nRowA,nColA), B(nRowB,nRowA)
real(kind=wp), intent(out) :: C(nColA,nRowB)
integer(kind=iwp) :: i, Incj, jj, k, mCache, nCache_, njVec

nCache_ = (64/8)*1024
mCache = (nCache_*3)/4-nRowA*nColA
Incj = mCache/(nRowA+nColA)

! Sectioning of long index

do jj=1,nRowB,Incj
  njVec = min(Incj,nRowB-jj+1)

  do i=1,nColA
    ! Set target to zero
    C(i,jj:jj+njVec-1) = Zero
    do k=1,nRowA
      if (A(k,i) /= Zero) C(i,jj:jj+njVec-1) = C(i,jj:jj+njVec-1)+A(k,i)*B(jj:jj+njVec-1,k)
    end do
  end do

end do    ! End of sectioning

return

end subroutine TTMul
