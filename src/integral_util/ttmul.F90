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

implicit none
integer nRowA, nColA, nRowB
real*8 A(nRowA,nColA), B(nRowB,nRowA), C(nColA,nRowB)
integer nCache_, mCache, Incj, jj, njVec, i, j, k

nCache_ = (64/8)*1024
mCache = (nCache_*3)/4-nRowA*nColA
Incj = mCache/(nRowA+nColA)

! Sectioning of long index

do jj=1,nRowB,Incj
  njVec = min(Incj,nRowB-jj+1)

  do i=1,nColA
    ! Set target to zero
    do j=jj,jj+njVec-1
      C(i,j) = Zero
    end do
    do k=1,nRowA
      if (A(k,i) /= Zero) then
        do j=jj,jj+njVec-1
          C(i,j) = C(i,j)+A(k,i)*B(j,k)
        end do
      end if
    end do
  end do

end do    ! End of sectioning

return

end subroutine TTMul
