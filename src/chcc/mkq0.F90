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

subroutine MkQ0(V)
! Q0(i,j,k,l) <- V(ij,kl)

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, Q0
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: V(nTri_Elem(no),nTri_Elem(no))
integer(kind=iwp) :: i, ij, k, kl

call mma_allocate(Q0,no,no,no,no,label='Q0')

kl = 0
do k=1,no
  ij = 0
  do i=1,no
    Q0(i,1:i-1,k,1:k-1) = V(ij+1:ij+i-1,kl+1:kl+k-1)
    Q0(i,1:i-1,1:k,k) = V(ij+1:ij+i-1,kl+1:kl+k)
    Q0(1:i,i,k,1:k-1) = V(ij+1:ij+i,kl+1:kl+k-1)
    Q0(1:i,i,1:k,k) = V(ij+1:ij+i,kl+1:kl+k)
    ij = ij+i
  end do
  kl = kl+k
end do

return

end subroutine MkQ0
