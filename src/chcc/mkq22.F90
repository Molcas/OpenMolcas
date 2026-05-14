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

subroutine MkQ22(V)
! Q22(a,b,k,l) <- V(ij,kl)

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, Q22
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: V(nTri_Elem(nv),nTri_Elem(no))
integer(kind=iwp) :: a, ab, k, kl

call mma_allocate(Q22,nv,nv,no,no,label='Q22')

kl = 0
do k=1,no
  ab = 0
  do a=1,nv
    Q22(a,1:a-1,k,1:k-1) = V(ab+1:ab+a-1,kl+1:kl+k-1)
    Q22(a,1:a-1,1:k,k) = V(ab+1:ab+a-1,kl+1:kl+k)
    Q22(1:a,a,k,1:k-1) = V(ab+1:ab+a,kl+1:kl+k-1)
    Q22(1:a,a,1:k,k) = V(ab+1:ab+a,kl+1:kl+k)
    ab = ab+a
  end do
  kl = kl+k
end do

return

end subroutine MkQ22
