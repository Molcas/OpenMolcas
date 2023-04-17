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
real(kind=wp) :: V(nTri_Elem(no),nTri_Elem(no))
integer(kind=iwp) :: i, ij, j, k, kl, l

call mma_allocate(Q0,no,no,no,no,label='Q0')

kl = 0
do k=1,no
  do l=1,k
    kl = kl+1
    ij = 0
    do i=1,no
      do j=1,i
        ij = ij+1
        Q0(i,j,k,l) = V(ij,kl)
        Q0(i,j,l,k) = V(ij,kl)
        Q0(j,i,k,l) = V(ij,kl)
        Q0(j,i,l,k) = V(ij,kl)
      end do
    end do
  end do
end do

return

end subroutine MkQ0
