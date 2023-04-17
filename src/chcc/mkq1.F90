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

subroutine MkQ1(V)
! Q1(a,j,k,l) <- V(aj,kl)

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, Q1
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: V(nv,no,nTri_Elem(no))
integer(kind=iwp) :: a, j, k, kl, l

call mma_allocate(Q1,nv,no,no,no,label='Q1')

kl = 0
do k=1,no
  do l=1,k
    kl = kl+1
    do j=1,no
      do a=1,nv
        Q1(a,j,k,l) = V(a,j,kl)
        Q1(a,j,l,k) = V(a,j,kl)
      end do
    end do
  end do
end do

return

end subroutine MkQ1
