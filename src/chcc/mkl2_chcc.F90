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

subroutine MkL2_chcc(V)
! L2(m,a,b) <- V(m,ab)

use Index_Functions, only: nTri_Elem
use chcc_global, only: L2k, nc, nv
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: V(nc,nTri_Elem(nv))
integer(kind=iwp) :: a, ab, b, m

call mma_allocate(L2k,nc,nv,nv,label='L2k')

ab = 0
do a=1,nv
  do b=1,a
    ab = ab+1
    do m=1,nc
      L2k(m,a,b) = V(m,ab)
      L2k(m,b,a) = V(m,ab)
    end do
  end do
end do

return

end subroutine MkL2_chcc
