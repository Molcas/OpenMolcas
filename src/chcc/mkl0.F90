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

subroutine MkL0(V)
! L0(m,i,j) <- V(m,ij)

use Index_Functions, only: nTri_Elem
use chcc_global, only: L0k, nc, no
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: V(nc,nTri_Elem(no))
integer(kind=iwp) :: i, ij

call mma_allocate(L0k,nc,no,no,label='L0k')

ij = 0
do i=1,no
  L0k(:,i,1:i-1) = V(:,ij+1:ij+i-1)
  L0k(:,1:i,i) = V(:,ij+1:ij+i)
  ij = ij+i
end do

return

end subroutine MkL0
