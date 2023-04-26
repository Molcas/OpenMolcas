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

subroutine MkQ3(V)
! Q1(a,b,c,l) <- V(ab,cl)

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, Q3
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: V(nTri_Elem(nv),no,nv)
integer(kind=iwp) :: a, ab, l

call mma_allocate(Q3,nv,nv,nv,no,label='Q3')

do l=1,no
  ab = 0
  do a=1,nv
    Q3(a,1:a-1,:,l) = V(ab+1:ab+a-1,l,:)
    Q3(1:a,a,:,l) = V(ab+1:ab+a,l,:)
    ab = ab+a
  end do
end do

return

end subroutine MkQ3
