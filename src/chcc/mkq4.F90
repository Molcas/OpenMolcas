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

subroutine MkQ4(V)
! Q4(a,b,c,d) <- V(ab,cd)

use Index_Functions, only: nTri_Elem
use chcc_global, only: nv, Q4
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: V(nTri_Elem(nv),nTri_Elem(nv))
integer(kind=iwp) :: a, ab, c, cd

call mma_allocate(Q4,nv,nv,nv,nv,label='Q4')

cd = 0
do c=1,nv
  ab = 0
  do a=1,nv
    Q4(a,1:a-1,c,1:c-1) = V(ab+1:ab+a-1,cd+1:cd+c-1)
    Q4(a,1:a-1,1:c,c) = V(ab+1:ab+a-1,cd+1:cd+c)
    Q4(1:a,a,c,1:c-1) = V(ab+1:ab+a,cd+1:cd+c-1)
    Q4(1:a,a,1:c,c) = V(ab+1:ab+a,cd+1:cd+c)
    ab = ab+a
  end do
  cd = cd+c
end do

return

end subroutine MkQ4
