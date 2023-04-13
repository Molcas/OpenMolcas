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
use Definitions, only: wp, iwp

implicit none
#include "chcc1.fh"
real(kind=wp) :: V(nTri_Elem(nv),nTri_Elem(no))
integer(kind=iwp) :: a, ab, b, k, kl, l

kl = 0
do k=1,no
  do l=1,k
    kl = kl+1
    ab = 0
    do a=1,nv
      do b=1,a
        ab = ab+1
        Q22(a,b,k,l) = V(ab,kl)
        Q22(a,b,l,k) = V(ab,kl)
        Q22(b,a,k,l) = V(ab,kl)
        Q22(b,a,l,k) = V(ab,kl)
      end do
    end do
  end do
end do

return

end subroutine MkQ22
