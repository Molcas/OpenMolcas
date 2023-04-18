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

subroutine MkV_Goo3(V,V2,dima,no)
! this routine does:
! Make AntiSymmetric integrals
! V2(a',j,i,u) <- 2 V(a,j|iu) - (a,i|ju)

use Index_Functions, only: nTri_Elem
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, no
real(kind=wp) :: V(dima,no,nTri_Elem(no)), V2(dima,no,no,no)
integer(kind=iwp) :: i, iu, j, ju, u

do u=1,no

  do i=1,no
    if (i > u) then
      iu = (i-1)*i/2+u
    else
      iu = (u-1)*u/2+i
    end if

    do j=1,no

      if (j > u) then
        ju = (j-1)*j/2+u
      else
        ju = (u-1)*u/2+j
      end if

      V2(:,j,i,u) = Two*V(:,j,iu)-V(:,i,ju)

    end do
  end do
end do

return

end subroutine MkV_Goo3
