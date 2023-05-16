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

subroutine MkT_red(T2red,T2full,dimbe,no)
! this routine produces reduced set of ampitudes for storing:
! T2red(be'ga',u,v) <- T2full(be',ga',u,v)
! for beGrp=gaGrp

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, no
real(kind=wp), intent(out) :: T2red(nTri_Elem(dimbe),no,no)
real(kind=wp), intent(in) :: T2full(dimbe,dimbe,no,no)
integer(kind=iwp) :: be, bega, i, j

do j=1,no
  do i=1,no

    bega = 0
    do be=1,dimbe

      T2red(bega+1:bega+be,i,j) = T2Full(be,1:be,i,j)

      bega = bega+be
    end do

  end do
end do

return

end subroutine MkT_red
