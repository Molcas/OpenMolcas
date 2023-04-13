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

subroutine ExV_X41(Vp,V,dimab,no)
! this routine does:
! Vp(a_b,ij) <- V(a_b,i,j) for i>=j

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimab, no
real(kind=wp) :: Vp(dimab,nTri_Elem(no)), V(dimab,no,no)
integer(kind=iwp) :: ab, i, ij, j

ij = 0
do i=1,no
  do j=1,i
    ij = ij+1
    do ab=1,dimab
      Vp(ab,ij) = V(ab,i,j)
    end do
  end do
end do

return

end subroutine ExV_X41
