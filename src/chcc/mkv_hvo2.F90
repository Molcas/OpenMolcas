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

subroutine MkV_Hvo2(V,V2,dimbe,dima,no)
! this routine does:
! Make AntiSymmetric integrals
! V2(a',i,be',j) <- [2 V(be',j|a'i) - V(be',i|a'j)]

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, dima, no
real(kind=wp), intent(in) :: V(dimbe,no,dima,no)
real(kind=wp), intent(out) :: V2(dima,no,dimbe,no)
integer(kind=iwp) :: be, i, j

do j=1,no
  do be=1,dimbe
    do i=1,no
      V2(:,i,be,j) = Two*V(be,j,:,i)-V(be,i,:,j)
    end do
  end do
end do

return

end subroutine MkV_Hvo2
