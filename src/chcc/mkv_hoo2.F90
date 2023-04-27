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

subroutine MkV_Hoo2(V2,V,dima,dimb,no)
! this routine do:
! Make AntiSymmetric integrals
! V2(i,a',b',j) <- 2 V(b'i|a'j) - V(b'j|a'i)

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimb, dima, no
real(kind=wp), intent(out) :: V2(no,dima,dimb,no)
real(kind=wp), intent(in) :: V(dimb,no,dima,no)
integer(kind=iwp) :: a, b, j

do j=1,no
  do b=1,dimb
    do a=1,dima
      V2(:,a,b,j) = Two*V(b,j,a,:)-V(b,:,a,j)
    end do
  end do
end do

return

end subroutine MkV_Hoo2
