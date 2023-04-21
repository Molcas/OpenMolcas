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

subroutine AdV_G2(G2,V,nv,dimbe,dima,no,addbe,adda,fact)
! this routine does:
! G2(be',a') <- fact sum(i) V(be',i,i,a')
!
! N.B. Kvajt odflaknute

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nv, dimbe, dima, no, addbe, adda
real(kind=wp), intent(inout) :: G2(nv,nv)
real(kind=wp), intent(in) :: V(dimbe,no,no,dima), fact
integer(kind=iwp) :: i

do i=1,no
  G2(addbe+1:addbe+dimbe,adda+1:adda+dima) = G2(addbe+1:addbe+dimbe,adda+1:adda+dima)+fact*V(:,i,i,:)
end do

return

end subroutine AdV_G2
