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

subroutine ExH_X2(Gvv,H,dima,dimbe,nv,adda,addbe)
! this routine does:
! H(a',be') <- Gvv(be,a)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimbe, nv, adda, addbe
real(kind=wp), intent(in) :: Gvv(nv,nv)
real(kind=wp), intent(out) :: H(dima,dimbe)
integer(kind=iwp) :: be

do be=1,dimbe
  H(:,be) = Gvv(addbe+be,adda+1:adda+dima)
end do

return

end subroutine ExH_X2
