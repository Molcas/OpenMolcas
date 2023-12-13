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

subroutine AdH_Hvv2(H,Hvv,dima,dimbe,adda,addbe,nv)
! this routine does:
! Hvv(a,be) <<- - H(be',a')

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimbe, adda, addbe, nv
real(kind=wp), intent(in) :: H(dimbe,dima)
real(kind=wp), intent(inout) :: Hvv(nv,nv)
integer(kind=iwp) :: be

do be=1,dimbe
  Hvv(adda+1:adda+dima,addbe+be) = Hvv(adda+1:adda+dima,addbe+be)-H(be,:)
end do

return

end subroutine AdH_Hvv2
