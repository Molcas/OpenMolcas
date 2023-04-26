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

subroutine ExH_T13(V,Hvv,dimbe,addbe,nv)
! this routine does:
! Extract V1(a,be') <- Hvv(a,be)
! Hvv(a,be') <- Fvv(a,be)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, addbe, nv
real(kind=wp), intent(out) :: V(nv,dimbe)
real(kind=wp), intent(in) :: Hvv(nv,nv)

V(:,:) = Hvv(:,addbe+1:addbe+dimbe)

return

end subroutine ExH_T13
