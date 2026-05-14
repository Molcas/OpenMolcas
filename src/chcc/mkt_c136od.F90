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

subroutine MkT_C136od(T2,X,Y,dimbe,dimga,no)
! this routine does:
! T2n(be',ga',u,v) <-
! C1                + 1/2 X(be',u,ga',v)
! C3                - 1/2 Y(be',u,ga',v)
! C6                - 1   Y(be',v,ga',u)
! for beGrp>gaGrp

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, dimga, no
real(kind=wp), intent(out) :: T2(dimbe,dimga,no,no)
real(kind=wp), intent(in) :: X(dimbe,no,dimga,no), Y(dimbe,no,dimga,no)
integer(kind=iwp) :: u, v

do v=1,no
  do u=1,no
    T2(:,:,u,v) = (X(:,u,:,v)-Y(:,u,:,v))*Half-Y(:,v,:,u)
  end do
end do

return

end subroutine MkT_C136od
