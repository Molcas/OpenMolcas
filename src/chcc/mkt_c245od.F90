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

subroutine MkT_C245od(T2,X,Y,dimbe,dimga,no)
! this routine does:
! T2n(be',ga',u,v) <<-
! C2                + 1/2 X(ga',v,be',u)
! C4                - 1/2 Y(ga',v,be',u)
! C5                - 1   Y(ga',u,be',v)
! for beGrp>gaGrp

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, dimga, no
real(kind=wp), intent(inout) :: T2(dimbe,dimga,no,no)
real(kind=wp), intent(in) :: X(dimga,no,dimbe,no), Y(dimga,no,dimbe,no)
integer(kind=iwp) :: ga, u, v

do v=1,no
  do u=1,no
    do ga=1,dimga
      T2(:,ga,u,v) = T2(:,ga,u,v)+(X(ga,v,:,u)-Y(ga,v,:,u))*Half-Y(ga,u,:,v)
    end do
  end do
end do

return

end subroutine MkT_C245od
