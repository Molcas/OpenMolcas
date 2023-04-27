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

subroutine MkT_CAlld(T2,X,Y,dimbe,no)
! this routine does:
! T2n(be',ga',u,v) <-
!    C1               + 1/2 X(be',u,ga',v)
!    C2               + 1/2 X(ga',v,be',u)
!    C3               - 1/2 Y(be',u,ga',v)
!    C4               - 1/2 Y(ga',v,be',u)
!    C5               - 1   Y(ga',u,be',v)
!    C6               - 1   Y(be',v,ga',u)
! for beGrp=gaGrp

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, no
real(kind=wp), intent(out) :: T2(dimbe,dimbe,no,no)
real(kind=wp), intent(in) :: X(dimbe,no,dimbe,no), Y(dimbe,no,dimbe,no)
integer(kind=iwp) :: be, ga, u, v

do v=1,no
  do u=1,no
    do ga=1,dimbe
      !do be=1,dimbe  - povodne, stacia iba cleny be>=ga
      T2(ga:dimbe,ga,u,v) = (X(ga:dimbe,u,ga,v)-Y(ga:dimbe,u,ga,v))*Half-Y(ga:dimbe,v,ga,u)
    end do
  end do
end do

do v=1,no
  do u=1,no
    do be=1,dimbe
      !do ga=1,dimbe  - povodne, stacia iba cleny be>=ga
      T2(be,1:be,u,v) = T2(be,1:be,u,v)+(X(1:be,v,be,u)-Y(1:be,v,be,u))*Half-Y(1:be,u,be,v)
    end do
  end do
end do

return

end subroutine MkT_CAlld
