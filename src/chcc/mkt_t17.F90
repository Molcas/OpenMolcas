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

subroutine MkT_T17(Ta,T,dimb,dimbe,no)
! this routine does:
! Ta(i,b',be',u) <- 2 T(be',b',u,i) - T(be',b',i,u)
!
! N.B. Qvajt odflaknute

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimb, dimbe, no
real(kind=wp), intent(out) :: Ta(no,dimb,dimbe,no)
real(kind=wp), intent(in) :: T(dimbe,dimb,no,no)
integer(kind=iwp) :: b, be, u

do u=1,no
  do be=1,dimbe
    do b=1,dimb
      Ta(:,b,be,u) = Two*T(be,b,u,:)-T(be,b,:,u)
    end do
  end do
end do

return

end subroutine MkT_T17
