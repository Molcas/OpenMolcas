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

subroutine MkT_T15(Tp,T2,T11,T12,dimbe,dima,no)
! this routine does:
! Tp(be',u,a',i) <- 2 t2(be,a,u,i)
!                 - t2(be,a,i,u)+t12(a,u).t11(be,i)/2
! N.B. mozno sa to da este vylepsit

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimbe, dima, no
real(kind=wp), intent(out) :: Tp(dimbe,no,dima,no)
real(kind=wp), intent(in) :: T2(dimbe,dima,no,no), T11(dimbe,no), T12(dima,no)
integer(kind=iwp) :: a, u

do a=1,dima
  do u=1,no
    Tp(:,u,a,:) = Two*T2(:,a,u,:)-T2(:,a,:,u)+T11(:,:)*T12(a,u)
  end do
end do

return

end subroutine MkT_T15
