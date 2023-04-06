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

subroutine unpckhelp5(a,b,dimp,dimj,dime,jadd,noj,eadd,noe)
! this routine does:
! b(j,e) = a(pj,qe)-a(qe,pj) for symp=symq

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimj, dime, jadd, noj, eadd, noe
real(kind=wp), intent(in) :: a(dimp,dimp)
real(kind=wp), intent(inout) :: b(dimj,dime)
integer(kind=iwp) :: e

do e=1,noe
  b(1:noj,e) = a(jadd+1:jadd+noj,eadd+e)-a(eadd+e,jadd+1:jadd+noj)
end do

return

end subroutine unpckhelp5
