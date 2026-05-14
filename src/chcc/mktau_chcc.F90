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

subroutine MkTau_chcc(T2,T11,T12,dima,dimb,no,f1,f2)
! this routine does:
! T2(a',b',i,j) = f1 . T2(a',b',i,j) + f2 . T11(a',i) . T12(b',j)
!
! N.B. Kvajt odflaknute

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, no
real(kind=wp), intent(inout) :: T2(dima,dimb,no,no)
real(kind=wp), intent(in) :: T11(dima,no), T12(dimb,no), f1, f2
integer(kind=iwp) :: b, j

do j=1,no
  do b=1,dimb
    T2(:,b,:,j) = f1*T2(:,b,:,j)+f2*T11(:,:)*T12(b,j)
  end do
end do

return

end subroutine MkTau_chcc
