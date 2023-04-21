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

subroutine MkT_QK42(T2,T1a,T1b,dima,dimb,no,f1,f2)
! this routine does:
! T2(a',b',i,j) <- f1 . T2(a',b',i,j) + f2 . T1a(a,i) . T1b(b,j)
!
! N.B. Kvajt odflaknute

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, no
real(kind=wp), intent(inout) :: T2(dima,dimb,no,no)
real(kind=wp), intent(in) :: T1a(dima,no), T1b(dimb,no), f1, f2
integer(kind=iwp) :: b, j

do j=1,no
  do b=1,dimb
    T2(:,b,:,j) = f1*T2(:,b,:,j)+f2*T1a(:,:)*T1b(j,b)
  end do
end do

return

end subroutine MkT_QK42
