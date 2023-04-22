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

subroutine MkT20u(T2,V,oe,dima,dimb,adda,addb,no)
! this routine does:
! T2(a',b',i,j)=(a'i|b'j)/Dija'b'
!  for aGrp/=bGrp
! N.B.2 qvajt odflaknute

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, adda, addb, no
real(kind=wp), intent(out) :: T2(dima,dimb,no,no)
real(kind=wp), intent(in) :: V(dima,no,dimb,no), oe(*)
integer(kind=iwp) :: b, i, j

do j=1,no
  do i=1,no
    do b=1,dimb
      T2(:,b,i,j) = V(:,i,b,j)/(oe(i)+oe(j)-oe(adda+1:adda+dima)-oe(addb+b))
    end do
  end do
end do

return

end subroutine MkT20u
