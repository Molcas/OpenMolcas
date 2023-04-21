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

subroutine MkD_Q46(D,T2,T1a,T1b,dima,dimb,no)
! this routine does:
! Create D(a',j,b',i) [NB. D(a,b,i,j) =  t2(a,b,j,i)-T2(a,b,i,j)]
! from following available arrays (permuted as given):
! T2(a',j,b',i) [NB. T2abij = 1/2 t2abij + tai . tjb ]
! T1a(a',p)
! T1b(b',p)
!
! N.B. Kvajt odflaknute

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, no
real(kind=wp), intent(out) :: D(dima,no,dimb,no)
real(kind=wp), intent(in) :: T2(dima,no,dimb,no), T1a(dima,no), T1b(dimb,no)
integer(kind=iwp) :: b, i

do i=1,no
  do b=1,dimb
    D(:,:,b,i) = Two*(T2(:,i,b,:)-T1a(:,:)*T1b(b,i))-T2(:,:,b,i)
  end do
end do

return

end subroutine MkD_Q46
