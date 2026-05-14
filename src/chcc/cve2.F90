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

subroutine CVE2(V,oe,dima,dimb,adda,addb,no,e2,e2os)
! this routine does:
! e2  =e2   + sum(a',b',i,j)[(a'i|b'j).{2(a'i|b'j)-(a'j|b'i)}]/Dija'b'
! e2os=e2os + sum(a',b',i,j)[(a'i|b'j).{ (a'i|b'j)          }]/Dija'b'
! N.B. Cierny Vypocet E2
! N.B.2 qvajt odflaknute

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, adda, addb, no
real(kind=wp), intent(in) :: V(dima,no,dimb,no), oe(*)
real(kind=wp), intent(inout) :: e2, e2os
integer(kind=iwp) :: a, b, i, j
real(kind=wp) :: dijab

do j=1,no
  do b=1,dimb
    do i=1,no
      do a=1,dima
        dijab = oe(i)+oe(j)-oe(adda+a)-oe(addb+b)
        e2 = e2+V(a,i,b,j)*(Two*V(a,i,b,j)-V(a,j,b,i))/dijab
        e2os = e2os+V(a,i,b,j)*(V(a,i,b,j))/dijab
      end do
    end do
  end do
end do

return

end subroutine CVE2
