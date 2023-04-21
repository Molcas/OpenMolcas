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

subroutine Chck_L1(L1,dima,adda)
! this routine tests L1

use chcc_global, only: L1k, nc, no
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, adda
real(kind=wp), intent(inout) :: L1(nc,dima,no)
integer(kind=iwp) :: a, bad, i, m, ntot
real(kind=wp) :: s

bad = 0
ntot = 0

do i=1,no
  do a=adda+1,adda+dima
    do m=1,nc

      s = L1k(m,i,a)

      if (abs(L1(m,a-adda,i)-s) > 1.0e-10_wp) then
        bad = bad+1
        L1(m,a-adda,i) = s
      end if
      ntot = ntot+1

    end do
  end do
end do

write(u6,*) ' L1   ',bad,ntot

return

end subroutine Chck_L1
