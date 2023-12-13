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

subroutine Chck_L2(L2,dima,dimb,adda,addb)
! this routine tests L2

use chcc_global, only: L1k, L2k, nc, no, T1c
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, adda, addb
real(kind=wp), intent(inout) :: L2(nc,dima,dimb)
integer(kind=iwp) :: a, b, bad, i, m, ntot
real(kind=wp) :: s

bad = 0
ntot = 0

do b=addb+1,addb+dimb
  do a=adda+1,adda+dima
    do m=1,nc

      s = L2k(m,a,b)
      do i=1,no
        s = s-L1k(m,i,a)*T1c(b,i)
      end do

      if (abs(L2(m,a-adda,b-addb)-s) > 1.0e-10_wp) then
        bad = bad+1
        L2(m,a-adda,b-addb) = s
      end if
      ntot = ntot+1

    end do
  end do
end do

write(u6,*) ' L2   ',bad,ntot

return

end subroutine Chck_L2
