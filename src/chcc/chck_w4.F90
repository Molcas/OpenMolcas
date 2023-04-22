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

subroutine Chck_W4(W4,dima,dimbe,dimb,dimga,adda,addbe,addb,addga)
! this routine tests W4

use chcc_global, only: no, Q3, Q4, T1c
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, dimbe, dimb, dimga, adda, addbe, addb, addga
real(kind=wp), intent(in) :: W4(dima,dimbe,dimb,dimga)
integer(kind=iwp) :: a, b, be, ga, i, bad, ntot
real(kind=wp) :: s

bad = 0
ntot = 0

do ga=1,dimga
  do b=1,dimb
    do be=1,dimbe
      do a=1,dima
        s = Q4(a+adda,be+addbe,b+addb,ga+addga)
        s = Zero
        do i=1,no
          s = s-Q3(a+adda,be+addbe,b+addb,i)*T1c(ga+addga,i)
          s = s-Q3(b+addb,ga+addga,a+adda,i)*T1c(be+addbe,i)
        end do
        !if (abs(W4(a,ga,b,be)-s) > 1.0e-10_wp) then
        if (abs(W4(a,be,b,ga)-s) > 1.0e-10_wp) then
          bad = bad+1
          !W4(a,be,b,ga) = s
        end if
        ntot = ntot+1
      end do
    end do
  end do
end do

write(u6,*) ' W4 test ',bad,ntot

return

end subroutine Chck_W4
