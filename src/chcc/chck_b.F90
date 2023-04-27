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

subroutine Chck_B(BB,dima,dimb,dimbe,dimga,adda,addb,addbe,addga)
! this routine tests L2

use chcc_global, only: Bc
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimbe, dimga, adda, addb, addbe, addga
real(kind=wp), intent(inout) :: BB(dima,dimbe,dimb,dimga)
integer(kind=iwp) :: a, b, bad, be, ga, ntot
real(kind=wp) :: s

bad = 0
ntot = 0

do ga=addga+1,addga+dimga
  do b=addb+1,addb+dimb
    do be=addbe+1,addbe+dimbe
      do a=adda+1,adda+dima

        s = Bc(a,b,be,ga)

        if (abs(BB(a-adda,be-addbe,b-addb,ga-addga)-s) > 1.0e-10_wp) then
          bad = bad+1
          BB(a-adda,be-addbe,b-addb,ga-addga) = s
        end if
        ntot = ntot+1

      end do
    end do
  end do
end do

write(u6,*) ' B test ',bad,ntot

return

end subroutine Chck_B
