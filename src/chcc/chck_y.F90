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

subroutine Chck_Y(Y,dimbe,addbe,dimga,addga)
! check Y(be,u,ga,v)

use chcc_global, only: Kc, no, nv, T2c
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dimbe, addbe, dimga, addga
real(kind=wp), intent(in) :: Y(dimbe,no,dimga,no)
integer(kind=iwp) :: a, bad, be, ga, i, u, v
real(kind=wp) :: s

bad = 0
do v=1,no
  do ga=addga+1,addga+dimga
    do u=1,no
      do be=addbe+1,addbe+dimbe

        s = Zero
        do i=1,no
          do a=1,nv
            s = s+Kc(i,be,u,a)*T2c(ga,a,i,v)
          end do
        end do

        if (abs(Y(be-addbe,u,ga-addga,v)-s) > 1.0e-10_wp) bad = bad+1

      end do
    end do
  end do
end do

write(u6,*) ' Chck Y :',bad

return

end subroutine Chck_Y
